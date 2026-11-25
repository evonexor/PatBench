"""Target protein 3D structures for the PatBench ranking tasks.

    python -I fetch_structures.py <release_dir> <stage> [--only ACC,ACC] [--outdir structures] [--jobs 4]

stages, run in this order (each one reads the cache of the earlier ones, so a stage can be re-run):

  meta       UniProt names for every target + the PatBench usage counts      -> target_name_map.tsv
  entries    PDBe best_structures (every PDB entry mapped to the accession) and the AlphaFold DB record
  gql        RCSB data API (GraphQL): method, resolution, date, entity sequences and UniProt alignments,
             ligands (CCD id, formula, InChIKey) for every PDB entry found
  select     per accession and per sequence variant (wild type, mutant, mature chain) choose the PDB entries
  download   the chosen mmCIF files, the SIFTS mapping of each, the AlphaFold models
  process    clean structures, ligands, pockets, residue maps, the four tables
  groups     the group -> structure table

Sources (public, read only): rest.uniprot.org, www.ebi.ac.uk/pdbe (best_structures and SIFTS mappings),
data.rcsb.org (GraphQL), files.rcsb.org (mmCIF), alphafold.ebi.ac.uk (API and PDB files).
Everything that is downloaded is treated as data: it is parsed with gemmi and never executed.
The date and the API answers behind a run are recorded in <outdir>/fetch_log.txt and <outdir>/api_cache/.
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import argparse
import collections
import concurrent.futures as cf
import datetime
import difflib
import gzip
import io
import json
import os
import random
import re
import shutil
import sys
import time
import urllib.error
import urllib.request

sys.path.insert(0, _ws + "/struct-deps-20261030")  # gemmi, installed with pip --target
sys.path.insert(0, _ws + "/outputs/uspto-sar-demo-20261003/runtime/deps")
import gemmi  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

UA = "PatBench-structure-fetch/1.0 (research use; public read-only APIs)"
RANK_TASKS = ["rank_ic50", "rank_ki_kd", "rank_ec50", "rank_highcontrast", "rank_complex_mutant"]
MAX_GENERIC_PICKS = 3        # per accession and variant: best_holo, best_coverage, best_resolution
MAX_MUTANT_PICKS = 3         # extra picks for a variant that carries substitutions
HOLO_RES_MAX = 3.0
CONTACT = 4.5                # a ligand counts as bound to the chain when an atom is this close (A)
POCKET = 6.0                 # pocket = chain residues with an atom this close to the ligand (A)
MATCH_IDENTITY = 0.98
GQL_BATCH = 40
TEST_REL = None
MAX_BYTES = int(1.8e9)
OLD_PDB_METHODS = ("X-RAY", "ELECTRON MICROSCOPY", "NMR")
# A target that has no PDB entry of its own but whose protein is a close homolog of one that has: the experimental
# structures of the homolog stand in for it. {target: (accession whose entries are used, position shift)};
# a residue at position p of the source is position p - shift of the target. The identity to the target is measured
# and reported, and such a structure never counts as matching the group's sequence unless it is >= 98% identical.
#   O91080 (HIV-1 group N isolate YBF30) reverse transcriptase 598-1157  <-  P03366 (group M, isolate BH10) 600-1159,
#   a gap-free alignment with 485 of 560 residues identical (86.6%).
PROXY = {"O91080": ("P03366", 2)}


def proxy(acc):
    return PROXY.get(acc, (acc, 0))


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# ------------------------------------------------------------------ http with cache
def http(url, post=None, retries=6, timeout=120):
    last = None
    for i in range(retries):
        try:
            hdr = {"User-Agent": UA}
            if post is not None:
                hdr["Content-Type"] = "application/json"
            with urllib.request.urlopen(urllib.request.Request(url, data=post, headers=hdr), timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if e.code == 400:
                raise
            last = e
        except Exception as e:  # noqa: BLE001 - network errors of every kind are retried
            last = e
        time.sleep(min(60, 2 ** i) + random.random())
    raise RuntimeError(f"{url}: {last}")


def atomic_write(path, data):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = f"{path}.part{os.getpid()}"
    with open(tmp, "wb") as f:
        f.write(data)
    os.replace(tmp, path)


def cached_json(path, url, post=None):
    if os.path.exists(path):
        return json.load(open(path))
    b = http(url, post)
    time.sleep(0.1)
    atomic_write(path, b"null" if b is None else b)
    return None if b is None else json.loads(b)


def pool_map(fn, items, jobs):
    out = {}
    with cf.ThreadPoolExecutor(max_workers=jobs) as ex:
        futs = {ex.submit(fn, it): it for it in items}
        for n, f in enumerate(cf.as_completed(futs), 1):
            it = futs[f]
            try:
                out[it] = f.result()
            except Exception as e:  # noqa: BLE001
                out[it] = e
            if n % 50 == 0:
                log(f"  {n}/{len(items)}")
    return out


# ------------------------------------------------------------------ ligand classes
ADDITIVES = set("""GOL EDO PEG PGE PG4 P6G 1PE 2PE 12P PE3 PE4 PE5 PE6 PE7 PE8 PE9 PEU PG0 PG5 PG6 PGO PGR MPD MRD MRY 15P 1BO
BU1 BU2 BU3 B3P BTB TRS MES EPE HEZ CXS CIT FLC TLA MLA MLI SO4 PO4 ACT ACY FMT EOH IPA MOH DMS BME DTT DTV DTD TAM IMD
NH4 TFA BCT GBL SPD SPM PEA PGA OCT NHE MPO POP PPV 2HP CAC ACN PIN SGM DIO P4C TBU 144 EGL PG3 PGT PE2 AZI SCN DOD
BCN CO3 TAR MAE MAL OXL SIN GLY2 TMA ETA NCO IPH CPS CHT CPT""".split())
COFACTORS = set("""ATP ADP AMP ANP ACP APC AGS ATG GTP GDP GMP GNP GCP GSP NAD NAP NDP NAI NAJ NADH FAD FMN HEM HEA HEC HEB SAM SAH
SFG COA ACO PLP PMP TPP THF FOL H4B BH4 MGD MTE MOS BCL CLA CHL PQQ F3S SF4 FES SRM ADN 5GP UDP UTP CTP CDP IMP XMP
NMN AP5 AP2 A3P PAP PPS TDP TPQ RBF 2AM""".split())
LIPIDS = set("""BOG LMT DDM OLC OLA OLB PLM MYR STE LDA C8E 6PL PC1 PCW PEF PGW POV DPP LPP HTG GLC12 CLR CHD CHO DGA DGD
12D D10 DD9 UND NNN HP6 LHG LMG LPC MC3 PGV PX4 PTY LAU CPL DMU UMQ MPG OCT CE9 DAO LI1 LIO TWT
NS3 D12 TRD FLC12 SDS CM5 BNG TCE""".split())
LIPID_WORDS = ("PHOSPHATIDYL", "OLEOYL", "PALMITOYL", "MYRISTOYL", "STEAROYL", "GLYCEROL", "POLYETHYLENE", "CHOLEST",
               "DODECYL", "OCTYL", "DECYL", "NONYL", "LAURYL", "GLUCOSIDE", "MALTOSIDE", "DIGITONIN", "TRITON",
               "DETERGENT", "SPHINGO", "CERAMIDE", "LIPID", "DIACYL", "ACYL-COA")


def formula_counts(formula):
    out = collections.Counter()
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula or ""):
        out[el] += int(n) if n else 1
    return out


def classify_ligand(cid, name, formula, ctype):
    """(class, heavy atoms). class: druglike, additive, cofactor, sugar, lipid_detergent, inorganic, small_organic."""
    cnt = formula_counts(formula)
    heavy = sum(v for k, v in cnt.items() if k not in ("H", "D"))
    nm = (name or "").upper()
    if cnt.get("C", 0) == 0:
        return "inorganic", heavy
    if cid in ADDITIVES:
        return "additive", heavy
    if cid in COFACTORS:
        return "cofactor", heavy
    if "saccharide" in (ctype or "").lower():
        return "sugar", heavy
    if cid in LIPIDS or any(w in nm for w in LIPID_WORDS):
        return "lipid_detergent", heavy
    if heavy < 10:
        return "small_organic", heavy
    return "druglike", heavy


# ------------------------------------------------------------------ variants of a target
def parse_variant(sid, acc, seq, kind, n_groups):
    off, subs, region = 0, [], None
    m = re.fullmatch(re.escape(acc) + r":(\d+)-(\d+)", sid)
    if m:
        region = (int(m.group(1)), int(m.group(2)))
        off = region[0] - 1
    m = re.fullmatch(re.escape(acc) + r":([A-Z]\d+[A-Z](?:\+[A-Z]\d+[A-Z])*)", sid)
    if m:
        subs = [(w, int(p), t) for w, p, t in re.findall(r"([A-Z])(\d+)([A-Z])", m.group(1))]
    return dict(sid=sid, kind=kind, seq=seq, off=off, subs=subs, region=region, n_groups=n_groups)


def load_variants(D, meta):
    ts = pd.read_parquet(f"{D}/target_sequences.parquet")
    out = collections.defaultdict(list)
    for (acc, sid), g in ts.groupby(["uniprot_accession", "sequence_id"], sort=True):
        out[acc].append(parse_variant(sid, acc, g.sequence.iloc[0], g.sequence_kind.iloc[0], len(g)))
    for acc in meta.index:
        if acc not in out:   # regression-only target: the UniProt sequence, no group
            out[acc].append(parse_variant(acc, acc, meta.sequence[acc], "regression_only", 0))
    for acc in out:
        out[acc].sort(key=lambda v: -v["n_groups"])
    return out


def variant_metrics(pos_aa, var):
    """pos_aa: UniProt position -> one-letter residue of a resolved structure residue."""
    n = len(var["seq"])
    cov = mis = 0
    for p, a in pos_aa.items():
        i = p - var["off"]
        if 1 <= i <= n:
            cov += 1
            mis += var["seq"][i - 1] != a
    detail, state = [], []
    for w, p, t in var["subs"]:
        a = pos_aa.get(p)
        s = "present" if a == t else ("absent:" + a if a else "unresolved")
        detail.append(f"{w}{p}{t}={s}")
        state.append(s)
    if not var["subs"]:
        mut = None
    elif all(s == "present" for s in state):
        mut = True
    elif any(s.startswith("absent") for s in state):
        mut = False
    else:
        mut = None
    return dict(n_mapped=cov, n_mismatch=int(mis), identity=(cov - mis) / cov if cov else None,
                coverage=cov / n if n else 0.0, mutation_in_structure=mut, mutation_detail=";".join(detail))


# ------------------------------------------------------------------ stage meta
def stage_meta(D, S, only):
    ts = pd.read_parquet(f"{D}/target_sequences.parquet")
    reg = pd.read_parquet(f"{D}/tasks/regression_affinity.parquet", columns=["uniprot_accession", "split_time"])
    accs = sorted(set(ts.uniprot_accession))
    extra = sorted(set(reg.uniprot_accession) - set(accs))
    allacc = accs + extra
    log(f"meta: {len(accs)} rank-task targets, {len(extra)} regression-only targets")
    rows = []
    for i in range(0, len(allacc), 100):
        chunk = allacc[i:i + 100]
        url = ("https://rest.uniprot.org/uniprotkb/accessions?accessions=" + ",".join(chunk) +
               "&fields=accession,id,reviewed,gene_primary,protein_name,organism_name,organism_id,length,sequence&format=tsv")
        b = http(url)
        t = pd.read_csv(io.BytesIO(b), sep="\t", dtype=str, keep_default_na=False)
        rows.append(t)
        time.sleep(0.3)
    u = pd.concat(rows).rename(columns={"Entry": "uniprot_accession", "Entry Name": "uniprot_entry",
                                        "Reviewed": "reviewed", "Gene Names (primary)": "gene_symbol",
                                        "Protein names": "protein_name", "Organism": "organism",
                                        "Organism (ID)": "taxid", "Length": "uniprot_length", "Sequence": "sequence"})
    missing = sorted(set(allacc) - set(u.uniprot_accession))
    log("meta: not returned by UniProt (merged or obsolete):", missing)
    u = u.set_index("uniprot_accession")
    audit = pd.read_parquet(f"{D}/target_audit.parquet", columns=["uniprot_accession", "target_class"]).drop_duplicates("uniprot_accession").set_index("uniprot_accession")

    def agg(g, col):
        vals = sorted({str(x) for x in g[col] if str(x) not in ("", "None", "nan")})
        return " | ".join(vals)

    nm = pd.DataFrame(index=allacc)
    nm.index.name = "uniprot_accession"
    for c in ["uniprot_entry", "reviewed", "gene_symbol", "protein_name", "organism", "taxid", "uniprot_length"]:
        nm[c] = [u[c].get(a, "") for a in allacc]
    g = ts.groupby("uniprot_accession")
    nm["patbench_target_symbols"] = [agg(g.get_group(a), "target_symbols") if a in g.groups else "" for a in allacc]
    nm["patbench_mapped_identity"] = [agg(g.get_group(a), "mapped_identity") if a in g.groups else "" for a in allacc]
    nm["patbench_target_class"] = [audit.target_class.get(a, "") for a in allacc]
    nm["sequence_kinds"] = [agg(g.get_group(a), "sequence_kind") if a in g.groups else "regression_only" for a in allacc]
    nm["sequence_ids"] = [agg(g.get_group(a), "sequence_id") if a in g.groups else "" for a in allacc]
    ng, npat = {}, {}
    for t in RANK_TASKS:
        d = pd.read_parquet(f"{D}/tasks/{t}.parquet", columns=["group_id", "patent_core", "uniprot_accession"])
        ng[t] = d.groupby("uniprot_accession").group_id.nunique()
        npat[t] = d.groupby("uniprot_accession").patent_core.nunique()
        nm[f"n_groups_{t}"] = [int(ng[t].get(a, 0)) for a in allacc]
        nm[f"n_patents_{t}"] = [int(npat[t].get(a, 0)) for a in allacc]
    nm["n_groups_all_rank_tasks"] = [int(g.get_group(a).group_id.nunique()) if a in g.groups else 0 for a in allacc]
    alld = pd.concat([pd.read_parquet(f"{D}/tasks/{t}.parquet", columns=["patent_core", "uniprot_accession", "group_id"]) for t in RANK_TASKS])
    pat = alld.groupby("uniprot_accession").patent_core.nunique()
    nm["n_patents_all_rank_tasks"] = [int(pat.get(a, 0)) for a in allacc]
    rc = reg.groupby("uniprot_accession").size()
    nm["in_regression"] = [a in rc.index for a in allacc]
    nm["n_rows_regression"] = [int(rc.get(a, 0)) for a in allacc]
    nm["scope"] = ["rank_tasks" if a in set(accs) else "regression_only" for a in allacc]
    nm["sequence"] = [u.sequence.get(a, "") for a in allacc]
    os.makedirs(S, exist_ok=True)
    nm.drop(columns=["sequence"]).to_csv(f"{S}/target_name_map.tsv", sep="\t", encoding="utf-8-sig")
    os.makedirs(f"{S}/api_cache", exist_ok=True)
    nm.to_pickle(f"{S}/api_cache/meta.pkl")
    log("meta: wrote target_name_map.tsv,", len(nm), "targets;", "isoform-style ids:", [a for a in allacc if "-" in a])


def load_meta(S):
    return pd.read_pickle(f"{S}/api_cache/meta.pkl")


# ------------------------------------------------------------------ stage entries
def stage_entries(D, S, only, jobs):
    meta = load_meta(S)
    accs = [a for a in meta.index if not only or a in only]

    def one(acc):
        b = cached_json(f"{S}/api_cache/pdbe_best/{acc}.json", f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{acc}")
        a = cached_json(f"{S}/api_cache/alphafold/{acc}.json", f"https://alphafold.ebi.ac.uk/api/prediction/{acc}")
        n_af = len([r for r in (a or []) if re.fullmatch(re.escape(f"AF-{acc}") + r"-F\d+", r["modelEntityId"])])
        return (0 if b is None else len(b.get(acc, [])), n_af)

    res = pool_map(one, accs, jobs)
    bad = {k: v for k, v in res.items() if isinstance(v, Exception)}
    log(f"entries: {len(accs)} accessions, errors {len(bad)}", list(bad.items())[:5])
    log("entries: accessions without any PDB row:", [k for k, v in res.items() if not isinstance(v, Exception) and v[0] == 0][:80])
    log("entries: accessions without AlphaFold:", [k for k, v in res.items() if not isinstance(v, Exception) and v[1] == 0])


def pdbe_rows(S, acc):
    b = json.load(open(f"{S}/api_cache/pdbe_best/{acc}.json"))
    return [] if b is None else b.get(acc, [])


# ------------------------------------------------------------------ stage gql
GQL = """query($ids:[String!]!){ entries(entry_ids:$ids){ rcsb_id
 exptl{method}
 rcsb_accession_info{initial_release_date}
 rcsb_entry_info{resolution_combined}
 struct{title}
 polymer_entities{ rcsb_polymer_entity_container_identifiers{entity_id auth_asym_ids asym_ids}
   entity_poly{pdbx_seq_one_letter_code_can rcsb_entity_polymer_type}
   rcsb_polymer_entity{pdbx_description}
   rcsb_polymer_entity_align{reference_database_name reference_database_accession aligned_regions{entity_beg_seq_id ref_beg_seq_id length}} }
 nonpolymer_entities{ rcsb_nonpolymer_entity_container_identifiers{auth_asym_ids entity_id}
   nonpolymer_comp{ chem_comp{id name formula type} rcsb_chem_comp_descriptor{InChIKey} } }
}}"""


def compact_entry(e):
    res = (e.get("rcsb_entry_info") or {}).get("resolution_combined") or [None]
    rec = dict(id=e["rcsb_id"], method=[x["method"] for x in e.get("exptl") or []], res=res[0],
               date=((e.get("rcsb_accession_info") or {}).get("initial_release_date") or "")[:10],
               title=(e.get("struct") or {}).get("title") or "", entities=[], ligands=[])
    for pe in e.get("polymer_entities") or []:
        ci = pe["rcsb_polymer_entity_container_identifiers"]
        poly = pe.get("entity_poly") or {}
        if (poly.get("rcsb_entity_polymer_type") or "") != "Protein":
            continue
        aligns = []
        for a in pe.get("rcsb_polymer_entity_align") or []:
            if a.get("reference_database_name") == "UniProt":
                aligns.append(dict(acc=a["reference_database_accession"],
                                   regions=[[r["entity_beg_seq_id"], r["ref_beg_seq_id"], r["length"]] for r in a.get("aligned_regions") or []]))
        rec["entities"].append(dict(entity_id=ci["entity_id"], chains=ci.get("auth_asym_ids") or [], asym=ci.get("asym_ids") or [],
                                    seq="".join((poly.get("pdbx_seq_one_letter_code_can") or "").split()),
                                    desc=(pe.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "", aligns=aligns))
    for ne in e.get("nonpolymer_entities") or []:
        comp = ne.get("nonpolymer_comp") or {}
        cc = comp.get("chem_comp") or {}
        cls, heavy = classify_ligand(cc.get("id"), cc.get("name"), cc.get("formula"), cc.get("type"))
        rec["ligands"].append(dict(comp=cc.get("id"), name=cc.get("name"), formula=cc.get("formula"), type=cc.get("type"),
                                   inchikey=(comp.get("rcsb_chem_comp_descriptor") or {}).get("InChIKey"),
                                   cls=cls, heavy=heavy,
                                   chains=(ne.get("rcsb_nonpolymer_entity_container_identifiers") or {}).get("auth_asym_ids") or []))
    return rec


def gql_fetch(ids):
    body = json.dumps({"query": GQL, "variables": {"ids": ids}}).encode()
    try:
        r = json.loads(http("https://data.rcsb.org/graphql", body))
    except urllib.error.HTTPError:
        if len(ids) == 1:
            return []
        h = len(ids) // 2
        return gql_fetch(ids[:h]) + gql_fetch(ids[h:])
    if r.get("errors") and not r.get("data"):
        if len(ids) == 1:
            return []
        h = len(ids) // 2
        return gql_fetch(ids[:h]) + gql_fetch(ids[h:])
    return [compact_entry(e) for e in (r["data"]["entries"] or []) if e]


def load_gql(S):
    out = {}
    p = f"{S}/api_cache/gql_entries.jsonl"
    if not os.path.exists(p) and os.path.exists(p + ".gz"):   # compressed once the run is finished
        lines = gzip.open(p + ".gz", "rt")
    elif os.path.exists(p):
        lines = open(p)
    else:
        return out
    for line in lines:
        r = json.loads(line)
        out[r["id"]] = r
    return out


def stage_gql(D, S, only, jobs):
    meta = load_meta(S)
    ids = set()
    for a in meta.index:
        if only and a not in only:
            continue
        ids |= {r["pdb_id"].upper() for r in pdbe_rows(S, a)}
    have = load_gql(S)
    todo = sorted(ids - set(have))
    log(f"gql: {len(ids)} unique PDB entries, {len(todo)} to fetch")
    batches = [tuple(todo[i:i + GQL_BATCH]) for i in range(0, len(todo), GQL_BATCH)]
    path = f"{S}/api_cache/gql_entries.jsonl"
    with open(path, "a") as out, cf.ThreadPoolExecutor(max_workers=jobs) as ex:
        futs = [ex.submit(lambda b: (time.sleep(0.1), gql_fetch(list(b)))[1], b) for b in batches]
        for n, f in enumerate(cf.as_completed(futs), 1):
            for rec in f.result():
                out.write(json.dumps(rec) + "\n")
            out.flush()
            if n % 20 == 0:
                log(f"  gql batches {n}/{len(batches)}")
    have = load_gql(S)
    log(f"gql: now {len(have)} entries; missing {len(ids - set(have))}")


# ------------------------------------------------------------------ stage select
def entry_candidates(S, acc, gq, variants):
    """One record per PDB entry mapped to acc, with SEQRES level metrics against each variant."""
    cands = []
    seen = {}
    src, shift = proxy(acc)
    for row in pdbe_rows(S, src):
        seen.setdefault(row["pdb_id"].upper(), row)
    for pid, row in seen.items():
        e = gq.get(pid)
        methods = " ".join(e["method"]).upper() if e else (row.get("experimental_method") or "").upper()
        if not any(m in methods for m in OLD_PDB_METHODS):
            continue
        res = e["res"] if e and e["res"] is not None else row.get("resolution")
        cand = dict(pdb_id=pid, method=methods, res=res, date=(e["date"] if e else ""),
                    has_gql=bool(e), ligands=(e["ligands"] if e else []), metrics={})
        cand["holo"] = any(l["cls"] == "druglike" for l in cand["ligands"])
        ents = [en for en in (e["entities"] if e else []) if any(a["acc"] == src for a in en["aligns"])]
        for v in variants:
            best = None
            for en in ents:
                pos_aa = {}
                for al in (a for a in en["aligns"] if a["acc"] == src):
                    for eb, rb, ln in al["regions"]:
                        for k in range(ln):
                            if eb - 1 + k < len(en["seq"]):
                                pos_aa[rb + k - shift] = en["seq"][eb - 1 + k]
                m = variant_metrics(pos_aa, v)
                if best is None or m["n_mapped"] > best["n_mapped"]:
                    best = m
            if best is None:   # no sequence level alignment: fall back on the PDBe coverage
                lo, hi = row["unp_start"] - shift, row["unp_end"] - shift
                n = len(v["seq"])
                ov = max(0, min(hi, v["off"] + n) - max(lo, v["off"] + 1) + 1)
                best = dict(n_mapped=ov, n_mismatch=0, identity=None, coverage=ov / n, mutation_in_structure=None, mutation_detail="")
            cand["metrics"][v["sid"]] = best
        cands.append(cand)
    return cands


def mut_rank(m, var):
    if not var["subs"]:
        return 0
    return {True: 0, None: 1, False: 2}[m["mutation_in_structure"]]


def mis_bucket(m):
    """0: identity >= 98% (the structure is the sequence), 1: 90-98% (a construct with changes), 2: below 90% (misaligned or a homolog)."""
    idt = m["identity"]
    if idt is None:
        return 1
    return 0 if idt >= MATCH_IDENTITY else (1 if idt >= 0.9 else 2)


def date_key(c):
    return -int((c["date"] or "0").replace("-", "") or 0)


def pick_for_variant(cands, var, generic):
    sid = var["sid"]
    pool = [c for c in cands if c["metrics"][sid]["coverage"] >= 0.05]
    if not pool:
        return []
    mr = lambda c: mut_rank(c["metrics"][sid], var)  # noqa: E731
    mb = lambda c: mis_bucket(c["metrics"][sid])     # noqa: E731
    cov = lambda c: c["metrics"][sid]["coverage"]    # noqa: E731
    rs = lambda c: c["res"] if c["res"] is not None else 99.0  # noqa: E731
    picks = []
    holo = [c for c in pool if c["holo"]]
    holo_ok = [c for c in holo if c["res"] is not None and c["res"] <= HOLO_RES_MAX] or holo
    if holo_ok:
        c = sorted(holo_ok, key=lambda c: (mr(c), mb(c), -round(cov(c) * 10), c["metrics"][sid]["n_mismatch"], rs(c), date_key(c)))[0]
        picks.append((c, "best_holo"))
    c = sorted(pool, key=lambda c: (mr(c), mb(c), -cov(c), rs(c), date_key(c)))[0]
    picks.append((c, "best_coverage"))
    maxcov = max(cov(c) for c in pool)
    rpool = [c for c in pool if c["res"] is not None and cov(c) >= 0.5 * maxcov]
    if rpool:
        c = sorted(rpool, key=lambda c: (mr(c), mb(c), rs(c), -cov(c)))[0]
        picks.append((c, "best_resolution"))
    picks = picks[:MAX_GENERIC_PICKS] if generic else []
    if var["subs"]:
        # the substitution may sit in a loop that the structure does not resolve: keep a few that carry it
        carriers = [c for c in pool if mr(c) == 0 and mb(c) <= 1]
        carriers.sort(key=lambda c: (mb(c), -int(c["holo"]), -round(cov(c) * 10), c["metrics"][sid]["n_mismatch"], rs(c), date_key(c)))
        picks += [(c, "mutant_match") for c in carriers[:MAX_MUTANT_PICKS]]
    return picks


def stage_select(D, S, only):
    meta = load_meta(S)
    gq = load_gql(S)
    variants = load_variants(D, meta)
    sel, summary = {}, []
    for acc in meta.index:
        if only and acc not in only:
            continue
        vs = variants[acc]
        cands = entry_candidates(S, acc, gq, vs)
        chosen = collections.OrderedDict()
        for i, v in enumerate(vs):
            generic = (i == 0) or not v["subs"]
            for c, role in pick_for_variant(cands, v, generic):
                e = chosen.setdefault(c["pdb_id"], dict(pdb_id=c["pdb_id"], roles=[], variants=[]))
                if role not in e["roles"]:
                    e["roles"].append(role)
                if v["sid"] not in e["variants"]:
                    e["variants"].append(v["sid"])
        sel[acc] = list(chosen.values())
        n_holo = sum(c["holo"] for c in cands)
        summary.append(dict(uniprot_accession=acc, n_pdb_total=len({r['pdb_id'] for r in pdbe_rows(S, proxy(acc)[0])}), n_pdb_allowed_method=len(cands),
                            n_pdb_with_druglike_ligand=n_holo, n_pdb_selected=len(chosen)))
    json.dump(sel, open(f"{S}/api_cache/selection.json", "w"), indent=1)
    pd.DataFrame(summary).to_csv(f"{S}/api_cache/selection_summary.tsv", sep="\t", index=False)
    s = pd.DataFrame(summary)
    log(f"select: {len(sel)} accessions; selected PDB entries {int(s.n_pdb_selected.sum())}; accessions without any: {int((s.n_pdb_selected == 0).sum())}")


# ------------------------------------------------------------------ stage download
def af_records(S, acc):
    """Canonical AlphaFold models of the accession (AF-<acc>-F<n>); isoform models (AF-<acc>-<k>-F<n>) are left out."""
    a = json.load(open(f"{S}/api_cache/alphafold/{acc}.json"))
    return [r for r in (a or []) if re.fullmatch(re.escape(f"AF-{acc}") + r"-F\d+", r["modelEntityId"])]


def af_isoform_count(S, acc):
    a = json.load(open(f"{S}/api_cache/alphafold/{acc}.json"))
    return len(a or []) - len(af_records(S, acc))


def stage_download(D, S, only, jobs):
    meta = load_meta(S)
    sel = json.load(open(f"{S}/api_cache/selection.json"))
    pdb_ids = sorted({e["pdb_id"] for a, es in sel.items() for e in es})
    af = []
    for acc in meta.index:
        if only and acc not in only:
            continue
        for r in af_records(S, acc):
            af.append((acc, r["modelEntityId"], r["pdbUrl"], r.get("latestVersion")))
    log(f"download: {len(pdb_ids)} PDB entries, {len(af)} AlphaFold models")

    def raw_bytes():
        return sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk(f"{S}/raw") for f in fs)

    def get_pdb(pid):
        p = f"{S}/raw/pdb/{pid}.cif.gz"
        if not os.path.exists(p):
            if raw_bytes() > MAX_BYTES:
                raise RuntimeError("size cap reached")
            b = http(f"https://files.rcsb.org/download/{pid}.cif.gz")
            if b is None:
                raise RuntimeError("404")
            atomic_write(p, b)
            time.sleep(0.15)
        return os.path.getsize(p)

    def get_sifts(pid):
        cached_json(f"{S}/raw/sifts/{pid}.json", f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pid.lower()}")
        return 0

    def get_af(item):
        acc, mid, url, ver = item
        p = f"{S}/raw/alphafold/{os.path.basename(url)}"
        if not os.path.exists(p):
            b = http(url)
            if b is None:
                raise RuntimeError("404")
            atomic_write(p, b)
            time.sleep(0.1)
        return os.path.getsize(p)

    with cf.ThreadPoolExecutor(max_workers=3) as top:   # three different hosts, at most `jobs` connections to each
        f1 = top.submit(pool_map, get_pdb, pdb_ids, jobs)
        f2 = top.submit(pool_map, get_sifts, pdb_ids, jobs)
        f3 = top.submit(pool_map, get_af, af, jobs)
        r1, r2, r3 = f1.result(), f2.result(), f3.result()
    for name, r in (("pdb", r1), ("sifts", r2), ("alphafold", r3)):
        bad = {k: str(v) for k, v in r.items() if isinstance(v, Exception)}
        log(f"download {name}: {len(r) - len(bad)} ok, {len(bad)} failed", list(bad.items())[:10])
    json.dump({"pdb_failed": {k: str(v) for k, v in r1.items() if isinstance(v, Exception)},
               "sifts_failed": {k: str(v) for k, v in r2.items() if isinstance(v, Exception)},
               "af_failed": {k[1]: str(v) for k, v in r3.items() if isinstance(v, Exception)}},
              open(f"{S}/api_cache/download_failures.json", "w"), indent=1)


# ------------------------------------------------------------------ structure helpers
def aa1(name):
    t = gemmi.find_tabulated_residue(name)
    if t is not None and t.is_amino_acid():
        c = (t.one_letter_code or "").strip().upper()
        return c if c.isalpha() else "X"
    return "X"


def heavy_coords(res):
    return np.array([[a.pos.x, a.pos.y, a.pos.z] for a in res if a.element.name not in ("H", "D")], dtype=float).reshape(-1, 3)


def segments(positions):
    ps = sorted(set(positions))
    out, start, prev = [], None, None
    for p in ps:
        if start is None:
            start = prev = p
        elif p == prev + 1:
            prev = p
        else:
            out.append((start, prev))
            start = prev = p
    if start is not None:
        out.append((start, prev))
    return ";".join(f"{a}-{b}" if a != b else str(a) for a, b in out)


def pdb_text(residues, chain_name, title, remarks):
    ns = gemmi.Structure()
    ns.name = title
    m = gemmi.Model("1")
    ch = gemmi.Chain(chain_name)
    for r in residues:
        ch.add_residue(r)
    m.add_chain(ch)
    ns.add_model(m)
    txt = "".join(l for l in ns.make_pdb_string().splitlines(True) if not l.startswith("CRYST1"))
    head = "".join(f"REMARK 999 {r}\n" for r in remarks)
    return head + txt


def write_text(path, txt):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(txt)


def residue_map_text(rows):
    cols = ["chain", "resnum", "icode", "resname", "aa", "label_seq", "uniprot_pos", "uniprot_aa", "group_pos_main", "source_uniprot_pos", "in_crop"]
    return "\t".join(cols) + "\n" + "".join("\t".join(str(x) for x in r) + "\n" for r in rows)


def main_variant(vs):
    return vs[0]


META = None   # the target table, set by stage_process; the UniProt sequences of stand-in sources are read from it


def segment_map(reslist, s0, s1, u0, u1, useq):
    """label_seq -> UniProt position for one SIFTS segment whose two ranges differ in length (an insertion or a
    deletion inside the construct): the residues of the chain are aligned to the UniProt slice and only the
    identical stretches (and substitutions of equal length) are mapped. A segment of equal lengths is linear."""
    rr = [(r.label_seq, aa1(r.name)) for r in reslist if r.label_seq is not None and s0 <= r.label_seq <= s1]
    a = "".join(x for _, x in rr)
    sm = difflib.SequenceMatcher(None, a, useq[u0 - 1:u1], autojunk=False)
    out = {}
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == "equal" or (tag == "replace" and i2 - i1 == j2 - j1):
            for k in range(i2 - i1):
                out[rr[i1 + k][0]] = u0 + j1 + k
    return out


KMER = 5
KMER_MIN_SUPPORT = 12


def kmer_remap(reslist, ref_seq):
    """Re-derive residue -> UniProt offsets from exact 5-mers that occur once in the reference sequence.

    SIFTS gives one offset per mapped segment; a construct with a deleted loop or a fused linker is mapped as a single
    segment although the offset changes at the junction, and every residue after it gets the wrong UniProt position.
    Here every residue takes the offset of the unique 5-mers that cover it, kept only if at least 12 residues share
    that offset; a residue between two anchors with the same offset inherits it. Returns {index in reslist: offset}."""
    idx = collections.defaultdict(list)
    for p in range(len(ref_seq) - KMER + 1):
        idx[ref_seq[p:p + KMER]].append(p + 1)
    items = sorted((r.label_seq, aa1(r.name), i) for i, r in enumerate(reslist) if r.label_seq is not None)
    votes = collections.defaultdict(collections.Counter)
    for j in range(len(items) - KMER + 1):
        w = items[j:j + KMER]
        if w[-1][0] - w[0][0] != KMER - 1:
            continue
        km = "".join(x[1] for x in w)
        hits = idx.get(km)
        if "X" not in km and hits and len(hits) == 1:
            for x in w:
                votes[x[2]][hits[0] - w[0][0]] += 1
    best = {i: c.most_common(1)[0][0] for i, c in votes.items()}
    support = collections.Counter(best.values())
    best = {i: d for i, d in best.items() if support[d] >= KMER_MIN_SUPPORT}
    order = [x[2] for x in items]
    prev, nxt, cur = {}, {}, None
    for i in order:
        cur = best.get(i, cur)
        prev[i] = cur
    cur = None
    for i in reversed(order):
        cur = best.get(i, cur)
        nxt[i] = cur
    out = {}
    for i in order:
        if i in best:
            out[i] = best[i]
        elif prev[i] is not None and prev[i] == nxt[i]:
            out[i] = prev[i]
    return out


def mapped_identity(reslist, upos, ref_seq, shift):
    n = m = 0
    for r, u in zip(reslist, upos):
        if u is None:
            continue
        i = u + shift
        if 1 <= i <= len(ref_seq):
            n += 1
            m += ref_seq[i - 1] == aa1(r.name)
    return m, n


def analyse_pdb(S, acc, entry, gq, sifts, vs, meta_row, mature_region, src_seq=None):
    pid = entry["pdb_id"]
    st = gemmi.read_structure(f"{S}/raw/pdb/{pid}.cif.gz")
    st.setup_entities()
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    if len(st) > 1:
        del st[1:]
    model = st[0]
    src, shift = proxy(acc)
    maps = ((sifts or {}).get(pid.lower()) or {}).get("UniProt", {}).get(src, {}).get("mappings", [])
    segs = collections.defaultdict(list)
    for m in maps:
        segs[m["struct_asym_id"]].append((m["start"]["residue_number"], m["end"]["residue_number"], m["unp_start"], m["unp_end"]))
    groups = collections.defaultdict(list)
    for ch in model:
        for res in ch:
            if res.entity_type == gemmi.EntityType.Polymer:
                groups[(ch.name, res.subchain)].append(res)
    gqe = gq.get(pid) or {"ligands": []}
    lig_info = {l["comp"]: l for l in gqe["ligands"]}
    cands = []
    for (cname, sub), reslist in groups.items():
        if sub not in segs:
            continue
        irregular = {}
        for s0, s1, u0, u1 in segs[sub]:
            if s1 - s0 != u1 - u0:
                irregular[(s0, s1)] = segment_map(reslist, s0, s1, u0, u1, META.sequence[src])
        mapped = []
        for r in reslist:
            ls = r.label_seq
            up = None
            if ls is not None:
                for s0, s1, u0, u1 in segs[sub]:
                    if s0 <= ls <= s1:
                        if (s0, s1) in irregular:
                            up = irregular[(s0, s1)].get(ls)
                            up = None if up is None else up - shift
                        elif u0 + (ls - s0) <= u1:
                            up = u0 + (ls - s0) - shift
                        break
            mapped.append(up)
        cands.append(dict(chain=cname, sub=sub, res=reslist, upos=mapped))
    if not cands:
        return None, "no chain of this entry is mapped to the accession by SIFTS"
    # ligands of the model
    ligs = []
    for ch in model:
        for res in ch:
            if res.entity_type == gemmi.EntityType.NonPolymer:
                info = lig_info.get(res.name)
                cls, heavy_f = (info["cls"], info["heavy"]) if info else classify_ligand(res.name, "", "", "")
                if cls == "druglike":
                    xyz = heavy_coords(res)
                    if len(xyz) >= 10:
                        ligs.append(dict(res=res, chain=ch.name, comp=res.name, xyz=xyz, info=info))
    for c in cands:
        keep = [(r, u) for r, u in zip(c["res"], c["upos"])]
        c["xyz_all"] = np.vstack([heavy_coords(r) for r, _ in keep]) if keep else np.zeros((0, 3))
        c["n_resolved_mapped"] = sum(u is not None for _, u in keep)
        best, bl = 0, None
        for l in ligs:
            if len(c["xyz_all"]) == 0:
                continue
            d = np.sqrt(((l["xyz"][:, None, :] - c["xyz_all"][None, :, :]) ** 2).sum(-1))
            n_c = int((d.min(axis=1) <= CONTACT).sum())
            if n_c > best:
                best, bl = n_c, l
        c["lig_contact"], c["lig"] = best, bl
    want_holo = "best_holo" in entry["roles"]
    key = (lambda c: (-c["lig_contact"], -c["n_resolved_mapped"], c["chain"])) if want_holo else (lambda c: (-c["n_resolved_mapped"], -c["lig_contact"], c["chain"]))
    cands.sort(key=key)
    ch0 = cands[0]
    if ch0["n_resolved_mapped"] == 0:
        return None, "no resolved residue of the chain is mapped to the accession"
    map_method = "SIFTS"
    if src_seq:
        m0, n0 = mapped_identity(ch0["res"], ch0["upos"], src_seq, shift)
        sifts_ident = m0 / n0 if n0 else 0.0
        if sifts_ident < MATCH_IDENTITY:
            rm = kmer_remap(ch0["res"], src_seq)
            new = [(r.label_seq + rm[i] - shift) if i in rm else None for i, r in enumerate(ch0["res"])]
            m1, n1 = mapped_identity(ch0["res"], new, src_seq, shift)
            if m1 > m0 + 5:
                ch0["upos"] = new
                map_method = f"kmer_realigned (SIFTS {sifts_ident:.0%} -> {m1 / n1:.0%} identical)"
    # crop to the mature chain when the target is a cleaved viral enzyme
    keep = []
    for r, u in zip(ch0["res"], ch0["upos"]):
        in_crop = True
        if mature_region is not None:
            in_crop = u is not None and mature_region[0] <= u <= mature_region[1]
        if in_crop:
            keep.append((r, u))
    pos_aa = {u: aa1(r.name) for r, u in keep if u is not None}
    res_rows = [(ch0["chain"], r.seqid.num, (r.seqid.icode or "").strip(), r.name, aa1(r.name), r.label_seq if r.label_seq is not None else "",
                 u if u is not None else "", meta_row.sequence[u - 1] if (u is not None and 1 <= u <= len(meta_row.sequence)) else "",
                 (u - vs[0]["off"]) if u is not None else "", (u + shift) if u is not None else "", True) for r, u in keep]
    metrics = {v["sid"]: variant_metrics(pos_aa, v) for v in vs}
    out = dict(pid=pid, chain=ch0["chain"], sub=ch0["sub"], map_method=map_method, keep=keep, pos_aa=pos_aa, metrics=metrics, res_rows=res_rows,
               n_chain_candidates=len(cands), n_resolved=len(keep))
    # primary ligand and pocket, measured on the chosen chain
    lig = ch0["lig"] if ch0["lig_contact"] > 0 else None
    out["ligands"] = []
    if lig is not None:
        # all drug-like ligand instances in contact with the chain, strongest contact first
        inst = []
        for l in ligs:
            d = np.sqrt(((l["xyz"][:, None, :] - ch0["xyz_all"][None, :, :]) ** 2).sum(-1))
            n_c = int((d.min(axis=1) <= CONTACT).sum())
            if n_c > 0:
                inst.append((n_c, l))
        inst.sort(key=lambda x: -x[0])
        out["ligands"] = inst
        pl = lig["xyz"]
        pocket, other = [], 0
        for r, u in keep:
            xyz = heavy_coords(r)
            if len(xyz) and np.sqrt(((xyz[:, None, :] - pl[None, :, :]) ** 2).sum(-1)).min() <= POCKET:
                pocket.append((r, u))
        for gkey, rl in groups.items():
            if gkey == (ch0["chain"], ch0["sub"]):
                continue
            for r in rl:
                xyz = heavy_coords(r)
                if len(xyz) and np.sqrt(((xyz[:, None, :] - pl[None, :, :]) ** 2).sum(-1)).min() <= POCKET:
                    other += 1
        out["pocket"] = pocket
        out["pocket_other_chain_residues"] = other
    return out, ""


def analyse_af(S, acc, rec, vs, mature_region, meta_row):
    path = f"{S}/raw/alphafold/{os.path.basename(rec['pdbUrl'])}"
    st = gemmi.read_structure(path)
    model = st[0]
    ch = model[0]
    reslist = list(ch)
    useq = rec.get("uniprotSequence") or rec.get("sequence") or ""
    start = rec.get("uniprotStart") or 1
    best_off, best_hit = 0, -1
    for off in {0, start - 1}:
        hit = sum(1 for r in reslist if 1 <= r.seqid.num + off <= len(useq) and useq[r.seqid.num + off - 1] == aa1(r.name))
        if hit > best_hit:
            best_off, best_hit = off, hit
    keep = []
    for r in reslist:
        u = r.seqid.num + best_off
        if mature_region is not None and not (mature_region[0] <= u <= mature_region[1]):
            continue
        keep.append((r, u))
    pos_aa = {u: aa1(r.name) for r, u in keep}
    pl = []
    for r, u in keep:
        for a in r:
            if a.name == "CA":
                pl.append(a.b_iso)
    pl = np.array(pl)
    res_rows = [(ch.name, r.seqid.num, "", r.name, aa1(r.name), r.seqid.num, u, meta_row.sequence[u - 1] if 1 <= u <= len(meta_row.sequence) else "",
                 u - vs[0]["off"], u, True) for r, u in keep]
    return dict(keep=keep, pos_aa=pos_aa, metrics={v["sid"]: variant_metrics(pos_aa, v) for v in vs}, res_rows=res_rows,
                plddt_mean=float(pl.mean()) if len(pl) else None, frac70=float((pl >= 70).mean()) if len(pl) else None,
                frac90=float((pl >= 90).mean()) if len(pl) else None, offset=best_off, seq_hits=best_hit, n=len(reslist)), ""


# ------------------------------------------------------------------ stage process
def stage_process(D, S, only):
    global META
    meta = load_meta(S)
    META = meta
    gq = load_gql(S)
    variants = load_variants(D, meta)
    sel = json.load(open(f"{S}/api_cache/selection.json"))
    fail = json.load(open(f"{S}/api_cache/download_failures.json"))
    rows, vrows, report = [], [], []
    for acc in meta.index:
        if only and acc not in only:
            continue
        vs = variants[acc]
        mrow = meta.loc[acc]
        mature = next((v["region"] for v in vs if v["region"]), None)
        src, shift = proxy(acc)
        n_own = len({r["pdb_id"] for r in pdbe_rows(S, acc)})
        n_total = len({r["pdb_id"] for r in pdbe_rows(S, src)})
        meths = collections.Counter()
        for pid in {r["pdb_id"].upper() for r in pdbe_rows(S, src)}:
            e = gq.get(pid)
            m = " ".join(e["method"]).upper() if e else ""
            meths["xray" if "X-RAY" in m else "em" if "ELECTRON MICROSCOPY" in m else "nmr" if "NMR" in m else "other"] += 1
        rep = dict(uniprot_accession=acc, scope=mrow.scope, structure_source_accession=src, n_pdb_own_entries=n_own, n_pdb_total=n_total, n_pdb_xray=meths["xray"], n_pdb_em=meths["em"],
                   n_pdb_nmr=meths["nmr"], n_pdb_selected=len(sel.get(acc, [])), n_pdb_downloaded=0, n_pdb_processed=0,
                   n_alphafold_models=0, n_alphafold_isoform_models_skipped=af_isoform_count(S, acc), has_experimental=False, has_holo=False, has_holo_le3A=False, fail_reasons=[])
        for entry in sel.get(acc, []):
            pid = entry["pdb_id"]
            if pid in fail["pdb_failed"] or pid in fail["sifts_failed"]:
                rep["fail_reasons"].append(f"{pid}: download failed")
                continue
            rep["n_pdb_downloaded"] += 1
            try:
                sifts = json.load(open(f"{S}/raw/sifts/{pid}.json"))
                a, why = analyse_pdb(S, acc, entry, gq, sifts, vs, mrow, mature, meta.sequence.get(proxy(acc)[0], ""))
            except Exception as e:  # noqa: BLE001
                a, why = None, f"{type(e).__name__}: {e}"
            if a is None:
                rep["fail_reasons"].append(f"{pid}: {why}")
                continue
            rep["n_pdb_processed"] += 1
            sid = f"{pid}_{a['chain']}"
            ent = gq.get(pid, {})
            ch_name = a["chain"] if len(a["chain"]) == 1 else "A"
            remarks = [f"PatBench r3 clean structure {sid}: source PDB {pid} chain {a['chain']} (label {a['sub']}), UniProt {acc}",
                       "protein only, model 1, first alternate location, no hydrogens, residue numbers are the author numbers of the entry",
                       "uniprot position of every residue: see the .residue_map.tsv next to this file" +
                       (f"; cropped to UniProt {mature[0]}-{mature[1]} (mature chain)" if mature else "")]
            write_text(f"{S}/clean/{acc}/{sid}.pdb", pdb_text([r for r, _ in a["keep"]], ch_name, sid, remarks))
            write_text(f"{S}/clean/{acc}/{sid}.residue_map.tsv", residue_map_text(a["res_rows"]))
            lig_files = []
            prim = None
            for k, (n_c, l) in enumerate(a["ligands"][:3]):
                lf = f"{pid}_{a['chain']}_{l['comp']}.pdb"
                if lf in lig_files:
                    lf = f"{pid}_{a['chain']}_{l['comp']}_{k + 1}.pdb"
                write_text(f"{S}/clean/{acc}/{lf}", pdb_text([l["res"]], l["chain"] if len(l["chain"]) == 1 else "L",
                                                            lf, [f"ligand {l['comp']} of {pid}, same coordinates as {sid}.pdb"]))
                lig_files.append(lf)
                if k == 0:
                    prim = l
            mm = a["metrics"][vs[0]["sid"]]
            info = prim["info"] if prim else None
            other = [x[1]["comp"] for x in a["ligands"][1:]]
            pocket = a.get("pocket", [])
            rows.append(dict(
                uniprot_accession=acc, gene_symbol=mrow.gene_symbol, protein_name=mrow.protein_name, organism=mrow.organism,
                structure_source="PDB", structure_id=sid, pdb_id=pid, chain_auth=a["chain"], chain_label=a["sub"],
                role="+".join(entry["roles"]), method=";".join(ent.get("method", [])), resolution=ent.get("res"), release_date=ent.get("date", ""),
                title=ent.get("title", ""), clean_path=f"clean/{acc}/{sid}.pdb", raw_path=f"raw/pdb/{pid}.cif.gz",
                residue_map_path=f"clean/{acc}/{sid}.residue_map.tsv", crop_region=f"{mature[0]}-{mature[1]}" if mature else "",
                uniprot_start=min(a["pos_aa"]) if a["pos_aa"] else None, uniprot_end=max(a["pos_aa"]) if a["pos_aa"] else None,
                resolved_segments_uniprot=segments(a["pos_aa"]), n_resolved_residues=a["n_resolved"],
                sequence_id_main=vs[0]["sid"], seq_identity_main=mm["identity"], n_mismatch_main=mm["n_mismatch"], seq_coverage_main=mm["coverage"],
                is_holo=bool(prim), ligand_comp_id=prim["comp"] if prim else "", ligand_name=(info or {}).get("name", ""),
                ligand_inchikey=(info or {}).get("inchikey") or "", ligand_heavy_atoms=len(prim["xyz"]) if prim else None,
                ligand_file=";".join(f"clean/{acc}/{x}" for x in lig_files), other_druglike_ligands=";".join(other),
                n_pocket_residues=len(pocket),
                pocket_residues_uniprot=";".join(f"{aa1(r.name)}{u}" for r, u in pocket if u is not None),
                pocket_residues_auth=";".join(f"{r.seqid.num}{(r.seqid.icode or '').strip()}" for r, _ in pocket),
                pocket_other_chain_residues=a.get("pocket_other_chain_residues"),
                plddt_mean=None, plddt_frac70=None, plddt_frac90=None, afdb_global_plddt_api=None, model_version="", model_date="",
                structure_source_accession=src, mapping_method=a["map_method"], selected_for=";".join(entry["variants"]), note="; ".join(filter(None, [
                    f"homolog stand-in: entry of {src} ({meta.gene_symbol.get(src, '')}), not of this target; sequence identity {(mm['identity'] or 0):.1%}" if src != acc else "",
                    "" if a["n_chain_candidates"] == 1 else f"{a['n_chain_candidates']} chains of the entry map to this target, one chosen",
                    "bound ligand sits in a dimer or cross-chain pocket: other chains contribute residues" if a.get("pocket_other_chain_residues") else "",
                    "resolution above 3.0 A" if (prim and ent.get("res") and ent["res"] > HOLO_RES_MAX) else "",
                    f"sequence identity to UniProt only {mm['identity']:.0%}: the entry's numbering or construct differs from UniProt, residue map unreliable" if (src == acc and mm["identity"] is not None and mm["identity"] < 0.9) else ""]))))
            for v in vs:
                m = a["metrics"][v["sid"]]
                vrows.append(dict(structure_id=sid, uniprot_accession=acc, sequence_id=v["sid"], **m))
        for rec in af_records(S, acc):
            mid = rec["modelEntityId"]
            if mid in fail["af_failed"]:
                rep["fail_reasons"].append(f"{mid}: download failed")
                continue
            try:
                a, why = analyse_af(S, acc, rec, vs, mature, mrow)
            except Exception as e:  # noqa: BLE001
                rep["fail_reasons"].append(f"{mid}: {type(e).__name__}: {e}")
                continue
            rep["n_alphafold_models"] += 1
            tag = f"_{mature[0]}-{mature[1]}" if mature else ""
            sid = mid
            remarks = [f"PatBench r3 clean structure {sid}: AlphaFold DB model v{rec.get('latestVersion')} of UniProt {acc}",
                       "B-factor column holds the per-residue pLDDT; residue numbers are UniProt positions" +
                       (f"; cropped to {mature[0]}-{mature[1]} (mature chain)" if mature else "")]
            write_text(f"{S}/clean/{acc}/{sid}{tag}.pdb", pdb_text([r for r, _ in a["keep"]], "A", sid, remarks))
            write_text(f"{S}/clean/{acc}/{sid}{tag}.residue_map.tsv", residue_map_text(a["res_rows"]))
            mm = a["metrics"][vs[0]["sid"]]
            rows.append(dict(
                uniprot_accession=acc, gene_symbol=mrow.gene_symbol, protein_name=mrow.protein_name, organism=mrow.organism,
                structure_source="AlphaFold", structure_id=sid, pdb_id="", chain_auth="A", chain_label="A", role="alphafold",
                method="AlphaFold prediction", resolution=None, release_date=(rec.get("modelCreatedDate") or "")[:10], title="",
                clean_path=f"clean/{acc}/{sid}{tag}.pdb", raw_path=f"raw/alphafold/{os.path.basename(rec['pdbUrl'])}",
                residue_map_path=f"clean/{acc}/{sid}{tag}.residue_map.tsv", crop_region=f"{mature[0]}-{mature[1]}" if mature else "",
                uniprot_start=min(a["pos_aa"]) if a["pos_aa"] else None, uniprot_end=max(a["pos_aa"]) if a["pos_aa"] else None,
                resolved_segments_uniprot=segments(a["pos_aa"]), n_resolved_residues=len(a["keep"]),
                sequence_id_main=vs[0]["sid"], seq_identity_main=mm["identity"], n_mismatch_main=mm["n_mismatch"], seq_coverage_main=mm["coverage"],
                is_holo=False, ligand_comp_id="", ligand_name="", ligand_inchikey="", ligand_heavy_atoms=None, ligand_file="",
                other_druglike_ligands="", n_pocket_residues=0, pocket_residues_uniprot="", pocket_residues_auth="",
                pocket_other_chain_residues=None, plddt_mean=a["plddt_mean"], plddt_frac70=a["frac70"], plddt_frac90=a["frac90"],
                afdb_global_plddt_api=rec.get("globalMetricValue"), model_version=f"v{rec.get('latestVersion')}",
                model_date=(rec.get("modelCreatedDate") or "")[:10], structure_source_accession=acc, mapping_method="AlphaFold numbering = UniProt", selected_for="", note=(
                    "" if len(af_records(S, acc)) == 1 else "AlphaFold fragment of a long protein")))
            for v in vs:
                vrows.append(dict(structure_id=sid, uniprot_accession=acc, sequence_id=v["sid"], **a["metrics"][v["sid"]]))
        t = [r for r in rows if r["uniprot_accession"] == acc and r["structure_source"] == "PDB"]
        rep["has_experimental"] = bool(t)
        rep["has_holo"] = any(r["is_holo"] for r in t)
        rep["has_holo_le3A"] = any(r["is_holo"] and r["resolution"] is not None and r["resolution"] <= HOLO_RES_MAX for r in t)
        rep["has_alphafold"] = rep["n_alphafold_models"] > 0
        rep["fail_reasons"] = " | ".join(rep["fail_reasons"])
        report.append(rep)
        log(f"process {acc}: {rep['n_pdb_processed']}/{rep['n_pdb_selected']} PDB, AF {rep['n_alphafold_models']}, holo {rep['has_holo']}")
    df = pd.DataFrame(rows)
    # leakage check: is the crystallised ligand a PatBench compound, and does it sit in a test group of the same target?
    # split_time comes from TEST_REL (the release whose test set is the one in use); group ids are identical across releases.
    sar = pd.read_parquet(f"{D}/sar_all.parquet", columns=["InChIKey", "uniprot_accession", "group_id"])
    split = pd.read_parquet(f"{TEST_REL}/sar_all.parquet", columns=["group_id", "split_time"]).drop_duplicates("group_id").set_index("group_id").split_time
    sar["split_time"] = sar.group_id.map(split).fillna("")
    sar["blk"] = sar.InChIKey.str[:14]
    rank_main = set()
    for tname in ("rank_ic50", "rank_ki_kd"):
        d = pd.read_parquet(f"{TEST_REL}/tasks/{tname}.parquet", columns=["group_id", "split_time"])
        rank_main |= set(d[d.split_time == "test"].group_id)
    blocks = set(df.ligand_inchikey.str[:14]) - {""}
    sub = sar[sar.blk.isin(blocks)]
    anyblk = set(sub.blk)
    same = sub.groupby(["blk", "uniprot_accession"]).group_id.nunique()
    test_same = sub[sub.split_time == "test"].groupby(["blk", "uniprot_accession"]).group_id.nunique()
    main_same = sub[sub.group_id.isin(rank_main)].groupby(["blk", "uniprot_accession"]).group_id.nunique()
    testblk = set(sub[sub.split_time == "test"].blk)
    mainblk = set(sub[sub.group_id.isin(rank_main)].blk)
    df["ligand_in_patbench"] = [bool(k) and k[:14] in anyblk for k in df.ligand_inchikey]
    df["ligand_in_patbench_same_target"] = [bool(k) and (k[:14], a) in same.index for k, a in zip(df.ligand_inchikey, df.uniprot_accession)]
    df["ligand_test_groups_same_target"] = [int(test_same.get((k[:14], a), 0)) if k else 0 for k, a in zip(df.ligand_inchikey, df.uniprot_accession)]
    df["ligand_in_patbench_test_any_target"] = [bool(k) and k[:14] in testblk for k in df.ligand_inchikey]
    df["ligand_test_groups_ic50_kikd_same_target"] = [int(main_same.get((k[:14], a), 0)) if k else 0 for k, a in zip(df.ligand_inchikey, df.uniprot_accession)]
    df["ligand_in_ic50_kikd_test_any_target"] = [bool(k) and k[:14] in mainblk for k in df.ligand_inchikey]
    df["leakage_split_source"] = os.path.basename(TEST_REL)
    df.to_pickle(f"{S}/api_cache/structures.pkl")
    pd.DataFrame(vrows).to_pickle(f"{S}/api_cache/variant_metrics.pkl")
    rp = pd.DataFrame(report)
    rp.to_pickle(f"{S}/api_cache/report.pkl")
    log(f"process: {len(df)} structure rows ({(df.structure_source == 'PDB').sum()} PDB, {(df.structure_source == 'AlphaFold').sum()} AlphaFold)")


# ------------------------------------------------------------------ stage groups
def stage_groups(D, S, only):
    meta = load_meta(S)
    variants = load_variants(D, meta)
    st = pd.read_pickle(f"{S}/api_cache/structures.pkl")
    vm = pd.read_pickle(f"{S}/api_cache/variant_metrics.pkl")
    rp = pd.read_pickle(f"{S}/api_cache/report.pkl")
    ts = pd.read_parquet(f"{D}/target_sequences.parquet")
    vmi = vm.set_index(["structure_id", "uniprot_accession", "sequence_id"])
    gblocks = collections.defaultdict(set)    # group -> InChIKey first blocks of its molecules
    for tname in RANK_TASKS:
        d = pd.read_parquet(f"{D}/tasks/{tname}.parquet", columns=["group_id", "InChIKey"])
        for gid, k in zip(d.group_id, d.InChIKey.str[:14]):
            gblocks[gid].add(k)
    out = []
    for _, g in ts.iterrows():
        acc, sid = g.uniprot_accession, g.sequence_id
        var = next(v for v in variants[acc] if v["sid"] == sid)
        cand = st[(st.uniprot_accession == acc)]
        pdbs = cand[cand.structure_source == "PDB"]
        afs = cand[cand.structure_source == "AlphaFold"]

        def mets(row):
            return vmi.loc[(row.structure_id, acc, sid)]

        def rank(row, coverage_first=False):
            m = mets(row)
            mr = 0 if not var["subs"] else {True: 0, False: 2}.get(m.mutation_in_structure, 1) if pd.notna(m.mutation_in_structure) else 1
            mb = 1 if pd.isna(m.identity) else (0 if m.identity >= MATCH_IDENTITY else (1 if m.identity >= 0.9 else 2))
            res = row.resolution if pd.notna(row.resolution) else 99.0
            holo = 2 if (row.is_holo and res <= HOLO_RES_MAX) else (1 if row.is_holo else 0)
            cov = -round(m.coverage * 10)
            date = -int((row.release_date or "0").replace("-", "") or 0)
            nmis = int(m.n_mismatch)
            return (mr, mb, cov, -holo, nmis, res, date) if coverage_first else (mr, mb, -holo, cov, nmis, res, date)
        ok = [r for _, r in pdbs.iterrows() if mets(r).n_mapped >= 20]
        ok.sort(key=rank)
        cov_first = sorted(ok, key=lambda r: rank(r, True))[0].structure_id if ok else ""
        best_af = None
        if len(afs):
            best_af = sorted([r for _, r in afs.iterrows()], key=lambda r: -mets(r).coverage)[0]
        if ok:
            rec = ok[0]
            m = mets(rec)
            src = "PDB"
        elif best_af is not None:
            rec, m, src = best_af, mets(best_af), "AlphaFold"
        else:
            rec, m, src = None, None, ""
        reason = []
        if rec is None:
            reason.append("no structure available")
        else:
            if src == "AlphaFold":
                reason.append("no experimental structure covers this target" if len(pdbs) == 0 else "no experimental structure with at least 20 mapped residues")
            if g.sequence_kind == "mutant_applied" and src == "PDB":
                reason.append("structure carries the substitutions" if m.mutation_in_structure is True else
                              ("structure is wild type or lacks the substitution site; " + m.mutation_detail if m.mutation_in_structure is False
                               else "substitution site not resolved; " + m.mutation_detail))
            if src == "PDB" and rec.structure_source_accession != acc:
                reason.append(f"homolog stand-in: experimental structure of {rec.structure_source_accession}, sequence identity {m.identity:.1%}")
            if src == "PDB" and rec.structure_source_accession == acc and m.identity is not None and m.identity < MATCH_IDENTITY \
                    and g.sequence_kind != "mutant_applied":
                reason.append(f"sequence identity {m.identity:.1%} to the group sequence: a construct with changes, or the entry's numbering differs")
            if g.sequence_kind == "complex_primary_subunit":
                reason.append(f"only the first subunit of the complex ({g.target_symbols}) is in the structure file")
            if g.sequence_kind == "mature_chain":
                reason.append(f"cropped to the mature chain {var['region'][0]}-{var['region'][1]}")
            if m.coverage < 0.5:
                reason.append(f"covers {m.coverage:.0%} of the sequence (domain structure)")
            if src == "PDB" and not rec.is_holo:
                reason.append("apo (no drug-like ligand in contact with the chain)")
        match = None
        if rec is not None:
            match = bool(m.identity is not None and m.identity >= MATCH_IDENTITY and m.n_mapped >= 20 and m.mutation_in_structure is not False)
        am = mets(best_af) if best_af is not None else None
        out.append(dict(
            group_id=g.group_id, uniprot_accession=acc, gene_symbol=meta.gene_symbol[acc], sequence_id=sid, sequence_kind=g.sequence_kind,
            mutation=g.mutation, recommended_structure_id=rec.structure_id if rec is not None else "", recommended_source=src,
            recommended_clean_path=rec.clean_path if rec is not None else "", recommended_role=rec.role if rec is not None else "",
            recommended_method=rec.method if rec is not None else "", recommended_resolution=rec.resolution if rec is not None else None,
            recommended_is_holo=bool(rec.is_holo) if rec is not None else None,
            seq_identity=m.identity if m is not None else None, n_mismatch=int(m.n_mismatch) if m is not None else None,
            seq_coverage=m.coverage if m is not None else None, structure_matches_sequence=match,
            mutation_in_structure=(m.mutation_in_structure if (m is not None and var["subs"]) else None),
            mutation_detail=m.mutation_detail if (m is not None and var["subs"]) else "",
            structure_is_homolog_standin=bool(rec is not None and rec.structure_source_accession != acc),
            coverage_first_structure_id=cov_first,
            alphafold_structure_id=best_af.structure_id if best_af is not None else "",
            alphafold_clean_path=best_af.clean_path if best_af is not None else "",
            af_identity=am.identity if am is not None else None, af_coverage=am.coverage if am is not None else None,
            af_mean_plddt=best_af.plddt_mean if best_af is not None else None,
            recommended_ligand_comp_id=rec.ligand_comp_id if src == "PDB" else "",
            recommended_ligand_in_group=bool(src == "PDB" and rec.ligand_inchikey and rec.ligand_inchikey[:14] in gblocks[g.group_id]),
            recommended_ligand_in_patbench_same_target=bool(src == "PDB" and rec.ligand_in_patbench_same_target),
            recommended_ligand_in_test_same_target=bool(src == "PDB" and rec.ligand_test_groups_same_target > 0),
            recommended_ligand_in_ic50_kikd_test_same_target=bool(src == "PDB" and rec.ligand_test_groups_ic50_kikd_same_target > 0),
            n_experimental_structures=len(pdbs), reason="; ".join(reason)))
    gm = pd.DataFrame(out)
    gm["mutation_in_structure"] = gm["mutation_in_structure"].astype("boolean")
    gm["structure_matches_sequence"] = gm["structure_matches_sequence"].astype("boolean")
    gm.to_pickle(f"{S}/api_cache/group_map.pkl")
    log(f"groups: {len(gm)} groups, unique {gm.group_id.nunique()}, matches {int(gm.structure_matches_sequence.fillna(False).sum())}")


def stage_write(D, S):
    st = pd.read_pickle(f"{S}/api_cache/structures.pkl")
    gm = pd.read_pickle(f"{S}/api_cache/group_map.pkl")
    rp = pd.read_pickle(f"{S}/api_cache/report.pkl")
    vm = pd.read_pickle(f"{S}/api_cache/variant_metrics.pkl")
    meta = load_meta(S)
    st = st.copy()
    st.insert(1, "scope", st.uniprot_accession.map(meta.scope))
    for c in ("ligand_heavy_atoms", "pocket_other_chain_residues", "uniprot_start", "uniprot_end"):
        st[c] = st[c].astype("Int64")
    for name, df in (("target_structures", st), ("group_structure_map", gm), ("structure_variant_metrics", vm)):
        df.to_parquet(f"{S}/{name}.parquet", index=False)
        df.to_csv(f"{S}/{name}.tsv", sep="\t", index=False, encoding="utf-8-sig")
    rp.to_csv(f"{S}/structure_report.tsv", sep="\t", index=False, encoding="utf-8-sig")
    log("write: tables written")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("release")
    ap.add_argument("stage")
    ap.add_argument("--only", default="")
    ap.add_argument("--outdir", default="structures")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--test-release", default="", help="release whose sar_all.parquet gives split_time for the ligand leakage count (default: the release itself)")
    ap.add_argument("--max-gb", type=float, default=1.8, help="stop downloading when raw/ is this large")
    a = ap.parse_args()
    D = a.release.rstrip("/")
    S = f"{D}/{a.outdir}"
    only = set(x for x in a.only.split(",") if x)
    global TEST_REL, MAX_BYTES
    TEST_REL = a.test_release.rstrip("/") or D
    MAX_BYTES = int(a.max_gb * 1e9)
    os.makedirs(S, exist_ok=True)
    with open(f"{S}/fetch_log.txt", "a") as f:
        f.write(f"{datetime.datetime.now().isoformat(timespec='seconds')} stage={a.stage} only={sorted(only) or 'all'} gemmi={gemmi.__version__}\n")
    jobs = min(a.jobs, 4)
    if a.stage == "meta":
        stage_meta(D, S, only)
    elif a.stage == "entries":
        stage_entries(D, S, only, jobs)
    elif a.stage == "gql":
        stage_gql(D, S, only, jobs)
    elif a.stage == "select":
        stage_select(D, S, only)
    elif a.stage == "download":
        stage_download(D, S, only, jobs)
    elif a.stage == "process":
        stage_process(D, S, only)
    elif a.stage == "groups":
        stage_groups(D, S, only)
    elif a.stage == "write":
        stage_write(D, S)
    else:
        raise SystemExit(f"unknown stage {a.stage}")


if __name__ == "__main__":
    main()
