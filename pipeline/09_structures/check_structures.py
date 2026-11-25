"""Independent checks of <release>/structures/: re-reads the written files and recomputes what the tables claim.

    python -I check_structures.py <release_dir> [seed] [outdir]

Prints one line per check, 'OK' or 'FAIL', and exits with 1 when any check failed.
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import os
import random
import sys

sys.path.insert(0, _ws + "/struct-deps-20261030")
sys.path.insert(0, _ws + "/outputs/uspto-sar-demo-20261003/runtime/deps")
import gemmi  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

D = sys.argv[1].rstrip("/")
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 20261030
S = f"{D}/{sys.argv[3] if len(sys.argv) > 3 else 'structures'}"
st = pd.read_parquet(f"{S}/target_structures.parquet")
gm = pd.read_parquet(f"{S}/group_structure_map.parquet")
rp = pd.read_csv(f"{S}/structure_report.tsv", sep="\t", encoding="utf-8-sig")
nm = pd.read_csv(f"{S}/target_name_map.tsv", sep="\t", encoding="utf-8-sig", keep_default_na=False)
ts = pd.read_parquet(f"{D}/target_sequences.parquet")
tg = pd.read_parquet(f"{D}/targets.parquet").set_index("uniprot_accession")
fails = []


def check(name, ok, detail=""):
    print(("OK   " if ok else "FAIL ") + name + (f" | {detail}" if detail else ""), flush=True)
    if not ok:
        fails.append(name)


def aa1(name):
    t = gemmi.find_tabulated_residue(name)
    if t is not None and t.is_amino_acid():
        c = (t.one_letter_code or "").strip().upper()
        return c if c.isalpha() else "X"
    return "X"


def read_clean(path):
    s = gemmi.read_structure(f"{S}/{path}")
    return s, [r for r in s[0][0]]


# 1. coverage of the group table
check("group map has one row per target-sequence group", len(gm) == len(ts), f"{len(gm)} vs {len(ts)}")
check("group ids unique", gm.group_id.is_unique)
check("group ids equal target_sequences", set(gm.group_id) == set(ts.group_id))
none_acc = set(rp[(rp.has_experimental.astype(str) == "False") & (rp.has_alphafold.astype(str) == "False")].uniprot_accession)
no_rec = gm[gm.recommended_structure_id == ""]
check("groups without a structure are exactly those of targets with neither a PDB nor an AlphaFold entry",
      set(no_rec.uniprot_accession) <= none_acc and set(gm[gm.uniprot_accession.isin(none_acc)].group_id) == set(no_rec.group_id),
      f"{len(no_rec)} groups, targets {sorted(set(no_rec.uniprot_accession))}")
paths = pd.concat([st.clean_path, st.raw_path, st.residue_map_path])
missing = [p for p in paths if not os.path.exists(f"{S}/{p}")]
check("all clean, raw and residue-map files of target_structures exist", not missing, f"missing {len(missing)}")
lig_missing = [p for x in st.ligand_file for p in str(x).split(";") if p and not os.path.exists(f"{S}/{p}")]
check("all ligand files exist", not lig_missing, f"missing {len(lig_missing)}")
check("recommended and AlphaFold paths exist", all(os.path.exists(f"{S}/{p}") for p in gm.recommended_clean_path if p) and
      all(os.path.exists(f"{S}/{p}") for p in gm.alphafold_clean_path if p))
_reg = set(pd.read_parquet(f"{D}/tasks/regression_affinity.parquet", columns=["uniprot_accession"]).uniprot_accession)
check("name map has one row per target", nm.uniprot_accession.is_unique and len(nm) == len(set(ts.uniprot_accession) | _reg), f"{len(nm)} vs {len(set(ts.uniprot_accession) | _reg)}")
check("every rank-task target is in the report", set(ts.uniprot_accession) <= set(rp.uniprot_accession))

# 2. random sample of experimental and AlphaFold structures: parse, count, compare with the residue map and the table
rng = random.Random(SEED)
pdbs = st[st.structure_source == "PDB"]
accs = sorted(pdbs.uniprot_accession.unique())
sample = rng.sample(accs, min(14, len(accs)))
for acc in sample:
    r = pdbs[pdbs.uniprot_accession == acc].iloc[0]
    s, res = read_clean(r.clean_path)
    rm = pd.read_csv(f"{S}/{r.residue_map_path}", sep="\t", keep_default_na=False)
    ok_n = len(res) == len(rm) == r.n_resolved_residues
    check(f"{acc} {r.structure_id}: residue count file = map = table", ok_n, f"{len(res)} / {len(rm)} / {r.n_resolved_residues}")
    ok_aa = all(aa1(x.name) == a for x, a in zip(res, rm.aa))
    check(f"{acc} {r.structure_id}: residue letters in file = map", ok_aa)
    mapped = rm[rm.uniprot_pos != ""]
    ref = tg.sequence[acc] if r.structure_source_accession == acc else tg.sequence[r.structure_source_accession]
    same = np.mean([ref[int(p) - 1] == a for p, a in zip(mapped.source_uniprot_pos, mapped.aa)]) if len(mapped) else 0
    check(f"{acc} {r.structure_id}: mapped residues equal the UniProt sequence of {r.structure_source_accession} (or the low identity is the one in the table)",
          same >= 0.9 or abs(same - r.seq_identity_main) < 0.01, f"{same:.3f} vs table {r.seq_identity_main:.3f}")
    ch = {x.name for x in s[0]}
    check(f"{acc} {r.structure_id}: one chain", len(ch) == 1, str(ch))
    has_h = any(a.element.name in ("H", "D") for x in res for a in x)
    check(f"{acc} {r.structure_id}: no hydrogens", not has_h)

afs = st[st.structure_source == "AlphaFold"]
for _, r in afs.sample(min(10, len(afs)), random_state=SEED).iterrows():
    s, res = read_clean(r.clean_path)
    b = np.array([a.b_iso for x in res for a in x if a.name == "CA"])
    ok = abs(b.mean() - r.plddt_mean) < 0.01 and abs((b >= 70).mean() - r.plddt_frac70) < 1e-6 and abs((b >= 90).mean() - r.plddt_frac90) < 1e-6
    check(f"{r.structure_id}: pLDDT from the B-factor column = table", ok, f"{b.mean():.2f} vs {r.plddt_mean:.2f}")
    check(f"{r.structure_id}: mean pLDDT within 0.5 of the AlphaFold DB global metric", r.crop_region != "" or abs(r.plddt_mean - r.afdb_global_plddt_api) < 0.5,
          f"{r.plddt_mean:.2f} vs {r.afdb_global_plddt_api}")

# 3. holo structures: ligand file and pocket recomputed from the files
holo = pdbs[pdbs.is_holo]
for _, r in holo.sample(min(12, len(holo)), random_state=SEED).iterrows():
    s, res = read_clean(r.clean_path)
    lf = r.ligand_file.split(";")[0]
    ls, lres = read_clean(lf)
    lat = np.array([[a.pos.x, a.pos.y, a.pos.z] for x in lres for a in x])
    check(f"{r.structure_id}: ligand file is one HETATM residue of {r.ligand_comp_id}", len(lres) == 1 and lres[0].name == r.ligand_comp_id[:3] and lres[0].het_flag == "H",
          "PDB format keeps 3 letters of the residue name; the full CCD id is in the file name" if len(r.ligand_comp_id) > 3 else "")
    check(f"{r.structure_id}: ligand heavy atoms = table", len(lat) == r.ligand_heavy_atoms, f"{len(lat)} vs {r.ligand_heavy_atoms}")
    pocket, dmin = [], 99.0
    for x in res:
        xyz = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in x])
        d = np.sqrt(((xyz[:, None, :] - lat[None, :, :]) ** 2).sum(-1)).min()
        dmin = min(dmin, d)
        if d <= 6.0:
            pocket.append(f"{x.seqid.num}{(x.seqid.icode or '').strip()}")
    check(f"{r.structure_id}: ligand within 4.5 A of the chain", dmin <= 4.5, f"{dmin:.2f}")
    check(f"{r.structure_id}: pocket residues recomputed = table", pocket == [p for p in r.pocket_residues_auth.split(";") if p], f"{len(pocket)} vs {r.n_pocket_residues}")

# 4. special sequences
seqs = ts.drop_duplicates("sequence_id").set_index("sequence_id")
for sid in seqs.index[seqs.sequence_kind == "mature_chain"]:
    g = gm[gm.sequence_id == sid].iloc[0]
    s, res = read_clean(g.recommended_clean_path)
    got = "".join(aa1(x.name) for x in res)
    want = seqs.sequence[sid]
    n_eq = sum(a == b for a, b in zip(got, want)) if len(got) == len(want) else -1
    check(f"mature chain {sid}: {g.recommended_structure_id} has {len(got)} residues (group sequence {len(want)})", len(got) >= 0.85 * len(want),
          f"identical positions {n_eq}; identity {g.seq_identity:.3f}; stand-in {bool(g.structure_is_homolog_standin)}")
    if sid == "P03366:501-599":
        n_same = sum(a == b for a, b in zip(got, want))
        check("HIV-1 protease crop is 99 residues and differs from the group sequence at most at 2 positions (engineered constructs)",
              len(got) == 99 and n_same >= 97, f"{len(got)} residues, {n_same} identical")
expect = {"P01116:G12C": [(12, "CYS")], "P15056:V600E": [(600, "GLU")], "P03372:Y537S": [(537, "SER")], "P03372:D538G": [(538, "GLY")],
          "Q5S007:G2019S": [(2019, "SER")], "O75874:R132H": [(132, "HIS")]}
for sid, subs in expect.items():
    rows = gm[gm.sequence_id == sid]
    if rows.empty:
        check(f"mutant {sid}: group present", False)
        continue
    g = rows.iloc[0]
    if g.recommended_source != "PDB" or pd.isna(g.mutation_in_structure) or not bool(g.mutation_in_structure):
        check(f"mutant {sid}: structure carries it", False, f"{g.recommended_structure_id} mutation_in_structure={g.mutation_in_structure} {g.mutation_detail}")
        continue
    s, res = read_clean(g.recommended_clean_path)
    rm = pd.read_csv(f"{S}/{g.recommended_clean_path[:-4]}.residue_map.tsv", sep="\t", keep_default_na=False)
    rm["uniprot_pos"] = pd.to_numeric(rm.uniprot_pos, errors="coerce")   # '' for tag residues that have no UniProt position
    rm = rm.dropna(subset=["uniprot_pos"]).astype({"uniprot_pos": int}).set_index("uniprot_pos")
    ok = True
    for pos, three in subs:
        row = rm.loc[pos] if pos in rm.index else None
        ok &= row is not None and row.resname == three
        if row is not None:   # the same residue read straight from the structure by its author number
            hit = [x for x in res if x.seqid.num == row.resnum and (x.seqid.icode or "").strip() == row.icode]
            ok &= bool(hit) and hit[0].name == three
    check(f"mutant {sid}: residue {subs} found in {g.recommended_structure_id} (file and map agree)", ok)
wt = gm[(gm.sequence_kind == "mutant_applied") & (gm.mutation_in_structure == False)]  # noqa: E712
print(f"INFO mutant groups whose recommended structure lacks the substitution: {len(wt)} of {int((gm.sequence_kind == 'mutant_applied').sum())}")

# 5. recompute structure_matches_sequence for a sample of groups from the files
for _, g in gm[gm.recommended_source == "PDB"].sample(20, random_state=SEED).iterrows():
    s, res = read_clean(g.recommended_clean_path)
    rm = pd.read_csv(f"{S}/{g.recommended_clean_path[:-4]}.residue_map.tsv", sep="\t", keep_default_na=False)
    seq = ts[ts.group_id == g.group_id].sequence.iloc[0]
    off = 0
    if ":" in g.sequence_id and "-" in g.sequence_id.split(":")[1] and g.sequence_id.split(":")[1].replace("-", "").isdigit():
        off = int(g.sequence_id.split(":")[1].split("-")[0]) - 1
    m = rm[rm.uniprot_pos != ""]
    pos = [(int(p) - off, a) for p, a in zip(m.uniprot_pos, m.aa) if 1 <= int(p) - off <= len(seq)]
    ident = np.mean([seq[i - 1] == a for i, a in pos]) if pos else 0
    cov = len(pos) / len(seq)
    check(f"group {g.group_id} ({g.recommended_structure_id}): identity and coverage recomputed from files", abs(ident - g.seq_identity) < 1e-6 and abs(cov - g.seq_coverage) < 1e-6,
          f"{ident:.4f}/{cov:.4f} vs {g.seq_identity:.4f}/{g.seq_coverage:.4f}")

print("\nFAILED CHECKS:", len(fails))
for f in fails:
    print("  -", f)
sys.exit(1 if fails else 0)
