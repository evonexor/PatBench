"""Per-group target sequences for the PatBench ranking tasks, with the audit behind them.

    python audit_targets.py <target_audit_groups.parquet> <targets.parquet> <out_dir>

target_audit_groups.parquet has one row per ranking-task group: the group joined to the target-map pair
it was mapped from (label, mapped identity, mutation, sequence caveat, cross-reference to BindingDB/ChEMBL,
patent title). The script writes

  target_sequences.parquet   one row per group: the sequence to use for that assay, how it was obtained
                             (sequence_kind) and the checks (label_check, xref_status)
  sequence_fixes.tsv         every sequence that differs from the UniProt entry, with its evidence

The sequence in targets.parquet is the UniProt canonical sequence of the accession. It is the wrong
protein for three kinds of group, and this script resolves the first two from UniProt itself:

  polyprotein   a viral enzyme (HIV-1 protease, reverse transcriptase, HCV NS5B) is mapped to the
                accession of the whole polyprotein; the mature chain is cut out with UniProt's own
                Chain feature, which is looked up by name and must match exactly one feature.
  mutant        the target is a named mutant, the stored sequence is wild type; the substitutions are
                applied, and only if the wild-type residue is found at every position.
  complex       the accession is the first subunit of a complex; flagged, not changed.

It needs the internet for the UniProt features (three requests).
"""
import json
import re
import sys
import time
import urllib.error
import urllib.request

import pandas as pd

AUD, TGT, OUT = sys.argv[1], sys.argv[2], sys.argv[3].rstrip("/")

# (accession, word in the mapped identity) -> exact Chain description in the UniProt entry
CHAINS = {("P03366", "protease"): "Protease",
          ("O91080", "reverse transcriptase"): "Reverse transcriptase/ribonuclease H",
          ("P26663", "NS5B"): "RNA-directed RNA polymerase"}
# the caveat also fires on very long ordinary proteins; these are whole proteins, not polyproteins
NOT_POLYPROTEIN = {"Q12721", "Q4WLT4", "O13423"}
# a haplotype is named, not written as substitutions
HAPLOTYPES = {"HAQ": "R71H;G230A;R293Q"}

a = pd.read_parquet(AUD)
tg = pd.read_parquet(TGT).set_index("uniprot_accession")
live = {}


def uniprot_json(acc):
    if acc not in live:
        for k in range(8):   # the UniProt REST service answers 500 now and then; the same request works seconds later
            try:
                with urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=60) as r:
                    live[acc] = json.load(r)
                break
            except urllib.error.HTTPError as e:
                if e.code < 500 or k == 7:
                    raise
                time.sleep(3 * (k + 1))
    return live[acc]


fixes, seq_of, kind_of, note_of = [], {}, {}, {}
for gid, r in a.set_index("group_id").iterrows():
    acc, cav = r.uniprot_accession, (r.sequence_caveat or "")
    base = tg.sequence[acc]
    seq, kind, note, sid = base, "full_length", "", acc
    if "polyprotein" in cav and acc not in NOT_POLYPROTEIN:
        ident = str(r.mapped_identity)
        hit = [d for (ac, w), d in CHAINS.items() if ac == acc and w.lower() in ident.lower()]
        if hit:
            feats = [f for f in uniprot_json(acc)["features"]
                     if f["type"] in ("Chain", "Peptide") and f.get("description") == hit[0]]
            assert len(feats) == 1, (acc, hit, len(feats))
            s, e = feats[0]["location"]["start"]["value"], feats[0]["location"]["end"]["value"]
            seq, kind, sid = base[s - 1:e], "mature_chain", f"{acc}:{s}-{e}"
            note = f"UniProt chain '{hit[0]}' {s}-{e} of the polyprotein ({len(base)} aa)"
            fixes.append(dict(sequence_id=sid, uniprot_accession=acc, kind=kind, description=hit[0],
                              start=s, end=e, evidence=f"UniProt {feats[0].get('featureId', 'Chain')}",
                              length=len(seq), sequence=seq))
        else:
            kind, note = "polyprotein_not_resolved", f"whole polyprotein kept; no chain rule for '{ident}'"
    elif "wild_type" in cav:
        mut = HAPLOTYPES.get(str(r.mutation).split(" ")[0], str(r.mutation))
        subs = re.findall(r"([A-Z])(\d+)([A-Z])", mut)
        if subs and all(base[int(p) - 1] == w for w, p, _ in subs):
            s = list(base)
            for w, p, m in subs:
                s[int(p) - 1] = m
            seq, kind, sid = "".join(s), "mutant_applied", f"{acc}:{'+'.join(w + p + m for w, p, m in subs)}"
            note = f"substitutions {', '.join(w + p + m for w, p, m in subs)} applied to the UniProt sequence; wild-type residues verified"
            fixes.append(dict(sequence_id=sid, uniprot_accession=acc, kind=kind, description=note,
                              start=0, end=len(seq), evidence="UniProt sequence + recorded mutation",
                              length=len(seq), sequence=seq))
        else:
            kind, note = "mutation_not_applied", f"mutation '{r.mutation}' could not be verified on the UniProt numbering"
    elif "complex_primary" in cav:
        kind, note = "complex_primary_subunit", f"first subunit only; complex = {r.target_symbols}"
    seq_of[gid], kind_of[gid], note_of[gid] = (sid, seq), kind, note

norm = lambda s: re.sub(r"[^a-z0-9]", "", str(s).casefold())


def label_check(r):
    lab = norm(f"{r.Target} {r.activity_header}")
    toks = [norm(x) for x in re.split(r"[;/| ,()]+", f"{r.mapped_identity} {r.gene}") if len(norm(x)) >= 3]
    toks += [norm(w) for w in re.findall(r"[A-Za-z0-9\-]{4,}", str(r.protein_name))]
    return "label_names_target" if any(t and t in lab for t in toks) else "label_does_not_name_target"


def xref_status(r):
    x = str(r.xref_accession or "")
    if not x:
        return "no_cross_reference"
    return "agrees" if r.uniprot_accession in x.split("|") else "differs"


a["gene"] = a.uniprot_accession.map(lambda k: str(tg.uniprot_entry.get(k, "")).split("_")[0])
out = pd.DataFrame({
    "group_id": a.group_id,
    "uniprot_accession": a.uniprot_accession,
    "sequence_id": a.group_id.map(lambda g: seq_of[g][0]),
    "sequence_kind": a.group_id.map(kind_of),
    "sequence_note": a.group_id.map(note_of),
    "sequence": a.group_id.map(lambda g: seq_of[g][1]),
    "mapped_identity": a.mapped_identity, "mutation": a.mutation, "target_symbols": a.target_symbols,
    "confidence": a.confidence,
    "label_check": a.apply(label_check, axis=1),
    "xref_status": a.apply(xref_status, axis=1),
})
out["seq_len"] = out.sequence.str.len()
out.to_parquet(f"{OUT}/target_sequences.parquet", index=False)
f = pd.DataFrame(fixes).drop_duplicates("sequence_id")
f.drop(columns=[]).to_csv(f"{OUT}/sequence_fixes.tsv", sep="\t", index=False)
print("groups:", len(out), "| sequence kinds:", out.sequence_kind.value_counts().to_dict())
print("fixes:", len(f))
print(f[["sequence_id", "kind", "description", "length"]].assign(description=f.description.str.slice(0, 60)).to_string(index=False))
