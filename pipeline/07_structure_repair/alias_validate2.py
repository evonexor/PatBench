"""Check a repair table against the molecules BindingDB / ChEMBL curated from the same patent (independent of the MOL files and of the rules).

    python alias_validate2.py <alias_dir> <records.tsv> <repairs.tsv> <public_keys.tsv> <report_prefix> [r5_dir]

hit = the connectivity block (first 14 characters of the InChIKey) of the molecule, or of its largest fragment, is among the public molecules of the same patent.
Reported for the repaired structure, for the structure stored now, and for a control: molecules of the same patents with no alias line at all (the rate at
which a correct structure of such a patent is found). Everything is computed only for patents that BindingDB/ChEMBL cover.
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
O = Path(sys.argv[1]).resolve()
REC, REP, PUB = (O / sys.argv[2]), (O / sys.argv[3]), (O / sys.argv[4])
PREFIX = sys.argv[5]
R5 = sys.argv[6] if len(sys.argv) > 6 else _ws + "/outputs/patbench-v1.0-20261030-r5"


def digits(s):
    return next(iter(sorted({str(int(x)) for x in re.findall(r"\d{6,8}", str(s))})), "")


def keys14(smiles):
    m = Chem.MolFromSmiles(smiles) if smiles else None
    if m is None:
        return set()
    out = set()
    try:
        out.add(Chem.MolToInchiKey(m)[:14])
        out.add(Chem.MolToInchiKey(max(Chem.GetMolFrags(m, asMols=True), key=lambda f: f.GetNumHeavyAtoms()))[:14])
    except Exception:
        pass
    return out


pub = pd.read_csv(PUB, sep="\t", dtype=str)
pubset = defaultdict(set)
for d, k in zip(pub.patent_digits, pub.inchikey14):
    pubset[d].add(k)
rec = pd.read_csv(REC, sep="\t", dtype=str, keep_default_na=False)
rep = pd.read_csv(REP, sep="\t", dtype=str, keep_default_na=False)
m = rec.merge(rep[["patent_core", "structure_file", "status", "smiles", "rules", "labels_expanded", "labels_unresolved"]], on=["patent_core", "structure_file"], how="left")
m["pd"] = m.patent_core.map(digits)
m = m[m.pd.isin(pubset)]
print(f"records {len(rec):,}; in patents with public data: {len(m):,}; status there: {m.status.value_counts().to_dict()}")

ex = m[m.status.isin(["expanded", "as_drawn", "counterion_only"]) & m.smiles.ne("")].copy()
ex["new_hit"] = [bool(keys14(s) & pubset[p]) for s, p in zip(ex.smiles, ex.pd)]
ex["old_hit"] = [bool(keys14(s) & pubset[p]) for s, p in zip(ex.SMILES_as_drawn, ex.pd)]

s = pq.read_table(R5 + "/sar_all.parquet", columns=["patent_core", "structure_file", "SMILES_as_drawn", "InChIKey_as_drawn"]).to_pandas()
s = s[~s.structure_file.isin(set(rec.structure_file))].drop_duplicates(["patent_core", "InChIKey_as_drawn"])
s["pd"] = s.patent_core.map(digits)
s = s[s.pd.isin(set(ex.pd))]
s = s.sample(min(len(s), 8000), random_state=1)
s["hit"] = [bool(keys14(x) & pubset[p]) for x, p in zip(s.SMILES_as_drawn, s.pd)]
base = s.groupby("pd").hit.mean()
ex["base_patent"] = ex.pd.map(base)
print(f"\ncontrol (no alias line), same patents, {len(s):,} sampled: found {s.hit.mean():.1%}")
b = ex[ex.base_patent.notna()]
print(f"repaired-or-as-drawn records: {len(ex):,} (with a control sample for their patent: {len(b):,})")
print(f"  stored structure found: {ex.old_hit.mean():.1%}   repaired structure found: {ex.new_hit.mean():.1%}   (control rate in those patents: {b.base_patent.mean():.1%}; repaired there: {b.new_hit.mean():.1%})")

def grp(df, key):
    g = df.groupby(key).agg(n=("new_hit", "size"), stored=("old_hit", "mean"), repaired=("new_hit", "mean"), control=("base_patent", "mean"))
    return g.sort_values("n", ascending=False)

pd.set_option("display.width", 200); pd.set_option("display.max_rows", 60)
print("\nby rule set:"); print(grp(ex, "rules").round(3).to_string())
ex["labs"] = ex.labels_expanded.map(lambda t: ";".join(sorted({x.split("->")[0].split(":", 1)[1] for x in t.split(";") if x})))
print("\nby label set (top 22):"); print(grp(ex[ex.labs.ne("")], "labs").head(22).round(3).to_string())
lost = ex[ex.old_hit & ~ex.new_hit]
print(f"\nrecords where the stored structure was found but the repaired one is not: {len(lost)}")
ex.to_csv(O / f"{PREFIX}_validated.tsv", sep="\t", index=False)
