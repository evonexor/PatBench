"""How much of this extraction is absent from the public patent-activity databases?

This is the novelty claim a data descriptor has to carry, so it is measured three ways:
  molecule level   InChIKeys not present in BindingDB / ChEMBL / HARVEST at all
  patent level     patent numbers those databases never cite
  pair level       (molecule, UniProt) pairs absent, over the rows where our target is mapped
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import json
import re
import sys

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
R = _ws + "/outputs/uspto-sar-recovery-20261006/results/"
M = _ws + "/outputs/target-uniprot-map-20261016/"
S = _ws + "/raw/patent-data/standardized/"
OUT = sys.argv[1].rstrip("/")
AFF = ["IC50", "Ki", "Kd", "EC50"]
KEY = ["patent_core", "Target", "Activity Type"]


def patent_numbers(s):
    """Every 6-8 digit run in an identifier, so granted numbers and publication numbers both hit."""
    return {str(int(x)) for x in re.findall(r"\d{6,8}", str(s))}


# ------------------------------------------------------------ public side
pub_mols, pub_patents, pub_pairs = set(), set(), set()
per_source = {}
for name in ["BindingDB", "ChEMBL", "HARVEST"]:
    t = pq.read_table(S + name + ".parquet", columns=["Patent", "SMILES", "UniProt"]).to_pandas()
    smiles = t.SMILES.dropna().unique()
    print(f"{name}: {len(t):,} rows, {len(smiles):,} unique SMILES", flush=True)
    key = {}
    for i, s in enumerate(smiles):
        if i and i % 200000 == 0:
            print(f"   {i:,}", flush=True)
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            continue
        try:
            key[s] = Chem.MolToInchiKey(mol)
        except Exception:
            pass
    iks = set(key.values())
    pats = set()
    for v in t.Patent.dropna().unique():
        pats |= patent_numbers(v)
    t = t.assign(ik=t.SMILES.map(key))
    pairs = set()
    for u, k in zip(t.UniProt.fillna(""), t.ik):
        if k and u:
            for acc in str(u).split("|"):
                acc = acc.strip()
                if acc:
                    pairs.add((k, acc))
    per_source[name] = {"rows": len(t), "molecules": len(iks), "patents": len(pats), "pairs": len(pairs)}
    print(f"   -> {len(iks):,} molecules, {len(pats):,} patent numbers, {len(pairs):,} pairs", flush=True)
    pub_mols |= iks
    pub_patents |= pats
    pub_pairs |= pairs
print(f"\nunion: {len(pub_mols):,} molecules, {len(pub_patents):,} patent numbers, "
      f"{len(pub_pairs):,} (molecule, target) pairs", flush=True)

# ------------------------------------------------------------- our side
m = ds.dataset(R + "SAR_merged.parquet").to_table(
    columns=["patent_core", "Target", "Activity Type", "Relation", "value_nM",
             "InChIKey", "SMILES"]).to_pandas()
m = m[m["Activity Type"].isin(AFF)].copy()
m["Target"] = m.Target.fillna("").str.strip()

ours_mols = set(m.InChIKey.dropna()) - {""}
new_mols = ours_mols - pub_mols
our_pats = {p for v in m.patent_core.unique() for p in patent_numbers(v)}
new_pats = our_pats - pub_patents
m["molecule_is_new"] = ~m.InChIKey.isin(pub_mols)
m["patent_is_new"] = [not (patent_numbers(p) & pub_patents) for p in m.patent_core]

tm = pq.read_table(M + "ranking_group_targets.parquet").to_pandas()
m = m.join(tm.set_index(KEY)[["uniprot_accession", "confidence"]], on=KEY)
mapped = m[m.uniprot_accession.fillna("").ne("") & m.confidence.isin(["high", "medium"])]
pair_new = ~pd.Series([(k, u) in pub_pairs for k, u in
                       zip(mapped.InChIKey, mapped.uniprot_accession)], index=mapped.index)

report = {
    "our_affinity_rows": len(m),
    "our_molecules": len(ours_mols),
    "our_patents": int(m.patent_core.nunique()),
    "public_sources": per_source,
    "public_union": {"molecules": len(pub_mols), "patents": len(pub_patents), "pairs": len(pub_pairs)},
    "molecule_level": {
        "new_molecules": len(new_mols),
        "new_molecule_share": round(len(new_mols) / len(ours_mols), 4),
        "rows_on_new_molecules": int(m.molecule_is_new.sum()),
    },
    "patent_level": {
        "new_patents": len(new_pats),
        "patents_checked": len(our_pats),
        "rows_in_new_patents": int(m.patent_is_new.sum()),
        "distinct_new_patents_in_rows": int(m.loc[m.patent_is_new, "patent_core"].nunique()),
    },
    "pair_level_mapped_targets_only": {
        "mapped_rows": len(mapped),
        "new_pairs_rows": int(pair_new.sum()),
        "new_pair_share": round(float(pair_new.mean()), 4),
    },
}
print("\n" + json.dumps(report, indent=1))
with open(OUT + "/novelty_vs_public.json", "w") as fh:
    json.dump(report, fh, indent=1)

m[["patent_core", "Activity Type", "InChIKey", "uniprot_accession",
   "molecule_is_new", "patent_is_new"]].to_parquet(OUT + "/work/novelty_flags.parquet", index=False)
print(f"\nwrote {OUT}/novelty_vs_public.json and work/novelty_flags.parquet")
