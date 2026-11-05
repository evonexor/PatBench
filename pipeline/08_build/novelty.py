"""How much of this extraction is absent from the public patent-activity databases?

Compared against BindingDB and ChEMBL. HARVEST is deliberately excluded: its public snapshot is
only the H-bench subset (48 targets, 1,617 patents), so treating it as a coverage baseline would
understate how much of our data is new.

Molecule and patent sets are cached per source, so changing the comparison set later is cheap.
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import json
import re
import os
import sys

import pandas as pd
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
R = os.environ.get("PATBENCH_SAR_DIR", _ws + "/outputs/uspto-sar-recovery-20261006/results/")
M = os.environ.get("PATBENCH_TARGET_MAP_DIR", _ws + "/outputs/target-uniprot-map-v2-20261019/")
S = _ws + "/raw/patent-data/standardized/"
OUT = sys.argv[1].rstrip("/")
SOURCES = sys.argv[2].split(",") if len(sys.argv) > 2 else ["BindingDB", "ChEMBL"]
AFF = ["IC50", "Ki", "Kd", "EC50"]
KEY = ["patent_core", "Target", "Activity Type"]


def patent_numbers(s):
    """Every 6-8 digit run in an identifier, so granted and publication numbers both hit."""
    return {str(int(x)) for x in re.findall(r"\d{6,8}", str(s))}


def source_sets(name):
    """InChIKeys, patent numbers and (molecule, accession) pairs for one public source."""
    cache = f"{OUT}/work/pubset_{name}.parquet"
    try:
        t = pq.read_table(cache).to_pandas()
        print(f"{name}: reusing cache {cache}", flush=True)
        mols = set(t.loc[t.kind == "molecule", "value"])
        pats = set(t.loc[t.kind == "patent", "value"])
        pairs = {tuple(v.split("\t")) for v in t.loc[t.kind == "pair", "value"]}
        return mols, pats, pairs
    except Exception:
        pass
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
    mols = set(key.values())
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
    print(f"   -> {len(mols):,} molecules, {len(pats):,} patents, {len(pairs):,} pairs", flush=True)
    pd.DataFrame(
        [("molecule", x) for x in mols] + [("patent", x) for x in pats]
        + [("pair", f"{a}\t{b}") for a, b in pairs],
        columns=["kind", "value"]).to_parquet(cache, index=False)
    return mols, pats, pairs


pub_mols, pub_patents, pub_pairs = set(), set(), set()
per_source = {}
for name in SOURCES:
    mols, pats, pairs = source_sets(name)
    per_source[name] = {"molecules": len(mols), "patents": len(pats), "pairs": len(pairs)}
    pub_mols |= mols
    pub_patents |= pats
    pub_pairs |= pairs
print(f"\nunion of {', '.join(SOURCES)}: {len(pub_mols):,} molecules, "
      f"{len(pub_patents):,} patent numbers, {len(pub_pairs):,} pairs", flush=True)

# ------------------------------------------------------------- our side
m = ds.dataset(R + "SAR_merged.parquet").to_table(
    columns=["patent_core", "Target", "Activity Type", "InChIKey"]).to_pandas()
m = m[m["Activity Type"].isin(AFF)].copy()
m["Target"] = m.Target.fillna("").str.strip()

ours_mols = set(m.InChIKey.dropna()) - {""}
new_mols = ours_mols - pub_mols
pat_new_map = {p: not (patent_numbers(p) & pub_patents) for p in m.patent_core.unique()}
m["molecule_is_new"] = ~m.InChIKey.isin(pub_mols)
m["patent_is_new"] = m.patent_core.map(pat_new_map)

tmap = pq.read_table(M + "target_map.parquet",
                     columns=["patent_core", "Target", "uniprot_accession", "confidence",
                              "target_class"]).to_pandas()
tmap = tmap[tmap.uniprot_accession.fillna("").ne("")].drop_duplicates(["patent_core", "Target"])
m = m.join(tmap.set_index(["patent_core", "Target"])[["uniprot_accession", "confidence", "target_class"]],
           on=["patent_core", "Target"])
mapped = m[m.uniprot_accession.fillna("").ne("") & m.confidence.isin(["high", "medium"])]
pair_new = pd.Series([(k, u) not in pub_pairs for k, u in
                      zip(mapped.InChIKey, mapped.uniprot_accession)], index=mapped.index)

report = {
    "compared_against": SOURCES,
    "note": "HARVEST excluded: its public snapshot is only the H-bench subset (48 targets).",
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
        "patents_checked": int(m.patent_core.nunique()),
        "new_patents": int(sum(pat_new_map.values())),
        "new_patent_share": round(sum(pat_new_map.values()) / len(pat_new_map), 4),
        "rows_in_new_patents": int(m.patent_is_new.sum()),
    },
    "pair_level_mapped_targets_only": {
        "mapped_rows": len(mapped),
        "new_pair_rows": int(pair_new.sum()),
        "new_pair_share": round(float(pair_new.mean()), 4),
    },
}
print("\n" + json.dumps(report, indent=1))
with open(OUT + "/novelty_vs_public.json", "w") as fh:
    json.dump(report, fh, indent=1)
m[["patent_core", "Activity Type", "InChIKey", "uniprot_accession",
   "molecule_is_new", "patent_is_new"]].to_parquet(OUT + "/work/novelty_flags.parquet", index=False)
print(f"\nwrote {OUT}/novelty_vs_public.json and work/novelty_flags.parquet")
