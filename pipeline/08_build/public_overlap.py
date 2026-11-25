"""Cache the (InChIKey, UniProt) pairs present in the public patent-activity datasets.

Reads the standardized BindingDB / ChEMBL / HARVEST tables that already sit on this host and
writes one row per (InChIKey, UniProt accession, source, activity type) with the median pAct.
Nothing here touches the USPTO extraction outputs.
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import json
import re
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

S = _ws + "/raw/patent-data/standardized/"
OUT = sys.argv[1]
SOURCES = ["BindingDB", "ChEMBL", "HARVEST"]
COLS = ["Patent", "SMILES", "UniProt", "Activity Type", "Relation", "Value", "Unit"]
TO_NM = {"nM": 1.0, "uM": 1e3, "µM": 1e3, "μM": 1e3, "mM": 1e6, "M": 1e9, "pM": 1e-3}


def patent_number(s):
    m = re.search(r"(\d{6,8})", str(s))
    return str(int(m.group(1))) if m else None


frames = []
patents = {}
for name in SOURCES:
    t = pq.read_table(S + name + ".parquet", columns=COLS).to_pandas()
    t["pnum"] = t.Patent.map(patent_number)
    patents[name] = sorted(set(t.pnum.dropna()))
    print(f"{name}: {len(t):,} rows, {len(patents[name]):,} patent numbers, "
          f"sample Patent values {t.Patent.dropna().head(3).tolist()}", flush=True)

    smiles = t.SMILES.dropna().unique()
    print(f"  computing InChIKeys for {len(smiles):,} unique SMILES", flush=True)
    key = {}
    for i, s in enumerate(smiles):
        if i and i % 100000 == 0:
            print(f"    {i:,} / {len(smiles):,}", flush=True)
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            continue
        try:
            key[s] = Chem.MolToInchiKey(mol)
        except Exception:
            pass
    print(f"  resolved {len(key):,} of {len(smiles):,}", flush=True)

    t["ik"] = t.SMILES.map(key)
    t = t[t.ik.notna() & t.UniProt.fillna("").ne("")]
    factor = t.Unit.map(TO_NM)
    nM = pd.to_numeric(t.Value, errors="coerce") * factor
    t["pAct"] = np.where(nM > 0, 9 - np.log10(nM.where(nM > 0)), np.nan)
    # a ChEMBL/BindingDB complex target lists several accessions joined by '|'
    t = t.assign(uniprot=t.UniProt.astype(str).str.split("|")).explode("uniprot")
    t["uniprot"] = t.uniprot.str.strip()
    t = t[t.uniprot.ne("")]
    t["source"] = name
    frames.append(t[["ik", "uniprot", "Activity Type", "source", "pAct", "pnum"]])

p = pd.concat(frames, ignore_index=True)
agg = (p.groupby(["ik", "uniprot", "Activity Type", "source"], dropna=False)
        .agg(pAct=("pAct", "median"), n=("pAct", "size")).reset_index())
agg.to_parquet(OUT + "/public_pairs.parquet", index=False)
with open(OUT + "/public_patents.json", "w") as fh:
    json.dump({k: v for k, v in patents.items()}, fh)
print(f"\nwrote {len(agg):,} (molecule, target, endpoint, source) rows")
print("distinct (molecule, target) pairs:", f"{agg.groupby(['ik', 'uniprot']).ngroups:,}")
print("per source:", agg.source.value_counts().to_dict())
