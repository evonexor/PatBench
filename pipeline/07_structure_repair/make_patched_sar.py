"""Write a copy of the upstream recovery results in which the structures are repaired. The original directory is not touched.

    python make_patched_sar.py <r5_dir> <repairs.tsv> <out_results_dir>

<repairs.tsv>  alias_repair.py output for all records with an alias line (status expanded / as_drawn / counterion_only / unresolved).
<out_results_dir>  gets SAR_merged.parquet (single file) with SMILES and InChIKey replaced where a repair was made, plus three columns:
    structure_repair         "" or the rules used: group, NH, nitrile (joined with ;)
    SMILES_before_repair     the SMILES that was stored before (only where repaired)
    InChIKey_before_repair
and symlinks to every other file of the original results directory.

Rows that have a MOL file with an alias line take the repaired structure of that file (status expanded; the nitrile step is already in it). Every other row
whose SMILES holds the reversed-nitrile pattern N(+)#C(-) is repaired from the SMILES (rule nitrile). Rows of unresolved records are left as stored; the build flags them.
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import os
import re
import sys
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
R5, REPAIRS, OUT = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve(), Path(sys.argv[3]).resolve()
SRC = Path(_ws + "/outputs/uspto-sar-recovery-20261006/results")
NITRILE = Chem.MolFromSmarts("[#6][N+]#[C-]")
ISO = re.compile(r"\[N\+\]#\[C-\]|\[C-\]#\[N\+\]")

sm = ds.dataset(str(SRC / "SAR_merged.parquet")).to_table().to_pandas()
assert sm.source_record_id.is_unique
print("SAR_merged rows:", len(sm), flush=True)
sar5 = pq.read_table(R5 / "sar_all.parquet", columns=["source_record_id", "structure_file"]).to_pandas()
sm = sm.merge(sar5, on="source_record_id", how="left")
rep0 = pd.read_csv(REPAIRS, sep="\t", dtype=str, keep_default_na=False)
unresolved = set(zip(rep0[rep0.status == "unresolved"].patent_core, rep0[rep0.status == "unresolved"].structure_file))
rep = rep0[rep0.status == "expanded"].set_index(["patent_core", "structure_file"])
key = list(zip(sm.patent_core, sm.structure_file.fillna("")))
sm["structure_repair"] = ""
sm["SMILES_before_repair"] = ""
sm["InChIKey_before_repair"] = ""
hit = [k in rep.index for k in key]
sm["_hit"] = hit
idx = sm.index[sm._hit]
new = rep.loc[[key[i] for i in idx]]
sm.loc[idx, "SMILES_before_repair"] = sm.loc[idx, "SMILES"].values
sm.loc[idx, "InChIKey_before_repair"] = sm.loc[idx, "InChIKey"].values
sm.loc[idx, "SMILES"] = new.smiles.values
sm.loc[idx, "InChIKey"] = new.inchikey.values
sm.loc[idx, "structure_repair"] = new.rules.values
print("rows repaired from a MOL file:", len(idx), " distinct files:", sm.loc[idx, "structure_file"].nunique(), flush=True)

# the reversed nitrile on rows that have no alias line
# rows of unresolved records stay as stored (the build flags them by their stored InChIKey)
sm["_unres"] = [k in unresolved for k in key]
rest = sm[~sm._hit & ~sm._unres & sm.SMILES.fillna("").map(lambda s: bool(ISO.search(s)))]
fix = {}
for smi in rest.SMILES.unique():
    m = Chem.MolFromSmiles(smi)
    if m is None:
        continue
    changed = False
    for _ in range(10):
        mt = m.GetSubstructMatch(NITRILE)
        if not mt:
            break
        _, n, c = mt
        m.GetAtomWithIdx(n).SetAtomicNum(6); m.GetAtomWithIdx(n).SetFormalCharge(0)
        m.GetAtomWithIdx(c).SetAtomicNum(7); m.GetAtomWithIdx(c).SetFormalCharge(0)
        Chem.SanitizeMol(m)
        changed = True
    if changed:
        fix[smi] = (Chem.MolToSmiles(m), Chem.MolToInchiKey(m))
j = rest.index[rest.SMILES.isin(fix)]
sm.loc[j, "SMILES_before_repair"] = sm.loc[j, "SMILES"].values
sm.loc[j, "InChIKey_before_repair"] = sm.loc[j, "InChIKey"].values
sm.loc[j, "InChIKey"] = [fix[s][1] for s in sm.loc[j, "SMILES"]]
sm.loc[j, "SMILES"] = [fix[s][0] for s in sm.loc[j, "SMILES"]]
sm.loc[j, "structure_repair"] = "nitrile"
print("rows repaired from the SMILES (nitrile only):", len(j), flush=True)
sm = sm.drop(columns=["_hit", "_unres", "structure_file"])
print("rows with a repair:", int(sm.structure_repair.ne("").sum()), " distinct drawn molecules before -> after:",
      sm.InChIKey_before_repair[sm.structure_repair.ne("")].nunique(), "->", sm.InChIKey[sm.structure_repair.ne("")].nunique())

OUT.mkdir(parents=True, exist_ok=True)
(OUT / "SAR_merged.parquet").mkdir(exist_ok=True)
pq.write_table(pa.Table.from_pandas(sm, preserve_index=False), OUT / "SAR_merged.parquet" / "part-00000.parquet")
for p in SRC.iterdir():
    if p.name == "SAR_merged.parquet":
        continue
    dst = OUT / p.name
    if not dst.exists():
        os.symlink(p, dst)
print("written", OUT)
