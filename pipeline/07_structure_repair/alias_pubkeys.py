"""For the patents that have records with unexpanded MOL aliases: the molecules BindingDB and ChEMBL curated from the same patent.

    python alias_pubkeys.py <alias_dir> [records.tsv] [out.tsv]

Writes <alias_dir>/public_keys_by_patent.tsv (patent_digits, source, inchikey14, inchikey) where inchikey14 is the connectivity block of the
InChIKey. A structure that really is the molecule of the patent should show up here; one with missing atoms should not.
This is independent of the MOL files and of any expansion rule.
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import re
import sys
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
O = Path(sys.argv[1]).resolve()
RECORDS = sys.argv[2] if len(sys.argv) > 2 else "records.tsv"
OUTNAME = sys.argv[3] if len(sys.argv) > 3 else "public_keys_by_patent.tsv"
S = _ws + "/raw/patent-data/standardized/"


def digits(s):
    return {str(int(x)) for x in re.findall(r"\d{6,8}", str(s))}


rec = pd.read_csv(O / RECORDS, sep="\t", dtype=str, keep_default_na=False)
want = {d for p in rec.patent_core for d in digits(p)}
print("patents wanted:", len(want), flush=True)
rows = []
for src in ["BindingDB", "ChEMBL"]:
    t = pq.read_table(S + src + ".parquet", columns=["Patent", "SMILES"]).to_pandas().dropna()
    t["pd"] = t.Patent.map(lambda v: sorted(digits(v) & want))
    t = t[t.pd.map(len) > 0]
    print(src, "rows in wanted patents:", len(t), "unique SMILES:", t.SMILES.nunique(), flush=True)
    key = {}
    for smi in t.SMILES.unique():
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        try:
            key[smi] = Chem.MolToInchiKey(m)
        except Exception:
            pass
    for smi, pds in zip(t.SMILES, t.pd):
        k = key.get(smi)
        if k:
            for d in pds:
                rows.append((d, src, k[:14], k))
out = pd.DataFrame(rows, columns=["patent_digits", "source", "inchikey14", "inchikey"]).drop_duplicates()
out.to_csv(O / OUTNAME, sep="\t", index=False)
print("public molecules (patent, source, key):", len(out), " patents covered:", out.patent_digits.nunique(), " of", len(want))
