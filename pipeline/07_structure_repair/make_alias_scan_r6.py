"""work/alias_scan.parquet for the r6 build: the (patent, molecule) pairs whose stored structure still rests on an alias nobody could expand.

    python make_alias_scan_r6.py <r5_dir> <r6_dir> <alias_dir> <all_records.tsv> <repairs.tsv>

The build flags a pair as structure_alias_suspect when it is listed here with n_bad == n_mol. A pair is listed when a record of that patent with an
unresolved alias has this stored InChIKey and no clean drawing (a file without alias lines, or a record that needed no change) gives the same key.
"""
import sys
import pandas as pd
import pyarrow.parquet as pq
R5, R6, A = [p.rstrip("/") + "/" for p in sys.argv[1:4]]
rec = pd.read_csv(A + sys.argv[4], sep="\t", dtype=str, keep_default_na=False)
rep = pd.read_csv(A + sys.argv[5], sep="\t", dtype=str, keep_default_na=False)
m = rec.merge(rep[["patent_core", "structure_file", "status", "inchikey", "labels_unresolved"]], on=["patent_core", "structure_file"])
s = pq.read_table(R5 + "sar_all.parquet", columns=["Patent", "structure_file", "InChIKey_as_drawn"]).to_pandas()
alias_files = set(rec.structure_file)
plain = s[~s.structure_file.isin(alias_files)]
clean = set(zip(plain.Patent, plain.InChIKey_as_drawn))
same = m[(m.status == "expanded") & (m.inchikey == m.InChIKey_as_drawn)]
clean |= set(zip(same.Patent, same.InChIKey_as_drawn))
un = m[m.status == "unresolved"].copy()
un = un[[(p, k) not in clean for p, k in zip(un.Patent, un.InChIKey_as_drawn)]]
out = pd.DataFrame({"Patent": un.Patent.values, "InChIKey": un.InChIKey_as_drawn.values, "n_mol": 1, "n_bad": 1,
                    "labels": un.labels_unresolved.str.slice(0, 200).values, "error": ""}).drop_duplicates(["Patent", "InChIKey"])
out.to_parquet(R6 + "work/alias_scan.parquet", index=False)
print("suspect (patent, molecule) pairs:", len(out), " unresolved records:", int((m.status == "unresolved").sum()))
