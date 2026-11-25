"""Two questions that must be settled before annotation starts.

1. What does legacy table markup actually look like, so the renderer can handle it?
2. Is the stored Value/Unit the pair printed in the patent, or a converted one? The annotation
   instruction is wrong if we guess this, and every annotator would repeat the error.
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds

pd.set_option("display.width", 220)
RAW = Path(_ws + "/raw/patent-data/")
R = _ws + "/outputs/uspto-sar-recovery-20261006/results/"

pt = ds.dataset(R + "patents.parquet").to_table(columns=["Patent", "archive", "source_format"]).to_pandas()
arch = dict(zip(pt.Patent, pt.archive))

print("=== 1. legacy table markup, US6225330B1")
with zipfile.ZipFile(RAW / arch["US6225330B1"]) as z:
    n = [x for x in z.namelist() if x.upper().endswith((".XML", ".SGM", ".SGML", ".TXT"))]
    raw = z.read(n[0]).decode("utf-8", "replace")
print("file:", n[0], len(raw), "chars")
m = re.search(r"<TABLE\b.*?</TABLE>", raw, re.S)
print("first <TABLE> block:\n", (m.group(0)[:1200] if m else "none"))
print("\ntags inside:", sorted(set(re.findall(r"<(/?[A-Za-z][\w-]*)", m.group(0))))[:30] if m else None)

print("\n=== 2. is Value/Unit the printed pair?")
TO_NM = {"nM": 1.0, "uM": 1e3, "mM": 1e6, "M": 1e9, "pM": 1e-3}
d = ds.dataset(R + "new_SAR.parquet").to_table(
    columns=["raw_value", "Value", "Unit", "value_nM", "activity_header",
             "unit_evidence"]).to_pandas()
num = pd.to_numeric(d.raw_value.str.extract(r"([\d.]+)")[0], errors="coerce")
d = d.assign(printed=num)
ok = d[d.printed.notna() & d.Value.notna()]
same = np.isclose(ok.printed, ok.Value, rtol=1e-6)
print(f"rows where stored Value equals the printed number: {int(same.sum()):,} of {len(ok):,} "
      f"= {same.mean():.1%}")
conv = np.isclose(ok.value_nM, ok.printed * ok.Unit.map(TO_NM).fillna(np.nan), rtol=1e-6)
print(f"rows where value_nM = printed x unit factor:        {int(conv.sum()):,} of {len(ok):,} "
      f"= {conv.mean():.1%}")
bad = ok[~same]
print(f"\n{len(bad):,} rows where Value differs from the printed number; a sample:")
print(bad.head(12)[["raw_value", "Value", "Unit", "value_nM", "activity_header"]].to_string(index=False))
print("\nunit stated in the header vs stored Unit, for those rows:")
hdr_unit = bad.activity_header.str.extract(r"[\(\[]\s*([numpM]{1,2}M|M)\s*[\)\]]", flags=re.I)[0]
print(pd.crosstab(hdr_unit.fillna("none in header"), bad.Unit, margins=True).to_string())
print("\n=== unit evidence kinds overall")
print(d.unit_evidence.value_counts().head(8).to_string())
