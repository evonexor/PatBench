"""Why does no table print the claimed value for some sampled items?

A failure here can be ours (markup we do not parse, number formatting, unit conversion) or the
data's. Telling them apart before annotation starts is the difference between measuring the
dataset and measuring our own evidence extractor.
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import html
import re
import unicodedata
import zipfile
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

pd.set_option("display.width", 220)
W = Path(_ws + "/outputs/sar-validation-20261019")
RAW = Path(_ws + "/raw/patent-data/")
R = _ws + "/outputs/uspto-sar-recovery-20261006/results/"

s = pq.read_table(W / "sample.parquet").to_pandas()
print("=== table located, by source format")
print(pd.crosstab(s.source_format, s.table_found, margins=True).to_string())
print("\n=== table located, by source")
print(pd.crosstab(s.source, s.table_found, margins=True).to_string())
print("\n=== table located, by endpoint")
print(pd.crosstab(s.endpoint, s.table_found, margins=True).to_string())

miss = s[~s.table_found]
print(f"\n{len(miss)} items with no table hit; {int(miss.has_recorded_locator.sum())} of them do "
      f"have a recorded cell locator from the pipeline")

# are the stored Value/Unit the printed pair, or the converted one?
pt = pq.read_table(R + "patents.parquet", columns=["Patent", "archive"]).to_pandas()
arch = dict(zip(pt.Patent, pt.archive))


def clean(t):
    t = re.sub(r"<sub>(.*?)</sub>", r"\1", t, flags=re.S)
    t = re.sub(r"<[^>]+>", " ", t)
    return re.sub(r"\s+", " ", unicodedata.normalize("NFKC", html.unescape(t))).strip()


def whole_text(patent):
    try:
        with zipfile.ZipFile(RAW / arch[patent]) as z:
            n = [x for x in z.namelist() if x.upper().endswith((".XML", ".SGM", ".SGML", ".TXT"))]
            return clean(z.read(n[0]).decode("utf-8", "replace")) if n else ""
    except Exception:
        return ""


def markup_kinds(patent):
    try:
        with zipfile.ZipFile(RAW / arch[patent]) as z:
            n = [x for x in z.namelist() if x.upper().endswith((".XML", ".SGM", ".SGML", ".TXT"))]
            raw = z.read(n[0]).decode("utf-8", "replace") if n else ""
    except Exception:
        return {}
    return {
        "<tables id=": len(re.findall(r"<tables\b[^>]*\bid=", raw)),
        "<tables no id": len(re.findall(r"<tables\b(?![^>]*\bid=)", raw)),
        "<table ": len(re.findall(r"<table\b", raw)),
        "<TABLE upper": len(re.findall(r"<TABLE\b", raw)),
        "CWU/maths": len(re.findall(r"<\?CWU|<maths", raw)),
    }


rows = []
for _, r in miss.head(40).iterrows():
    txt = whole_text(r.Patent)
    forms = {f"{r.Value:g}", f"{r.Value:.2f}", f"{r.Value:.3f}", f"{r.Value:.1f}"}
    forms |= {f.rstrip("0").rstrip(".") for f in list(forms) if "." in f}
    in_text = any(f in txt for f in forms if f)
    # the value as it would read if printed in another unit
    alt = {}
    for name, factor in [("nM", 1.0), ("uM", 1e-3), ("mM", 1e-6), ("M", 1e-9), ("pM", 1e3)]:
        v = r.value_nM * factor
        alt[name] = any(x in txt for x in {f"{v:g}", f"{v:.2f}", f"{v:.1f}"} if x)
    rows.append({"item": r.item_id, "patent": r.Patent, "fmt": r.source_format,
                 "endpoint": r.endpoint, "Value": r.Value, "Unit": r.Unit,
                 "value_nM": r.value_nM, "printed_value_in_text": in_text,
                 "value_found_as": ",".join(k for k, v in alt.items() if v) or "-",
                 "markup": markup_kinds(r.Patent)})
d = pd.DataFrame(rows)
print("\n=== first 40 misses: is the number anywhere in the patent text?")
print(d[["item", "patent", "fmt", "endpoint", "Value", "Unit", "printed_value_in_text",
         "value_found_as"]].to_string(index=False))
print("\nprinted value present in full text:", int(d.printed_value_in_text.sum()), "of", len(d))
print("\n=== table markup counts in those patents")
mk = pd.DataFrame(list(d.markup))
print(mk.describe().loc[["mean", "min", "max"]].round(1).to_string())
print("\npatents with zero '<tables id=' but some '<table':",
      int(((mk["<tables id="] == 0) & (mk["<table "] > 0)).sum()), "of", len(mk))
