"""Compare two release directories table by table.

    python tools/compare_releases.py <reference_release> <rebuilt_release> [--skip name,name]

Used to check that rebuilding a release from the same inputs with this repository's code gives the
release that was published. For every parquet file present in both directories it compares shape,
column names and order, dtypes, and every value (floats to a relative tolerance of 1e-9, missing
values equal to missing values, row order significant). TSV/CSV/JSON files are compared byte for
byte, except release_manifest.json, which is compared key by key without the fields that record where
and when it was built. Exit status 1 when anything differs.
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

MANIFEST_IGNORED = {"built", "inputs", "log"}      # build date, absolute paths, console log with timings


def compare_frames(a, b, rtol=1e-9):
    """[] when equal, else a list of reasons."""
    if list(a.columns) != list(b.columns):
        return [f"columns differ: only in reference {sorted(set(a.columns) - set(b.columns))}, "
                f"only in rebuilt {sorted(set(b.columns) - set(a.columns))}" if set(a.columns) != set(b.columns)
                else "column order differs"]
    if len(a) != len(b):
        return [f"{len(a):,} rows vs {len(b):,}"]
    out = []
    for c in a.columns:
        x, y = a[c], b[c]
        if str(x.dtype) != str(y.dtype):
            out.append(f"{c}: dtype {x.dtype} vs {y.dtype}")
            continue
        if x.dtype.kind in "fc":
            bad = ~(((x - y).abs() <= rtol * x.abs().clip(lower=1e-300)) | (x == y) | (x.isna() & y.isna()))
        else:
            try:
                bad = ~((x == y) | (x.isna() & y.isna()))
            except Exception:       # arrays, lists in a cell
                bad = pd.Series([not _same(u, v) for u, v in zip(x, y)], index=x.index)
        n = int(bad.sum())
        if n:
            i = bad.idxmax()
            out.append(f"{c}: {n:,} of {len(x):,} values differ (first at row {i}: {x[i]!r} vs {y[i]!r})")
    return out


def _same(u, v):
    try:
        if u is None or v is None:
            return u is v
        if hasattr(u, "tolist"):
            u = u.tolist()
        if hasattr(v, "tolist"):
            v = v.tolist()
        return u == v or (u != u and v != v)
    except Exception:
        return False


def compare_manifests(a, b):
    out = []
    for k in sorted((set(a) | set(b)) - MANIFEST_IGNORED):
        if a.get(k) != b.get(k):
            out.append(f"{k}")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("reference")
    ap.add_argument("rebuilt")
    ap.add_argument("--skip", default="", help="comma separated path prefixes to leave out")
    ap.add_argument("--rtol", type=float, default=1e-9)
    a = ap.parse_args(argv)
    ref, new = Path(a.reference), Path(a.rebuilt)
    skip = tuple(s for s in a.skip.split(",") if s)
    files = sorted(p.relative_to(new) for p in new.rglob("*")
                   if p.is_file() and p.suffix in (".parquet", ".tsv", ".csv", ".json")
                   and not str(p.relative_to(new)).startswith(("work/",) + skip)
                   and not p.name.endswith((".log", ".report.json")))
    bad = ok = 0
    for rel in files:
        r = ref / rel
        if not r.exists():
            print(f"ONLY IN REBUILT  {rel}")
            bad += 1
            continue
        if rel.suffix == ".parquet":
            why = compare_frames(pd.read_parquet(r), pd.read_parquet(new / rel), a.rtol)
        elif rel.name == "release_manifest.json":
            why = [f"manifest key differs: {k}" for k in compare_manifests(json.loads(r.read_text()), json.loads((new / rel).read_text()))]
        else:
            why = [] if r.read_bytes() == (new / rel).read_bytes() else ["bytes differ"]
        print(f"{'OK  ' if not why else 'DIFF'} {rel}" + ("" if not why else "  <- " + "; ".join(why[:4])))
        ok += not why
        bad += bool(why)
    missing = [p for p in sorted(p.relative_to(ref) for p in ref.rglob("*.parquet"))
               if not str(p).startswith(("work/", "structures/", "eval_support/") + skip) and not (new / p).exists()]
    for p in missing:
        print(f"NOT REBUILT      {p}")
    print(f"\n{ok} identical, {bad} different, {len(missing)} reference tables not rebuilt")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
