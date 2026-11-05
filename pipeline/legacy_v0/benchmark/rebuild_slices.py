"""Rebuild baseline_slices.parquet from the saved per-group metric files.

Only needed because the first write mixed string and boolean slice labels in one column.
"""
import glob
import os
import sys

import pandas as pd

import evaluate as ev

OUT = sys.argv[1].rstrip("/")
frames = []
for path in sorted(glob.glob(f"{OUT}/work/groups_*.parquet")):
    name = os.path.basename(path)[len("groups_"):-len(".parquet")]
    gm = pd.read_parquet(path)
    slices = pd.concat([
        ev.aggregate(gm, "split").assign(slice_kind="split"),
        ev.aggregate(gm[gm.split == "test"], "test_tier").assign(slice_kind="test tier"),
        ev.aggregate(gm, "Activity Type").assign(slice_kind="endpoint"),
        ev.aggregate(gm[gm.split == "test"], "unseen_target").assign(slice_kind="unseen target"),
        ev.aggregate(gm).assign(slice_kind="all", slice="all"),
    ], ignore_index=True)
    slices["slice"] = slices["slice"].astype(str)
    frames.append(slices.assign(model=name))
out = pd.concat(frames, ignore_index=True)
out.to_parquet(OUT + "/baseline_slices.parquet", index=False)
print(f"wrote {OUT}/baseline_slices.parquet  ({len(out)} rows, {out.model.nunique()} models)")

head = out[(out.slice_kind == "split") & (out["slice"] == "test")]
print("\n=== test split, median per-group Spearman (95% CI over groups)")
for _, r in head.sort_values("spearman_median", ascending=False).iterrows():
    print(f"  {r.model:20s} {r.spearman_median:+.3f}  [{r.spearman_ci_lo:+.3f}, {r.spearman_ci_hi:+.3f}]"
          f"   pair_acc {r.pair_accuracy_median:.3f}  ndcg10 {r.ndcg10_median:.3f}")
nov = out[(out.slice_kind == "test tier") & (out["slice"] == "novel")]
print("\n=== novel tier only (57 groups absent from public data)")
for _, r in nov.sort_values("spearman_median", ascending=False).iterrows():
    print(f"  {r.model:20s} {r.spearman_median:+.3f}  [{r.spearman_ci_lo:+.3f}, {r.spearman_ci_hi:+.3f}]"
          f"   pair_acc {r.pair_accuracy_median:.3f}")
