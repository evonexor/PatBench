"""target_audit.parquet: the join of join_target_map.py (one row per ranking-task group, with the target-map pair it was
mapped from) plus, for each of the five ranking tasks, whether the group is in it, its time split there and its test tier.

    python make_target_audit.py <release_dir> <join_target_map output> <out target_audit.parquet>
"""
import sys

import pandas as pd

R, JOIN, OUT = sys.argv[1].rstrip("/"), sys.argv[2], sys.argv[3]
TASKS = ["rank_ic50", "rank_ki_kd", "rank_ec50", "rank_highcontrast", "rank_complex_mutant"]
a = pd.read_parquet(JOIN)
for t in TASKS:
    d = pd.read_parquet(f"{R}/tasks/{t}.parquet", columns=["group_id", "split_time", "test_tier"]).drop_duplicates("group_id")
    d = d.set_index("group_id")
    a[f"in_{t}"] = a.group_id.isin(d.index)
    a[f"split_{t}"] = a.group_id.map(d.split_time)
    a[f"tier_{t}"] = a.group_id.map(d.test_tier)
a.to_parquet(OUT, index=False)
print(len(a), "groups,", a.shape[1], "columns ->", OUT)
