"""Score every baseline prediction file and print the leaderboard.

    python leaderboard.py <release dir> [task ...]

Reads baselines/<task>_<model>[_k<K>].parquet, scores each against the task's time-split test set,
and writes baselines/leaderboard.csv. Zero-shot rows also carry cliff accuracy.
"""
import glob
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evaluate as E

D = sys.argv[1].rstrip("/")
TASKS = sys.argv[2:] or ["rank_ic50", "rank_ki_kd", "rank_highcontrast", "rank_complex_mutant"]
SPLIT, PART = "split_time", "test"
rows = []
for task in TASKS:
    full = pd.read_parquet(f"{D}/tasks/{task}.parquet")
    te = full[full[SPLIT] == PART]
    cliffs = pd.read_parquet(f"{D}/protocols/{task}_cliffs.parquet")
    for f in sorted(glob.glob(f"{D}/baselines/{task}_*.parquet")):
        name = os.path.basename(f)[len(task) + 1:-len(".parquet")]
        m = re.match(r"(.+?)_k(\d+)$", name)
        model, k = (m.group(1), int(m.group(2))) if m else (name, 0)
        pred = pd.read_parquet(f)
        if k:
            proto = pd.read_parquet(f"{D}/protocols/{task}_{SPLIT.replace('split_', '')}_k{k}.parquet")
            gm = E.per_group_fewshot(te, pred, proto)
        else:
            gm = E.per_group(te, pred)
        s = E.summary(gm)
        rec = {"task": task, "model": model, "shots": k, "groups": s["groups"],
               "spearman": round(s["spearman_median"], 3),
               "ci": f"[{s['spearman_ci_lo']:.3f}, {s['spearman_ci_hi']:.3f}]",
               "pair_acc": round(s["pair_accuracy_median"], 3),
               "ndcg10": round(s["ndcg10_median"], 3),
               "top10pct": round(s["top10pct_hit_mean"], 3)}
        if not k:
            c = E.cliff_accuracy(te, pred, cliffs, set(te.group_id))
            rec["cliff_acc"] = round(c["accuracy"], 3) if c else None
            rec["cliff_pairs"] = c["pairs"] if c else 0
        rows.append(rec)
        print(f"  {task:20s} {model:14s} k={k:<3d} rho {rec['spearman']:+.3f}", flush=True)

lb = pd.DataFrame(rows)
lb.to_csv(f"{D}/baselines/leaderboard.csv", index=False)
for task, g in lb.groupby("task"):
    print(f"\n=== {task} | {SPLIT}={PART}")
    print(g.drop(columns="task").sort_values(["shots", "spearman"], ascending=[True, False])
          .to_string(index=False, na_rep="-"))
