"""Evaluation sets for models that were NOT trained on PatBench.

The time split keeps a small test set (rank_ic50: 84 groups, rank_ki_kd: 10) because the rest of the
task must stay available for training. A model that is not trained on PatBench has no use for that
reserve: every group of the task can serve as test data. This script writes eval_all/<task>.parquet,
a copy of each ranking task in which

  split_all       is "test" on every row, so evaluate.py --split split_all --part test scores all groups
  public_share    is the share of the group's rows whose (molecule, target) pair is in BindingDB/ChEMBL
  test_tier       is recomputed over ALL groups: "novel" when public_share is 0, "reference" otherwise
  test_tier_time  keeps the tier of the time split ("" outside the time-split test)
  priority_year   is the group's priority year, for stratifying by age

A model trained on PatBench, or tuned on any split of it, must not use these files: they contain the
training rows. A model pretrained on BindingDB or ChEMBL should be judged on the novel tier.

    python make_eval_sets.py <release_dir>
    python evaluate.py eval_all/rank_ic50.parquet predictions.parquet --split split_all --part test --tier novel
"""
import os
import sys

import pandas as pd

D = sys.argv[1].rstrip("/")
os.makedirs(D + "/eval_all", exist_ok=True)

for name in ["rank_ic50", "rank_ki_kd", "rank_ec50", "rank_highcontrast", "rank_complex_mutant"]:
    d = pd.read_parquet(f"{D}/tasks/{name}.parquet")
    share = d.groupby("group_id").in_public_data.mean()
    d["public_share"] = d.group_id.map(share)
    d["test_tier_time"] = d.test_tier
    d["test_tier"] = (d.public_share == 0).map({True: "novel", False: "reference"})
    d["split_all"] = "test"
    d["priority_year"] = (pd.to_datetime(d.priority_date, errors="coerce").dt.year
                          .fillna(d.year).astype(int))
    d.to_parquet(f"{D}/eval_all/{name}.parquet", index=False)
    g = d.drop_duplicates("group_id")
    nov = d[d.test_tier == "novel"]
    print(f"{name:20s} all: {len(g):4d} groups {d.patent_core.nunique():4d} patents "
          f"{d.uniprot_accession.nunique():4d} targets {d.InChIKey.nunique():6d} structures | "
          f"novel: {nov.group_id.nunique():4d} groups {nov.patent_core.nunique():4d} patents "
          f"{nov.uniprot_accession.nunique():4d} targets {nov.InChIKey.nunique():6d} structures")
