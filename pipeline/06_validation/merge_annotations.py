"""Merge the two annotation passes into one per-item table, ready for adjudication."""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

W = sys.argv[1].rstrip("/")
REL = _ws + "/outputs/affinity-sar-release-20261019/"
KEY = ["patent_core", "Target", "Activity Type"]

s = pq.read_table(W + "/sample.parquet").to_pandas()
A = pd.read_csv(W + "/annotations/annotator_A.tsv", sep="\t", dtype=str).fillna("")
B = pd.read_csv(W + "/annotations/annotator_B.tsv", sep="\t", dtype=str).fillna("")
m = s.merge(A.add_suffix("_A").rename(columns={"item_id_A": "item_id"}), on="item_id").merge(
    B.add_suffix("_B").rename(columns={"item_id_B": "item_id"}), on="item_id")
for q in ["q1", "q2", "q3"]:
    m[f"{q}_consensus"] = np.where(m[f"{q}_A"] == m[f"{q}_B"], m[f"{q}_A"], "DISAGREE")
m["needs_adjudication"] = (m[["q1_consensus", "q2_consensus", "q3_consensus"]] == "DISAGREE").any(axis=1)

sar = pq.read_table(REL + "sar_all.parquet", columns=KEY + ["InChIKey", "heavy_atoms",
                                                           "group_id"]).to_pandas()
gmed = sar.groupby("group_id").heavy_atoms.median().rename("group_median_heavy_atoms")
sar = sar.join(gmed, on="group_id").drop_duplicates(KEY + ["InChIKey"])
m = m.rename(columns={"endpoint": "Activity Type"}).merge(
    sar[KEY + ["InChIKey", "group_median_heavy_atoms"]], on=KEY + ["InChIKey"], how="left")
m["group_is_fragment_suspect"] = m.group_median_heavy_atoms < 16

cols = ["item_id", "stratum", "source", "source_format", "confidence", "Patent", "Target",
        "Activity Type", "Value", "Unit", "InChIKey", "uniprot_accession", "mapped_identity",
        "mol_inchikey_matches", "group_median_heavy_atoms", "group_is_fragment_suspect",
        "q1_A", "q1_B", "q1_consensus", "q1_reason_A", "q1_reason_B",
        "q2_A", "q2_B", "q2_consensus", "q2_reason_A", "q2_reason_B",
        "q3_A", "q3_B", "q3_consensus", "q3_reason_A", "q3_reason_B",
        "needs_adjudication"]
m[cols].sort_values("item_id").to_csv(W + "/merged_annotations.tsv", sep="\t", index=False)
adj = m[m.needs_adjudication]
adj[cols].sort_values("item_id").to_csv(W + "/adjudication_queue.tsv", sep="\t", index=False)
print(f"merged {len(m)} items -> merged_annotations.tsv")
print(f"{len(adj)} items need adjudication -> adjudication_queue.tsv")
print("by question:", {q: int((m[f'{q}_consensus'] == 'DISAGREE').sum()) for q in ["q1", "q2", "q3"]})
print("adjudication items by source format:", adj.source_format.value_counts().to_dict())
