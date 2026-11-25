"""Size the two failure modes the annotators found, and test detectors for them.

Mode 1, Markush fragments: a table laid out as "Example | R | IC50" draws only the substituent
in the structure column, and the pipeline stored that substituent as the compound. The tell is a
molecule far smaller than the rest of its group.

Mode 2, alias atoms: a MOL file labels a group with an MDL alias ("HN", "Boc"), RDKit reads the
aliased atom as a placeholder, and the stored structure loses atoms. This one is invisible to the
InChIKey check, because the stored and recomputed keys come from the same faulty parse.

Each detector is scored against the annotators' consensus before it is applied dataset-wide, so
the flag that ships has a measured precision and recall rather than an assumed one.
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import json
import multiprocessing as mp
import sys
import warnings
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")
pd.set_option("display.width", 220)
pd.set_option("display.max_colwidth", 80)
W = sys.argv[1].rstrip("/")
REL = _ws + "/outputs/affinity-sar-release-20261019/"
RAW = Path(_ws + "/raw/patent-data/")
KEY = ["patent_core", "Target", "Activity Type"]
FRAG_RATIO = 0.5

# ------------------------------------------------------------ annotations
s = pq.read_table(W + "/sample.parquet").to_pandas()
A = pd.read_csv(W + "/annotations/annotator_A.tsv", sep="\t", dtype=str).fillna("")
B = pd.read_csv(W + "/annotations/annotator_B.tsv", sep="\t", dtype=str).fillna("")
s = s.merge(A[["item_id", "q2", "q2_reason"]].rename(columns={"q2": "q2A", "q2_reason": "rA"}),
            on="item_id").merge(
    B[["item_id", "q2", "q2_reason"]].rename(columns={"q2": "q2B", "q2_reason": "rB"}), on="item_id")
s["q2_consensus"] = np.where(s.q2A == s.q2B, s.q2A, "disagree")
reasons = (s.rA + " " + s.rB).str.lower()
s["said_fragment"] = reasons.str.contains(
    r"markush|substituent|r-group|r group|r column|r1 |r3 |reagent|side chain|fragment")
s["said_alias"] = reasons.str.contains(r"alias|placeholder|drops the|aliased")

print("=== Q2 consensus by stratum (judgeable items only)")
j = s[s.q2_consensus.isin(["correct", "incorrect"])]
t = j.groupby("stratum").q2_consensus.agg(
    n="size", correct=lambda x: int((x == "correct").sum()))
t["accuracy"] = (t.correct / t.n).round(3)
print(t.sort_values("accuracy").to_string())
print("\nby source format:")
f = j.groupby("source_format").q2_consensus.agg(
    n="size", correct=lambda x: int((x == "correct").sum()))
f["accuracy"] = (f.correct / f.n).round(3)
print(f.to_string())

print(f"\nitems whose reasons cite a Markush / substituent / reagent fragment: "
      f"{int(s.said_fragment.sum())}")
print(f"items whose reasons cite an alias / placeholder atom: {int(s.said_alias.sum())}")

# ------------------------------------------------------ mode 1 detector
sar = pq.read_table(REL + "sar_all.parquet",
                    columns=KEY + ["InChIKey", "heavy_atoms", "group_id", "source",
                                   "is_exact", "is_affinity"]).to_pandas()
gmed = sar.groupby("group_id").heavy_atoms.median().rename("group_median_heavy")
sar = sar.join(gmed, on="group_id")
sar["frag_ratio"] = sar.heavy_atoms / sar.group_median_heavy
sar["suspect_fragment"] = sar.frag_ratio < FRAG_RATIO

s2 = s.rename(columns={"endpoint": "Activity Type"}).merge(
    sar[KEY + ["InChIKey", "heavy_atoms", "group_median_heavy", "frag_ratio",
               "suspect_fragment"]].drop_duplicates(KEY + ["InChIKey"]),
    on=KEY + ["InChIKey"], how="left")
truth = s2[s2.q2_consensus.isin(["correct", "incorrect"])].copy()
truth["is_fragment_error"] = truth.q2_consensus.eq("incorrect") & truth.said_fragment
tp = int((truth.suspect_fragment & truth.is_fragment_error).sum())
fp = int((truth.suspect_fragment & truth.q2_consensus.eq("correct")).sum())
fn = int((~truth.suspect_fragment.fillna(False) & truth.is_fragment_error).sum())
prec = tp / (tp + fp) if tp + fp else float("nan")
rec = tp / (tp + fn) if tp + fn else float("nan")
print(f"\n=== mode 1 detector: heavy atoms < {FRAG_RATIO} x group median")
print(f"  against annotated fragment errors: TP {tp}  FP {fp}  FN {fn}  "
      f"precision {prec:.2f}  recall {rec:.2f}")
print("  heavy-atom ratio of annotated fragment errors:",
      truth.loc[truth.is_fragment_error, "frag_ratio"].round(2).describe()[
          ["min", "25%", "50%", "75%", "max"]].to_dict())
print("  heavy-atom ratio of annotated correct items:  ",
      truth.loc[truth.q2_consensus.eq("correct"), "frag_ratio"].round(2).describe()[
          ["min", "25%", "50%", "75%", "max"]].to_dict())
for thr in [0.3, 0.4, 0.5, 0.6, 0.7]:
    flag = truth.frag_ratio < thr
    tp_ = int((flag & truth.is_fragment_error).sum())
    fp_ = int((flag & truth.q2_consensus.eq("correct")).sum())
    fn_ = int((~flag & truth.is_fragment_error).sum())
    print(f"    ratio < {thr}: TP {tp_:>2} FP {fp_:>2} FN {fn_:>2}  "
          f"precision {tp_ / (tp_ + fp_) if tp_ + fp_ else float('nan'):.2f}  "
          f"recall {tp_ / (tp_ + fn_) if tp_ + fn_ else float('nan'):.2f}")

print("\n  dataset-wide prevalence of the flag:")
grp = sar[sar.is_affinity & sar.is_exact]
print(f"    affinity exact rows flagged: {int(grp.suspect_fragment.sum()):,} of {len(grp):,} "
      f"= {grp.suspect_fragment.mean():.2%}")
print("    by extraction pass:", grp.groupby("source").suspect_fragment.mean().round(4).to_dict())


# ------------------------------------------------------ mode 2 prevalence
def alias_job(args):
    archive, items = args
    out = []
    try:
        with zipfile.ZipFile(RAW / archive) as z:
            names = {n.split("/")[-1].upper(): n for n in z.namelist() if n.upper().endswith(".MOL")}
            for it in items:
                fn = it["mol"]
                if fn not in names:
                    out.append(dict(it, has_alias=None, has_dummy=None))
                    continue
                block = z.read(names[fn]).decode("utf-8", "replace")
                mol = Chem.MolFromMolBlock(block, sanitize=False)
                if mol is None:
                    out.append(dict(it, has_alias=None, has_dummy=None))
                    continue
                alias = any(a.HasProp("molFileAlias") for a in mol.GetAtoms())
                dummy = any(a.GetAtomicNum() == 0 for a in mol.GetAtoms())
                # a V2000 alias block is "A  nnn" followed by the label on the next line
                alias = alias or any(l.startswith("A  ") for l in block.splitlines())
                out.append(dict(it, has_alias=alias, has_dummy=dummy))
    except Exception:
        out.extend(dict(it, has_alias=None, has_dummy=None) for it in items)
    return out


sc = pq.read_table(W + "/structure_check_records.parquet").to_pandas()
sc = sc[sc.verdict.eq("exact") & sc.detail.fillna("").ne("")]
sc["mol"] = sc.detail.str.upper()
jobs = [(a, g[["source_record_id", "stratum", "source_format", "pass_name", "mol"]].to_dict("records"))
        for a, g in sc.groupby("archive", sort=False)]
res = []
with mp.Pool(18) as pool:
    for b in pool.imap_unordered(alias_job, jobs, chunksize=4):
        res.extend(b)
al = pd.DataFrame(res)
al = al[al.has_alias.notna()]
print(f"\n=== mode 2: alias atoms in the matched MOL file ({len(al):,} of the 20k sample)")
print(f"  MOL carries an alias label: {int(al.has_alias.sum()):,} = {al.has_alias.mean():.2%}")
print(f"  MOL carries a dummy atom:   {int(al.has_dummy.sum()):,} = {al.has_dummy.mean():.2%}")
print("  alias rate by stratum:")
print(al.groupby("stratum").has_alias.agg(["size", "mean"]).round(4).to_string())

out = {
    "q2_by_source_format": {k: {"n": int(v.n), "accuracy": float(v.accuracy)}
                            for k, v in f.iterrows()},
    "fragment_detector": {"rule": f"heavy_atoms < {FRAG_RATIO} x group median heavy_atoms",
                          "tp": tp, "fp": fp, "fn": fn,
                          "precision": round(prec, 3), "recall": round(rec, 3),
                          "dataset_rows_flagged": int(grp.suspect_fragment.sum()),
                          "dataset_share_flagged": round(float(grp.suspect_fragment.mean()), 5)},
    "alias_atoms": {"records_checked": len(al),
                    "with_alias": int(al.has_alias.sum()),
                    "alias_rate": round(float(al.has_alias.mean()), 5),
                    "with_dummy_atom": int(al.has_dummy.sum()),
                    "by_stratum": {k: round(float(v), 5)
                                   for k, v in al.groupby("stratum").has_alias.mean().items()}},
}
with open(W + "/failure_modes.json", "w") as fh:
    json.dump(out, fh, indent=1, default=float)
print(f"\nwrote {W}/failure_modes.json")
