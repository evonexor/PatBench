"""Second attempt at detecting fragment groups, and a direct test of the alias flag.

The row-level ratio failed because a Markush table extracted wholesale makes every molecule in
the group a fragment, so the group median is a fragment size too. If that is right, fragment
groups should stand out on the group median itself, as an absolute size.
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
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
W = sys.argv[1].rstrip("/")
REL = _ws + "/outputs/affinity-sar-release-20261019/"
RAW = Path(_ws + "/raw/patent-data/")
KEY = ["patent_core", "Target", "Activity Type"]

s = pq.read_table(W + "/sample.parquet").to_pandas()
A = pd.read_csv(W + "/annotations/annotator_A.tsv", sep="\t", dtype=str).fillna("")
B = pd.read_csv(W + "/annotations/annotator_B.tsv", sep="\t", dtype=str).fillna("")
s = s.merge(A[["item_id", "q2", "q2_reason"]].rename(columns={"q2": "q2A", "q2_reason": "rA"}),
            on="item_id").merge(
    B[["item_id", "q2", "q2_reason"]].rename(columns={"q2": "q2B", "q2_reason": "rB"}), on="item_id")
s["q2c"] = np.where(s.q2A == s.q2B, s.q2A, "disagree")
txt = (s.rA + " " + s.rB).str.lower()
s["frag"] = txt.str.contains(r"markush|substituent|r-group|r group|r column|reagent|side chain|fragment")
s["alias_reason"] = txt.str.contains(r"alias|placeholder|drops the|aliased")

sar = pq.read_table(REL + "sar_all.parquet",
                    columns=KEY + ["InChIKey", "heavy_atoms", "group_id"]).to_pandas()
g = sar.groupby("group_id").agg(gmed=("heavy_atoms", "median"), gn=("InChIKey", "nunique"))
sar = sar.join(g, on="group_id")
s = s.rename(columns={"endpoint": "Activity Type"}).merge(
    sar[KEY + ["InChIKey", "heavy_atoms", "gmed", "gn"]].drop_duplicates(KEY + ["InChIKey"]),
    on=KEY + ["InChIKey"], how="left")

judg = s[s.q2c.isin(["correct", "incorrect"])].copy()
judg["frag_err"] = judg.q2c.eq("incorrect") & judg.frag
print("=== absolute size: group-median heavy atoms")
for label, m in [("fragment errors", judg.frag_err),
                 ("other incorrect", judg.q2c.eq("incorrect") & ~judg.frag),
                 ("correct", judg.q2c.eq("correct"))]:
    v = judg.loc[m, "gmed"]
    print(f"  {label:16s} n={len(v):>3}  group-median HA  min {v.min():>5.1f}  "
          f"q25 {v.quantile(.25):>5.1f}  median {v.median():>5.1f}  q75 {v.quantile(.75):>5.1f}")
print("\n  threshold scan on the group median (flag the whole group):")
best = None
for thr in [12, 14, 16, 18, 20, 22, 24]:
    flag = judg.gmed < thr
    tp = int((flag & judg.frag_err).sum())
    fp = int((flag & judg.q2c.eq("correct")).sum())
    fn = int((~flag & judg.frag_err).sum())
    p = tp / (tp + fp) if tp + fp else float("nan")
    r = tp / (tp + fn) if tp + fn else float("nan")
    print(f"    group median HA < {thr:>2}: TP {tp:>2} FP {fp:>2} FN {fn:>2}  "
          f"precision {p:.2f}  recall {r:.2f}")
    if tp + fp and (best is None or (p >= 0.9 and r > best[2])):
        best = (thr, p, r)

grp = sar.drop_duplicates("group_id")
print("\n  dataset-wide, groups by median heavy atoms:")
for thr in [14, 18, 22]:
    sm = grp[grp.gmed < thr]
    rows = sar[sar.group_id.isin(set(sm.group_id))]
    print(f"    median HA < {thr}: {len(sm):,} groups, {len(rows):,} rows "
          f"({len(rows) / len(sar):.2%} of all formal SAR)")


# ------------------------------------------------------ does an alias predict a wrong structure?
def job(args):
    archive, items = args
    out = []
    try:
        with zipfile.ZipFile(RAW / archive) as z:
            names = {n.split("/")[-1].upper(): n for n in z.namelist() if n.upper().endswith(".MOL")}
            for it in items:
                fn = (it["mol_file"] or "").upper()
                if fn not in names:
                    out.append(dict(it, has_alias=None, alias_on_carbon=None, aliases=""))
                    continue
                block = z.read(names[fn]).decode("utf-8", "replace")
                lines = block.splitlines()
                labels, on_c = [], False
                mol = Chem.MolFromMolBlock(block, sanitize=False)
                for k, line in enumerate(lines):
                    if line.startswith("A  ") and k + 1 < len(lines):
                        idx = int(line[3:].strip() or 0) - 1
                        lab = lines[k + 1].strip()
                        labels.append(lab)
                        if mol is not None and 0 <= idx < mol.GetNumAtoms():
                            sym = mol.GetAtomWithIdx(idx).GetSymbol()
                            # a label naming a heteroatom or group sitting on a carbon changes the
                            # chemistry; a label that is just "C" or a number does not
                            if sym == "C" and any(ch in lab for ch in "NOSPFIB") and lab not in ("C",):
                                on_c = True
                out.append(dict(it, has_alias=bool(labels), alias_on_carbon=on_c,
                                aliases="|".join(labels[:6])))
    except Exception:
        out.extend(dict(it, has_alias=None, alias_on_carbon=None, aliases="") for it in items)
    return out


jobs = [(a, g_[["item_id", "mol_file"]].to_dict("records"))
        for a, g_ in s.merge(
            pq.read_table(W + "/sample.parquet", columns=["item_id"]).to_pandas(), on="item_id"
        ).assign(archive=lambda d: d.Patent).groupby("Patent")]
# archive path comes from the patents table
pt = pq.read_table(_ws + "/outputs/uspto-sar-recovery-20261006/results/patents.parquet",
                   columns=["Patent", "archive"]).to_pandas()
amap = dict(zip(pt.Patent, pt.archive))
jobs = [(amap[p], items) for p, items in jobs if p in amap]
res = []
with mp.Pool(12) as pool:
    for b in pool.imap_unordered(job, jobs):
        res.extend(b)
al = pd.DataFrame(res)
s = s.merge(al[["item_id", "has_alias", "alias_on_carbon", "aliases"]], on="item_id", how="left")
judg = s[s.q2c.isin(["correct", "incorrect"])]
print("\n=== alias flags against Q2 consensus (judgeable items)")
print(pd.crosstab(judg.alias_on_carbon.fillna(False), judg.q2c, margins=True).to_string())
tp = int((judg.alias_on_carbon.fillna(False) & judg.q2c.eq("incorrect") & judg.alias_reason).sum())
fp = int((judg.alias_on_carbon.fillna(False) & judg.q2c.eq("correct")).sum())
fn = int((~judg.alias_on_carbon.fillna(False) & judg.q2c.eq("incorrect") & judg.alias_reason).sum())
print(f"\n  alias-on-carbon vs annotated alias errors: TP {tp} FP {fp} FN {fn}  "
      f"precision {tp / (tp + fp) if tp + fp else float('nan'):.2f}  "
      f"recall {tp / (tp + fn) if tp + fn else float('nan'):.2f}")
print("\n  example alias labels on items annotators called alias errors:")
for _, r in judg[judg.alias_reason & judg.q2c.eq("incorrect")].head(8).iterrows():
    print(f"    {r.item_id}  {str(r.aliases)[:70]}")
