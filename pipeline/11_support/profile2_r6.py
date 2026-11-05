"""Second part of the r6 profile: molecular weight against activity, over time, extremes, and the histogram for the figure.

    python profile2_r6.py <release dir> <out dir>
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from scipy.stats import spearmanr

RDLogger.DisableLog("rdApp.*")
D = Path(sys.argv[1]).resolve()
OUT = Path(sys.argv[2]).resolve()
R = {}
cache = {}


def mw(s):
    if s not in cache:
        m = Chem.MolFromSmiles(s)
        cache[s] = Descriptors.MolWt(max(Chem.GetMolFrags(m, asMols=True), key=lambda f: f.GetNumHeavyAtoms())) if m else np.nan
    return cache[s]


for t in ("rank_ic50", "rank_ki_kd", "rank_highcontrast", "rank_complex_mutant", "rank_ec50"):
    df = pd.read_parquet(D / "tasks" / f"{t}.parquet")
    df["mw"] = [mw(s) for s in df.SMILES]
    out = {}
    out["pooled_spearman_mw_pAct"] = round(float(spearmanr(df.mw, df.pAct)[0]), 3)
    rho = []
    for gid, g in df.groupby("group_id"):
        if g.mw.nunique() > 2 and g.pAct.nunique() > 2:
            rho.append((gid, spearmanr(g.mw, g.pAct)[0], g.split_time.iloc[0], g.test_tier.iloc[0], len(g)))
    r = pd.DataFrame(rho, columns=["gid", "rho", "split", "tier", "n"])
    out["within_group_rho_median"] = round(float(r.rho.median()), 3)
    out["within_group_rho_p10_p90"] = [round(float(r.rho.quantile(.1)), 3), round(float(r.rho.quantile(.9)), 3)]
    out["share_groups_rho_gt_0.3"] = round(float((r.rho > 0.3).mean()), 3)
    out["share_groups_rho_lt_-0.3"] = round(float((r.rho < -0.3).mean()), 3)
    for sp in ("train", "valid", "test"):
        x = r[r.split == sp].rho
        out[f"rho_median_{sp}"] = round(float(x.median()), 3) if len(x) else None
    for ti in ("novel", "reference"):
        x = r[(r.split == "test") & (r.tier == ti)].rho
        out[f"rho_median_test_{ti}"] = round(float(x.median()), 3) if len(x) else None
    out["pAct_gt_10.5_rows"] = int((df.pAct > 10.5).sum())
    out["pAct_lt_3.5_rows"] = int((df.pAct < 3.5).sum())
    out["mw_gt_800_rows"] = int((df.mw > 800).sum())
    out["mw_gt_800_molecules"] = int(df[df.mw > 800].InChIKey.nunique())
    out["pAct_mean_by_split"] = df.groupby("split_time").pAct.mean().round(2).to_dict()
    out["mw_median_by_year_groups"] = None
    if t == "rank_ic50":
        out["mw_median_by_year"] = df.groupby("year").mw.median().round(0).to_dict()
        out["rows_by_year"] = df.groupby("year").size().to_dict()
        out["test_years"] = df[df.split_time == "test"].year.value_counts().sort_index().to_dict()
        out["train_years"] = df[df.split_time == "train"].year.value_counts().sort_index().to_dict()
        out["valid_years"] = df[df.split_time == "valid"].year.value_counts().sort_index().to_dict()
        big = df[df.mw > 800].drop_duplicates("InChIKey").head(6)
        out["examples_mw_gt_800"] = [(round(m, 0), s[:90]) for m, s in zip(big.mw, big.SMILES)]
        tiny = df[df.mw < 200].drop_duplicates("InChIKey").head(6)
        out["examples_mw_lt_200"] = [(round(m, 0), s[:90]) for m, s in zip(tiny.mw, tiny.SMILES)]
    R[t] = out
    print(t, "done", flush=True)
json.dump(R, open(OUT / "profile2.json", "w"), indent=1, ensure_ascii=False, default=str)
h = pd.read_csv(OUT / "mw_hist.tsv", sep="\t")
print(h.groupby("set")[["count_all", "count_test"]].sum())
