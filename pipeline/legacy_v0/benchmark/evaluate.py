"""Scoring for the affinity-ranking benchmark.

A prediction file is a parquet with columns patent_core, Target, Activity Type, InChIKey, score.
Higher score = predicted more potent. Every metric is computed inside one ranking group and then
aggregated across groups, so a model cannot win by getting the between-assay offsets right.
"""
import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr

KEY = ["patent_core", "Target", "Activity Type"]
PAIR_MARGIN = 0.5      # log units; pairs closer than this are inside assay noise
NDCG_K = 10
TOP_FRACTION = 0.1
BOOTSTRAP = 2000
RNG = np.random.default_rng(7)


def _ndcg(truth, score, k=NDCG_K):
    rng = truth.max() - truth.min()
    rel = (truth - truth.min()) / rng if rng > 0 else np.zeros_like(truth)
    k = min(k, len(truth))
    disc = 1.0 / np.log2(np.arange(2, k + 2))
    got = rel[np.argsort(-score)[:k]]
    best = rel[np.argsort(-truth)[:k]]
    ideal = float((best * disc).sum())
    return float((got * disc).sum() / ideal) if ideal > 0 else np.nan


def _pair_accuracy(truth, score, margin=PAIR_MARGIN):
    dt = truth[:, None] - truth[None, :]
    dsc = score[:, None] - score[None, :]
    use = np.triu(np.abs(dt) >= margin, k=1)
    n = int(use.sum())
    if n == 0:
        return np.nan, 0
    agree = (np.sign(dt) == np.sign(dsc)) & use
    ties = (dsc == 0) & use
    return float((agree.sum() + 0.5 * ties.sum()) / n), n


def _top_hit(truth, score, frac=TOP_FRACTION):
    k = max(1, int(round(len(truth) * frac)))
    true_top = set(np.argsort(-truth)[:k])
    pred_top = set(np.argsort(-score)[:k])
    return len(true_top & pred_top) / k


def group_metrics(truth, score):
    """Per-group metrics. truth = pAct (higher is more potent), score = model output."""
    out = {"n_molecules": len(truth)}
    if len(truth) < 3 or np.ptp(truth) == 0:
        return None
    out["spearman"] = float(spearmanr(score, truth).statistic)
    out["kendall"] = float(kendalltau(score, truth).statistic)
    acc, npairs = _pair_accuracy(truth, score)
    out["pair_accuracy"] = acc
    out["n_pairs"] = npairs
    out["ndcg10"] = _ndcg(truth, score)
    out["top10pct_hit"] = _top_hit(truth, score)
    centred = truth - truth.mean()
    s = score - score.mean()
    denom = s.std()
    if denom > 0:
        s = s / denom * centred.std()
    out["centred_rmse"] = float(np.sqrt(np.mean((centred - s) ** 2)))
    return out


def per_group(pred, rows):
    """Join predictions to the benchmark and compute metrics for every group."""
    d = rows.merge(pred[KEY + ["InChIKey", "score"]], on=KEY + ["InChIKey"], how="left")
    missing = d.score.isna().sum()
    if missing:
        raise ValueError(f"{missing:,} benchmark rows have no prediction")
    recs = []
    for key, g in d.groupby(KEY, sort=False):
        met = group_metrics(g.pAct.to_numpy(float), g.score.to_numpy(float))
        if met is None:
            continue
        met.update(dict(zip(KEY, key)))
        for col in ["split", "test_tier", "Activity Type", "uniprot_accession", "unseen_target"]:
            if col in g:
                met[col] = g[col].iloc[0]
        met["public_share"] = float(g.in_public_data.mean()) if "in_public_data" in g else np.nan
        recs.append(met)
    return pd.DataFrame(recs)


def _ci(values):
    v = np.asarray([x for x in values if np.isfinite(x)], float)
    if len(v) < 2:
        return (np.nan, np.nan)
    idx = RNG.integers(0, len(v), size=(BOOTSTRAP, len(v)))
    meds = np.median(v[idx], axis=1)
    return (float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5)))


METRICS = ["spearman", "kendall", "pair_accuracy", "ndcg10", "top10pct_hit", "centred_rmse"]


def aggregate(gm, by=None):
    """Median (with bootstrap CI over groups) and mean for each metric."""
    rows = []
    groups = [("all", gm)] if by is None else list(gm.groupby(by, sort=False))
    for label, g in groups:
        rec = {"slice": label if not isinstance(label, tuple) else " / ".join(map(str, label)),
               "groups": len(g), "molecules": int(g.n_molecules.sum())}
        for met in METRICS:
            if met not in g:
                continue
            v = g[met].to_numpy(float)
            lo, hi = _ci(v)
            rec[f"{met}_median"] = float(np.nanmedian(v))
            rec[f"{met}_mean"] = float(np.nanmean(v))
            rec[f"{met}_ci_lo"], rec[f"{met}_ci_hi"] = lo, hi
        rows.append(rec)
    return pd.DataFrame(rows)


def report(gm, name):
    """Print the headline table for one model."""
    head = aggregate(gm)
    print(f"\n### {name}")
    # top-10% hit is reported as a mean: with a median group of ~36 molecules k is about 4,
    # so its per-group median sits at 0 even for a model that is clearly working
    cols = ["groups", "spearman_median", "spearman_ci_lo", "spearman_ci_hi",
            "pair_accuracy_median", "ndcg10_median", "top10pct_hit_mean", "centred_rmse_median"]
    print(head[[c for c in cols if c in head]].round(3).to_string(index=False))
    return head
