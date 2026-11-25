"""Score a model on a PatBench ranking task.

    python evaluate.py tasks/rank_ic50.parquet predictions.parquet [--split split_time]
                       [--part test] [--tier novel] [--collapse-stereo] [--out per_group.csv]

A prediction file (parquet or csv) has columns group_id, InChIKey, score; higher score means
predicted more potent. Every metric is computed inside one group (one patent, one assay) and then
summarised across groups, so a model cannot score by getting the offsets between assays right.

--collapse-stereo merges each stereo set (molecules of a group that differ only in stereo or
isotopes, column stereo_set_id) into one entry, taking the median of the truth and of the score.
Report both views: a model reading 2D graphs gives every member of a set the same score, and the
uncollapsed view charges it for orderings it has no input to make.

--cliffs protocols/<task>_cliffs.parquet reports cliff accuracy over the frozen pairs of
near-identical molecules (ECFP4 Tanimoto >= 0.85) that differ by >= 1 log in potency: the share
the model orders correctly, ties counting half. Chance is 0.500, and a model that has learned only
the scaffold's average potency lands there while its Spearman still looks healthy.

--fewshot protocols/<task>_<split>_k<K>.parquet scores the few-shot protocol: in each repeat the
"support" molecules of a group are revealed to the model with their values and only the "query"
molecules are scored. The prediction file then needs a `repeat` column alongside group_id,
InChIKey and score; support rows are dropped before scoring, and the metrics are averaged over
repeats per group before the usual across-group summary.
"""
import argparse

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr

PAIR_MARGIN = 0.5      # log units; pairs closer than this are inside assay noise
NDCG_K = 10
TOP_FRACTION = 0.1
BOOTSTRAP = 2000
RNG = np.random.default_rng(7)
METRICS = ["spearman", "kendall", "pair_accuracy", "ndcg10", "top10pct_hit"]


def _ndcg(truth, score, k=NDCG_K):
    span = truth.max() - truth.min()
    rel = (truth - truth.min()) / span if span > 0 else np.zeros_like(truth)
    k = min(k, len(truth))
    disc = 1.0 / np.log2(np.arange(2, k + 2))
    ideal = float((rel[np.argsort(-truth)[:k]] * disc).sum())
    return float((rel[np.argsort(-score)[:k]] * disc).sum() / ideal) if ideal > 0 else np.nan


def _pair_accuracy(truth, score, margin=PAIR_MARGIN):
    dt = truth[:, None] - truth[None, :]
    ds = score[:, None] - score[None, :]
    use = np.triu(np.abs(dt) >= margin, k=1)
    n = int(use.sum())
    if n == 0:
        return np.nan, 0
    agree = ((np.sign(dt) == np.sign(ds)) & use).sum()
    ties = ((ds == 0) & use).sum()
    return float((agree + 0.5 * ties) / n), n


def _top_hit(truth, score, frac=TOP_FRACTION):
    k = max(1, int(round(len(truth) * frac)))
    return len(set(np.argsort(-truth)[:k]) & set(np.argsort(-score)[:k])) / k


def group_metrics(truth, score):
    if len(truth) < 3 or np.ptp(truth) == 0:
        return None
    acc, npairs = _pair_accuracy(truth, score)
    return {"n_molecules": len(truth),
            "spearman": float(spearmanr(score, truth).statistic) if np.ptp(score) > 0 else 0.0,
            "kendall": float(kendalltau(score, truth).statistic) if np.ptp(score) > 0 else 0.0,
            "pair_accuracy": acc, "n_pairs": npairs,
            "ndcg10": _ndcg(truth, score), "top10pct_hit": _top_hit(truth, score)}


def collapse_stereo(d):
    unit = np.where(d.stereo_set_id.ne(""), d.stereo_set_id, d.group_id.astype(str) + ":" + d.InChIKey)
    keep = [c for c in d.columns if c not in ("pAct", "score", "InChIKey")]
    return (d.assign(_u=unit).groupby("_u", sort=False)
            .agg({**{c: "first" for c in keep}, "pAct": "median", "score": "median"}).reset_index(drop=True))


def per_group(task, pred, collapse=False):
    d = task.merge(pred[["group_id", "InChIKey", "score"]], on=["group_id", "InChIKey"], how="left")
    missing = int(d.score.isna().sum())
    if missing:
        raise SystemExit(f"{missing:,} task rows have no prediction")
    if collapse:
        d = collapse_stereo(d)
    recs = []
    for gid, g in d.groupby("group_id", sort=False):
        m = group_metrics(g.pAct.to_numpy(float), g.score.to_numpy(float))
        if m is None:
            continue
        m.update(group_id=gid, patent_core=g.patent_core.iloc[0], uniprot_accession=g.uniprot_accession.iloc[0],
                 endpoint=g["Activity Type"].iloc[0])
        recs.append(m)
    return pd.DataFrame(recs)


def per_group_fewshot(task, pred, protocol, collapse=False):
    """One row per (group, repeat): only the query molecules of that repeat are scored."""
    if "repeat" not in pred.columns:
        raise SystemExit("few-shot predictions need a `repeat` column (one block per repeat)")
    p = protocol[protocol.role.eq("query")].merge(
        pred[["group_id", "repeat", "InChIKey", "score"]], on=["group_id", "repeat", "InChIKey"],
        how="left")
    missing = int(p.score.isna().sum())
    if missing:
        raise SystemExit(f"{missing:,} query rows have no prediction")
    d = task.merge(p, on=["group_id", "InChIKey"], how="inner")
    recs = []
    for (gid, rep), g in d.groupby(["group_id", "repeat"], sort=False):
        gg = collapse_stereo(g) if collapse else g
        m = group_metrics(gg.pAct.to_numpy(float), gg.score.to_numpy(float))
        if m is None:
            continue
        m.update(group_id=gid, repeat=rep, patent_core=g.patent_core.iloc[0],
                 uniprot_accession=g.uniprot_accession.iloc[0], endpoint=g["Activity Type"].iloc[0])
        recs.append(m)
    gm = pd.DataFrame(recs)
    # average the repeats of one group first, so a group counts once however many repeats it has
    agg = {c: "mean" for c in METRICS + ["n_molecules", "n_pairs"] if c in gm}
    first = {c: "first" for c in ["patent_core", "uniprot_accession", "endpoint"]}
    return gm.groupby("group_id", as_index=False).agg({**agg, **first})


def cliff_accuracy(task, pred, cliffs, part_groups):
    """Share of cliff pairs ordered correctly (ties half), over the pairs inside `part_groups`."""
    score = pred.set_index(["group_id", "InChIKey"]).score
    truth = task.set_index(["group_id", "InChIKey"]).pAct
    c = cliffs[cliffs.group_id.isin(part_groups)]
    keys_a = list(zip(c.group_id, c.InChIKey_a))
    keys_b = list(zip(c.group_id, c.InChIKey_b))
    use = [k in score.index and j in score.index for k, j in zip(keys_a, keys_b)]
    c, keys_a, keys_b = c[use], [k for k, u in zip(keys_a, use) if u], [k for k, u in zip(keys_b, use) if u]
    if not len(c):
        return None
    ds = score.loc[keys_a].to_numpy() - score.loc[keys_b].to_numpy()
    dt = truth.loc[keys_a].to_numpy() - truth.loc[keys_b].to_numpy()
    correct = float(((np.sign(ds) == np.sign(dt)).sum() + 0.5 * (ds == 0).sum()) / len(c))
    z = 1.959963985
    lo, hi = _wilson(correct * len(c), len(c), z)
    return {"pairs": len(c), "groups": int(c.group_id.nunique()), "accuracy": correct,
            "ci_lo": lo, "ci_hi": hi}


def _wilson(k, n, z=1.959963985):
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return float(centre - half), float(centre + half)


def _ci(v):
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    if len(v) < 2:
        return np.nan, np.nan
    meds = np.median(v[RNG.integers(0, len(v), size=(BOOTSTRAP, len(v)))], axis=1)
    return float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))


def summary(gm):
    rec = {"groups": len(gm), "molecules": int(gm.n_molecules.sum())}
    for m in METRICS:
        v = gm[m].to_numpy(float)
        rec[f"{m}_median"] = float(np.nanmedian(v))
        rec[f"{m}_mean"] = float(np.nanmean(v))
        rec[f"{m}_ci_lo"], rec[f"{m}_ci_hi"] = _ci(v)
    return rec


def read(path):
    return pd.read_parquet(path) if str(path).endswith(".parquet") else pd.read_csv(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("task")
    ap.add_argument("predictions")
    ap.add_argument("--split", default="split_time")
    ap.add_argument("--part", default="test")
    ap.add_argument("--tier", default=None, help="novel or reference (time split only)")
    ap.add_argument("--collapse-stereo", action="store_true")
    ap.add_argument("--cliffs", default=None, help="protocols/<task>_cliffs.parquet")
    ap.add_argument("--fewshot", default=None, help="protocols/<task>_<split>_k<K>.parquet")
    ap.add_argument("--out", default=None, help="write per-group metrics here")
    a = ap.parse_args()
    task = read(a.task)
    task = task[task[a.split] == a.part]
    if a.tier:
        task = task[task.test_tier == a.tier]
    pred = read(a.predictions)
    if a.fewshot:
        gm = per_group_fewshot(task, pred, read(a.fewshot), a.collapse_stereo)
    else:
        gm = per_group(task, pred, a.collapse_stereo)
    if a.out:
        gm.to_csv(a.out, index=False)
    s = summary(gm)
    print(f"{a.task} | {a.split}={a.part}{' tier=' + a.tier if a.tier else ''}"
          f"{' | stereo sets collapsed' if a.collapse_stereo else ''}"
          f"{' | few-shot ' + a.fewshot.split('/')[-1] if a.fewshot else ''}")
    print(f"groups {s['groups']:,}  molecules {s['molecules']:,}")
    for m in METRICS:
        print(f"  {m:14s} median {s[m + '_median']:.3f} [{s[m + '_ci_lo']:.3f}, {s[m + '_ci_hi']:.3f}]"
              f"  mean {s[m + '_mean']:.3f}")
    if a.cliffs:
        # cliffs are scored zero-shot: a few-shot prediction file holds one score per repeat
        p = pred[pred.repeat.eq(pred.repeat.min())] if "repeat" in pred.columns else pred
        c = cliff_accuracy(task, p, read(a.cliffs), set(task.group_id))
        if c is None:
            print("  cliff_accuracy  no pairs in this slice")
        else:
            print(f"  cliff_accuracy  {c['accuracy']:.3f} [{c['ci_lo']:.3f}, {c['ci_hi']:.3f}]"
                  f"  over {c['pairs']:,} pairs in {c['groups']} groups (chance 0.500)")


if __name__ == "__main__":
    main()
