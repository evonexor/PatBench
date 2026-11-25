"""Metric definitions of pipeline/10_evaluation/evaluate.py, on small synthetic score vectors."""
import importlib.util

import numpy as np
import pandas as pd
import pytest

from patbench import PIPELINE

spec = importlib.util.spec_from_file_location("pb_evaluate", PIPELINE / "10_evaluation/evaluate.py")
E = importlib.util.module_from_spec(spec)
spec.loader.exec_module(E)

TRUTH = np.array([5.0, 6.0, 7.0, 8.0, 9.0, 4.0, 5.5, 6.5, 7.5, 8.5, 9.5, 3.0])


def test_perfect_ranking():
    m = E.group_metrics(TRUTH, TRUTH * 2 + 1)
    assert m["spearman"] == pytest.approx(1.0) and m["kendall"] == pytest.approx(1.0)
    assert m["pair_accuracy"] == pytest.approx(1.0) and m["ndcg10"] == pytest.approx(1.0)
    assert m["top10pct_hit"] == 1.0


def test_reversed_ranking():
    m = E.group_metrics(TRUTH, -TRUTH)
    assert m["spearman"] == pytest.approx(-1.0) and m["pair_accuracy"] == pytest.approx(0.0)


def test_constant_scores_score_zero_correlation_and_half_pair_accuracy():
    m = E.group_metrics(TRUTH, np.zeros_like(TRUTH))
    assert m["spearman"] == 0.0 and m["kendall"] == 0.0
    assert m["pair_accuracy"] == pytest.approx(0.5)


def test_groups_that_cannot_be_ranked_are_skipped():
    assert E.group_metrics(np.array([1.0, 2.0]), np.array([1.0, 2.0])) is None     # fewer than 3 molecules
    assert E.group_metrics(np.array([5.0, 5.0, 5.0]), np.array([1.0, 2.0, 3.0])) is None   # no spread


def test_pair_accuracy_ignores_pairs_inside_the_noise_margin():
    truth = np.array([5.0, 5.2, 7.0])
    acc, n = E._pair_accuracy(truth, np.array([2.0, 1.0, 3.0]))
    assert n == 2 and acc == pytest.approx(1.0)       # (5.0, 5.2) is closer than 0.5 log and not counted


def test_pair_accuracy_counts_ties_half():
    acc, n = E._pair_accuracy(np.array([5.0, 6.0]), np.array([1.0, 1.0]))
    assert n == 1 and acc == pytest.approx(0.5)


def test_top_hit_share_of_the_true_top_fraction():
    truth = np.arange(20.0)
    score = truth.copy()
    score[[18, 19]] = [0.0, 1.0]                      # the two truly best molecules are ranked last
    assert E._top_hit(truth, truth) == 1.0
    assert E._top_hit(truth, score) == 0.0


def test_ndcg_of_the_ideal_order_is_one_and_of_the_worst_is_lower():
    assert E._ndcg(TRUTH, TRUTH) == pytest.approx(1.0)
    assert E._ndcg(TRUTH, -TRUTH) < 0.7


def test_wilson_interval_contains_the_estimate():
    lo, hi = E._wilson(60, 100)
    assert lo < 0.6 < hi and 0 < lo and hi < 1


def test_bootstrap_interval_brackets_the_median():
    lo, hi = E._ci(np.linspace(-0.1, 0.4, 50))
    assert lo < np.median(np.linspace(-0.1, 0.4, 50)) < hi


def _task(n_groups=3, n=12, seed=0):
    r = np.random.default_rng(seed)
    rows = []
    for g in range(n_groups):
        for i in range(n):
            rows.append(dict(group_id=g, InChIKey=f"K{g}_{i}", pAct=float(r.normal(7, 1)), patent_core=f"US{g}",
                             uniprot_accession="P00001", stereo_set_id="", **{"Activity Type": "IC50"}))
    return pd.DataFrame(rows)


def test_per_group_end_to_end_with_a_perfect_model():
    task = _task()
    pred = task[["group_id", "InChIKey"]].assign(score=task.pAct)
    gm = E.per_group(task, pred)
    assert len(gm) == 3 and (gm.spearman > 0.999).all()
    assert E.summary(gm)["spearman_median"] == pytest.approx(1.0)


def test_missing_predictions_abort():
    task = _task()
    pred = task[["group_id", "InChIKey"]].assign(score=task.pAct).iloc[:-1]
    with pytest.raises(SystemExit):
        E.per_group(task, pred)


def test_collapse_stereo_merges_a_stereo_set_into_one_entry():
    task = _task(1, 6)
    task["stereo_set_id"] = ["s1", "s1", "", "", "", ""]
    pred = task[["group_id", "InChIKey"]].assign(score=task.pAct)
    d = task.merge(pred, on=["group_id", "InChIKey"])
    assert len(E.collapse_stereo(d)) == 5


def test_cliff_accuracy_counts_correctly_ordered_pairs():
    task = pd.DataFrame(dict(group_id=[1, 1, 1, 1], InChIKey=list("abcd"), pAct=[5.0, 7.0, 5.0, 8.0]))
    pred = task[["group_id", "InChIKey"]].assign(score=[1.0, 2.0, 2.0, 1.0])     # first pair right, second wrong
    cliffs = pd.DataFrame(dict(group_id=[1, 1], InChIKey_a=["b", "d"], InChIKey_b=["a", "c"]))
    r = E.cliff_accuracy(task, pred, cliffs, {1})
    assert r["pairs"] == 2 and r["accuracy"] == pytest.approx(0.5)
