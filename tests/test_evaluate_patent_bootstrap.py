"""The patent-cluster evaluator named in the r6 freeze record, on small synthetic data (not PatBench data)."""
import importlib.util
import hashlib

import numpy as np
import pandas as pd
import pytest

from patbench import PIPELINE


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, PIPELINE / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


B = _load("pb_evaluate_patent", "10_evaluation/evaluate_patent_bootstrap.py")
G = _load("pb_evaluate_groups", "10_evaluation/evaluate.py")

FROZEN_SHA256 = "35b70dd019dc8da9fdb43a587dc8182e6d7f4c6c30a9268b1baefb1a13fc4f09"      # FROZEN.json of r6


def test_file_is_the_evaluator_named_in_the_freeze_record():
    assert hashlib.sha256((PIPELINE / "10_evaluation/evaluate_patent_bootstrap.py").read_bytes()).hexdigest() == FROZEN_SHA256


def test_group_metrics_agree_with_the_build_time_evaluator():
    rng = np.random.default_rng(3)
    for _ in range(5):
        truth = rng.normal(7, 1, 30)
        score = truth + rng.normal(0, 1, 30)
        a, b = B.group_metrics(truth, score), G.group_metrics(truth, score)
        for k in ("spearman", "kendall", "pair_accuracy", "ndcg10", "top10pct_hit"):
            assert a[k] == pytest.approx(b[k]), k


def test_summary_resamples_patents():
    gm = pd.DataFrame({"group_id": range(8), "patent_core": [f"US{i // 2}" for i in range(8)],
                       "n_molecules": 20, "n_pairs": 100,
                       **{m: np.linspace(0.0, 0.7, 8) for m in B.METRICS}})
    s = B.summary(gm)
    assert s["ci_unit"] == "patent" and s["patents"] == 4 and s["groups"] == 8
    assert s["spearman_ci_lo"] <= s["spearman_median"] <= s["spearman_ci_hi"]


def test_group_only_intervals_are_refused():
    with pytest.raises(ValueError):
        B._ci(np.array([0.1, 0.2, 0.3]))


def test_duplicate_prediction_keys_are_refused():
    pred = pd.DataFrame({"group_id": [1, 1], "InChIKey": ["A", "A"], "score": [1.0, 2.0]})
    with pytest.raises(ValueError):
        B._check_pred(pred, ["group_id", "InChIKey"])
