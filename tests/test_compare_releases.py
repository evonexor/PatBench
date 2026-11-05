"""tools/compare_releases.py: the comparison behind the reproduction claim must catch what it claims to catch.

Frames are small fixtures built here, not PatBench data."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tools"))
import compare_releases as C  # noqa: E402


def frame():
    return pd.DataFrame({"id": [1, 2, 3], "x": [0.5, np.nan, 2.0], "s": ["a", None, "c"], "b": [True, False, True]})


def test_equal_frames():
    assert C.compare_frames(frame(), frame()) == []


def test_missing_equals_missing_but_not_a_value():
    b = frame()
    b.loc[1, "x"] = 1.0
    assert any("x:" in w for w in C.compare_frames(frame(), b))


def test_float_noise_below_tolerance_passes_and_above_fails():
    b = frame()
    b.loc[0, "x"] = 0.5 * (1 + 1e-12)
    assert C.compare_frames(frame(), b) == []
    b.loc[0, "x"] = 0.5 * (1 + 1e-6)
    assert any("x:" in w for w in C.compare_frames(frame(), b))


def test_changed_string_row_count_dtype_and_columns_are_caught():
    b = frame()
    b.loc[2, "s"] = "z"
    assert any("s:" in w for w in C.compare_frames(frame(), b))
    assert any("rows" in w for w in C.compare_frames(frame(), frame().iloc[:2]))
    assert any("dtype" in w for w in C.compare_frames(frame(), frame().astype({"id": "float64"})))
    assert any("columns differ" in w for w in C.compare_frames(frame(), frame().drop(columns="b")))
    assert any("order" in w for w in C.compare_frames(frame(), frame()[["x", "id", "s", "b"]]))


def test_row_order_matters():
    assert C.compare_frames(frame(), frame().iloc[::-1].reset_index(drop=True))


def test_array_cells():
    a = pd.DataFrame({"v": [np.array([1, 2]), np.array([3])]})
    assert C.compare_frames(a, a.copy()) == []
    b = pd.DataFrame({"v": [np.array([1, 2]), np.array([4])]})
    assert C.compare_frames(a, b)


def test_manifest_comparison_ignores_build_time_fields_only():
    a = {"config": {"A": 1}, "tasks": {"t": 1}, "built": "x", "inputs": {"p": "/a"}, "log": ["1"]}
    b = {"config": {"A": 1}, "tasks": {"t": 1}, "built": "y", "inputs": {"p": "/b"}, "log": ["2"]}
    assert C.compare_manifests(a, b) == []
    b["config"]["WS"] = "/w"
    assert C.compare_manifests(a, b) == ["config"]


def test_directory_comparison_end_to_end(tmp_path, capsys):
    ref, new = tmp_path / "ref", tmp_path / "new"
    for d in (ref, new):
        (d / "tasks").mkdir(parents=True)
        frame().to_parquet(d / "tasks/t.parquet", index=False)
        (d / "release_manifest.json").write_text(json.dumps({"config": {"A": 1}}))
    assert C.main([str(ref), str(new)]) == 0
    other = frame()
    other.loc[0, "id"] = 9
    other.to_parquet(new / "tasks/t.parquet", index=False)
    assert C.main([str(ref), str(new)]) == 1
    out = capsys.readouterr().out
    assert "DIFF tasks/t.parquet" in out and "1 different" in out
