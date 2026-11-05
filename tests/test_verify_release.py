"""verify_release on small synthetic release directories.

The tables below are test fixtures built in memory for this file. They are not PatBench data and carry
no chemistry: identifiers are placeholders chosen to exercise each rule."""
import copy
import json

import pandas as pd
import pytest

from patbench import CONFIGS
from patbench.verify import REQUIRED, TASKS, verify_release

CFG = json.loads((CONFIGS / "release_v1.0.json").read_text())


def _rows():
    r = []
    # five patents; two in test (priority 2018-2019, granted 2021-2022: above every task's cutoff and the grant floor)
    plan = [("US1", "train", 2012, 2015), ("US2", "train", 2013, 2016), ("US3", "valid", 2014, 2017),
            ("US4", "test", 2018, 2021), ("US5", "test", 2019, 2022)]
    gid = 0
    for pat, part, prio, year in plan:
        gid += 1
        for i in range(4):
            r.append(dict(group_id=gid, patent_core=pat, family_id=f"F{pat}", InChIKey=f"{pat}-MOL{i}",
                          uniprot_accession=f"P{gid:05d}", year=year, priority_date=f"{prio}-03-01",
                          split_time=part, split_time_grant=part, split_random=part, split_scaffold_group=part))
    return pd.DataFrame(r)


def _release(tmp_path, mutate=None):
    root = tmp_path / "rel"
    (root / "tasks").mkdir(parents=True)
    base = _rows()
    tasks = {}
    for t in TASKS:
        d = base.copy()
        if mutate:
            d = mutate(t, d)
        d.to_parquet(root / f"tasks/{t}.parquet", index=False)
        tasks[t] = d
    for f in REQUIRED:
        p = root / f
        if not p.exists() and f.endswith(".parquet"):
            pd.DataFrame({"x": [1]}).to_parquet(p, index=False)
    stats = {t: dict(rows=len(d), groups=int(d.group_id.nunique()), patents=int(d.patent_core.nunique()),
                     molecules=int(d.InChIKey.nunique())) for t, d in tasks.items()}
    (root / "release_manifest.json").write_text(json.dumps({"config": {**CFG, "VAL": "x", "BUILT": "y"}, "tasks": stats}))
    return root


def _failed(checks):
    return {c.name for c in checks if not c.ok}


def test_a_consistent_release_passes(tmp_path):
    checks = verify_release(_release(tmp_path))
    assert _failed(checks) == set(), [c for c in checks if not c.ok]
    assert any(c.name.startswith("no_pair_leak") for c in checks)


def test_a_patent_in_train_and_test_is_caught(tmp_path):
    def m(t, d):
        d = d.copy()
        d.loc[d.index[-1], "patent_core"] = "US1"           # a test row now carries a train patent
        return d
    assert any(n.startswith("no_patent_leak[rank_ic50:split_time]") for n in _failed(verify_release(_release(tmp_path, m))))


def test_a_family_across_parts_is_caught(tmp_path):
    def m(t, d):
        d = d.copy()
        d.loc[d.split_time == "test", "family_id"] = "FUS1"
        return d
    assert any(n.startswith("no_family_leak") for n in _failed(verify_release(_release(tmp_path, m))))


def test_a_shared_molecule_target_pair_is_caught(tmp_path):
    def m(t, d):
        d = d.copy()
        i = d.index[d.split_time == "test"][0]
        j = d.index[d.split_time == "train"][0]
        d.loc[i, ["InChIKey", "uniprot_accession"]] = d.loc[j, ["InChIKey", "uniprot_accession"]].tolist()
        return d
    assert any(n.startswith("no_pair_leak") for n in _failed(verify_release(_release(tmp_path, m))))


def test_a_test_row_older_than_the_cutoff_is_caught(tmp_path):
    def m(t, d):
        d = d.copy()
        d.loc[d.split_time == "test", "priority_date"] = "2011-01-01"
        return d
    assert any(n.startswith("time_rule") for n in _failed(verify_release(_release(tmp_path, m))))


def test_a_test_row_granted_before_the_floor_is_caught(tmp_path):
    def m(t, d):
        d = d.copy()
        d.loc[d.split_time == "test", "year"] = 2018
        return d
    assert any(n.startswith("time_rule") for n in _failed(verify_release(_release(tmp_path, m))))


def test_an_undated_test_row_is_caught(tmp_path):
    def m(t, d):
        d = d.copy()
        d.loc[d.index[-1], "priority_date"] = None
        return d
    assert any(n.startswith("time_rule") for n in _failed(verify_release(_release(tmp_path, m))))


def test_a_repeated_group_molecule_is_caught(tmp_path):
    def m(t, d):
        return pd.concat([d, d.iloc[[0]]], ignore_index=True)
    assert any(n.startswith("unique_rows") for n in _failed(verify_release(_release(tmp_path, m))))


def test_a_group_in_two_parts_is_caught(tmp_path):
    def m(t, d):
        d = d.copy()
        d.loc[d.index[0], "split_random"] = "test"
        return d
    assert any(n.startswith("group_one_split") for n in _failed(verify_release(_release(tmp_path, m))))


def test_config_drift_is_caught(tmp_path):
    root = _release(tmp_path)
    man = json.loads((root / "release_manifest.json").read_text())
    man["config"]["MIN_GROUP"] = 5
    (root / "release_manifest.json").write_text(json.dumps(man))
    assert "config" in _failed(verify_release(root))


def test_count_drift_is_caught(tmp_path):
    root = _release(tmp_path)
    man = json.loads((root / "release_manifest.json").read_text())
    man["tasks"]["rank_ic50"]["rows"] += 1
    (root / "release_manifest.json").write_text(json.dumps(man))
    assert "manifest_counts[rank_ic50]" in _failed(verify_release(root))


def test_missing_files_are_reported(tmp_path):
    root = _release(tmp_path)
    (root / "groups.parquet").unlink()
    assert "files" in _failed(verify_release(root))
