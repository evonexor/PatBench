"""Integrity checks on a built PatBench release directory.

These are the checks the benchmark's design promises, run against the files rather than the code:

  files            the release has every table the README lists
  unique_rows      one value per (group, molecule) in every ranking task
  group_one_split  a group is in one split only, for every split column
  no_patent_leak   no patent in two of train / valid / test (time split and random split)
  no_family_leak   no patent family in two of them
  no_pair_leak     no (molecule, target) pair shared between train and test in a task file
  time_rule        every test row: earliest priority year >= the task's cutoff, grant year >= the floor
  config           the parameters recorded in release_manifest.json equal configs/release_v1.0.json
  manifest_counts  rows / groups / patents / molecules per task equal the manifest's

`split_scaffold_group` is reported but not required to be patent-disjoint (the README says so).
"""
import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from patbench import CONFIGS

TASKS = ["rank_ic50", "rank_ki_kd", "rank_ec50", "rank_highcontrast", "rank_complex_mutant"]
DISJOINT_SPLITS = ["split_time", "split_random"]
ALL_SPLITS = ["split_time", "split_time_grant", "split_random", "split_scaffold_group"]
PARTS = {"train", "valid", "test"}
NOT_COMPARED = {"VAL", "BUILT"}          # a server path and a date, not parameters
REQUIRED = ["sar_all.parquet", "groups.parquet", "targets.parquet", "release_manifest.json",
            *[f"tasks/{t}.parquet" for t in TASKS], "tasks/regression_affinity.parquet"]


@dataclass
class Check:
    name: str
    ok: bool
    detail: str = ""


def _crossing(df, key, split):
    """Values of `key` that occur in more than one part of `split`."""
    d = df[[key, split]].dropna()
    d = d[d[split].isin(PARTS)]
    n = d.groupby(key)[split].nunique()
    return n[n > 1].index.tolist()


def verify_release(path, config_file=None):
    root = Path(path)
    out = []
    missing = [f for f in REQUIRED if not (root / f).exists()]
    out.append(Check("files", not missing, "all present" if not missing else "missing: " + ", ".join(missing)))
    man = json.loads((root / "release_manifest.json").read_text()) if (root / "release_manifest.json").exists() else {}
    cfg = man.get("config", {})
    cutoffs = {t: cfg.get("TASK_TEST_FROM", {}).get(t, cfg.get("TEST_FROM_PRIORITY_YEAR")) for t in TASKS}
    floor = cfg.get("TEST_MIN_GRANT_YEAR")

    ref = json.loads(Path(config_file or CONFIGS / "release_v1.0.json").read_text())
    drift = {k: (cfg.get(k), v) for k, v in ref.items() if k not in NOT_COMPARED and cfg.get(k) != v}
    drift |= {k: (cfg[k], None) for k in cfg if k not in ref and k not in NOT_COMPARED}
    out.append(Check("config", not drift and bool(cfg), "equal to configs/release_v1.0.json" if not drift and cfg
                     else f"differs: {drift}" if drift else "no config in release_manifest.json"))

    cols = ["group_id", "patent_core", "family_id", "InChIKey", "uniprot_accession", "year", "priority_date",
            *ALL_SPLITS]
    for t in TASKS:
        f = root / f"tasks/{t}.parquet"
        if not f.exists():
            continue
        d = pd.read_parquet(f)
        d = d[[c for c in cols if c in d.columns]]
        dup = int(d.duplicated(["group_id", "InChIKey"]).sum())
        out.append(Check(f"unique_rows[{t}]", dup == 0, f"{len(d):,} rows, {dup} repeated (group, molecule)"))
        for s in [c for c in ALL_SPLITS if c in d.columns]:
            vals = set(d[s].dropna().unique()) - PARTS - {""}
            n_multi = int((d.groupby("group_id")[s].nunique() > 1).sum())
            out.append(Check(f"group_one_split[{t}:{s}]", n_multi == 0 and not vals,
                             f"{n_multi} groups in more than one part" + (f"; unexpected values {sorted(vals)}" if vals else "")))
        for s in [c for c in DISJOINT_SPLITS if c in d.columns]:
            pat = _crossing(d, "patent_core", s)
            out.append(Check(f"no_patent_leak[{t}:{s}]", not pat, f"{len(pat)} patents cross parts"))
            if "family_id" in d.columns:
                fam = _crossing(d[d.family_id.notna() & (d.family_id.astype(str) != "")], "family_id", s)
                out.append(Check(f"no_family_leak[{t}:{s}]", not fam, f"{len(fam)} families cross parts"))
        if "split_time" in d.columns and {"InChIKey", "uniprot_accession"} <= set(d.columns):
            p = d[d.uniprot_accession.notna() & (d.uniprot_accession.astype(str) != "")]
            tr = set(zip(p.InChIKey[p.split_time == "train"], p.uniprot_accession[p.split_time == "train"]))
            te = set(zip(p.InChIKey[p.split_time == "test"], p.uniprot_accession[p.split_time == "test"]))
            out.append(Check(f"no_pair_leak[{t}]", not (tr & te), f"{len(tr & te)} (molecule, target) pairs in train and test"))
        if "split_time" in d.columns and cutoffs[t] is not None:
            te = d[d.split_time == "test"]
            py = pd.to_datetime(te.priority_date, errors="coerce").dt.year
            bad_prio = int((py.isna() | (py < cutoffs[t])).sum())
            bad_year = int((te.year < floor).sum()) if floor is not None and "year" in te else 0
            out.append(Check(f"time_rule[{t}]", bad_prio == 0 and bad_year == 0,
                             f"{len(te):,} test rows; {bad_prio} with priority year < {cutoffs[t]} or undated, "
                             f"{bad_year} granted before {floor}"))
        stats = man.get("tasks", {}).get(t, {})
        have = {"rows": len(d), "groups": int(d.group_id.nunique()), "patents": int(d.patent_core.nunique()),
                "molecules": int(d.InChIKey.nunique())}
        diff = {k: (have[k], stats[k]) for k in have if k in stats and stats[k] != have[k]}
        out.append(Check(f"manifest_counts[{t}]", bool(stats) and not diff,
                         "equal" if stats and not diff else f"differs (file, manifest): {diff}" if diff else "task not in manifest"))
    return out


def report(checks):
    w = max(len(c.name) for c in checks)
    for c in checks:
        print(f"{'OK  ' if c.ok else 'FAIL'} {c.name:<{w}}  {c.detail}")
    bad = [c for c in checks if not c.ok]
    print(f"\n{len(checks) - len(bad)} of {len(checks)} checks passed")
    return not bad
