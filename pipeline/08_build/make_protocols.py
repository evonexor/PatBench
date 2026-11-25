"""Freeze the few-shot supports and the activity-cliff pairs, so every model sees the same ones.

    python make_protocols.py <release dir>

Few-shot (industrial use: a project has measured a handful of compounds of a new series and wants
the rest ranked): in every test group, K molecules are revealed with their values (role
"support") and the model ranks the others (role "query"). K = 5 and 10, five draws per group;
a group takes part only if at least 10 query molecules remain.
    protocols/<task>_<split>_k<K>.parquet: group_id, repeat, InChIKey, role

Activity cliffs: pairs inside one group that are structurally near-identical (ECFP4 Tanimoto
>= 0.85) yet differ by >= 1 log (10-fold) in potency. Pairs inside one stereo set are left out:
their 2D fingerprints are identical by construction, and they are reported by --collapse-stereo.
    protocols/<task>_cliffs.parquet: group_id, InChIKey_a, InChIKey_b, tanimoto, delta
"""
import os
import sys
import warnings

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")
D = sys.argv[1].rstrip("/")
OUT = D + "/protocols"
TASKS = ["rank_ic50", "rank_ki_kd", "rank_ec50", "rank_highcontrast", "rank_complex_mutant"]
SPLITS = ["split_time", "split_random"]
KS = (5, 10)
REPEATS = 5
MIN_QUERY = 10
CLIFF_SIM, CLIFF_DELTA = 0.85, 1.0
SEED = 20261020
os.makedirs(OUT, exist_ok=True)
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

for task in TASKS:
    d = pd.read_parquet(f"{D}/tasks/{task}.parquet")
    for split in SPLITS:
        te = d[d[split] == "test"]
        for k in KS:
            rng = np.random.default_rng([SEED, k, TASKS.index(task), SPLITS.index(split)])
            rows = []
            for gid, g in te.groupby("group_id", sort=True):
                keys = g.InChIKey.to_numpy()
                if len(keys) < k + MIN_QUERY:
                    continue
                for r in range(REPEATS):
                    sup = set(rng.choice(keys, size=k, replace=False))
                    rows += [(gid, r, key, "support" if key in sup else "query") for key in keys]
            f = pd.DataFrame(rows, columns=["group_id", "repeat", "InChIKey", "role"])
            f.to_parquet(f"{OUT}/{task}_{split.replace('split_', '')}_k{k}.parquet", index=False)
            print(f"{task:20s} {split:12s} k={k:<2d} groups {f.group_id.nunique():>4} of {te.group_id.nunique()}")
    fps = {k: gen.GetFingerprint(Chem.MolFromSmiles(s)) for k, s in zip(d.InChIKey, d.SMILES)}
    pairs = []
    for gid, g in d.groupby("group_id", sort=True):
        keys, vals, st = g.InChIKey.to_list(), g.pAct.to_numpy(), g.stereo_set_id.to_list()
        fl = [fps[k] for k in keys]
        for i in range(len(keys) - 1):
            sims = DataStructs.BulkTanimotoSimilarity(fl[i], fl[i + 1:])
            for j, s in enumerate(sims, start=i + 1):
                if s >= CLIFF_SIM and abs(vals[i] - vals[j]) >= CLIFF_DELTA and not (st[i] and st[i] == st[j]):
                    pairs.append((gid, keys[i], keys[j], round(s, 3), round(vals[i] - vals[j], 3)))
    c = pd.DataFrame(pairs, columns=["group_id", "InChIKey_a", "InChIKey_b", "tanimoto", "delta"])
    c.to_parquet(f"{OUT}/{task}_cliffs.parquet", index=False)
    sp = d.drop_duplicates("group_id").set_index("group_id").split_time
    print(f"{task:20s} cliff pairs {len(c):,} in {c.group_id.nunique()} groups; "
          f"time-split test {int(c.group_id.map(sp).eq('test').sum()):,}")
