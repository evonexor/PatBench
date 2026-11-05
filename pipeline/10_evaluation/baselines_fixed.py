"""Reference baselines for the PatBench ranking tasks.

    python baselines_fixed.py <release dir> <out dir> [task ...]

Writes <out dir>/baselines/<task>_<model>[_k<K>].parquet prediction files and prints the leaderboard.

Fixed copy of r6 scripts/baselines.py (2026-11-01). The only changes: predictions go to <out dir>, because the
frozen release is read-only, and the frozen zero-shot forest scored on the few-shot query molecules is written as
`ecfp_rf_zero_k<K>`. In r6 it was written to `ecfp_rf_k<K>` and overwrote the few-shot forest trained on the k
support molecules, so r6's ecfp_rf k5/k10 numbers were zero-shot numbers.

The trivial ones exist to say what a metric is worth: a benchmark where molecular weight already
scores 0.1 Spearman tells you that 0.1 is noise, not signal. A model worth reporting has to beat
`heavy_atoms` by more than its confidence interval.

  random         a seeded random score            - the zero line
  heavy_atoms    heavier is more potent           - the size trend inside a series
  logp           crippen logP                     - the lipophilicity trend
  ecfp_rf        ECFP4 + random forest            - a real model, trained on the task's own
                 train split (zero-shot) or on the k revealed molecules (few-shot)

ecfp_rf zero-shot trains ONE model over every training group on centred pAct (each group's values
minus its own median), because a model that predicts raw pAct across groups is mostly learning
which assays are easy. Few-shot trains one model per (group, repeat) on the k support molecules.
"""
import os
import sys
import warnings

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, rdFingerprintGenerator
from sklearn.ensemble import RandomForestRegressor

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")
D = sys.argv[1].rstrip("/")
O = sys.argv[2].rstrip("/")
TASKS = sys.argv[3:] or ["rank_ic50", "rank_ki_kd", "rank_highcontrast", "rank_complex_mutant"]
SPLIT, PART = "split_time", "test"
KS = (5, 10)
SEED = 20261020
# cores for the one large forest; the many small few-shot forests run on one core each, because
# n_jobs=-1 on hundreds of tiny forests only oversubscribes the machine
JOBS = int(os.environ.get("PATBENCH_JOBS", "4"))
OUT = O + "/baselines"
os.makedirs(OUT, exist_ok=True)
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def fingerprints(smiles):
    X = np.zeros((len(smiles), 2048), dtype=np.int8)
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s)
        if m is not None:
            X[i] = np.frombuffer(gen.GetFingerprintAsNumPy(m).tobytes(), dtype=np.uint8)[:2048]
    return X


def logp(smiles):
    out = []
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        out.append(Crippen.MolLogP(m) if m is not None else 0.0)
    return np.array(out)


for task in TASKS:
    d = pd.read_parquet(f"{D}/tasks/{task}.parquet")
    te = d[d[SPLIT] == PART].copy()
    tr = d[d[SPLIT] == "train"].copy()
    rng = np.random.default_rng([SEED, TASKS.index(task)])
    te["random"] = rng.random(len(te))
    te["heavy_atoms"] = te.heavy_atoms.astype(float)
    te["logp"] = logp(te.SMILES)

    Xtr = fingerprints(tr.SMILES.to_list())
    ytr = (tr.pAct - tr.groupby("group_id").pAct.transform("median")).to_numpy()
    rf = RandomForestRegressor(n_estimators=300, min_samples_leaf=2, n_jobs=JOBS, random_state=SEED)
    rf.fit(Xtr, ytr)
    Xte = fingerprints(te.SMILES.to_list())
    te["ecfp_rf"] = rf.predict(Xte)
    for name in ["random", "heavy_atoms", "logp", "ecfp_rf"]:
        te.rename(columns={name: "score"})[["group_id", "InChIKey", "score"]].to_parquet(
            f"{OUT}/{task}_{name}.parquet", index=False)
    print(f"{task}: trained ecfp_rf on {len(tr):,} rows / {tr.group_id.nunique()} groups, "
          f"scored {len(te):,} test rows", flush=True)

    fp_te = {k: x for k, x in zip(te.InChIKey, Xte)}
    truth = te.set_index(["group_id", "InChIKey"]).pAct
    for k in KS:
        f = f"{D}/protocols/{task}_{SPLIT.replace('split_', '')}_k{k}.parquet"
        if not os.path.exists(f):
            continue
        proto = pd.read_parquet(f)
        recs = []
        for (gid, rep), g in proto.groupby(["group_id", "repeat"], sort=False):
            sup = g[g.role.eq("support")].InChIKey.to_list()
            qry = g[g.role.eq("query")].InChIKey.to_list()
            Xs = np.array([fp_te[i] for i in sup])
            ys = truth.loc[[(gid, i) for i in sup]].to_numpy()
            m = RandomForestRegressor(n_estimators=200, min_samples_leaf=1, n_jobs=1,
                                      random_state=SEED).fit(Xs, ys)
            pred = m.predict(np.array([fp_te[i] for i in qry]))
            recs += [(gid, rep, i, float(p)) for i, p in zip(qry, pred)]
        pd.DataFrame(recs, columns=["group_id", "repeat", "InChIKey", "score"]).to_parquet(
            f"{OUT}/{task}_ecfp_rf_k{k}.parquet", index=False)
        # the zero-shot models need a repeat column too, to be scored on the same query molecules;
        # the zero-shot forest is written as ecfp_rf_zero so that it no longer overwrites the few-shot forest
        for name in ["random", "heavy_atoms", "logp", "ecfp_rf"]:
            s = te.set_index(["group_id", "InChIKey"])[name]
            q = proto[proto.role.eq("query")].copy()
            q["score"] = [float(s.loc[(gi, ik)]) for gi, ik in zip(q.group_id, q.InChIKey)]
            q[["group_id", "repeat", "InChIKey", "score"]].to_parquet(
                f"{OUT}/{task}_{'ecfp_rf_zero' if name == 'ecfp_rf' else name}_k{k}.parquet", index=False)
        print(f"  few-shot k={k}: {proto.group_id.nunique()} groups x {proto.repeat.nunique()} repeats",
              flush=True)
