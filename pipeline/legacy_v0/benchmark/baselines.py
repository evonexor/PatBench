"""Reference baselines for the affinity-ranking benchmark.

  random          scores drawn at random - confirms the metrics sit at 0 when nothing is learned
  heavy_atoms     bigger molecule = more potent; the trivial chemistry prior a real model must beat
  ecfp_knn        nearest training molecule on the same target by Tanimoto, take its pAct
  ecfp_hgb        global histogram gradient boosting on ECFP4 counts + target as a categorical

Any new model only has to emit (patent_core, Target, Activity Type, InChIKey, score) and can be
scored with the same evaluate.py.
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import json
import sys
import warnings

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import HistGradientBoostingRegressor

import evaluate as ev

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")
KEY = ev.KEY
OUT = sys.argv[1].rstrip("/")
NBITS = 1024
RADIUS = 2
SEED = 20261019

bench = pd.read_parquet(OUT + "/benchmark_rows.parquet")
pool = pd.read_parquet(OUT + "/train_pool_rows.parquet")
held = bench[bench.split.isin(["valid", "test"])].copy()
print(f"benchmark {len(bench):,} rows; held-out {len(held):,} rows "
      f"({held.groupby(KEY).ngroups:,} groups); training pool {len(pool):,} rows", flush=True)

gen = rdFingerprintGenerator.GetMorganGenerator(radius=RADIUS, fpSize=NBITS)
smiles = pd.concat([bench[["InChIKey", "SMILES"]], pool[["InChIKey", "SMILES"]]]).drop_duplicates("InChIKey")
print(f"featurising {len(smiles):,} molecules", flush=True)
fps, order = {}, []
for k, s in zip(smiles.InChIKey, smiles.SMILES):
    mol = Chem.MolFromSmiles(s)
    if mol is None:
        continue
    fps[k] = gen.GetCountFingerprint(mol)
    order.append(k)
print(f"  {len(fps):,} fingerprints", flush=True)

mat = np.zeros((len(order), NBITS), dtype=np.uint8)
idx = {k: i for i, k in enumerate(order)}
for k in order:
    for bit, count in fps[k].GetNonzeroElements().items():
        mat[idx[k], bit] = min(count, 255)

results, tables = {}, {}


def score_model(name, pred):
    gm = ev.per_group(pred, held)
    head = ev.report(gm, name)
    gm.to_parquet(f"{OUT}/work/groups_{name}.parquet", index=False)
    slices = pd.concat([
        ev.aggregate(gm, "split").assign(slice_kind="split"),
        ev.aggregate(gm[gm.split == "test"], "test_tier").assign(slice_kind="test tier"),
        ev.aggregate(gm, "Activity Type").assign(slice_kind="endpoint"),
        ev.aggregate(gm[gm.split == "test"], "unseen_target").assign(slice_kind="unseen target"),
    ], ignore_index=True)
    tables[name] = slices
    results[name] = {
        "overall": head.to_dict("records")[0],
        "slices": slices.to_dict("records"),
    }
    print(slices[["slice_kind", "slice", "groups", "spearman_median",
                  "pair_accuracy_median", "ndcg10_median"]].round(3).to_string(index=False))


# ----------------------------------------------------------------- random
rng = np.random.default_rng(SEED)
score_model("random", held.assign(score=rng.normal(size=len(held))))

# ------------------------------------------------------------ heavy atoms
score_model("heavy_atoms", held.assign(score=held.heavy_atoms.astype(float)))

# ------------------------------------------------- nearest training neighbour
train = pool[pool.uniprot_accession.fillna("").ne("") & pool.InChIKey.isin(fps)]
global_mean = float(train.pAct.mean())
# one fingerprint list and one pAct per (target, molecule), built once
by_target = {}
for acc, d in train.groupby("uniprot_accession"):
    p = d.groupby("InChIKey").pAct.median()
    keys = list(p.index)
    by_target[acc] = (keys, [fps[k] for k in keys], p.to_numpy(float),
                      {k: i for i, k in enumerate(keys)})
scores = np.full(len(held), np.nan)
hits = 0
for i, (k, acc) in enumerate(zip(held.InChIKey, held.uniprot_accession)):
    entry = by_target.get(acc)
    if entry is None or k not in fps:
        continue
    keys, keyfps, pacts, pos = entry
    sims = np.array(DataStructs.BulkTanimotoSimilarity(fps[k], keyfps))
    self_i = pos.get(k)
    if self_i is not None:
        sims[self_i] = -1.0          # never score a molecule with its own measurement
    if sims.max() < 0:
        continue
    scores[i] = pacts[int(sims.argmax())]
    hits += 1
print(f"\nknn: found a same-target training neighbour for {hits:,} of {len(held):,} held-out rows",
      flush=True)
score_model("ecfp_knn", held.assign(score=np.where(np.isnan(scores), global_mean, scores)))

# ---------------------------------------------------- global gradient boosting
fit = pool[pool.uniprot_accession.fillna("").ne("") & pool.InChIKey.isin(idx)]
fit = fit[~fit.patent_core.isin(set(held.patent_core))]
ep = sorted(set(fit["Activity Type"]) | set(held["Activity Type"]))
ep_code = {e: i for i, e in enumerate(ep)}

# the target is described by its amino-acid composition, not by its identity: a one-hot would
# cap out at sklearn's 255-category limit and would say nothing about an unseen protein
AA = list("ACDEFGHIKLMNPQRSTVWY")
tmap = pd.read_parquet(_ws + "/outputs/target-uniprot-map-20261016/target_map.parquet",
                       columns=["uniprot_accession", "sequence"])
tmap = tmap[tmap.uniprot_accession.fillna("").ne("") & tmap.sequence.fillna("").ne("")]
seq_of = dict(zip(tmap.uniprot_accession, tmap.sequence))
comp = {}
for acc, s in seq_of.items():
    n = len(s)
    comp[acc] = (np.array([s.count(a) for a in AA], dtype=np.float32) / n if n else
                 np.zeros(len(AA), dtype=np.float32))
zero = np.zeros(len(AA), dtype=np.float32)


def design(frame):
    rows = np.array([idx[k] for k in frame.InChIKey])
    x = mat[rows].astype(np.float32)
    aa = np.vstack([comp.get(a, zero) for a in frame.uniprot_accession])
    length = np.array([[len(seq_of.get(a, ""))] for a in frame.uniprot_accession], dtype=np.float32)
    code = frame["Activity Type"].map(ep_code).fillna(-1).to_numpy(np.float32)[:, None]
    return np.hstack([x, aa, length, code])


print(f"fitting HistGradientBoosting on {len(fit):,} rows, "
      f"{fit.uniprot_accession.nunique():,} targets", flush=True)
X = design(fit)
cat = np.zeros(X.shape[1], dtype=bool)
cat[-1] = True          # endpoint only; 4 categories


def hgb():
    return HistGradientBoostingRegressor(
        max_iter=250, learning_rate=0.1, max_leaf_nodes=63, min_samples_leaf=20,
        l2_regularization=1.0, max_features=0.3, categorical_features=cat,
        random_state=SEED, early_stopping=False)


Xh = design(held)
model = hgb().fit(X, fit.pAct.to_numpy(np.float32))
score_model("ecfp_hgb", held.assign(score=model.predict(Xh)))

# same features, but fitted to the within-group deviation, which is what the metric asks for
centred = (fit.pAct - fit.groupby(KEY).pAct.transform("mean")).to_numpy(np.float32)
model_c = hgb().fit(X, centred)
score_model("ecfp_hgb_centred", held.assign(score=model_c.predict(Xh)))

# ----------------------------------------------------------------- write out
manifest = json.load(open(OUT + "/splits_manifest.json"))
summary = {"assay_noise_floor_log": manifest.get("assay_noise_floor_log"),
           "pair_margin_log": ev.PAIR_MARGIN, "bootstrap_samples": ev.BOOTSTRAP,
           "fingerprint": f"Morgan radius {RADIUS}, {NBITS} count bits", "baselines": results}
with open(OUT + "/baseline_results.json", "w") as fh:
    json.dump(summary, fh, indent=1, default=float)
pd.concat([t.assign(model=n) for n, t in tables.items()], ignore_index=True).to_parquet(
    OUT + "/baseline_slices.parquet", index=False)

print("\n=== headline: median per-group Spearman on the test split")
for name in results:
    gm = pd.read_parquet(f"{OUT}/work/groups_{name}.parquet")
    t = gm[gm.split == "test"]
    nov = t[t.test_tier == "novel"]
    print(f"  {name:14s} all test {np.nanmedian(t.spearman):+.3f}   "
          f"novel tier {np.nanmedian(nov.spearman):+.3f}  ({len(nov):,} groups)")
print(f"\nwrote {OUT}/baseline_results.json and baseline_slices.parquet")
