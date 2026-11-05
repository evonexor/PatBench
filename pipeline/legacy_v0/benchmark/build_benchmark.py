"""Build the USPTO patent affinity-ranking benchmark.

Three products, written to the output directory:

  benchmark/   strict, for scoring models      groups with >= MIN_CORE exact-value molecules,
                                               a single mapped protein target, biochemical assay
  small/       strict but 5-9 molecules        same filters, smaller groups, kept separate
  train_pool/  permissive, for fitting models  exact values, >= MIN_POOL molecules, any target,
                                               with every group that leaks into valid/test removed

Reads only from the extraction results and the target map; writes nothing back to them.
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import json
import sys
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")

R = _ws + "/outputs/uspto-sar-recovery-20261006/results/"
M = _ws + "/outputs/target-uniprot-map-20261016/"
OUT = sys.argv[1].rstrip("/")
SEED = 20261019

AFF = ["IC50", "Ki", "Kd", "EC50"]
KEY = ["patent_core", "Target", "Activity Type"]
MIN_CORE = 10          # molecules per benchmark group
MIN_SMALL = 5          # molecules per small-slice group
MIN_POOL = 5           # molecules per training-pool group
MAX_REPLICATE_SPREAD = 1.0   # log units between repeat measurements of one molecule
MIN_DYNAMIC_RANGE = 1.0      # log units between the weakest and strongest molecule in a group
MIN_DISTINCT_FRAC = 0.3      # distinct values / molecules, to reject disguised grade columns
MIN_DISTINCT_VALUES = 4
HEAVY_MIN, HEAVY_MAX = 10, 70
FAMILY_OVERLAP = 0.5         # shared molecules / smaller group, above which two groups are one family
TEST_FROM_YEAR = 2021        # patents filed/issued from this year form the test set
VALID_PATENT_FRACTION = 0.15 # share of pre-2021 patents held out for validation
ORGANIC = set("C H N O S P F Cl Br I B Si Se".split())
cfg = {k: v for k, v in sorted(globals().items())
       if k.isupper() and k not in {"R", "M", "OUT", "AFF", "KEY", "ORGANIC"}}
log = []


def note(msg):
    print(msg, flush=True)
    log.append(msg)


def describe(label, d):
    s = {"stage": label, "rows": len(d), "groups": int(d.groupby(KEY).ngroups),
         "patents": int(d.patent_core.nunique()), "molecules": int(d.InChIKey.nunique())}
    note(f"{label:54s} rows {s['rows']:>8,d}  groups {s['groups']:>6,d}  "
         f"patents {s['patents']:>5,d}  molecules {s['molecules']:>7,d}")
    return s


funnel = []

# ---------------------------------------------------------------- source rows
m = ds.dataset(R + "SAR_merged.parquet").to_table(
    columns=["source", "source_record_id", "patent_core", "Target", "Activity Type",
             "Relation", "value_nM", "SMILES", "InChIKey"]).to_pandas()
m = m[m["Activity Type"].isin(AFF)].copy()
m["Target"] = m.Target.fillna("").str.strip()
funnel.append(describe("verified affinity rows", m))

nk = ds.dataset(R + "new_SAR.parquet").to_table(
    columns=["source_record_id", "measurement_kind"]).to_pandas()
graded = set(nk.loc[nk.measurement_kind.isin(["legend_threshold", "legend_interval"]),
                    "source_record_id"])
m = m[~m.source_record_id.isin(graded)]
funnel.append(describe("drop rows derived from A/B/C grades", m))

m = m[m.Relation.eq("=") & m.value_nM.gt(0)].copy()
m["pAct"] = 9 - np.log10(m.value_nM)
m = m[(m.pAct >= 2) & (m.pAct <= 13)]
funnel.append(describe("exact '=' values, pAct in 2..13", m))

# -------------------------------------------- one value per (group, molecule)
dup = m.groupby(KEY + ["InChIKey"]).pAct.agg(["median", "max", "min", "size"])
spread = (dup["max"] - dup["min"])[lambda s: s > 0]
note(f"repeat measurements: {int((dup['size'] > 1).sum()):,} (group, molecule) pairs measured more "
     f"than once; median spread {spread.median():.3f} log, {int((spread > MAX_REPLICATE_SPREAD).sum()):,} "
     f"above {MAX_REPLICATE_SPREAD} log")
note(f"ASSAY NOISE FLOOR: median absolute difference between repeat measurements is "
     f"{spread.median():.3f} log units - no ranking metric can do better than this")
agg = m.groupby(KEY + ["InChIKey"]).agg(
    pAct=("pAct", "median"), pAct_spread=("pAct", lambda s: s.max() - s.min()),
    n_measurements=("pAct", "size"), SMILES=("SMILES", "first"),
    value_nM=("value_nM", "median"), source=("source", "first")).reset_index()
agg = agg[agg.pAct_spread <= MAX_REPLICATE_SPREAD]
funnel.append(describe("one value per (group, molecule)", agg))

# -------------------------------------------------------------- molecule sanity
uniq = agg.drop_duplicates("InChIKey")[["InChIKey", "SMILES"]]
heavy, odd = {}, {}
for k, s in zip(uniq.InChIKey, uniq.SMILES):
    mol = Chem.MolFromSmiles(s) if isinstance(s, str) else None
    if mol is None:
        heavy[k], odd[k] = -1, True
        continue
    heavy[k] = mol.GetNumHeavyAtoms()
    odd[k] = any(at.GetSymbol() not in ORGANIC for at in mol.GetAtoms())
agg["heavy_atoms"] = agg.InChIKey.map(heavy)
agg["unusual_elements"] = agg.InChIKey.map(odd)
note(f"molecules: {int((agg.heavy_atoms < 0).sum()):,} unparseable rows, "
     f"{int((agg.heavy_atoms.between(0, HEAVY_MIN - 1)).sum()):,} under {HEAVY_MIN} heavy atoms, "
     f"{int((agg.heavy_atoms > HEAVY_MAX).sum()):,} over {HEAVY_MAX}, "
     f"{int(agg.unusual_elements.sum()):,} with elements outside {sorted(ORGANIC)}")
agg = agg[agg.heavy_atoms.between(HEAVY_MIN, HEAVY_MAX) & ~agg.unusual_elements]
funnel.append(describe(f"drug-like molecules ({HEAVY_MIN}-{HEAVY_MAX} heavy atoms)", agg))

# ------------------------------------------------------------------- patent year
pt = ds.dataset(R + "patents.parquet").to_table(columns=["Patent", "archive"]).to_pandas()
pt["patent_core"] = pt.Patent.str.replace(r"(B\d|A\d|E|S|H|P\d)$", "", regex=True)
pt["year"] = pt.archive.str.extract(r"/data/(\d{4})/")[0].astype(int)
year = pt.drop_duplicates("patent_core").set_index("patent_core").year
agg["year"] = agg.patent_core.map(year)
note(f"patent year missing for {int(agg.year.isna().sum()):,} rows - dropped")
agg = agg[agg.year.notna()].copy()
agg["year"] = agg.year.astype(int)

# ------------------------------------------------------------------ target map
tm = pq.read_table(M + "ranking_group_targets.parquet").to_pandas()
agg = agg.join(tm.set_index(KEY)[["target_class", "confidence", "uniprot_accession",
                                  "mapped_identity", "measured_in_cells", "seq_len"]], on=KEY)
seq = pq.read_table(M + "target_map.parquet",
                    columns=["uniprot_accession", "sequence", "protein_name", "organism"]).to_pandas()
seq = seq[seq.uniprot_accession.fillna("").ne("")].drop_duplicates("uniprot_accession")

# ------------------------------------------------------- strict benchmark filters
strict = agg[agg.confidence.isin(["high", "medium"])
             & agg.uniprot_accession.fillna("").ne("")
             & agg.seq_len.gt(0)]
funnel.append(describe("target mapped to a sequence (high/medium)", strict))
strict = strict[strict.target_class.eq("protein")]
funnel.append(describe("single protein target", strict))
strict = strict[strict.measured_in_cells.ne("Y")]
funnel.append(describe("biochemical assay (not measured in cells)", strict))

gstat = strict.groupby(KEY).agg(
    n_molecules=("InChIKey", "size"), pAct_range=("pAct", lambda s: s.max() - s.min()),
    n_distinct=("pAct", lambda s: s.round(2).nunique())).reset_index()
usable = gstat[(gstat.pAct_range >= MIN_DYNAMIC_RANGE)
               & (gstat.n_distinct >= MIN_DISTINCT_VALUES)
               & (gstat.n_distinct / gstat.n_molecules >= MIN_DISTINCT_FRAC)]
note(f"groups rejected for a flat or grade-like value column: "
     f"{len(gstat) - len(usable):,} of {len(gstat):,}")
strict = strict.merge(usable[KEY], on=KEY)
funnel.append(describe(f"dynamic range >= {MIN_DYNAMIC_RANGE} log, real numbers", strict))


# ----------------------------------------------------------------- family dedup
def family_clusters(frame, overlap):
    """Union groups that measure the same endpoint on a largely shared molecule set."""
    g = frame.groupby(KEY).InChIKey.agg(frozenset).reset_index(name="mols")
    g["n"] = g.mols.map(len)
    index = defaultdict(list)
    for i, (et, mols) in enumerate(zip(g["Activity Type"], g.mols)):
        for k in mols:
            index[(et, k)].append(i)
    shared = defaultdict(int)
    for members in index.values():
        if 1 < len(members) < 200:
            for i in range(len(members)):
                for j in range(i + 1, len(members)):
                    shared[(members[i], members[j])] += 1
    parent = list(range(len(g)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    size = g.n.values
    for (i, j), c in shared.items():
        if c / min(size[i], size[j]) >= overlap:
            ri, rj = find(i), find(j)
            if ri != rj:
                parent[ri] = rj
    g["family"] = [find(i) for i in range(len(g))]
    return g


fam = family_clusters(strict, FAMILY_OVERLAP)
fam = fam.merge(strict.groupby(KEY).year.min().reset_index(), on=KEY)
note(f"family clustering: {len(fam):,} groups collapse into {fam.family.nunique():,} families "
     f"({int((fam.groupby('family').size() > 1).sum()):,} families hold more than one group)")
# keep the earliest, then largest, group in each family
keep = fam.sort_values(["family", "year", "n"], ascending=[True, True, False]).drop_duplicates("family")
strict = strict.merge(keep[KEY], on=KEY)
funnel.append(describe("one group per patent family", strict))

# ------------------------------------------------------------- core vs small
size = strict.groupby(KEY).InChIKey.size()
core = strict.merge(size[size >= MIN_CORE].reset_index()[KEY], on=KEY)
small = strict.merge(size[size.between(MIN_SMALL, MIN_CORE - 1)].reset_index()[KEY], on=KEY)
funnel.append(describe(f"BENCHMARK: groups with >= {MIN_CORE} molecules", core))
funnel.append(describe(f"SMALL SLICE: groups with {MIN_SMALL}-{MIN_CORE - 1} molecules", small))

# ------------------------------------------------------------- public overlap
pub = pq.read_table(OUT + "/work/public_pairs.parquet").to_pandas()
pub_pairs = set(zip(pub.ik, pub.uniprot))
pub_patents = json.load(open(OUT + "/work/public_patents.json"))
all_pub_patents = {p for v in pub_patents.values() for p in v}


def mark_public(frame):
    frame = frame.copy()
    frame["in_public_data"] = [(k, u) in pub_pairs
                               for k, u in zip(frame.InChIKey, frame.uniprot_accession)]
    return frame


core, small = mark_public(core), mark_public(small)
note(f"rows whose (molecule, UniProt) pair already exists in BindingDB/ChEMBL/HARVEST: "
     f"{core.in_public_data.mean():.1%} of the benchmark")

# ------------------------------------------------------------------- splits
rng = np.random.default_rng(SEED)
core["split"] = np.where(core.year >= TEST_FROM_YEAR, "test", "train")
train_patents = sorted(core.loc[core.split == "train", "patent_core"].unique())
n_valid = int(round(len(train_patents) * VALID_PATENT_FRACTION))
valid_patents = set(rng.choice(train_patents, size=n_valid, replace=False))
core.loc[core.patent_core.isin(valid_patents), "split"] = "valid"
note(f"validation holds out {len(valid_patents):,} of {len(train_patents):,} pre-{TEST_FROM_YEAR} patents")

# test tiering: a group is "novel" when none of its rows appear in the public datasets
pub_share = core[core.split == "test"].groupby(KEY).in_public_data.mean()
novel = set(pub_share[pub_share == 0].reset_index()[KEY].apply(tuple, axis=1))
core["test_tier"] = np.where(
    core.split != "test", "",
    np.where(core[KEY].apply(tuple, axis=1).isin(novel), "novel", "reference"))
note(f"test tiers: {len(novel):,} novel groups (absent from public data), "
     f"{int(pub_share.shape[0] - len(novel)):,} reference groups")

held_out = core[core.split.isin(["valid", "test"])]
unseen_targets = set(held_out.uniprot_accession) - set(core.loc[core.split == "train", "uniprot_accession"])
core["unseen_target"] = core.uniprot_accession.isin(unseen_targets) & core.split.ne("train")
note(f"held-out targets never seen in training: {len(unseen_targets):,}")

train_pairs = set(zip(core.loc[core.split == "train", "InChIKey"],
                      core.loc[core.split == "train", "uniprot_accession"]))
leak = [(k, u) in train_pairs for k, u in zip(held_out.InChIKey, held_out.uniprot_accession)]
note(f"held-out rows sharing a (molecule, target) pair with training: {int(np.sum(leak)):,} "
     f"of {len(held_out):,} - these carry leaks_from_train = True")
core["leaks_from_train"] = False
core.loc[held_out.index[np.array(leak)], "leaks_from_train"] = True

# ------------------------------------------------- permissive training pool
pool_size = agg.groupby(KEY).InChIKey.size()
pool = agg.merge(pool_size[pool_size >= MIN_POOL].reset_index()[KEY], on=KEY)
funnel.append(describe(f"training pool, groups with >= {MIN_POOL} molecules", pool))
# drop any pool group that overlaps a held-out group on the same endpoint
held_sets = {}
for (pc, tg, et), d in held_out.groupby(KEY):
    held_sets.setdefault(et, []).append(frozenset(d.InChIKey))
drop_keys = []
for (pc, tg, et), d in pool.groupby(KEY):
    mine = frozenset(d.InChIKey)
    for other in held_sets.get(et, []):
        if len(mine & other) / min(len(mine), len(other)) >= FAMILY_OVERLAP:
            drop_keys.append((pc, tg, et))
            break
if drop_keys:
    dk = pd.DataFrame(drop_keys, columns=KEY)
    pool = pool.merge(dk, on=KEY, how="left", indicator=True)
    pool = pool[pool._merge == "left_only"].drop(columns="_merge")
note(f"training pool: removed {len(drop_keys):,} groups that overlap a held-out group")
pool = pool[~pool.patent_core.isin(set(held_out.patent_core))]
funnel.append(describe("training pool after leakage removal", pool))

# ----------------------------------------------------------------- write out
seq_cols = seq.rename(columns={"sequence": "target_sequence"})
core = core.merge(seq_cols[["uniprot_accession", "target_sequence", "protein_name", "organism"]],
                  on="uniprot_accession", how="left")
small = small.merge(seq_cols[["uniprot_accession", "target_sequence", "protein_name", "organism"]],
                    on="uniprot_accession", how="left")

keep_cols = ["patent_core", "year", "Target", "Activity Type", "uniprot_accession", "mapped_identity",
             "protein_name", "organism", "target_sequence", "InChIKey", "SMILES", "pAct", "value_nM",
             "heavy_atoms", "n_measurements", "pAct_spread", "in_public_data", "split", "test_tier",
             "unseen_target", "leaks_from_train", "source"]
core[keep_cols].to_parquet(OUT + "/benchmark_rows.parquet", index=False)
small[[c for c in keep_cols if c in small]].to_parquet(OUT + "/small_groups_rows.parquet", index=False)
pool_cols = ["patent_core", "year", "Target", "Activity Type", "uniprot_accession", "confidence",
             "target_class", "measured_in_cells", "InChIKey", "SMILES", "pAct", "value_nM",
             "heavy_atoms", "source"]
pool[[c for c in pool_cols if c in pool]].to_parquet(OUT + "/train_pool_rows.parquet", index=False)

groups = core.groupby(KEY + ["split", "test_tier", "uniprot_accession", "mapped_identity", "year"]).agg(
    n_molecules=("InChIKey", "size"), pAct_min=("pAct", "min"), pAct_max=("pAct", "max"),
    pAct_std=("pAct", "std"), public_share=("in_public_data", "mean")).reset_index()
groups.to_parquet(OUT + "/benchmark_groups.parquet", index=False)

for name, d in [("train", core[core.split == "train"]), ("valid", core[core.split == "valid"]),
                ("test", core[core.split == "test"])]:
    d[keep_cols].to_parquet(f"{OUT}/splits/{name}.parquet", index=False)

summary = {
    "built": "2026-10-19", "seed": SEED, "config": cfg,
    "inputs": {"sar": R + "SAR_merged.parquet", "target_map": M + "ranking_group_targets.parquet",
               "public": {k: len(v) for k, v in pub_patents.items()}},
    "assay_noise_floor_log": round(float(spread.median()), 3),
    "funnel": funnel,
    "benchmark": {},
    "small_slice": {"rows": len(small), "groups": int(small.groupby(KEY).ngroups),
                    "patents": int(small.patent_core.nunique())},
    "train_pool": {"rows": len(pool), "groups": int(pool.groupby(KEY).ngroups),
                   "patents": int(pool.patent_core.nunique()),
                   "molecules": int(pool.InChIKey.nunique())},
}
for sp in ["train", "valid", "test"]:
    d = core[core.split == sp]
    summary["benchmark"][sp] = {
        "rows": len(d), "groups": int(d.groupby(KEY).ngroups), "patents": int(d.patent_core.nunique()),
        "molecules": int(d.InChIKey.nunique()), "targets": int(d.uniprot_accession.nunique()),
        "by_endpoint": {k: int(v) for k, v in d.groupby("Activity Type").apply(
            lambda x: x.groupby(KEY).ngroups).items()},
        "public_share": round(float(d.in_public_data.mean()), 3),
    }
summary["benchmark"]["test"]["tiers"] = {
    t: int(core[(core.split == "test") & (core.test_tier == t)].groupby(KEY).ngroups)
    for t in ["novel", "reference"]}
summary["log"] = log
with open(OUT + "/splits_manifest.json", "w") as fh:
    json.dump(summary, fh, indent=1, ensure_ascii=False)

note("\n=== benchmark groups by endpoint and split")
tab = core.groupby(["Activity Type", "split"]).apply(lambda d: d.groupby(KEY).ngroups).unstack(fill_value=0)
note(tab.to_string())
note(f"\nwrote {OUT}/benchmark_rows.parquet, small_groups_rows.parquet, train_pool_rows.parquet, "
     f"benchmark_groups.parquet, splits/, splits_manifest.json")
