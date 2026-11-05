"""Assemble PatBench: the full formal SAR table plus quality and split flags, then cut the task
subsets from it.

Nothing is deleted. Every selection the paper reports is reproducible from flag columns on
`sar_all.parquet`, which is the point of publishing it this way.

v1.0 adds the industrial-grade cleaning on top of the 2026-10-19 release:
  - molecules are keyed by their parent (largest organic fragment, neutralised); the drawn form is kept
  - one value per (group, compound) in every task file: the median, with the count and spread
  - group-level unit and column checks, within-patent duplicate tables, uninformative labels
  - target labels normalised before grouping, so "BTK" and "Btk" in one patent are one assay
  - every task row traces back to its patent cell (compound label, table locator, header, raw value)
  - the time split follows the priority date, when the chemistry was made, not the grant year
  - stereo sets, undefined stereocentres, warheads, species and out-of-range values are flagged
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import json
import os
import re
import sys
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")

NAME, VERSION, BUILT = "PatBench", "1.0", "2026-10-30"
R = os.environ.get("PATBENCH_SAR_DIR", _ws + "/outputs/uspto-sar-recovery-20261006/results/")
FULL = _ws + "/outputs/uspto-sar-full-20261004/results/"
M = os.environ.get("PATBENCH_TARGET_MAP_DIR", _ws + "/outputs/target-uniprot-map-v2-20261019/")
VAL = _ws + "/outputs/sar-validation-20261019"
OUT = sys.argv[1].rstrip("/")
SEED = 20261019

AFFINITY = ["IC50", "Ki", "Kd", "EC50"]
# labels are compared case- and punctuation-blind, and a group never mixes two mapped proteins
KEY = ["patent_core", "Target_norm", "Activity Type", "_uniprot_key"]
# HARVEST is excluded from the public-coverage comparison: its public snapshot is only the
# H-bench subset (48 targets, 1,617 patents), so it would understate what is new here.
PUBLIC_SOURCES = ["BindingDB", "ChEMBL"]
ORGANIC = set("C H N O S P F Cl Br I B Si Se".split())
HEAVY_MIN, HEAVY_MAX = 10, 70
MAX_REPLICATE_SPREAD = 1.0
MIN_GROUP = 10
# a patent contributes to a task only with at least this many distinct molecules in that task;
# counted per task so the rule holds in every task file on its own
MIN_PATENT_MOLECULES = 20
MIN_RANGE = 1.0
HIGH_CONTRAST_RANGE = 2.0
# A group whose median molecule has fewer heavy atoms than this is a substituent table, not a
# compound series: "Example | R | IC50" layouts draw only the R group, and the pipeline stored the
# R group as the compound. The manual validation (sar-validation-20261019) separated the two at
# group-median 10 vs 31 heavy atoms; at 16 the flag caught 36 of 37 annotated fragment errors
# with no false positives. Applied per group, because every molecule in such a group is a
# fragment and no within-group comparison can see it.
FRAGMENT_GROUP_MEDIAN_HEAVY = 16
FAMILY_OVERLAP = 0.5
# the time split follows the earliest priority date of a split unit: chemistry filed from this
# year on is the test set. The grant-year split of the earlier release is kept as split_time_grant.
TEST_FROM_PRIORITY_YEAR = 2016
# A task whose patents are mostly recent would have a test set larger than its training set under
# the common cutoff (rank_complex_mutant: 24 train groups against 40 test at 2016). It keeps a later
# cutoff of its own; split_time and test_tier in that task file carry it. sar_all keeps the common one.
TASK_TEST_FROM = {"rank_complex_mutant": 2018}
# The test set is also limited to patents granted from this year on. A unit that is test by priority
# year but holds a patent granted earlier leaves the test set and joins train/valid, so that every
# test patent is a 2020-2025 grant and no test unit shares a patent or family with train/valid.
TEST_MIN_GRANT_YEAR = 2020
TEST_FROM_YEAR = 2021
VALID_FRACTION = 0.15
RANDOM_FRACTIONS = (0.70, 0.15, 0.15)
SCAFFOLD_FRACTIONS = (0.70, 0.15, 0.15)
# two groups of one patent and endpoint are the same table read twice when they share this much of
# their molecules and agree to within DUPLICATE_MAX_DIFF log on the shared ones
DUPLICATE_OVERLAP = 0.8
DUPLICATE_MAX_DIFF = 0.05
# a group whose values sit this far from BindingDB/ChEMBL on >= PUBLIC_MIN_MATCHES shared
# measurements was read on the wrong scale
PUBLIC_MAX_OFFSET = 1.0
PUBLIC_MIN_MATCHES = 3
# flags only: values outside a usual assay window, and molecules beyond rule-of-5 size
PACT_FLAG_HIGH, PACT_FLAG_LOW = 11.0, 3.0
BRO5_HEAVY = 50
cfg = {k: v for k, v in sorted(globals().items())
       if k.isupper() and k not in {"R", "FULL", "M", "OUT", "AFFINITY", "KEY", "ORGANIC"}}
log = []


def note(msg):
    print(msg, flush=True)
    log.append(msg)


# ----------------------------------------------------------------- full SAR
REPAIR_COLS = ["structure_repair", "SMILES_before_repair", "InChIKey_before_repair"]   # written by make_patched_sar.py (r6)
_names = ds.dataset(R + "SAR_merged.parquet").schema.names
sar = ds.dataset(R + "SAR_merged.parquet").to_table(
    columns=["source", "source_record_id", "patent_core", "Patent", "Target", "Activity Type",
             "Relation", "Value", "Unit", "value_nM", "SMILES", "InChIKey"]
    + [c for c in REPAIR_COLS if c in _names]).to_pandas()
for _c in REPAIR_COLS:
    sar[_c] = sar[_c].fillna("").astype(str) if _c in sar else ""
sar["Target"] = sar.Target.fillna("").str.strip()
note(f"formal SAR rows: {len(sar):,}  endpoints: {sar['Activity Type'].nunique()}  "
     f"patents: {sar.patent_core.nunique():,}  molecules: {sar.InChIKey.nunique():,}")

# ------------------------------------------------------------------ provenance
# Each run kept, per record, the cell it was read from. Carrying it lets a user open the patent at
# the right table row: compound label, table/row/column locator, column header, the raw text.
EVIDENCE = ["compound_alias", "activity_locator", "activity_header", "raw_value", "structure_file"]
ev = pd.concat([
    ds.dataset(R + "new_SAR.parquet").to_table(columns=["source_record_id"] + EVIDENCE).to_pandas(),
    ds.dataset(FULL + "SAR_with_evidence.parquet").to_table(
        columns=["source_record_id"] + EVIDENCE).to_pandas(),
]).drop_duplicates("source_record_id")
sar = sar.merge(ev, on="source_record_id", how="left")
for c in EVIDENCE:
    sar[c] = sar[c].fillna("").astype(str)
has_ev = sar.activity_locator.ne("") | sar.activity_header.ne("")
note(f"provenance found for {int(has_ev.sum()):,} of {len(sar):,} rows ({has_ev.mean():.2%}); "
     f"by run: {sar.assign(e=has_ev).groupby('source').e.mean().round(4).to_dict()}")

# rows the upstream run derived from an A/B/C legend rather than a printed number
mk = ds.dataset(R + "new_SAR.parquet").to_table(
    columns=["source_record_id", "measurement_kind"]).to_pandas()
graded = set(mk.loc[mk.measurement_kind.isin(["legend_threshold", "legend_interval"]),
                    "source_record_id"])
sar["is_categorical_derived"] = sar.source_record_id.isin(graded)
sar["is_affinity"] = sar["Activity Type"].isin(AFFINITY)
sar["is_exact"] = sar.Relation.eq("=") & sar.value_nM.gt(0) & ~sar.is_categorical_derived
sar["pAct"] = np.where(sar.is_exact, 9 - np.log10(sar.value_nM.where(sar.value_nM > 0)), np.nan)
sar.loc[~sar.pAct.between(2, 13), "pAct"] = np.nan
sar["is_exact"] = sar.is_exact & sar.pAct.notna()
note(f"is_affinity {int(sar.is_affinity.sum()):,} | is_exact {int(sar.is_exact.sum()):,} | "
     f"is_categorical_derived {int(sar.is_categorical_derived.sum()):,}")

# --------------------------------------------- molecules: parent form, size, scaffold
# standardize.py keyed every drawn molecule to its parent (largest organic fragment, neutralised).
# From here on InChIKey and SMILES mean the parent; the drawn form stays in *_as_drawn and is what
# the alias scan and the public-coverage cache were keyed on.
props = pq.read_table(OUT + "/work/mol_std.parquet").to_pandas().set_index("InChIKey")
missing = set(sar.InChIKey.dropna()) - set(props.index)
if missing:
    raise SystemExit(f"mol_std.parquet lacks {len(missing):,} molecules; rerun scripts/standardize.py")
sar["InChIKey_as_drawn"] = sar.InChIKey
sar["SMILES_as_drawn"] = sar.SMILES
sar["InChIKey"] = sar.InChIKey_as_drawn.map(props.parent_inchikey)
sar["SMILES"] = sar.InChIKey_as_drawn.map(props.parent_smiles)
for c in ["heavy_atoms", "scaffold", "n_stereocentres", "has_undefined_stereo",
          "is_isotope_labelled", "covalent_warhead"]:
    sar[c] = sar.InChIKey_as_drawn.map(props[c])
sar["is_druglike"] = (sar.heavy_atoms.between(HEAVY_MIN, HEAVY_MAX)
                      & ~sar.InChIKey_as_drawn.map(props.unusual_elements).fillna(True))
sar["is_bro5_size"] = sar.heavy_atoms > BRO5_HEAVY
n_salt = int(sar.InChIKey.ne(sar.InChIKey_as_drawn).sum())
note(f"parent standardisation: {sar.InChIKey_as_drawn.nunique():,} drawn molecules -> "
     f"{sar.InChIKey.nunique():,} parents; {n_salt:,} rows carry a salt, charge or extra component")
note(f"is_druglike {int(sar.is_druglike.sum()):,} rows; "
     f"{int((sar.scaffold == '').sum()):,} rows have no ring scaffold")

# ---------------------------------------------------------------- patent year
pt = ds.dataset(R + "patents.parquet").to_table(columns=["Patent", "archive"]).to_pandas()
pt["patent_core"] = pt.Patent.str.replace(r"(B\d|A\d|E|S|H|P\d)$", "", regex=True)
pt["year"] = pt.archive.str.extract(r"/data/(\d{4})/")[0].astype(int)
sar["year"] = sar.patent_core.map(pt.drop_duplicates("patent_core").set_index("patent_core").year)
note(f"patent year missing for {int(sar.year.isna().sum()):,} rows")

# priority date: when the chemistry was filed. The earliest claimed priority, else the PCT filing,
# else the US filing. A patent with none of the three keeps an upper bound - its grant year.
bib = ds.dataset(FULL + "patents.parquet").to_table(
    columns=["patent_core", "priority_date", "pct_filing_date", "application_date"]).to_pandas()
bib["priority"] = pd.NaT
bib["priority_source"] = ""
for col, label in [("priority_date", "priority_claim"), ("pct_filing_date", "pct_filing"),
                   ("application_date", "us_filing")]:
    d = pd.to_datetime(bib[col], errors="coerce")
    take = bib.priority.isna() & d.notna()
    bib.loc[take, "priority"] = d[take]
    bib.loc[take, "priority_source"] = label
bib = bib.dropna(subset=["priority"]).sort_values("priority").drop_duplicates("patent_core")
sar["priority_date"] = sar.patent_core.map(bib.set_index("patent_core").priority.dt.strftime("%Y-%m-%d")).fillna("")
sar["priority_source"] = sar.patent_core.map(bib.set_index("patent_core").priority_source).fillna("none")
sar["priority_year"] = sar.patent_core.map(bib.set_index("patent_core").priority.dt.year)
pp = sar.drop_duplicates("patent_core")
note(f"priority date known for {int(pp.priority_year.notna().sum()):,} of {len(pp):,} patents; "
     f"source {pp.priority_source.value_counts().to_dict()}")

# first public date: the earliest of the US application publication, the PCT publication and the grant.
# It is when the chemistry became readable by anyone. A patent with none of the three has no value
# (empty), never an estimate.
pub = ds.dataset(FULL + "patents.parquet").to_table(
    columns=["patent_core", "publication_date", "application_publication_date",
             "wo_publication_date"]).to_pandas()
_pub_cols = {"publication_date": "grant", "application_publication_date": "us_application",
             "wo_publication_date": "wo"}
for _c in _pub_cols:
    _raw = pub[_c].fillna("").astype(str).str.strip().str[:10]
    pub[_c + "_dt"] = pd.to_datetime(_raw.where(_raw.ne(""), None), errors="coerce")
    note(f"{_c}: {int(_raw.ne('').sum()):,} archives with a value, "
         f"{int((_raw.ne('') & pub[_c + '_dt'].isna()).sum()):,} unparseable")
_order = ["application_publication_date", "wo_publication_date", "publication_date"]
pub = pub.groupby("patent_core")[[c + "_dt" for c in _order]].min()
pub["first_public"] = pub.min(axis=1)
pub["first_public_source"] = ""
for _c in reversed(_order):
    pub.loc[pub[_c + "_dt"].eq(pub.first_public), "first_public_source"] = _pub_cols[_c]
sar["first_public_date"] = sar.patent_core.map(pub.first_public.dt.strftime("%Y-%m-%d")).fillna("")
sar["first_public_year"] = sar.patent_core.map(pub.first_public.dt.year).astype("Int64")
sar["first_public_source"] = sar.patent_core.map(pub.first_public_source).fillna("")
sar.loc[sar.first_public_date.eq(""), "first_public_source"] = "none"
pp = sar.drop_duplicates("patent_core")
note(f"first public date known for {int(pp.first_public_date.ne('').sum()):,} of {len(pp):,} patents; "
     f"source {pp.first_public_source.value_counts().to_dict()}")

# ------------------------------------------------------------- target mapping
tmap = pq.read_table(M + "target_map.parquet").to_pandas()
tmap = tmap[tmap.uniprot_accession.fillna("").ne("")]
pair_map = tmap.drop_duplicates(["patent_core", "Target"]).set_index(["patent_core", "Target"])
sar = sar.join(pair_map[["uniprot_accession", "mapped_identity", "target_class",
                         "measured_in_cells", "confidence", "mappability", "source_run"]],
               on=["patent_core", "Target"])
sar["target_mapped"] = (sar.uniprot_accession.fillna("").ne("")
                        & sar.confidence.isin(["high", "medium"])
                        & sar.target_class.eq("protein"))
sar["is_biochemical"] = sar.measured_in_cells.eq("N")
sar["mapping_run"] = sar.source_run
note("mappability of the (patent, Target) pairs behind these rows: "
     + str(sar.drop_duplicates(["patent_core", "Target"]).mappability.value_counts(dropna=False).to_dict()))
all_pairs = sar.drop_duplicates(["patent_core", "Target"])
note(f"target mapping covers {len(pair_map):,} of {len(all_pairs):,} distinct "
     f"(patent, Target) pairs = {len(pair_map) / len(all_pairs):.1%}; "
     f"target_mapped is True on {int(sar.target_mapped.sum()):,} rows "
     f"({sar.target_mapped.mean():.1%})")
note("target_mapped requires a single wild-type protein at high or medium confidence. Every "
     "pair has now been attempted, so a False value means the mapping genuinely failed or the "
     "assay has no single protein target - not that it was skipped. See mappability and "
     "target_class for which of the two.")
organism = tmap.drop_duplicates("uniprot_accession").set_index("uniprot_accession").organism
sar["target_organism"] = sar.uniprot_accession.map(organism).fillna("")
sar["is_human_target"] = sar.target_organism.str.contains("Homo sapiens", regex=False)
# drug-metabolism and safety panels (CYP inhibition, hERG, Nav1.5, Cav1.2, efflux and uptake
# transporters): measured in many patents, but on an anti-target, not the series' own target.
# MAO, CYP17 and CYP11B stay unflagged: they are drug targets in their own right.
ADME_SAFETY = r"^(CP3A4|CP3A5|CP2C9|CP2CJ|CP2D6|CP1A2|CP2B6|CP2C8|KCNH2|SCN5A|CAC1C|MDR1|ABCG2|SO1B1|SO1B3)_"
entry = tmap.drop_duplicates("uniprot_accession").set_index("uniprot_accession").uniprot_entry.fillna("")
sar["is_adme_safety_target"] = sar.uniprot_accession.map(entry).fillna("").str.match(ADME_SAFETY)

# ------------------------------------------------------------ assay labels
# The label is whatever text sat in the column header next to the endpoint, so the same assay can
# arrive as "BTK" and "Btk", or with stray separators ("SHP2 ,"). Grouping on a case- and
# punctuation-blind form joins those; Unicode letters stay, so PI3Kα and PI3Kβ remain apart.
sar["Target_norm"] = sar.Target.map(lambda t: re.sub(r"[\W_]+", "", t.casefold()))
sar["_uniprot_key"] = sar.uniprot_accession.fillna("")


def clean_label(t):
    t = re.sub(r"\s+", " ", t)
    t = re.sub(r"(?:\s*[;,]\s*)+(?=[;,]|$)", "", t)            # empty list items: "a; 200 ; ;"
    return re.sub(r"^[\s,;:*=\-_/|.]+|[\s,;:*=\-_/|.]+$", "", t)


sar["target_label"] = sar.Target.map(clean_label)

# A label that is only a number ("9.7", "<1", "±1") is a value from a neighbouring column that the
# extraction took for the assay label; every distinct value then became a group of its own (one
# patent was cut into 109 groups of 1-3 rows). The activity header names the assay, so those rows
# take their label from it and are regrouped. Labels that are names, such as cell lines, are real
# assays and are left alone. Target keeps the original text.
NUMERIC_LABEL = re.compile(r"[\s\d.,+\-±<>=%/]+")
HEADER_ENDPOINT = re.compile(r"\s*[(\[]?\s*\b(?:IC50|EC50|Ki|Kd|GI50|CC50|AC50|MIC)\b.*$", re.I)


def label_from_header(h):
    h = clean_label(h)
    stem = clean_label(HEADER_ENDPOINT.sub("", h))
    return stem if len(re.sub(r"[\W_]", "", stem)) >= 2 else h


_num = sar.target_label.map(lambda t: bool(NUMERIC_LABEL.fullmatch(t))) & sar.activity_header.ne("")
sar["label_repaired"] = _num
sar.loc[_num, "target_label"] = sar.loc[_num, "activity_header"].map(label_from_header)
sar.loc[_num, "Target_norm"] = sar.loc[_num, "target_label"].map(
    lambda t: re.sub(r"[\W_]+", "", t.casefold()))
note(f"label_repaired: {int(_num.sum()):,} rows whose label was a bare number now take the label "
     f"from their activity header ({sar.loc[_num, 'Target'].nunique():,} distinct numbers)")

# --------------------------------------------- one value per (group, molecule)
n_raw_groups = sar.groupby(["patent_core", "Target", "Activity Type"], sort=False).ngroup().nunique()
sar["group_id"] = sar.groupby(KEY, sort=False).ngroup()
note(f"groups: {n_raw_groups:,} on the raw label -> {sar.group_id.nunique():,} after joining labels "
     f"that differ only by case or punctuation (and never across two mapped proteins)")
glabel = sar.groupby("group_id").target_label.agg(lambda s: s.value_counts().index[0])
sar["target_label"] = sar.group_id.map(glabel)
# a bare-number label with no header to repair it from stays uninformative
sar["group_label_uninformative"] = sar.target_label.map(
    lambda t: len(re.sub(r"[\W_]", "", t)) < 2 or bool(NUMERIC_LABEL.fullmatch(t)))
note(f"group_label_uninformative: {int(sar.loc[sar.group_label_uninformative, 'group_id'].nunique()):,} "
     f"groups whose label holds fewer than two letters or digits, or is only a number")
_rep = sar[sar.label_repaired].groupby(["patent_core", "target_label"]).agg(
    activity_header=("activity_header", "first"), n_groups_before=("Target", "nunique"),
    n_rows=("InChIKey", "size"), n_molecules=("InChIKey", "nunique"),
    activity_type=("Activity Type", lambda s: "|".join(sorted(set(s))))).reset_index()
_rep.to_csv(OUT + "/label_repaired.tsv", sep="\t", index=False)
note(f"label_repaired: {len(_rep):,} (patent, label) combinations; the groups they now form are "
     f"unmapped (the target map never saw these labels), see label_repaired.tsv")
ex = sar[sar.is_exact]
spread = ex.groupby(["group_id", "InChIKey"]).pAct.agg(lambda s: s.max() - s.min())
noise = float(spread[spread > 0].median())
note(f"ASSAY NOISE FLOOR: median spread between repeat measurements of one molecule in one group "
     f"is {noise:.3f} log units")
wide = set(spread[spread > MAX_REPLICATE_SPREAD].index)
sar["replicate_disagrees"] = [
    (g, k) in wide for g, k in zip(sar.group_id, sar.InChIKey)]
note(f"replicate_disagrees on {int(sar.replicate_disagrees.sum()):,} rows "
     f"(spread over {MAX_REPLICATE_SPREAD} log)")

# ------------------------------------------------------------- group statistics
usable = sar[sar.is_exact & sar.is_druglike & ~sar.replicate_disagrees]
gstat = usable.groupby("group_id").agg(
    group_n_molecules=("InChIKey", "nunique"),
    group_pAct_min=("pAct", "min"), group_pAct_max=("pAct", "max"),
    group_pAct_std=("pAct", "std"))
gstat["group_pAct_range"] = gstat.group_pAct_max - gstat.group_pAct_min
gstat["group_n_distinct"] = usable.groupby("group_id").pAct.agg(lambda s: s.round(2).nunique())
sar = sar.join(gstat, on="group_id")
sar["group_is_gradelike"] = ((sar.group_n_distinct < 4)
                             | (sar.group_n_distinct / sar.group_n_molecules < 0.3)).fillna(False)
note(f"groups with >= {MIN_GROUP} usable molecules: "
     f"{int((gstat.group_n_molecules >= MIN_GROUP).sum()):,} of {len(gstat):,}")


# ------------------------------------------------------------------- families
def family_ids(frame):
    """Union groups measuring the same endpoint on a largely shared molecule set."""
    g = frame.groupby(["group_id", "Activity Type"]).InChIKey.agg(frozenset).reset_index(name="mols")
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
        if c / min(size[i], size[j]) >= FAMILY_OVERLAP:
            ri, rj = find(i), find(j)
            if ri != rj:
                parent[ri] = rj
    g["family_id"] = [find(i) for i in range(len(g))]
    # the union-find root depends on the order molecules come out of the frozensets, which changes
    # with the hash seed; the partition does not. Name each family by its smallest group id so two
    # builds of the same data give the same family_id and series_id.
    g["family_id"] = g.groupby("family_id").group_id.transform("min")
    return g[["group_id", "family_id", "n"]]


# ---------------------------------------------------- substituent-table groups
# median over every row of the group, exactly as the validation measured it
group_median_heavy = sar.groupby("group_id").heavy_atoms.median()
sar["group_median_heavy_atoms"] = sar.group_id.map(group_median_heavy)
sar["group_is_fragment_suspect"] = sar.group_median_heavy_atoms < FRAGMENT_GROUP_MEDIAN_HEAVY
n_frag = int((group_median_heavy < FRAGMENT_GROUP_MEDIAN_HEAVY).sum())
note(f"group_is_fragment_suspect: {n_frag:,} groups, "
     f"{int(sar.group_is_fragment_suspect.sum()):,} rows "
     f"({sar.group_is_fragment_suspect.mean():.2%}) with group-median heavy atoms "
     f"< {FRAGMENT_GROUP_MEDIAN_HEAVY}")

# ------------------------------------------------------- units and table layout
# A unit is read once per column, so when it is wrong the whole table shifts by 3 or 6 log and
# the within-group ranking survives while every absolute value is wrong. The audit against
# BindingDB/ChEMBL traced the shifted tables to four layouts, each checked here per group:
#   mixed_units          one assay column read in two units
#   molar_scale          millimolar, or molar values of 1 mM and above: a lost micro sign
#                        ("IC50, -M") or a typed mM where the patent means uM
#   header_unit_conflict the header names a unit other than the one applied, or names only an
#                        assay concentration ("MEK1 1 mM ATP IC50") and that was applied
#   percent_in_column    the column also holds % inhibition ("Ki [nM] or % inhibition [1 uM]")
#   public_offset        >= 3 shared measurements sit > 1 log from BindingDB/ChEMBL on median
#   tables_disagree      two source tables joined into one group disagree: shared compounds
#                        differ by > 1 log on median, or the tables are in different units and
#                        their medians sit > 2.5 log apart (one of them carries a wrong unit;
#                        BindingDB often copied the same reading, so the public check misses it)
UNIT_TOKEN = re.compile(r"(?<![A-Za-z])(picomolar|nanomolar|micromolar|millimolar|molar|"
                        r"pM|nM|nm|[uμµ]M|[uμµ]m|mM|M)(?![A-Za-z])")
UNIT_NAME = {"picomolar": "pM", "nanomolar": "nM", "micromolar": "uM", "millimolar": "mM",
             "molar": "M", "pM": "pM", "nM": "nM", "nm": "nM", "mM": "mM", "M": "M"}


def norm_unit(u):
    u = str(u).strip().replace("μ", "u").replace("µ", "u")
    return UNIT_NAME.get(u, "uM" if u.lower() == "um" else u)


def header_units(h):
    value, conc = set(), set()
    h = re.sub(r"(?i)\b(IC|EC|GI|ED|LD)\s+50\b", r"\g<1>50", h)      # "IC 50 uM" is IC50 in uM
    for m in UNIT_TOKEN.finditer(h):
        u = norm_unit(m.group(1))
        # a concentration is a number standing alone before the unit ("1 mM ATP", "[1 uM]"); the
        # digits ending a word ("IC50 uM", "NaV1.7") are not
        standalone_number = re.search(r"(?<![\w.])\d+(?:\.\d+)?\s*$", h[:m.start()])
        (conc if standalone_number else value).add(u)
    return value, conc


exact_rows = sar[sar.is_exact].copy()
exact_rows["unit_n"] = exact_rows.Unit.map(norm_unit)
hdr = {h: header_units(h) for h in exact_rows.activity_header.unique() if h}


def unit_conflict(h, u):
    if not h:
        return False
    value, conc = hdr[h]
    return (u not in value) if value else (u in conc)


exact_rows["unit_conflict"] = [unit_conflict(h, u) for h, u in
                               zip(exact_rows.activity_header, exact_rows.unit_n)]
exact_rows["molar_scale"] = exact_rows.unit_n.eq("mM") | (
    exact_rows.unit_n.eq("M") & exact_rows.value_nM.ge(1e6))
exact_rows["percent_col"] = exact_rows.activity_header.str.contains(
    r"%|(?i:inhibition\s+at)|@", regex=True)
pub = pq.read_table(OUT + "/work/public_pairs.parquet").to_pandas()
pub = pub[pub.source.isin(PUBLIC_SOURCES)]
pub_val = pub.groupby(["ik", "uniprot", "Activity Type"]).pAct.median()
look = list(zip(exact_rows.InChIKey_as_drawn, exact_rows.uniprot_accession.fillna(""),
                exact_rows["Activity Type"]))
look_parent = list(zip(exact_rows.InChIKey, exact_rows.uniprot_accession.fillna(""),
                       exact_rows["Activity Type"]))
pv = pub_val.to_dict()
exact_rows["public_pAct"] = [pv.get(a, pv.get(b, np.nan)) for a, b in zip(look, look_parent)]
# mixed units are judged per header: joined label variants may legitimately be two tables in
# two units, but one column read in two units is a parsing error
col_units = exact_rows.groupby(["group_id", "activity_header"]).unit_n.nunique().groupby("group_id").max()
gu = exact_rows.groupby("group_id").agg(
    molar_scale=("molar_scale", "any"),
    header_unit_conflict=("unit_conflict", "any"), percent_in_column=("percent_col", "any"))
off = (exact_rows.dropna(subset=["public_pAct"])
       .assign(d=lambda x: x.pAct - x.public_pAct).groupby("group_id").d.agg(["size", "median"]))
gu["public_offset"] = gu.index.map(
    lambda g: bool(g in off.index and off.at[g, "size"] >= PUBLIC_MIN_MATCHES
                   and abs(off.at[g, "median"]) > PUBLIC_MAX_OFFSET))
gu["mixed_units"] = gu.index.map(col_units).fillna(1) > 1
exact_rows["table"] = exact_rows.activity_locator.str.split("/").str[0]
multi_t = exact_rows[exact_rows.groupby("group_id").table.transform("nunique") > 1]
disagree = set()
for gid, x in multi_t.groupby("group_id"):
    per = x.groupby(["table", "InChIKey"]).pAct.median().unstack(0)
    tabs = list(per.columns)
    tmed = x.groupby("table").pAct.median()
    tunit = x.groupby("table").unit_n.agg(lambda v: v.mode().iloc[0])
    for i in range(len(tabs)):
        for j in range(i + 1, len(tabs)):
            a, b = tabs[i], tabs[j]
            both = per[[a, b]].dropna()
            if len(both) >= 2 and abs((both[a] - both[b]).median()) > 1.0:
                disagree.add(gid)
            if tunit[a] != tunit[b] and abs(tmed[a] - tmed[b]) > 2.5:
                disagree.add(gid)
gu["tables_disagree"] = gu.index.isin(disagree)
UNIT_CHECKS = ["mixed_units", "molar_scale", "header_unit_conflict", "percent_in_column",
               "public_offset", "tables_disagree"]
gu["issue"] = gu[UNIT_CHECKS].apply(lambda r: ";".join(c for c in UNIT_CHECKS if r[c]), axis=1)
sar["group_unit_issue"] = sar.group_id.map(gu.issue).fillna("")
sar["group_unit_suspect"] = sar.group_unit_issue.ne("")
note("group_unit_suspect: " + ", ".join(f"{c} {int(gu[c].sum()):,}" for c in UNIT_CHECKS)
     + f" -> {int(gu.issue.ne('').sum()):,} groups, {int(sar.group_unit_suspect.sum()):,} rows")
note(f"public offset check: {int((off['size'] >= PUBLIC_MIN_MATCHES).sum()):,} groups have >= "
     f"{PUBLIC_MIN_MATCHES} BindingDB/ChEMBL matches")

# ------------------------------------------------- one table read twice in a patent
# Two groups of one patent and endpoint that share most molecules and agree to the second decimal
# are one column read under two labels; scoring both counts that series twice. The group kept is
# the mapped one, then the larger, then the first.
mv = usable.groupby(["group_id", "InChIKey"]).pAct.median()
gm_pe = usable.drop_duplicates("group_id").set_index("group_id")[["patent_core", "Activity Type"]]
mapped_g = sar.drop_duplicates("group_id").set_index("group_id").target_mapped
dup_of = {}
for _, gids in gm_pe.groupby(["patent_core", "Activity Type"]).groups.items():
    if len(gids) < 2:
        continue
    vals = {g: mv.loc[g] for g in gids}
    order = sorted(gids, key=lambda g: (not mapped_g.get(g, False), -len(vals[g]), g))
    for i, a in enumerate(order):
        if a in dup_of:
            continue
        for b in order[i + 1:]:
            if b in dup_of:
                continue
            A, Bv = vals[a], vals[b]
            common = A.index.intersection(Bv.index)
            if len(common) < 2 or len(common) / min(len(A), len(Bv)) < DUPLICATE_OVERLAP:
                continue
            if ((A[common] - Bv[common]).abs().median() < DUPLICATE_MAX_DIFF
                    and Bv[common].max() - Bv[common].min() >= MIN_RANGE):
                dup_of[b] = a
sar["group_duplicate_of"] = sar.group_id.map(dup_of).fillna(-1).astype(int)
sar["group_is_duplicate"] = sar.group_duplicate_of.ge(0)
note(f"group_is_duplicate: {len(dup_of):,} groups repeat another group of the same patent and "
     f"endpoint, {int(sar.group_is_duplicate.sum()):,} rows")

# ------------------------------------- validation sample: disputed and confirmed-wrong
# A record the two annotators disagreed on stays out of every task until a person adjudicates
# it; a record both annotators called wrong stays out for good. Both are held out per
# (patent, Target) pair, because the target mapping and the table layout behind a record are
# shared by every record of that pair: dropping the one sampled row would keep its neighbours,
# which carry the same doubt.
QS = ["q1", "q2", "q3"]
val = pd.read_csv(VAL + "/merged_annotations.tsv", sep="\t", dtype=str).fillna("")
val = val[["item_id"] + [f"{q}_consensus" for q in QS]].merge(
    pq.read_table(VAL + "/sample.parquet", columns=["item_id", "patent_core", "Target"]).to_pandas(),
    on="item_id")
val["disputed"] = val[[f"{q}_consensus" for q in QS]].eq("DISAGREE").any(axis=1)
val["wrong"] = val[[f"{q}_consensus" for q in QS]].eq("incorrect").any(axis=1)
# a person's adjudication, once given, settles a disputed item: every open question judged
# correct (or not applicable) puts it back, any judged incorrect makes it a confirmed error, and
# an open question still unjudged or judged uncertain keeps it held out
adj_file = VAL + "/adjudication_results.tsv"
n_adjudicated = 0
if os.path.exists(adj_file):
    adj = pd.read_csv(adj_file, sep="\t", dtype=str).fillna("")
    adj = adj[adj.adjudicated.ne("")]
    for item, a in adj.groupby("item_id"):
        row = val.index[val.item_id == item]
        if not len(row) or not val.loc[row[0], "disputed"]:
            continue
        open_qs = [q for q in QS if val.loc[row[0], f"{q}_consensus"] == "DISAGREE"]
        got = dict(zip(a.question, a.adjudicated))
        if not all(q in got for q in open_qs):
            continue
        verdicts = [got[q] for q in open_qs]
        if "incorrect" in verdicts:
            val.loc[row[0], ["disputed", "wrong"]] = [False, True]
            n_adjudicated += 1
        elif all(v in ("correct", "not_applicable") for v in verdicts):
            val.loc[row[0], "disputed"] = False
            n_adjudicated += 1
pair_key = list(zip(sar.patent_core, sar.Target))
disputed_pairs = set(zip(val.loc[val.disputed, "patent_core"], val.loc[val.disputed, "Target"]))
wrong_pairs = set(zip(val.loc[val.wrong, "patent_core"], val.loc[val.wrong, "Target"]))
sar["pair_under_adjudication"] = [k in disputed_pairs for k in pair_key]
sar["pair_validated_error"] = [k in wrong_pairs for k in pair_key]
# A patent where the validation found one doubtful or wrong record is held out whole: one bad
# extraction says the patent's tables were read in a way that can go wrong, and the other
# assays of the same patent were read the same way.
doubt_patents = ({p for p, _ in disputed_pairs} | {p for p, _ in wrong_pairs})
sar["patent_has_validation_issue"] = sar.patent_core.isin(doubt_patents)
note(f"patents held out whole for a disputed or wrong record: {len(doubt_patents):,} patents, "
     f"{int(sar.patent_has_validation_issue.sum()):,} rows")

# ------------------------------------------------ structures resting on an abbreviation
# ChemDraw writes an abbreviated group ("F3C", "MeO", "Ph") as one atom with an alias; RDKit
# keeps the atom and drops the group, so the stored structure lacks atoms. alias_scan.py checked
# every drawing of every stored molecule: a structure is suspect when every drawing of it in its
# patent carries such an alias (one clean drawing of the same molecule clears it).
alias = pq.read_table(OUT + "/work/alias_scan.parquet").to_pandas()
alias = alias[alias.error.eq("") & alias.n_bad.gt(0) & alias.n_bad.eq(alias.n_mol)]
suspect_pairs = set(zip(alias.Patent, alias.InChIKey))
sar["structure_alias_suspect"] = [(p, k) in suspect_pairs
                                  for p, k in zip(sar.Patent, sar.InChIKey_as_drawn)]
note(f"structure_alias_suspect: {len(suspect_pairs):,} (patent, molecule) pairs, "
     f"{int(sar.structure_alias_suspect.sum()):,} rows")
note(f"validation holdout: {int(val.disputed.sum())} disputed items -> {len(disputed_pairs)} pairs, "
     f"{int(sar.pair_under_adjudication.sum()):,} rows; {int(val.wrong.sum())} confirmed-wrong items "
     f"-> {len(wrong_pairs)} pairs, {int(sar.pair_validated_error.sum()):,} rows; "
     f"{n_adjudicated} items settled by adjudication")

fam = family_ids(usable)
fam = fam.merge(usable.groupby("group_id").year.min().reset_index(), on="group_id", how="left")
sar = sar.join(fam.set_index("group_id")[["family_id"]], on="group_id")
# the representative is the group a strict user would keep: one that passes the task filters
# first, then the earliest, then the largest. Choosing it blind to the filters would let an
# unusable group represent a family and silently delete the usable ones. A substituent-table
# group never qualifies, for the same reason.
bad_group = sar.groupby("group_id")[["group_unit_suspect", "group_is_duplicate",
                                      "group_label_uninformative"]].first().any(axis=1)
qualifies = (gstat.group_n_molecules.ge(MIN_GROUP) & gstat.group_pAct_range.ge(MIN_RANGE)
             & (group_median_heavy.reindex(gstat.index) >= FRAGMENT_GROUP_MEDIAN_HEAVY)
             & ~bad_group.reindex(gstat.index).fillna(False).astype(bool))
fam["qualifies"] = fam.group_id.map(qualifies).fillna(False)
rep = (fam.sort_values(["family_id", "qualifies", "year", "n"],
                       ascending=[True, False, True, False])
       .drop_duplicates("family_id"))
sar["is_family_representative"] = sar.group_id.isin(set(rep.group_id))
note(f"family clustering: {len(fam):,} groups -> {fam.family_id.nunique():,} families; "
     f"{int((fam.groupby('family_id').size() > 1).sum()):,} families hold more than one group")

# -------------------------------------------------------------- public overlap
# the public cache is keyed on drawn structures, so a row counts as covered when either its drawn
# form or its parent is there
pub_pairs = set(zip(pub.ik, pub.uniprot))
note(f"public comparison set: {', '.join(PUBLIC_SOURCES)} "
     f"({len(pub_pairs):,} distinct (molecule, target) pairs)")
sar["in_public_data"] = [((k, u) in pub_pairs or (d, u) in pub_pairs) if isinstance(u, str) and u
                         else False
                         for k, d, u in zip(sar.InChIKey, sar.InChIKey_as_drawn, sar.uniprot_accession)]
try:
    nov = pq.read_table(OUT + "/work/novelty_flags.parquet").to_pandas()
    mol_new = dict(zip(nov.InChIKey, nov.molecule_is_new))
    pat_new = dict(zip(nov.patent_core, nov.patent_is_new))
    sar["molecule_is_new"] = sar.InChIKey_as_drawn.map(mol_new)
    sar["patent_is_new"] = sar.patent_core.map(pat_new)
    note(f"molecule_is_new known for {int(sar.molecule_is_new.notna().sum()):,} rows; "
         f"True on {int(sar.molecule_is_new.fillna(False).sum()):,}")
except Exception as exc:                                    # novelty run not finished
    sar["molecule_is_new"] = pd.NA
    sar["patent_is_new"] = pd.NA
    note(f"novelty flags unavailable ({exc}); molecule_is_new / patent_is_new left null")
note(f"in_public_data True on {int(sar.in_public_data.sum()):,} rows "
     f"({sar.in_public_data.mean():.1%} of all formal SAR, "
     f"{sar.loc[sar.target_mapped, 'in_public_data'].mean():.1%} of mapped-target rows)")

# ------------------------------------------------------------------- splits
rng = np.random.default_rng(SEED)


def label_by_unit(units, fractions, order=None):
    """Assign train/valid/test over whole units (patents or scaffold sets)."""
    counts = units.value_counts()
    keys = list(counts.index) if order == "size" else list(rng.permutation(counts.index))
    total, acc, out = counts.sum(), 0, {}
    tr, va = fractions[0] * total, (fractions[0] + fractions[1]) * total
    for k in keys:
        out[k] = "train" if acc < tr else ("valid" if acc < va else "test")
        acc += counts[k]
    return units.map(out)


# Splits are assigned per CONNECTED COMPONENT of patents and families. Two reasons to go wider
# than the group: continuation and divisional patents reuse the same SAR table, so one family
# must stay inside one split; and one patent can hold several families (different targets or
# endpoints over the same molecules), so a patent must not straddle splits either. Linking both
# and splitting the components is the only assignment that closes both leaks.
grp_family = fam.set_index("group_id").family_id
sar["_family"] = sar.group_id.map(grp_family)

parent = {}


def find(x):
    parent.setdefault(x, x)
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[ra] = rb


for pc, fid in sar[["patent_core", "_family"]].drop_duplicates().itertuples(index=False):
    pnode = f"pat:{pc}"
    find(pnode)
    if pd.notna(fid):
        union(f"fam:{int(fid)}", pnode)
sar["_unit"] = ["cmp:" + find(f"pat:{pc}") for pc in sar.patent_core]
n_units = sar._unit.nunique()
note(f"split units: {n_units:,} connected components over "
     f"{sar.patent_core.nunique():,} patents and {fam.family_id.nunique():,} families")
unit_year = sar.groupby("_unit").year.min()


def time_split(before, after, seed):
    """before/after: boolean Series over units; the rest stays unassigned. 15% of `before` -> valid."""
    out = pd.Series("", index=before.index, dtype=object)
    out[before] = "train"
    out[after] = "test"
    early = sorted(out[out == "train"].index)
    pick = np.random.default_rng([SEED, seed]).choice(
        early, size=int(round(len(early) * VALID_FRACTION)), replace=False)
    out[out.index.isin(set(pick))] = "valid"
    return out


# time split: by the unit's earliest PRIORITY year, when its chemistry was filed. A patent without
# any filing date has its grant year as an upper bound. A unit is before the cutoff if any patent
# in it is (known priority or upper bound below the cutoff); it is test only if every patent in it
# has a known priority on or after the cutoff; otherwise its date is unknowable and it is left out.
pat = sar.drop_duplicates("patent_core").set_index("patent_core")
bound = pat.priority_year.fillna(pat.year)
unit_of = sar.drop_duplicates("patent_core").set_index("patent_core")._unit
u_bound = bound.groupby(unit_of).min()
u_dated = pat.priority_year.notna().groupby(unit_of).all()


def unit_time_split(cutoff):
    """Test = every patent has a known priority year >= cutoff AND the unit's earliest grant year is
    >= TEST_MIN_GRANT_YEAR. The train/valid draw is made over the units before the cutoff exactly as
    before, so their assignment does not move; units that lose the test place by grant year (their
    chemistry is not older than the cutoff) get their own 15% valid draw."""
    out = time_split(u_bound < cutoff, (u_bound >= cutoff) & u_dated, 1)
    early_grant = (out == "test") & ~(unit_year.reindex(out.index) >= TEST_MIN_GRANT_YEAR)
    moved = sorted(out[early_grant].index)
    pick = set(np.random.default_rng([SEED, 7]).choice(
        moved, size=int(round(len(moved) * VALID_FRACTION)), replace=False)) if moved else set()
    out[early_grant] = ["valid" if u in pick else "train" for u in out[early_grant].index]
    return out


unit_split = unit_time_split(TEST_FROM_PRIORITY_YEAR)
sar["split_time"] = sar._unit.map(unit_split).fillna("")
note(f"split_time (priority >= {TEST_FROM_PRIORITY_YEAR} and every patent granted >= {TEST_MIN_GRANT_YEAR} is test): "
     + str(sar.split_time.value_counts().to_dict())
     + f"; units left out for an unknowable date: {int((unit_split == '').sum()):,}")
_by_priority = ((u_bound >= TEST_FROM_PRIORITY_YEAR) & u_dated)
note(f"units test by priority year alone: {int(_by_priority.sum()):,}; moved to train/valid because a patent in "
     f"the unit was granted before {TEST_MIN_GRANT_YEAR}: {int((_by_priority & (unit_split != 'test')).sum()):,}")
sar["unit_priority_year"] = sar._unit.map(u_bound)

# the grant-year split of the 2026-10-19 release, kept for comparison
grant_split = time_split(unit_year < TEST_FROM_YEAR, unit_year >= TEST_FROM_YEAR, 2)
sar["split_time_grant"] = sar._unit.map(grant_split).fillna("")
note("split_time_grant (grant >= %d is test): %s" % (TEST_FROM_YEAR, sar.split_time_grant.value_counts().to_dict()))

# random split: by the same unit
sar["split_random"] = label_by_unit(sar._unit, RANDOM_FRACTIONS).fillna("")
note("split_random: " + str(sar.split_random.value_counts().to_dict()))

# scaffold split, molecule level - for the regression task
mol_scaf = sar.drop_duplicates("InChIKey").set_index("InChIKey").scaffold
sm = label_by_unit(mol_scaf.replace("", "<none>"), SCAFFOLD_FRACTIONS, order="size")
sar["split_scaffold_mol"] = sar.InChIKey.map(sm)
note("split_scaffold_mol: " + str(sar.split_scaffold_mol.value_counts().to_dict()))

# scaffold split, group level - for the ranking tasks, so a group never spans splits
modal = (usable.groupby(["group_id", "scaffold"]).size().reset_index(name="n")
         .sort_values(["group_id", "n"], ascending=[True, False]).drop_duplicates("group_id"))
gs = modal.set_index("group_id").scaffold.replace("", "<none>")
gsplit = label_by_unit(gs, SCAFFOLD_FRACTIONS, order="size")
sar["split_scaffold_group"] = sar.group_id.map(gsplit).fillna("")
note("split_scaffold_group: " + str(sar.split_scaffold_group.value_counts().to_dict()))

# Novel tier on the time-split test rows. A group is scored whenever it has an assigned protein
# at high or medium confidence, whatever its class: restricting this to target_mapped (single
# wild-type protein) left every complex and mutant group unscored, and the old code then labelled
# all of them "reference". A test group that cannot be scored is now "unscored", not "reference".
has_target = (sar.uniprot_accession.fillna("").ne("") & sar.confidence.isin(["high", "medium"]))
test_groups = sar[(sar.split_time == "test") & has_target & sar.is_exact]
share = test_groups.groupby("group_id").in_public_data.mean()
novel = set(share[share == 0].index)
scored = set(share.index)
sar["test_tier"] = np.select(
    [sar.split_time != "test", sar.group_id.isin(novel), sar.group_id.isin(scored)],
    ["", "novel", "reference"], default="unscored")
note(f"test_tier: {len(novel):,} novel groups, {len(scored) - len(novel):,} reference, "
     f"{int(sar.loc[sar.test_tier == 'unscored', 'group_id'].nunique()):,} unscored")

# --------------------------------------------------------------- write records
sar["value_beyond_assay_range"] = sar.pAct.gt(PACT_FLAG_HIGH) | sar.pAct.lt(PACT_FLAG_LOW)
cols = ["source", "source_record_id", "Patent", "patent_core", "year", "priority_date",
        "priority_source", "first_public_date", "first_public_year", "first_public_source", "Target", "target_label", "Activity Type",
        "Relation", "Value", "Unit", "value_nM", "pAct", "SMILES", "InChIKey", "SMILES_as_drawn",
        "InChIKey_as_drawn", "scaffold", "heavy_atoms",
        "compound_alias", "activity_locator", "activity_header", "raw_value", "structure_file",
        "group_id", "family_id", "uniprot_accession", "mapped_identity",
        "target_class", "target_organism", "is_human_target", "confidence", "measured_in_cells",
        "mappability", "mapping_run",
        "is_affinity", "is_exact", "is_categorical_derived", "is_druglike", "target_mapped",
        "is_biochemical", "replicate_disagrees", "is_family_representative",
        "group_n_molecules", "group_pAct_min", "group_pAct_max", "group_pAct_range",
        "group_pAct_std", "group_n_distinct", "group_is_gradelike",
        "group_median_heavy_atoms", "group_is_fragment_suspect",
        "group_unit_suspect", "group_unit_issue", "group_is_duplicate", "group_duplicate_of",
        "group_label_uninformative", "label_repaired",
        "pair_under_adjudication", "pair_validated_error",
        "patent_has_validation_issue", "structure_alias_suspect",
        "n_stereocentres", "has_undefined_stereo", "is_isotope_labelled", "covalent_warhead",
        "is_bro5_size", "value_beyond_assay_range", "is_adme_safety_target",
        "in_public_data", "molecule_is_new", "patent_is_new",
        "split_time", "split_time_grant", "split_random", "split_scaffold_mol",
        "split_scaffold_group", "test_tier", "structure_repair", "SMILES_before_repair", "InChIKey_before_repair"]
sar[cols].to_parquet(OUT + "/sar_all.parquet", index=False)

seq = tmap.drop_duplicates("uniprot_accession")[
    ["uniprot_accession", "uniprot_entry", "protein_name", "organism", "tax_id", "seq_len", "sequence"]]
seq.to_parquet(OUT + "/targets.parquet", index=False)

# every group is listed; those with no exact, drug-like, non-disagreeing measurement have no
# statistics (NaN) and has_usable_rows False. The earlier release dropped them silently.
groups = sar.groupby("group_id").agg(
    n_rows=("InChIKey", "size"),
    patent_core=("patent_core", "first"), year=("year", "first"),
    priority_date=("priority_date", "first"), first_public_year=("first_public_year", "first"),
    target_label=("target_label", "first"),
    endpoint=("Activity Type", "first"), uniprot_accession=("uniprot_accession", "first"),
    target_organism=("target_organism", "first"),
    is_adme_safety_target=("is_adme_safety_target", "first"),
    family_id=("family_id", "first"), n_molecules=("group_n_molecules", "first"),
    pAct_range=("group_pAct_range", "first"), pAct_std=("group_pAct_std", "first"),
    target_mapped=("target_mapped", "first"), is_biochemical=("is_biochemical", "first"),
    is_family_representative=("is_family_representative", "first"),
    median_heavy_atoms=("group_median_heavy_atoms", "first"),
    is_fragment_suspect=("group_is_fragment_suspect", "first"),
    unit_issue=("group_unit_issue", "first"), is_duplicate=("group_is_duplicate", "first"),
    duplicate_of=("group_duplicate_of", "first"),
    label_uninformative=("group_label_uninformative", "first"),
    label_repaired=("label_repaired", "max"),
    public_share=("in_public_data", "mean"), split_time=("split_time", "first"),
    split_time_grant=("split_time_grant", "first"),
    split_scaffold_group=("split_scaffold_group", "first"), test_tier=("test_tier", "first"),
).reset_index()
groups["has_usable_rows"] = groups.group_id.isin(gstat.index)
groups["n_molecules"] = groups.n_molecules.astype("Int64")
groups.to_parquet(OUT + "/groups.parquet", index=False)

# ------------------------------------------------------------------- tasks
base = (sar.is_exact & sar.is_druglike & ~sar.replicate_disagrees)
validation_out = (sar.pair_under_adjudication | sar.pair_validated_error
                  | sar.patent_has_validation_issue | sar.structure_alias_suspect)
# v1.0 industrial checks: a table on the wrong unit scale, a table read twice, a label naming nothing
industrial_out = sar.group_unit_suspect | sar.group_is_duplicate | sar.group_label_uninformative
held_out = validation_out | industrial_out
# every task also drops substituent-table groups; base stays as it was for the cost table
clean = base & ~sar.group_is_fragment_suspect & ~held_out
# no is_family_representative filter here: splits already keep each family inside one split, so
# family duplicates cannot leak. The flag stays on the rows for anyone wanting strict dedup.
rank_pre = (base & sar.target_mapped & sar.is_biochemical
            & sar.group_n_molecules.ge(MIN_GROUP) & sar.group_pAct_range.ge(MIN_RANGE)
            & ~sar.group_is_gradelike)
rank_nofrag = rank_pre & ~sar.group_is_fragment_suspect
rank_base = rank_nofrag & ~held_out

# columns constant within one (group, molecule), and columns listing every measurement behind it.
# Provenance lists keep measurement order, so the n-th locator, raw value and record id belong
# together; the other lists are de-duplicated.
FIRST = ["patent_core", "year", "priority_date", "first_public_date", "first_public_year",
         "first_public_source", "target_label", "Activity Type",
         "uniprot_accession", "mapped_identity", "target_organism", "is_human_target", "family_id",
         "is_family_representative", "SMILES", "scaffold", "heavy_atoms", "n_stereocentres",
         "has_undefined_stereo", "is_isotope_labelled", "covalent_warhead", "is_bro5_size",
         "is_adme_safety_target",
         "split_time", "split_time_grant", "split_random", "split_scaffold_group", "test_tier", "structure_repair"]
ALIGNED = ["source_record_id", "compound_alias", "activity_locator", "raw_value", "structure_file"]
DISTINCT = ["Patent", "Target", "InChIKey_as_drawn", "SMILES_as_drawn", "activity_header", "Unit"]
TASK_COLS = ["group_id", "patent_core", "Patent", "year", "priority_date", "first_public_date",
             "first_public_year", "first_public_source", "target_label", "Target",
             "Activity Type", "uniprot_accession", "mapped_identity", "target_organism",
             "is_human_target", "is_adme_safety_target", "family_id", "is_family_representative",
             "series_id", "series_n_copies",
             "InChIKey", "SMILES", "InChIKey_as_drawn", "SMILES_as_drawn", "scaffold", "heavy_atoms",
             "pAct", "value_nM", "n_measurements", "pAct_replicate_spread",
             "compound_alias", "activity_locator", "activity_header", "raw_value", "Unit",
             "structure_file", "source_record_id",
             "group_n_molecules", "group_pAct_range",
             "stereo_set_id", "has_undefined_stereo", "n_stereocentres", "is_isotope_labelled",
             "covalent_warhead", "is_bro5_size", "value_beyond_assay_range", "at_assay_ceiling",
             "at_assay_floor", "in_public_data", "molecule_is_new",
             "split_time", "split_time_grant", "split_random", "split_scaffold_group", "test_tier", "structure_repair"]
# round concentrations a series is often capped at; a value there shared by >= 3 molecules at the
# weak end of a group is the assay's top concentration written as a number
ASSAY_LIMITS_NM = (1e3, 3e3, 5e3, 1e4, 2e4, 3e4, 5e4, 1e5)
# the same artefact at the strong end: a value shared by >= 3 molecules at the strongest end of a
# group, at a round concentration, is the lower limit of what the assay resolves (a Ki written as
# 0.001 nM for every compound below it), not a measured rank
ASSAY_FLOORS_NM = (1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5, 1.0)
tasks = {}


def _distinct(s):
    return " | ".join(dict.fromkeys(str(x) for x in s if isinstance(x, str) and x))


def _aligned(s):
    return " | ".join(str(x) for x in s)


def collapse(rows):
    """One row per (group, molecule): the median of its measurements, their count and spread, and
    every source cell behind it."""
    key = ["group_id", "InChIKey"]
    rows = rows.sort_values(key + ["source_record_id"])
    rep = rows.duplicated(key, keep=False)
    one = rows[~rep].copy()
    one["n_measurements"] = 1
    one["pAct_replicate_spread"] = 0.0
    many = rows[rep]
    if len(many):
        g = many.groupby(key, sort=False)
        agg = g[FIRST].first()
        agg["pAct"] = g.pAct.median()
        agg["n_measurements"] = g.size()
        agg["pAct_replicate_spread"] = g.pAct.max() - g.pAct.min()
        agg["in_public_data"] = g.in_public_data.max()
        agg["molecule_is_new"] = g.molecule_is_new.agg(lambda s: bool(s.eq(True).all()))
        for c in ALIGNED:
            agg[c] = g[c].agg(_aligned)
        for c in DISTINCT:
            agg[c] = g[c].agg(_distinct)
        one = pd.concat([one, agg.reset_index()], ignore_index=True)
    one["value_nM"] = 10 ** (9 - one.pAct)
    one["value_beyond_assay_range"] = one.pAct.gt(PACT_FLAG_HIGH) | one.pAct.lt(PACT_FLAG_LOW)
    return one


series_removed = {}


def finalize(mask, min_range, min_patent_molecules=MIN_PATENT_MOLECULES, one_per_series=True):
    """Collapse to one value per (group, molecule), then recount each group and each patent.

    Group size and range were measured before row-level exclusions such as the alias check, so a
    group can fall under the bar once its suspect rows are gone; it then leaves the task whole.
    A patent left with fewer than `min_patent_molecules` distinct molecules in the task leaves it
    whole too. Dropping whole patents cannot shrink another patent's groups, so one pass suffices.
    """
    d = collapse(sar.loc[mask])
    g = d.groupby("group_id").agg(n=("InChIKey", "nunique"), lo=("pAct", "min"), hi=("pAct", "max"))
    ok = g[(g.n >= MIN_GROUP) & ((g.hi - g.lo) >= min_range)].index
    d = d[d.group_id.isin(ok)].copy()
    # one group per series: continuation and divisional patents of one family repeat the same table
    # on the same protein (96% shared molecules, 90% identical values), and every copy would count
    # that series again in any per-group average. The copy kept is the largest, then the earliest.
    d["series_id"] = d.family_id.fillna(-1).astype(int).astype(str) + ":" + d.uniprot_accession.fillna("")
    gi = d.groupby("group_id").agg(series=("series_id", "first"), n=("InChIKey", "nunique"),
                                   prio=("priority_date", "first"), yr=("year", "first")).reset_index()
    gi["prio"] = gi.prio.replace("", "9999")
    copies = gi.groupby("series").size()
    d["series_n_copies"] = d.series_id.map(copies)
    if one_per_series:
        keep = (gi.sort_values(["series", "n", "prio", "yr", "group_id"],
                               ascending=[True, False, True, True, True]).drop_duplicates("series").group_id)
        series_removed["last"] = len(gi) - len(keep)
        d = d[d.group_id.isin(set(keep))].copy()
    if min_patent_molecules:
        per_patent = d.groupby("patent_core").InChIKey.nunique()
        d = d[d.patent_core.isin(per_patent[per_patent >= min_patent_molecules].index)].copy()
    d["group_n_molecules"] = d.group_id.map(g.n)
    d["group_pAct_range"] = d.group_id.map(g.hi - g.lo)
    # stereo set: molecules of one group sharing a connectivity layer (first InChIKey block) differ
    # only in stereo or isotopes, which a model reading 2D graphs cannot tell apart
    skel = d.InChIKey.str[:14]
    size = d.groupby([d.group_id, skel]).InChIKey.transform("nunique")
    d["stereo_set_id"] = np.where(size > 1, d.group_id.astype(str) + ":" + skel, "")
    top = d.groupby("group_id").value_nM.transform("max")
    at_top = np.isclose(d.value_nM, top, rtol=1e-6)
    n_top = pd.Series(at_top, index=d.index).groupby(d.group_id).transform("sum")
    round_top = top.map(lambda v: any(np.isclose(v, x, rtol=1e-6) for x in ASSAY_LIMITS_NM))
    d["at_assay_ceiling"] = at_top & round_top & (n_top >= 3)
    bottom = d.groupby("group_id").value_nM.transform("min")
    at_bottom = np.isclose(d.value_nM, bottom, rtol=1e-6)
    n_bottom = pd.Series(at_bottom, index=d.index).groupby(d.group_id).transform("sum")
    round_bottom = bottom.map(lambda v: any(np.isclose(v, x, rtol=1e-6) for x in ASSAY_FLOORS_NM))
    d["at_assay_floor"] = at_bottom & round_bottom & (n_bottom >= 3)
    return d[TASK_COLS].sort_values(["group_id", "pAct"], ascending=[True, False]).reset_index(drop=True)


for name, mask in [
    ("rank_ic50", rank_base & sar["Activity Type"].eq("IC50")),
    ("rank_ki_kd", rank_base & sar["Activity Type"].isin(["Ki", "Kd"])),
    ("rank_ec50", rank_base & sar["Activity Type"].eq("EC50")),
    ("rank_highcontrast", rank_base & sar.group_pAct_range.ge(HIGH_CONTRAST_RANGE)),
    # complexes and named mutants carry a real sequence, they are simply not one wild-type
    # protein, so they ship as their own track rather than being discarded
    ("rank_complex_mutant", clean & sar.is_affinity & sar.is_biochemical
     & sar.uniprot_accession.fillna("").ne("") & sar.confidence.isin(["high", "medium"])
     & sar.target_class.isin(["complex", "mutant"])
     & sar.group_n_molecules.ge(MIN_GROUP) & sar.group_pAct_range.ge(MIN_RANGE)
     & ~sar.group_is_gradelike),
]:
    d = finalize(mask, HIGH_CONTRAST_RANGE if name == "rank_highcontrast" else MIN_RANGE)
    if name in TASK_TEST_FROM:
        d["split_time"] = d.patent_core.map(unit_of).map(unit_time_split(TASK_TEST_FROM[name])).fillna("")
        d["test_tier"] = d.test_tier.where(d.split_time.eq("test"), "")
        note(f"{name}: split_time uses its own cutoff, priority >= {TASK_TEST_FROM[name]}")
    n_series_copies = series_removed["last"]
    d.to_parquet(f"{OUT}/tasks/{name}.parquet", index=False)
    tasks[name] = {"rows": len(d), "groups": int(d.group_id.nunique()),
                   "patents": int(d.patent_core.nunique()),
                   "molecules": int(d.InChIKey.nunique()),
                   "targets": int(d.uniprot_accession.nunique()),
                   "measurements_behind_rows": int(d.n_measurements.sum()),
                   "by_split_time": d.split_time.value_counts().to_dict(),
                   "groups_by_split_time": d.groupby("split_time").group_id.nunique().to_dict(),
                   "groups_by_split_time_grant": d.groupby("split_time_grant").group_id.nunique().to_dict(),
                   "test_groups_by_tier": d[d.split_time == "test"].groupby(
                       "test_tier").group_id.nunique().to_dict(),
                   "series_copies_removed": n_series_copies,
                   "adme_safety_groups": int(d.loc[d.is_adme_safety_target, "group_id"].nunique()),
                   "rows_at_assay_floor": int(d.at_assay_floor.sum()),
                   "groups_at_assay_floor": int(d.loc[d.at_assay_floor, "group_id"].nunique()),
                   "rows_in_stereo_sets": int(d.stereo_set_id.ne("").sum()),
                   "groups_with_stereo_sets": int(d.loc[d.stereo_set_id.ne(""), "group_id"].nunique())}
    note(f"task {name:20s} rows {len(d):>7,d}  groups {d.group_id.nunique():>5,d}  "
         f"patents {d.patent_core.nunique():>5,d}  targets {d.uniprot_accession.nunique():>4,d}  "
         f"(from {int(d.n_measurements.sum()):,} measurements)")

# regression task: one row per (molecule, target, endpoint), median pAct
# biochemical assays only, as in the ranking tasks: a median over a cell EC50 and an enzyme IC50 of
# one compound on one protein mixes two different quantities
reg_src = sar[clean & sar.target_mapped & sar.is_affinity & sar.is_biochemical]
# the same per-patent floor: a patent with fewer distinct molecules contributes no rows
_reg_n = reg_src.groupby("patent_core").InChIKey.nunique()
reg_src = reg_src[reg_src.patent_core.isin(_reg_n[_reg_n >= MIN_PATENT_MOLECULES].index)].copy()
# a pair measured in patents of different splits belongs to the earliest one: it was known then.
# (The 2026-10-19 release took the alphabetical minimum, which put such pairs in test.)
ORDER = {"train": 0, "valid": 1, "test": 2}
INV = {v: k for k, v in ORDER.items()}
rkey = ["InChIKey", "uniprot_accession", "Activity Type"]
for c in ["split_time", "split_time_grant", "split_random"]:
    reg_src[c + "_o"] = reg_src[c].map(ORDER)
rg = reg_src.groupby(rkey)
reg = rg.agg(pAct=("pAct", "median"), n_records=("pAct", "size"),
             pAct_spread=("pAct", lambda s: s.max() - s.min()), SMILES=("SMILES", "first"),
             scaffold=("scaffold", "first"), heavy_atoms=("heavy_atoms", "first"),
             n_patents=("patent_core", "nunique"), year=("year", "min"),
             in_public_data=("in_public_data", "max"),
             target_organism=("target_organism", "first"), is_human_target=("is_human_target", "first"),
             has_undefined_stereo=("has_undefined_stereo", "first"),
             covalent_warhead=("covalent_warhead", "first"),
             is_isotope_labelled=("is_isotope_labelled", "first"), is_bro5_size=("is_bro5_size", "first"),
             is_adme_safety_target=("is_adme_safety_target", "first"),
             split_scaffold_mol=("split_scaffold_mol", "first"),
             split_time=("split_time_o", "min"), split_time_grant=("split_time_grant_o", "min"),
             split_random=("split_random_o", "min"))
for c in ["split_time", "split_time_grant", "split_random"]:
    reg[c] = reg[c].map(INV).fillna("")
reg["patents"] = rg.patent_core.agg(_distinct)
reg["source_record_id"] = rg.source_record_id.agg(_aligned)
reg = reg.reset_index()
reg.to_parquet(OUT + "/tasks/regression_affinity.parquet", index=False)
tasks["regression_affinity"] = {
    "rows": len(reg), "molecules": int(reg.InChIKey.nunique()),
    "targets": int(reg.uniprot_accession.nunique()),
    "source_patents": int(reg_src.patent_core.nunique()),
    "measured_more_than_once": int((reg.n_records > 1).sum()),
    "median_spread_when_repeated": round(float(reg.loc[reg.n_records > 1, "pAct_spread"].median()), 3),
    "by_split_scaffold": reg.split_scaffold_mol.value_counts().to_dict(),
    "by_split_time": reg.split_time.value_counts().to_dict()}
note(f"task regression_affinity  rows {len(reg):,}  molecules {reg.InChIKey.nunique():,}  "
     f"targets {reg.uniprot_accession.nunique():,}  patents {reg_src.patent_core.nunique():,}")

# ----------------------------------------- novelty, with and without fragment groups
# Substituent fragments are not real compounds, so public databases rarely hold them, and they
# read as "new". The headline novelty figures must therefore be quoted without them.
novelty = {}
aff = sar[sar.is_affinity]
for label, x in [("all_affinity_rows", aff),
                 ("excluding_fragment_groups", aff[~aff.group_is_fragment_suspect]),
                 ("fragment_groups_only", aff[aff.group_is_fragment_suspect])]:
    mols = x.drop_duplicates("InChIKey")
    pats = x.drop_duplicates("patent_core")
    mapped = x[x.uniprot_accession.fillna("").ne("") & x.confidence.isin(["high", "medium"])]
    novelty[label] = {
        "molecules": len(mols),
        "new_molecules": int(mols.molecule_is_new.astype("boolean").fillna(False).sum()),
        "new_molecule_share": round(float(mols.molecule_is_new.astype("boolean").fillna(False).mean()), 4),
        "patents": len(pats),
        "new_patents": int(pats.patent_is_new.astype("boolean").fillna(False).sum()),
        "new_patent_share": round(float(pats.patent_is_new.astype("boolean").fillna(False).mean()), 4),
        "mapped_rows": len(mapped),
        "new_pair_rows": int((~mapped.in_public_data).sum()),
        "new_pair_share": round(float(1 - mapped.in_public_data.mean()), 4) if len(mapped) else None}
    v = novelty[label]
    note(f"novelty {label:26s} molecules {v['molecules']:,} new {v['new_molecule_share']:.1%} | "
         f"patents {v['patents']:,} new {v['new_patent_share']:.1%} | "
         f"pairs new {v['new_pair_share']:.1%}")

# ------------------------------------------------- what each filter costs in patents
note("\n=== patents retained as each ranking filter is added (cumulative)")
steps = [("all formal SAR", pd.Series(True, index=sar.index)),
         ("affinity endpoints", sar.is_affinity),
         ("exact values", sar.is_affinity & sar.is_exact),
         ("drug-like molecules", base & sar.is_affinity),
         ("group has >= 10 molecules", base & sar.is_affinity & sar.group_n_molecules.ge(MIN_GROUP)),
         ("dynamic range >= 1 log", base & sar.is_affinity & sar.group_n_molecules.ge(MIN_GROUP)
          & sar.group_pAct_range.ge(MIN_RANGE) & ~sar.group_is_gradelike),
         ("target mapped to a sequence", base & sar.is_affinity & sar.group_n_molecules.ge(MIN_GROUP)
          & sar.group_pAct_range.ge(MIN_RANGE) & ~sar.group_is_gradelike & sar.target_mapped),
         ("biochemical assays only", rank_pre),
         ("drop substituent-table groups", rank_nofrag),
         ("hold out disputed and confirmed-wrong pairs", rank_nofrag
          & ~sar.pair_under_adjudication & ~sar.pair_validated_error),
         ("hold out whole patents with a validation issue", rank_nofrag
          & ~sar.pair_under_adjudication & ~sar.pair_validated_error
          & ~sar.patent_has_validation_issue),
         ("drop structures resting on an unexpanded abbreviation", rank_nofrag & ~validation_out),
         ("drop tables on a wrong unit scale", rank_nofrag & ~validation_out & ~sar.group_unit_suspect),
         ("drop tables read twice in one patent", rank_nofrag & ~validation_out
          & ~sar.group_unit_suspect & ~sar.group_is_duplicate),
         ("drop labels that name no assay", rank_base)]
regrouped = finalize(rank_base, MIN_RANGE, 0, one_per_series=False)
steps_extra = [("one value per compound; regroup >= 10 molecules", regrouped),
               ("one group per series (family copies removed)", finalize(rank_base, MIN_RANGE, 0)),
               (f"patents with >= {MIN_PATENT_MOLECULES} molecules",
                finalize(rank_base, MIN_RANGE, MIN_PATENT_MOLECULES))]
cost = []
for label, mask in steps:
    d = sar[mask]
    cost.append({"filter": label, "rows": len(d), "patents": int(d.patent_core.nunique()),
                 "groups": int(d.group_id.nunique()), "molecules": int(d.InChIKey.nunique())})
    note(f"  {label:32s} rows {len(d):>7,d}  patents {d.patent_core.nunique():>5,d}  "
         f"groups {d.group_id.nunique():>6,d}  molecules {d.InChIKey.nunique():>7,d}")
for label, d in steps_extra:
    cost.append({"filter": label, "rows": len(d), "patents": int(d.patent_core.nunique()),
                 "groups": int(d.group_id.nunique()), "molecules": int(d.InChIKey.nunique())})
    note(f"  {label:32s} rows {len(d):>7,d}  patents {d.patent_core.nunique():>5,d}  "
         f"groups {d.group_id.nunique():>6,d}  molecules {d.InChIKey.nunique():>7,d}")
relaxed = {
    "ranking incl. cell-based assays": int(sar.loc[
        clean & sar.target_mapped & sar.group_n_molecules.ge(MIN_GROUP)
        & sar.group_pAct_range.ge(MIN_RANGE) & ~sar.group_is_gradelike, "patent_core"].nunique()),
    "ranking without the target-mapping requirement": int(sar.loc[
        clean & sar.is_affinity & sar.group_n_molecules.ge(MIN_GROUP)
        & sar.group_pAct_range.ge(MIN_RANGE) & ~sar.group_is_gradelike, "patent_core"].nunique()),
    "ranking with groups >= 5 molecules": int(sar.loc[
        clean & sar.target_mapped & sar.is_biochemical & sar.group_n_molecules.ge(5)
        & sar.group_pAct_range.ge(MIN_RANGE) & ~sar.group_is_gradelike, "patent_core"].nunique()),
}
note("\n=== patents recoverable by relaxing one filter at a time")
for k, v in relaxed.items():
    note(f"  {k:48s} {v:>5,d} patents")

manifest = {"name": NAME, "version": VERSION, "built": BUILT, "seed": SEED, "config": cfg,
            "patent_cost_per_filter": cost, "patents_if_relaxed": relaxed,
            "assay_noise_floor_log": round(noise, 3),
            "inputs": {"sar": R, "target_map": M, "public_cache": OUT + "/work/"},
            "full_record": {"rows": len(sar), "patents": int(sar.patent_core.nunique()),
                            "molecules": int(sar.InChIKey.nunique()),
                            "endpoints": int(sar["Activity Type"].nunique()),
                            "groups": int(sar.group_id.nunique())},
            "flag_counts": {c: int(sar[c].fillna(False).sum()) for c in
                            ["is_affinity", "is_exact", "is_categorical_derived", "is_druglike",
                             "target_mapped", "is_biochemical", "replicate_disagrees",
                             "is_family_representative", "in_public_data",
                             "group_is_fragment_suspect", "pair_under_adjudication",
                             "pair_validated_error", "patent_has_validation_issue",
                             "structure_alias_suspect", "group_unit_suspect", "group_is_duplicate",
                             "group_label_uninformative", "label_repaired",
                             "has_undefined_stereo",
                             "is_isotope_labelled", "is_bro5_size", "value_beyond_assay_range",
                             "is_human_target", "is_adme_safety_target"]},
            "industrial_checks": {
                "provenance_coverage": round(float(has_ev.mean()), 4),
                "rows_with_salt_or_charge": n_salt,
                "drawn_molecules": int(sar.InChIKey_as_drawn.nunique()),
                "parent_molecules": int(sar.InChIKey.nunique()),
                "groups_on_raw_label": int(n_raw_groups), "groups_after_label_join": int(sar.group_id.nunique()),
                "unit_suspect_groups_by_check": {c: int(gu[c].sum()) for c in UNIT_CHECKS},
                "unit_suspect_groups": int(gu.issue.ne("").sum()),
                "duplicate_groups": len(dup_of),
                "uninformative_label_groups": int(sar.loc[sar.group_label_uninformative, "group_id"].nunique()),
                "label_repaired_rows": int(sar.label_repaired.sum()),
                "label_repaired_groups_now": int(sar.loc[sar.label_repaired, "group_id"].nunique()),
                "groups_without_usable_rows": int((~groups.has_usable_rows).sum()),
                "covalent_warhead_molecules": int(sar.loc[sar.covalent_warhead.fillna("").ne(""), "InChIKey"].nunique()),
                "priority_date_source": pp.priority_source.value_counts().to_dict(),
                "time_split_units_left_out": int((unit_split == "").sum())},
            "validation_holdout": {"disputed_items": int(val.disputed.sum()),
                                   "disputed_pairs": len(disputed_pairs),
                                   "confirmed_wrong_items": int(val.wrong.sum()),
                                   "confirmed_wrong_pairs": len(wrong_pairs),
                                   "settled_by_adjudication": n_adjudicated,
                                   "patents_held_out_whole": len(doubt_patents),
                                   "alias_suspect_pairs": len(suspect_pairs)},
            "test_tier_groups": {t: int(sar.loc[sar.test_tier == t, "group_id"].nunique())
                                 for t in ["novel", "reference", "unscored"]},
            "target_mapping_pair_coverage": round(len(pair_map) / len(all_pairs), 4),
            "novelty_vs_bindingdb_chembl": novelty,
            "tasks": tasks, "log": log}
with open(OUT + "/release_manifest.json", "w") as fh:
    json.dump(manifest, fh, indent=1, ensure_ascii=False, default=str)
note(f"\nwrote {OUT}/sar_all.parquet, targets.parquet, groups.parquet, tasks/, release_manifest.json")
