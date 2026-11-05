"""Assemble the data-descriptor release: the full formal SAR table plus quality and split flags,
then cut the task subsets from it.

Nothing is deleted. Every selection the paper reports is reproducible from flag columns on
`sar_all.parquet`, which is the point of publishing it this way.
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import json
import os
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

R = _ws + "/outputs/uspto-sar-recovery-20261006/results/"
M = _ws + "/outputs/target-uniprot-map-v2-20261019/"
B = _ws + "/outputs/affinity-rank-benchmark-20261019/"
VAL = _ws + "/outputs/sar-validation-20261019"
OUT = sys.argv[1].rstrip("/")
SEED = 20261019

AFFINITY = ["IC50", "Ki", "Kd", "EC50"]
KEY = ["patent_core", "Target", "Activity Type"]
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
TEST_FROM_YEAR = 2021
VALID_FRACTION = 0.15
RANDOM_FRACTIONS = (0.70, 0.15, 0.15)
SCAFFOLD_FRACTIONS = (0.70, 0.15, 0.15)
cfg = {k: v for k, v in sorted(globals().items())
       if k.isupper() and k not in {"R", "M", "B", "OUT", "AFFINITY", "KEY", "ORGANIC"}}
log = []


def note(msg):
    print(msg, flush=True)
    log.append(msg)


# ----------------------------------------------------------------- full SAR
sar = ds.dataset(R + "SAR_merged.parquet").to_table(
    columns=["source", "source_record_id", "patent_core", "Patent", "Target", "Activity Type",
             "Relation", "Value", "Unit", "value_nM", "SMILES", "InChIKey"]).to_pandas()
sar["Target"] = sar.Target.fillna("").str.strip()
note(f"formal SAR rows: {len(sar):,}  endpoints: {sar['Activity Type'].nunique()}  "
     f"patents: {sar.patent_core.nunique():,}  molecules: {sar.InChIKey.nunique():,}")

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

# ------------------------------------------------------- molecules: size, scaffold
uniq = sar.drop_duplicates("InChIKey")[["InChIKey", "SMILES"]]
cache = OUT + "/work/mol_props.parquet"
try:
    props = pq.read_table(cache).to_pandas()
    if set(props.InChIKey) >= set(uniq.InChIKey):
        note(f"reusing molecule properties for {len(props):,} molecules from {cache}")
    else:
        raise ValueError("cache does not cover every molecule")
except Exception:
    note(f"parsing {len(uniq):,} unique molecules for heavy-atom count and Bemis-Murcko scaffold")
    recs = []
    for k, s in zip(uniq.InChIKey, uniq.SMILES):
        mol = Chem.MolFromSmiles(s) if isinstance(s, str) else None
        if mol is None:
            recs.append((k, -1, True, ""))
            continue
        try:
            core = MurckoScaffold.GetScaffoldForMol(mol)
            sc = Chem.MolToSmiles(core) if core.GetNumAtoms() else ""
        except Exception:
            sc = ""
        recs.append((k, mol.GetNumHeavyAtoms(),
                     any(a.GetSymbol() not in ORGANIC for a in mol.GetAtoms()), sc))
    props = pd.DataFrame(recs, columns=["InChIKey", "heavy_atoms", "unusual_elements", "scaffold"])
    props.to_parquet(cache, index=False)
props = props.drop_duplicates("InChIKey").set_index("InChIKey")
sar["heavy_atoms"] = sar.InChIKey.map(props.heavy_atoms)
sar["scaffold"] = sar.InChIKey.map(props.scaffold)
sar["is_druglike"] = (sar.heavy_atoms.between(HEAVY_MIN, HEAVY_MAX)
                      & ~sar.InChIKey.map(props.unusual_elements).fillna(True))
note(f"is_druglike {int(sar.is_druglike.sum()):,} rows; "
     f"{int((sar.scaffold == '').sum()):,} rows have no ring scaffold")

# ---------------------------------------------------------------- patent year
pt = ds.dataset(R + "patents.parquet").to_table(columns=["Patent", "archive"]).to_pandas()
pt["patent_core"] = pt.Patent.str.replace(r"(B\d|A\d|E|S|H|P\d)$", "", regex=True)
pt["year"] = pt.archive.str.extract(r"/data/(\d{4})/")[0].astype(int)
sar["year"] = sar.patent_core.map(pt.drop_duplicates("patent_core").set_index("patent_core").year)
note(f"patent year missing for {int(sar.year.isna().sum()):,} rows")

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

# --------------------------------------------- one value per (group, molecule)
sar["group_id"] = sar.groupby(KEY, sort=False).ngroup()
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
sar["structure_alias_suspect"] = [(p, k) in suspect_pairs for p, k in zip(sar.Patent, sar.InChIKey)]
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
qualifies = (gstat.group_n_molecules.ge(MIN_GROUP) & gstat.group_pAct_range.ge(MIN_RANGE)
             & (group_median_heavy.reindex(gstat.index) >= FRAGMENT_GROUP_MEDIAN_HEAVY))
fam["qualifies"] = fam.group_id.map(qualifies).fillna(False)
rep = (fam.sort_values(["family_id", "qualifies", "year", "n"],
                       ascending=[True, False, True, False])
       .drop_duplicates("family_id"))
sar["is_family_representative"] = sar.group_id.isin(set(rep.group_id))
note(f"family clustering: {len(fam):,} groups -> {fam.family_id.nunique():,} families; "
     f"{int((fam.groupby('family_id').size() > 1).sum()):,} families hold more than one group")

# -------------------------------------------------------------- public overlap
pub = pq.read_table(OUT + "/work/public_pairs.parquet").to_pandas()
pub = pub[pub.source.isin(PUBLIC_SOURCES)]
pub_pairs = set(zip(pub.ik, pub.uniprot))
note(f"public comparison set: {', '.join(PUBLIC_SOURCES)} "
     f"({len(pub_pairs):,} distinct (molecule, target) pairs)")
sar["in_public_data"] = [(k, u) in pub_pairs if isinstance(u, str) and u else False
                         for k, u in zip(sar.InChIKey, sar.uniprot_accession)]
try:
    nov = pq.read_table(OUT + "/work/novelty_flags.parquet").to_pandas()
    mol_new = dict(zip(nov.InChIKey, nov.molecule_is_new))
    pat_new = dict(zip(nov.patent_core, nov.patent_is_new))
    sar["molecule_is_new"] = sar.InChIKey.map(mol_new)
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

# time split: by the unit's earliest year
unit_split = pd.Series(np.where(unit_year >= TEST_FROM_YEAR, "test", "train"), index=unit_year.index)
early = sorted(unit_split[unit_split == "train"].index)
valid_units = set(rng.choice(early, size=int(round(len(early) * VALID_FRACTION)), replace=False))
unit_split[unit_split.index.isin(valid_units)] = "valid"
unit_split[unit_year.isna()] = ""
sar["split_time"] = sar._unit.map(unit_split).fillna("")
note("split_time: " + str(sar.split_time.value_counts().to_dict()))

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
cols = ["source", "source_record_id", "Patent", "patent_core", "year", "Target", "Activity Type",
        "Relation", "Value", "Unit", "value_nM", "pAct", "SMILES", "InChIKey", "scaffold",
        "heavy_atoms", "group_id", "family_id", "uniprot_accession", "mapped_identity",
        "target_class", "confidence", "measured_in_cells", "mappability",
        "mapping_run",
        "is_affinity", "is_exact", "is_categorical_derived", "is_druglike", "target_mapped",
        "is_biochemical", "replicate_disagrees", "is_family_representative",
        "group_n_molecules", "group_pAct_min", "group_pAct_max", "group_pAct_range",
        "group_pAct_std", "group_n_distinct", "group_is_gradelike",
        "group_median_heavy_atoms", "group_is_fragment_suspect",
        "pair_under_adjudication", "pair_validated_error",
        "patent_has_validation_issue", "structure_alias_suspect",
        "in_public_data", "molecule_is_new", "patent_is_new",
        "split_time", "split_random", "split_scaffold_mol", "split_scaffold_group", "test_tier"]
sar[cols].to_parquet(OUT + "/sar_all.parquet", index=False)

seq = tmap.drop_duplicates("uniprot_accession")[
    ["uniprot_accession", "uniprot_entry", "protein_name", "organism", "tax_id", "seq_len", "sequence"]]
seq.to_parquet(OUT + "/targets.parquet", index=False)

groups = sar[sar.group_id.isin(gstat.index)].groupby("group_id").agg(
    patent_core=("patent_core", "first"), year=("year", "first"), Target=("Target", "first"),
    endpoint=("Activity Type", "first"), uniprot_accession=("uniprot_accession", "first"),
    family_id=("family_id", "first"), n_molecules=("group_n_molecules", "first"),
    pAct_range=("group_pAct_range", "first"), pAct_std=("group_pAct_std", "first"),
    target_mapped=("target_mapped", "first"), is_biochemical=("is_biochemical", "first"),
    is_family_representative=("is_family_representative", "first"),
    median_heavy_atoms=("group_median_heavy_atoms", "first"),
    is_fragment_suspect=("group_is_fragment_suspect", "first"),
    public_share=("in_public_data", "mean"), split_time=("split_time", "first"),
    split_scaffold_group=("split_scaffold_group", "first"), test_tier=("test_tier", "first"),
).reset_index()
groups.to_parquet(OUT + "/groups.parquet", index=False)

# ------------------------------------------------------------------- tasks
base = (sar.is_exact & sar.is_druglike & ~sar.replicate_disagrees)
# every task also drops substituent-table groups; base stays as it was for the cost table
held_out = (sar.pair_under_adjudication | sar.pair_validated_error
            | sar.patent_has_validation_issue | sar.structure_alias_suspect)
clean = base & ~sar.group_is_fragment_suspect & ~held_out
# no is_family_representative filter here: splits already keep each family inside one split, so
# family duplicates cannot leak. The flag stays on the rows for anyone wanting strict dedup.
rank_pre = (base & sar.target_mapped & sar.is_biochemical
            & sar.group_n_molecules.ge(MIN_GROUP) & sar.group_pAct_range.ge(MIN_RANGE)
            & ~sar.group_is_gradelike)
rank_nofrag = rank_pre & ~sar.group_is_fragment_suspect
rank_base = rank_nofrag & ~held_out
task_cols = ["patent_core", "year", "Target", "Activity Type", "group_id", "family_id",
             "is_family_representative", "uniprot_accession", "mapped_identity", "InChIKey",
             "SMILES", "scaffold", "heavy_atoms", "pAct", "value_nM", "group_n_molecules",
             "group_pAct_range", "in_public_data", "molecule_is_new", "split_time",
             "split_random", "split_scaffold_group", "test_tier"]
tasks = {}


def finalize(mask, min_range, min_patent_molecules=MIN_PATENT_MOLECULES):
    """Recount each group, then each patent, on the rows that survived every exclusion.

    Group size and range were measured before row-level exclusions such as the alias check, so a
    group can fall under the bar once its suspect rows are gone; it then leaves the task whole.
    A patent left with fewer than `min_patent_molecules` distinct molecules in the task leaves it
    whole too. Dropping whole patents cannot shrink another patent's groups, so one pass suffices.
    """
    d = sar.loc[mask, task_cols].copy()
    g = d.groupby("group_id").agg(n=("InChIKey", "nunique"), lo=("pAct", "min"), hi=("pAct", "max"))
    ok = g[(g.n >= MIN_GROUP) & ((g.hi - g.lo) >= min_range)].index
    d = d[d.group_id.isin(ok)].copy()
    if min_patent_molecules:
        per_patent = d.groupby("patent_core").InChIKey.nunique()
        d = d[d.patent_core.isin(per_patent[per_patent >= min_patent_molecules].index)].copy()
    d["group_n_molecules"] = d.group_id.map(g.n)
    d["group_pAct_range"] = d.group_id.map(g.hi - g.lo)
    return d


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
    d.to_parquet(f"{OUT}/tasks/{name}.parquet", index=False)
    tasks[name] = {"rows": len(d), "groups": int(d.group_id.nunique()),
                   "patents": int(d.patent_core.nunique()),
                   "molecules": int(d.InChIKey.nunique()),
                   "targets": int(d.uniprot_accession.nunique()),
                   "by_split_time": d.split_time.value_counts().to_dict(),
                   "groups_by_split_time": d.groupby("split_time").group_id.nunique().to_dict(),
                   "test_groups_by_tier": d[d.split_time == "test"].groupby(
                       "test_tier").group_id.nunique().to_dict()}
    note(f"task {name:20s} rows {len(d):>7,d}  groups {d.group_id.nunique():>5,d}  "
         f"patents {d.patent_core.nunique():>5,d}  targets {d.uniprot_accession.nunique():>4,d}")

# regression task: one row per (molecule, target, endpoint), median pAct
reg_src = sar[clean & sar.target_mapped & sar.is_affinity]
# the same per-patent floor: a patent with fewer distinct molecules contributes no rows
_reg_n = reg_src.groupby("patent_core").InChIKey.nunique()
reg_src = reg_src[reg_src.patent_core.isin(_reg_n[_reg_n >= MIN_PATENT_MOLECULES].index)]
reg = (reg_src.groupby(["InChIKey", "uniprot_accession", "Activity Type"])
       .agg(pAct=("pAct", "median"), n_records=("pAct", "size"),
            pAct_spread=("pAct", lambda s: s.max() - s.min()), SMILES=("SMILES", "first"),
            scaffold=("scaffold", "first"), heavy_atoms=("heavy_atoms", "first"),
            n_patents=("patent_core", "nunique"), year=("year", "min"),
            in_public_data=("in_public_data", "max"),
            split_scaffold_mol=("split_scaffold_mol", "first"),
            split_time=("split_time", "min"), split_random=("split_random", "min"))
       .reset_index())
reg.to_parquet(OUT + "/tasks/regression_affinity.parquet", index=False)
tasks["regression_affinity"] = {
    "rows": len(reg), "molecules": int(reg.InChIKey.nunique()),
    "targets": int(reg.uniprot_accession.nunique()),
    "source_patents": int(reg_src.patent_core.nunique()),
    "measured_more_than_once": int((reg.n_records > 1).sum()),
    "median_spread_when_repeated": round(float(reg.loc[reg.n_records > 1, "pAct_spread"].median()), 3),
    "by_split_scaffold": reg.split_scaffold_mol.value_counts().to_dict()}
note(f"task regression_affinity  rows {len(reg):,}  molecules {reg.InChIKey.nunique():,}  "
     f"targets {reg.uniprot_accession.nunique():,}")

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
         ("drop structures resting on an unexpanded abbreviation", rank_base)]
regrouped = finalize(rank_base, MIN_RANGE, 0)
steps_extra = [("regroup: >= 10 molecules after row exclusions", regrouped),
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

manifest = {"built": "2026-10-19", "seed": SEED, "config": cfg,
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
                             "structure_alias_suspect"]},
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
