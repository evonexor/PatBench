"""Write docs/DATA_SCHEMA.md from a schema dump (tools/dump_schemas.py) and the descriptions below.

    python tools/make_schema_doc.py docs/data/schemas_r6.json > docs/DATA_SCHEMA.md

Every column of every table is listed with its type. A description is given only where the pipeline's
own code or its release documentation states what the column is; the rest are listed bare rather than
guessed at.
"""
import json
import sys
from collections import OrderedDict

D = {
    # identity and provenance
    "source": "Extraction run: `USPTO_MOL` (original extraction) or `USPTO_MOL_RECOVERY` (recovery run)",
    "source_record_id": "Extraction record id; unique per row of `sar_all`. In task files, rows merged from repeated measurements join the ids with ` | ` in measurement order",
    "Patent": "Patent number with kind code (e.g. US10899738B2)",
    "patent_core": "Patent number without kind code; the patent key across tables",
    "year": "Grant year (from the year directory of the archive)",
    "priority_date": "Priority date: earliest priority claim, else PCT filing date, else US filing date; the grant year is used as an upper bound when none of the three exists",
    "priority_source": "Where `priority_date` comes from",
    "first_public_date": "First public date: the earliest of US application publication, PCT publication and grant date; left empty (not estimated) when none exists",
    "first_public_year": "Year of the first public date",
    "first_public_source": "`us_application` / `wo` / `grant` / `none`",
    "compound_alias": "Compound number printed in the patent",
    "activity_locator": "Cell holding the activity: table / tgroup / row / column, e.g. `TABLE-US-00007/tgroup-2/row-152/column-8`",
    "activity_header": "Verbatim header of that column",
    "raw_value": "Verbatim text of the cell",
    "structure_file": "Original MOL file of the compound (path inside the archive)",
    # assay and value
    "Target": "Assay label as printed in the patent (verbatim)",
    "target_label": "Cleaned label (case and punctuation normalised), used for grouping and display",
    "Activity Type": "Endpoint: IC50, Ki, Kd, EC50, ...",
    "Relation": "Relation sign (`=`, `>`, `<`, ...)",
    "Value": "Value in the original unit",
    "Unit": "Unit",
    "value_nM": "Value converted to nM",
    "pAct": "9 - log10(value_nM), only for exact values within 2-13; in task files the median over repeated measurements",
    "n_measurements": "Number of measurements merged into this row",
    "pAct_replicate_spread": "Largest difference between repeated measurements (log units)",
    # molecule
    "SMILES": "Parent molecule (largest organic fragment, neutralised)",
    "InChIKey": "InChIKey of the parent molecule",
    "SMILES_as_drawn": "Structure as drawn in the patent",
    "InChIKey_as_drawn": "InChIKey of the structure as drawn",
    "scaffold": "Murcko scaffold",
    "heavy_atoms": "Number of heavy atoms",
    "stereo_set_id": "Id shared by molecules that differ only in stereochemistry or isotopes (same first InChIKey block); empty otherwise",
    "n_stereocentres": "Number of stereocentres",
    "has_undefined_stereo": "At least one stereocentre is undefined",
    "is_isotope_labelled": "Isotope-labelled (mostly deuterated)",
    "covalent_warhead": "Contains a covalent warhead (acrylamide, chloroacetamide, vinyl sulfone, boronic acid, ...)",
    "is_bro5_size": "More than 50 heavy atoms (macrocycles, PROTACs, ...)",
    "structure_repair": "Structure-repair rules applied to the stored structure: `NH`, `group`, `nitrile`, `text`, joined by `;`; empty when unchanged",
    "SMILES_before_repair": "SMILES stored before the structure repair",
    "InChIKey_before_repair": "InChIKey stored before the structure repair",
    # groups and series
    "group_id": "Group = one patent x one assay label x one endpoint",
    "family_id": "Patent family (clustered by shared molecules, not INPADOC); numbered by the smallest `group_id` in the family",
    "is_family_representative": "Whether this group represents its patent family",
    "series_id": "Series = one patent family x one protein; only one group per series is kept in a task",
    "series_n_copies": "How many copies the series originally had (continuations and divisionals reuse the same table)",
    "n_rows": "Number of records in the group",
    "n_molecules": "Number of distinct molecules in the group",
    "group_n_molecules": "Number of distinct molecules in the group (recounted after row hold-outs)",
    "pAct_range": "Largest pAct difference within the group",
    "pAct_std": "Standard deviation of pAct within the group",
    "group_pAct_range": "Largest pAct difference within the group",
    "group_pAct_min": "Smallest pAct in the group",
    "group_pAct_max": "Largest pAct in the group",
    "group_pAct_std": "Standard deviation of pAct within the group",
    "group_n_distinct": "Number of distinct pAct values in the group",
    "group_is_gradelike": "Values in the group look like grades (A/B/C, ...)",
    "group_median_heavy_atoms": "Median heavy-atom count of the molecules in the group",
    "median_heavy_atoms": "Median heavy-atom count of the molecules in the group",
    "endpoint": "Endpoint",
    "has_usable_rows": "The group has at least one exact, drug-like record whose repeated measurements do not contradict each other",
    # target
    "uniprot_accession": "UniProt accession of the target (from the target mapping)",
    "mapped_identity": "Mapped target identity (gene name; non-protein assays are written `ORG:species|name`)",
    "target_class": "Target class",
    "target_organism": "Target organism",
    "is_human_target": "Whether the target is human",
    "confidence": "Confidence of the target mapping: high / medium / low / none",
    "measured_in_cells": "The assay was run in cells",
    "mappability": "Judgement of whether the label can be mapped to a target",
    "mapping_run": "Target-mapping run (v1 / v2)",
    "target_mapped": "Target mapped at high or medium confidence",
    "is_adme_safety_target": "Off-target panels: CYP inhibition, hERG, Nav1.5, Cav1.2, transporters, ...",
    # quality flags
    "is_affinity": "Endpoint is an affinity type (IC50, Ki, Kd, EC50)",
    "is_exact": "Exact value (relation `=`, and not derived from an A/B/C legend)",
    "is_categorical_derived": "Threshold or interval derived from an A/B/C legend, not a printed number",
    "is_druglike": "Drug-like (heavy-atom count, element composition, ...)",
    "is_biochemical": "Biochemical (cell-free) assay",
    "replicate_disagrees": "Repeated measurements of one compound differ by more than 1 log; the pair is not used",
    "group_is_fragment_suspect": "Median heavy-atom count of the group below 16: a substituent table",
    "group_unit_suspect": "Unit or table-layout problem",
    "group_unit_issue": "Reason for `group_unit_suspect` (one of the six unit checks)",
    "group_is_duplicate": "Duplicate group: the same table read twice",
    "group_duplicate_of": "The kept group a duplicate group points to",
    "group_label_uninformative": "The label does not point to any assay",
    "label_repaired": "The label was a bare number; the group was rebuilt from the header text",
    "pair_under_adjudication": "(patent, target) pair on which the two annotators of the validation sample disagreed",
    "pair_validated_error": "(patent, target) pair judged wrong by both annotators",
    "patent_has_validation_issue": "The patent has any disputed or wrong validation record; the whole patent is held out",
    "structure_alias_suspect": "Structure rests on an unexpanded MOL abbreviation; the stored structure is missing atoms",
    "value_beyond_assay_range": "pAct above 11 or below 3, mostly extrapolations or tight-binding limits",
    "at_assay_ceiling": "The weakest integer concentration of the group is shared by at least 3 molecules (assay limit written as an exact value)",
    "at_assay_floor": "The strongest integer concentration of the group is shared by at least 3 molecules (assay floor written as an exact value)",
    # novelty
    "in_public_data": "The (molecule, target) pair occurs in BindingDB or ChEMBL",
    "molecule_is_new": "The molecule is not in BindingDB / ChEMBL",
    "patent_is_new": "The patent is not cited by BindingDB / ChEMBL",
    "public_share": "Share of the group's records that fall in the public databases",
    # splits
    "split_time": "Main split: test when the earliest priority year of the component is at or after the cutoff and every patent was granted in 2020 or later; 15% of the rest random valid",
    "split_time_grant": "Comparison: test when the earliest grant year of the component is 2021 or later (rule of an earlier version)",
    "split_random": "Random 70/15/15 over components",
    "split_scaffold_mol": "By molecular scaffold (regression)",
    "split_scaffold_group": "By the dominant scaffold of the group (ranking)",
    "test_tier": "Test groups only: `novel` = no (molecule, target) pair of the group is in the public databases; `reference` = the rest; `unscored` = no high/medium-confidence target",
    # targets / sequences
    "uniprot_entry": "UniProt entry name",
    "protein_name": "Protein name",
    "organism": "Organism",
    "tax_id": "NCBI taxonomy id",
    "seq_len": "Sequence length",
    "sequence": "UniProt canonical sequence (`targets.parquet`); in `target_sequences.parquet` the sequence to use for the group",
    "sequence_id": "Sequence id: `ACC`, `ACC:start-end` (mature chain) or `ACC:mutation`",
    "sequence_kind": "`full_length` / `mature_chain` / `mutant_applied` / `complex_primary_subunit`, ...",
    "label_check": "Whether the label names the target (crude string comparison)",
    "xref_status": "Relation to the public-database cross-reference: `agrees` / `differs` / `no_cross_reference`",
    # regression
    "n_records": "Regression: number of records merged into this row",
    "pAct_spread": "Regression: largest difference between the merged records",
    "n_patents": "Regression: number of patents this (molecule, target, endpoint) comes from",
    "patents": "Regression: source patents",
    # eval support
    "domain_start": "First residue modelled by the recommended structure (position in the group's sequence)",
    "domain_end": "Last residue modelled by the recommended structure",
    "domain_len": "Length of the structural domain",
    "domain_sequence": "Slice of the group's sequence between [domain_start, domain_end]",
    "target_seen_in_train_same_task": "The same UniProt accession occurs in a training group of the same task",
    "mapping_strict": "High confidence, and the label names the target, and the public cross-reference is not `differs`",
    "docking_eligible": "Has a structure with a ligand, the sequence matches, and the co-crystallised ligand is not among the target's test molecules",
    "structure_ligand_leak": "The co-crystallised ligand of the recommended structure is also a test molecule of this target",
    "split_all": "All `test` in `eval_all`: a model not trained on PatBench can use every group as test data",
    "test_tier_time": "Tier kept from the time split in `eval_all` (`test_tier` is recomputed over all groups)",
    "priority_year": "Priority year, for stratification by year",
    "is_fragment_suspect": "Median heavy-atom count of the group below 16: a substituent table",
    "unit_issue": "Reason for the unit flag of the group (one of the six unit checks)",
    "is_duplicate": "Duplicate group: the same table read twice",
    "duplicate_of": "The kept group a duplicate group points to",
    "label_uninformative": "The label does not point to any assay",
    "mutation": "Variant or mutation recorded for the target of the group",
    "organism_stated": "Species stated in the patent",
    "organism_detected": "Species detected for the assay",
    "organism_resolution": "How the species of the mapping was resolved",
    "sequence_caveat": "Caveat recorded for the sequence",
    "xref_db": "Public database that cross-references the patent (BindingDB or ChEMBL)",
    "xref_accession": "UniProt accession given by the public cross-reference",
    "xref_symbol": "Gene symbol given by the public cross-reference",
    # protocol / baseline files
    "score": "Predicted score; higher means predicted stronger",
    "role": "Few-shot protocol: `support` (value revealed) or `query` (scored)",
    "repeat": "Which draw of the few-shot protocol",
    "tanimoto": "ECFP4 Tanimoto similarity",
    "delta": "Difference in pAct between the two molecules",
}

TITLES = OrderedDict([
    ("sar_all.parquet", "All formal SAR records with provenance, quality flags and splits (one measurement per row)"),
    ("groups.parquet", "Metadata of all groups"),
    ("targets.parquet", "UniProt accession -> sequence, protein name, organism"),
    ("tasks/rank_ic50.parquet", "Ranking task (the five ranking files have the same columns)"),
    ("tasks/regression_affinity.parquet", "Regression task: one row = one (parent molecule, target, endpoint)"),
    ("target_sequences.parquet", "The target sequence to use for each ranking-task group, with the checks"),
    ("target_audit.parquet", "Group-level table behind the target checks"),
    ("structures/target_structures.parquet", "PDB / AlphaFold structures of every target"),
    ("structures/group_structure_map.parquet", "Recommended structure of every group"),
    ("eval_support/target_domain_sequences.parquet", "Structural-domain sequence of every group"),
    ("eval_support/test_group_profile.parquet", "Profile of the test groups of the five ranking tasks"),
    ("eval_all/rank_ic50.parquet", "Full copy for models not trained on PatBench (four more columns than the task file)"),
    ("protocols/rank_ic50_time_k5.parquet", "Few-shot protocol"),
    ("protocols/rank_ic50_cliffs.parquet", "Activity-cliff pairs"),
    ("baselines/rank_ic50_heavy_atoms.parquet", "Prediction file format"),
])


def main():
    d = json.load(open(sys.argv[1]))
    out = ["# Data schema", "",
           "Generated by `tools/make_schema_doc.py` from the column headers of the r6 release (`docs/data/schemas_r6.json`,",
           "written by `tools/dump_schemas.py`). A description is given only for columns whose meaning the pipeline code or the",
           "release documentation states; the remaining columns are listed with their type and nothing else.", ""]
    for f, title in TITLES.items():
        if f not in d:
            continue
        v = d[f]
        out += [f"## `{f}`", "", f"{title}. {v['rows']:,} rows, {len(v['cols'])} columns.", "",
                "| Column | Type | Meaning |", "|---|---|---|"]
        bare = []
        for c, t in v["cols"]:
            if c in D:
                out.append(f"| `{c}` | {t} | {D[c]} |")
            else:
                bare.append(f"`{c}` ({t})")
        out.append("")
        if bare:
            out += ["Columns without a description: " + ", ".join(bare), ""]
    print("\n".join(out))


if __name__ == "__main__":
    main()
