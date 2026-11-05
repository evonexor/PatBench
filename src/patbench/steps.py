"""The processing steps, in the order a PatBench release is built.

Each step is one script under pipeline/. The registry records where the script is, how it is called,
what it reads and what it writes. `RELEASE_CHAIN` is the sequence that turns the upstream tables into a
release directory (the sequence of run_r6.sh and run_baselines_r6.sh of the original run). The steps
outside the chain are the structure repair, the target audit, the structure download, the support
tables and the verification; they are run one at a time, with the arguments given in `notes`.

The upstream stages (extraction, recovery, target mapping, validation) are long jobs with their own
command lines and, in two places, manual curation. They are listed in `STAGES` so `patbench stages`
documents them, and can be started with `patbench exec`, but they are not part of the release chain.
"""
from dataclasses import dataclass


@dataclass
class Ctx:
    out: str = ""          # the release directory being built
    workers: int = 20


@dataclass
class Step:
    name: str
    script: str                       # relative to pipeline/
    summary: str
    argv: object = None               # callable(Ctx) -> list[str]
    reads: tuple = ()
    writes: tuple = ()
    omp_threads: int = 1
    heavy: bool = False
    network: bool = False
    needs_out: bool = True
    notes: str = ""

    def args(self, ctx):
        return list(self.argv(ctx)) if self.argv else []


STEPS = {s.name: s for s in [
    # ---- release chain ------------------------------------------------------------------------------
    Step("public-pairs", "08_build/public_overlap.py",
         "cache the (InChIKey, UniProt) pairs of BindingDB / ChEMBL / HARVEST with their median pAct",
         lambda c: [f"{c.out}/work"], reads=("standardized/*.parquet",), writes=("work/public_pairs.parquet",)),
    Step("standardize", "08_build/standardize.py",
         "key every drawn molecule to its parent (largest organic fragment, neutralised); size, scaffold, flags",
         lambda c: [c.out, str(c.workers)], reads=("SAR_merged.parquet",), writes=("work/mol_std.parquet",)),
    Step("novelty", "08_build/novelty.py",
         "which molecules / patents of the extraction are absent from BindingDB and ChEMBL",
         lambda c: [c.out], reads=("SAR_merged.parquet", "target_map.parquet", "standardized/*.parquet"),
         writes=("work/novelty_flags.parquet", "novelty_vs_public.json")),
    Step("build", "08_build/build_patbench.py",
         "assemble sar_all.parquet, groups, targets, splits and the task files with all quality flags",
         lambda c: [c.out],
         reads=("SAR_merged", "new_SAR", "SAR_with_evidence", "patents", "target_map", "validation annotations",
                "work/{mol_std,public_pairs,alias_scan,novelty_flags}.parquet"),
         writes=("sar_all.parquet", "groups.parquet", "targets.parquet", "tasks/*.parquet", "release_manifest.json"),
         omp_threads=2, notes="2-10 minutes"),
    Step("protocols", "08_build/make_protocols.py",
         "freeze the few-shot supports and the activity-cliff pairs",
         lambda c: [c.out], reads=("tasks/*.parquet",), writes=("protocols/*.parquet",)),
    Step("eval-sets", "08_build/make_eval_sets.py",
         "full copies of the ranking tasks for models not trained on PatBench (every group is test)",
         lambda c: [c.out], reads=("tasks/*.parquet",), writes=("eval_all/*.parquet",)),
    Step("baselines", "10_evaluation/baselines_fixed.py",
         "reference baselines: random, heavy atoms, logP, ECFP4 + random forest (zero-shot and few-shot)",
         lambda c: [c.out, c.out], reads=("tasks/*.parquet", "protocols/*.parquet"), writes=("baselines/*.parquet",),
         notes="the corrected copy of baselines.py: the zero-shot forest on the few-shot queries is written as "
               "ecfp_rf_zero_k<K> instead of overwriting the few-shot forest"),
    Step("baselines-r6-asrun", "10_evaluation/baselines.py",
         "the baselines script exactly as it ran for r6 (KNOWN DEFECT: ecfp_rf k5/k10 files hold zero-shot predictions)",
         lambda c: [c.out], reads=("tasks/*.parquet", "protocols/*.parquet"), writes=("baselines/*.parquet",),
         notes="kept for provenance; use `baselines`"),
    Step("leaderboard", "10_evaluation/leaderboard.py",
         "score the reference baselines and write baselines/leaderboard.csv",
         lambda c: [c.out], reads=("baselines/*.parquet",), writes=("baselines/leaderboard.csv",)),

    # ---- structure repair (stage 07); input is the release built BEFORE the repair ----------------------
    Step("alias-scan", "07_structure_repair/alias_scan.py",
         "first scan: stored structures that rest on a MOL abbreviation RDKit could not expand",
         lambda c: [str(c.workers)], reads=("raw archives", "first-release sar_all.parquet"),
         writes=("<first release>/work/alias_scan.parquet",), heavy=True, needs_out=False,
         notes="opens the raw archives; writes into the first release's work/ directory"),
    Step("alias-extract", "07_structure_repair/alias_extract.py",
         "collect the MOL files and an alias inventory of every record flagged structure_alias_suspect",
         None, needs_out=False, heavy=True,
         notes="arguments after `--`: <pre-repair release dir> <alias_dir> [workers]"),
    Step("alias-scan-all", "07_structure_repair/alias_scan_all.py",
         "alias inventory of EVERY structure file behind the pre-repair release",
         None, needs_out=False, heavy=True,
         notes="arguments after `--`: <pre-repair release dir> <alias_dir> [workers]"),
    Step("alias-pubkeys", "07_structure_repair/alias_pubkeys.py",
         "molecules BindingDB / ChEMBL curated from the same patents (reference for the repair validation)",
         None, needs_out=False,
         notes="arguments after `--`: <alias_dir> [records.tsv] [out.tsv]"),
    Step("alias-repair", "07_structure_repair/alias_repair.py",
         "first implementation of the alias expansion (label catalogue, placeholder atom rewritten in place)",
         None, needs_out=False,
         notes="arguments after `--`: <alias_dir> <records.tsv> <mols_dir> <out.tsv> [workers]"),
    Step("alias-expand", "07_structure_repair/alias_expand_independent.py",
         "second, independent implementation of the alias expansion (the one the repair uses)",
         None, needs_out=False,
         notes="reads records.tsv and mols/ beside the script: copy it into <alias_dir> as expand.py and run it "
               "there; see docs/REPRODUCE.md"),
    Step("alias-validate", "07_structure_repair/alias_validate2.py",
         "check a repair table against the molecules BindingDB / ChEMBL curated from the same patent",
         None, needs_out=False,
         notes="arguments after `--`: <alias_dir> <records.tsv> <repairs.tsv> <public_keys.tsv> <report_prefix> [pre-repair release dir]"),
    Step("alias-merge", "07_structure_repair/merge_final.py",
         "turn the expander output into the repair table (adds the nitrile step)",
         None, needs_out=False,
         notes="arguments after `--`: <alias_dir> <expansions.tsv> <records.tsv> <out.tsv>"),
    Step("patch-sar", "07_structure_repair/make_patched_sar.py",
         "write a copy of the recovery results with the repaired structures (SAR_merged.parquet + 3 repair columns)",
         None, needs_out=False, heavy=True,
         notes="arguments after `--`: <pre-repair release dir> <repairs.tsv> <out_results_dir>"),
    Step("alias-scan-final", "07_structure_repair/make_alias_scan_r6.py",
         "work/alias_scan.parquet of the repaired build: pairs whose structure still rests on an unexpanded alias",
         None, needs_out=False,
         notes="arguments after `--`: <pre-repair release dir> <new release dir> <alias_dir> <all_records.tsv> <repairs.tsv>"),
    Step("list-unresolved", "07_structure_repair/list_unresolved_labels.py",
         "every unresolved alias label with its context (structure_repair/unresolved_labels.tsv)",
         None, needs_out=False,
         notes="arguments after `--`: <alias_dir> <unresolved_records.tsv> <out.tsv>"),

    Step("nitrile-inventory", "07_structure_repair/nitrile_audit/nitrile_inventory.py",
         "audit of the reversed-nitrile conversion, 1/6: the structure files and archives behind every record with the `nitrile` rule",
         None, needs_out=False,
         notes="no arguments; paths are fixed inside the script (workspace-relative); writes inventory.parquet to outputs/nitrile-evidence-audit-20261031"),
    Step("nitrile-extract", "07_structure_repair/nitrile_audit/nitrile_extract.py",
         "audit 2/6: the N(+)#C(-) motifs of every file with their labels, bond geometry and the original image",
         None, needs_out=False, heavy=True, notes="no arguments; writes motif_evidence.parquet and originals/"),
    Step("nitrile-text", "07_structure_repair/nitrile_audit/nitrile_text.py",
         "audit 3/6: the patent text around each structure image (headings, chemical-name flags)",
         None, needs_out=False,
         notes="no arguments; reads motif_evidence_geometry.parquet, a table whose producing step is not among the scripts"),
    Step("nitrile-panels", "07_structure_repair/nitrile_audit/nitrile_panels.py",
         "audit 4/6: image panels with the motif marked, for the visual review; writes visual_review_index.parquet",
         None, needs_out=False, notes="no arguments; the visual review itself is not scripted"),
    Step("nitrile-finalize", "07_structure_repair/nitrile_audit/nitrile_finalize.py",
         "audit 5/6: decisions per motif, record and row; the derived sar_all_reviewed.parquet; the quarantined records",
         None, needs_out=False, notes="no arguments; writes into <r6 release>/structure_repair/nitrile_audit"),
    Step("nitrile-verify", "07_structure_repair/nitrile_audit/nitrile_verify.py",
         "audit 6/6: checks the decision tables with the evidence guard (verified_nitrile_guard.py)",
         None, needs_out=False, notes="no arguments"),

    # ---- target audit (after the build) ------------------------------------------------------------
    Step("join-target-map", "08_build/join_target_map.py",
         "join each ranking-task group to the target-map pair it came from (input of target-audit)",
         lambda c: [c.out], reads=("sar_all.parquet", "tasks/*.parquet", "target_map.parquet"),
         notes="pass the map and output paths after `--`: <target_map.parquet> <join output .parquet>"),
    Step("target-audit", "08_build/make_target_audit.py",
         "target_audit.parquet: the join plus, per ranking task, membership, time split and test tier",
         lambda c: [c.out], writes=("target_audit.parquet",),
         notes="pass `<join output .parquet> <release>/target_audit.parquet` after `--`"),
    Step("audit-targets", "08_build/audit_targets.py",
         "per-group target sequences (polyprotein chains, mutants) and the label / cross-reference checks",
         lambda c: [], network=True, needs_out=False,
         writes=("target_sequences.parquet", "sequence_fixes.tsv"),
         notes="pass `<release>/target_audit.parquet <release>/targets.parquet <release>` after `--`; needs UniProt"),

    # ---- structures ---------------------------------------------------------------------------------
    Step("fetch-structures", "09_structures/fetch_structures.py",
         "download and clean the PDB / AlphaFold structure of every target; build the group -> structure table",
         lambda c: [c.out], network=True, heavy=True,
         writes=("structures/*",), notes="then a stage name after `--`: meta entries gql select download process groups write"),
    Step("check-structures", "09_structures/check_structures.py",
         "independent re-check of structures/: re-reads the files and recomputes what the tables claim",
         lambda c: [c.out], reads=("structures/*",), notes="optional arguments after `--`: [seed] [outdir]"),

    # ---- support tables, review and profile -------------------------------------------------------------
    Step("eval-support", "11_support/build_eval_support.py",
         "domain sequences per group, test-group profile (seen target, docking eligibility, mapping status)",
         lambda c: [c.out], reads=("tasks/*.parquet", "structures/*", "target_sequences.parquet"),
         writes=("eval_support/*",)),
    Step("target-names", "11_support/build_target_names.py",
         "what each group's target IS (accession, gene, protein, organism, mutation) next to what the patent printed",
         lambda c: [c.out], writes=("eval_support/group_target_names.*",)),
    Step("review-table", "11_support/make_review_table.py",
         "manual-review tables of the time-split test groups (reviewer columns left blank)",
         lambda c: [c.out], reads=("target_audit.parquet", "target_sequences.parquet", "tasks/*.parquet"),
         notes="pass the output directory after `--`: <out_dir>"),
    Step("review-page", "11_support/review_page/make_review_page.py",
         "self-contained HTML page for the manual review of the target mappings of every test group",
         lambda c: [c.out], reads=("target_review_test.tsv", "target_audit.parquet", "eval_support/*"),
         writes=("review/target_review.html",)),
    Step("profile", "11_support/profile_r6.py",
         "profile of a release: size, tasks, splits, tiers, molecular weight, remaining problems",
         lambda c: [c.out], notes="pass the output directory after `--`: <out_dir>"),
    Step("profile-extra", "11_support/profile2_r6.py",
         "second part of the profile: molecular weight against activity, over time, extremes",
         lambda c: [c.out], notes="pass the output directory after `--`: <out_dir>"),

    # ---- evaluation and verification ------------------------------------------------------------------
    Step("evaluate", "10_evaluation/evaluate_patent_bootstrap.py",
         "score a prediction file on a ranking task; intervals resample whole patents, paired comparison with --reference",
         None, needs_out=False,
         notes="the evaluator named as recommended in the r6 freeze record; arguments after `--`: <task.parquet> "
               "<predictions> [--split ... --part ... --tier ... --reference ... --profile ... --json-out ...]"),
    Step("evaluate-groups", "10_evaluation/evaluate.py",
         "build-time evaluator used by `leaderboard`: Spearman, Kendall, pair accuracy, NDCG@10, top-10% hit; intervals resample groups",
         None, needs_out=False, notes="arguments after `--`: <task.parquet> <predictions> [--split ... --tier ...]"),
    Step("verify-provenance", "12_verify/provenance_check.py",
         "trace a random sample of benchmark rows back to the original patent archive",
         lambda c: [c.out], reads=("sar_all.parquet", "raw archives"), heavy=True,
         notes="pass the sample sizes after `--`: <n_test_rows> <n_train_rows>"),
]}

RELEASE_CHAIN = ["standardize", "novelty", "build", "protocols", "eval-sets", "baselines", "leaderboard"]

# (number, directory, what it does, how to run it)
STAGES = [
    ("00", "pipeline/00_ops", "one-off operations: copy the patent archives; launch the recovery service",
     "site specific, kept as run"),
    ("01", "pipeline/01_public_data", "standardise BindingDB / ChEMBL / HARVEST / SureChEMBL and merge them",
     "python pipeline/01_public_data/build_standardized.py ... (see its README.md)"),
    ("02", "pipeline/02_extract", "SAR tables from the USPTO MOL archives: 345,940 ZIPs -> formal SAR + evidence",
     "python pipeline/02_extract/run_full.py --root ... --manifest ... --output ... --full"),
    ("03", "pipeline/03_recovery", "recover activity links the first pass skipped; merge into SAR_merged",
     "python pipeline/03_recovery/run_recovery.py ... (code frozen; four path constants, see docs/REPRODUCE.md)"),
    ("04", "pipeline/04_catalog", "catalogue every MOL of every archive and the patent metadata (all domains)",
     "python pipeline/04_catalog/run_structure_catalog.py ... (read by the recovery stage; not by the build)"),
    ("05", "pipeline/05_target_map", "assay label -> UniProt accession for every (patent, label) pair (v1, v2, v2.1)",
     "v2_universe -> v2_candidates -> v2_make_batches -> (curation) -> v2_read_curation -> v2_assemble -> v2_xref -> "
     "v2_finalize -> v2_1/apply_target_overrides"),
    ("06", "pipeline/06_validation", "the manual validation sample, annotation pack, agreement and failure analysis",
     "build_sample -> (annotation) -> merge_annotations -> analyze_agreement ..."),
    ("07", "pipeline/07_structure_repair", "repair structures stored wrongly from MOL alias lines (NH, groups, nitriles)",
     "patbench run alias-extract ... (see docs/PIPELINE.md)"),
    ("08", "pipeline/08_build", "molecule standardisation, public novelty, the PatBench build, protocols, target audit",
     "patbench release --out <dir>"),
    ("09", "pipeline/09_structures", "protein 3D structures for every target, and their independent check",
     "patbench run fetch-structures --out <dir> -- <stage>"),
    ("10", "pipeline/10_evaluation", "evaluation script, reference baselines, leaderboard", "patbench run evaluate -- ..."),
    ("11", "pipeline/11_support", "evaluation-support tables; the target-mapping review page; release profile",
     "patbench run eval-support --out <dir>"),
    ("12", "pipeline/12_verify", "provenance check against the original patent archives", "patbench run verify-provenance --out <dir>"),
    ("--", "pipeline/legacy_v0", "the first release and benchmark builders PatBench grew out of", "history only"),
]
