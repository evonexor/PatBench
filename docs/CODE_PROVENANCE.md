# Code provenance and integrity

The processing code in this repository was not rewritten. It was **copied from the server that produced the data**, in
several directories, one per stage. The reason is simple: the data are the output of that code, and a rewrite could not
show that the repository contains the code that produced them.

## Three kinds of files

`provenance/CODE_MANIFEST.tsv` has one row for every file under `pipeline/`: the repository path, the source path on the
server, the status, the number of path edits, the SHA-256 of the original and the SHA-256 of the copy. The manifest lists
153 files: 83 verbatim, 65 patched and 5 authored. `patbench verify-code` requires every imported file to match its
recorded hash and every file under `pipeline/` to be registered; for `authored` files it only requires that they exist and are
registered, because they are edited.

| Status | Meaning | Files |
|---|---|---|
| `verbatim` | copied byte for byte; original hash = repository hash (83 files) | all of `00_ops`, `02_extract`, `03_recovery` and `04_catalog`, and the scripts of the other directories that never named a server path |
| `patched-paths` | only the way a path is written was changed, see below (65 files) | scripts in `01_public_data`, `05_target_map`, `06_validation`, `07_structure_repair` (including `nitrile_audit/`), `08_build`, `09_structures`, `11_support` and `legacy_v0` that named a server path |
| `authored` | written in this repository (5 files) | `11_support/review_page/make_review_page.py` and `template.html` (the target-mapping review page), `12_verify/provenance_check.py`, `01_public_data/README.md` and `pipeline/README.md` |

## How paths were rewritten

The scripts on the server contain 148 string literals that begin with the absolute path of the server's workspace or, in `01_public_data`, with the absolute path where the
public-database scripts were first written. The two prefixes are given to the import tool through the environment variables `PATBENCH_IMPORT_SERVER_PREFIX` and `PATBENCH_IMPORT_PUBLIC_PREFIX`; their values are not recorded in this repository. Each literal was rewritten to `_ws + "<relative path>"`, where `_ws` is
defined once near the top of the file:

```python
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
```

If the variable is not set the script stops immediately instead of reading from somewhere else. The variable is
deliberately lower case: `build_patbench.py` writes every upper-case global into the `config` of `release_manifest.json`,
so an upper-case name would add a key to a rebuilt manifest. Two special cases: in the target-mapping scripts, the path that
put "the directory of the script itself" on `sys.path` became relative to the file (`v1` and `v2` sit side by side), and in
`v2_assemble.py` the expression `V1 + 'scripts'` (data directory plus `scripts`) now points to the `v1` directory of the repository.

The rewriting is done by `tools/import_pipeline.py`, which processes string literals with Python's tokenizer and **does
not touch comments or docstrings**. `tests/test_import_tool.py` tests its behaviour, and `patbench verify-code` checks that
no server path is left in the result. Two files also had a comment or docstring changed to remove an absolute path (`07_structure_repair/alias_expand_independent.py`, `07_structure_repair/nitrile_audit/nitrile_panels.py`); they are marked `patched-paths`. To repeat the import:

```bash
python tools/import_pipeline.py <snapshot of the server directories>
```

## Code hashes of the three frozen stages

When the extraction, the recovery and the structure catalogue ran, each wrote the SHA-256 of its code directory into a
status file (`job_config.json`, `completion.json`). `provenance/STAGE_HASHES.json` stores these three values, and
`patbench verify-code` **recomputes** them from the files in this repository with the formula used at that time:

| Directory | Formula | Recorded hash (first 12 digits) | Original run |
|---|---|---|---|
| `02_extract` | the bytes of the `*.py` files of the directory, sorted by file name, concatenated | `148c103fb46f` | uspto-sar-full-20261004, 0.8.0-full |
| `03_recovery` | all `*.py` files recursively; for each, relative path + `\0` + bytes, sorted by path, concatenated | `81fd68912e6f` | uspto-sar-recovery-20261006 |
| `04_catalog` | as `02_extract` | `67b554fd4599` | uspto-all-structures-20261004, 0.9.1-sharded-aggregation |

All three match. The code in these directories is therefore byte-identical to the code that produced `SAR_merged`.
**These three directories must not be edited**, including their hard-coded paths (see [REPRODUCE](REPRODUCE.md)).
Metadata files that macOS leaves behind (`._name`) are ignored by the check; the original directories contained none.

## Language of the imported code

The imported files are kept as they were. Chinese-language text therefore remains in a few of them: in the frozen stages
(`patent_metadata.py`, `patent_documents.py`, `build_review_html.py`, `prepare_review.py`, `catalog_metadata.py` and the
unit test `test_metadata.py`), where it is part of the code whose hash was recorded, in two note strings of
`11_support/build_target_names.py`, and in the text that `07_structure_repair/nitrile_audit/nitrile_finalize.py` writes into its audit tables. In the frozen stages, some descriptive metadata values written to the `patents` tables
(for example the `indication` and `technical_field_source` columns) are Chinese labels. The PatBench build reads only the date columns of
those tables. Everything written for this repository (documentation, tests, tools, the review page) is in English.

## Scripts taken from the evaluation runs

Two scripts in `10_evaluation/` do not come from the build directory: `evaluate_patent_bootstrap.py`, the evaluator that the freeze record of r6 names as recommended (patent-cluster bootstrap, paired comparison; SHA-256 `35b70dd019dc8da9...`, recorded in `FROZEN.json` and checked by `tests/test_evaluate_patent_bootstrap.py`), comes from the frozen evaluation run; `baselines_fixed.py` comes from the corrected few-shot baseline run (`patbench-r6-fewshot-baselines-20261101`). The other files of those directories (the model-run harness and the scoring helpers) are listed in `provenance/NOT_INCLUDED.tsv`.

## Scripts that have no runtime hash

The scripts from `07_structure_repair` onwards (repair, build, evaluation, structures) did not record a code hash when they
ran. Their evidence is the original hash in `CODE_MANIFEST.tsv` (for the build scripts, taken from the `scripts/` directory of the r6
release on the server, which is read-only and covered by that release's checksum list) and the checks listed in [VALIDATION](VALIDATION.md).

## Scripts deliberately not included

`provenance/NOT_INCLUDED.tsv` lists the server scripts that were not imported, with their SHA-256 and the reason:

- release housekeeping that writes text into the data release README or manifest, or freezes the release directory (`freeze_r6.py`, `freeze_r6.sh`, `note_unresolved_r6.py`, `stamp_baselines_r6.py`, `stamp_revision.py`);
- superseded or exploratory scripts of the alias work (`alias_expand.py`, `alias_validate.py`, `make_expand_v7.py`, `draw_labels.py`); the rule extension tried with `make_expand_v7.py` was not merged into r6;
- `carry_over_target_files.py`, which reused the target files of an earlier release and is not used by r6.

## Code of uncertain origin

About 1,200 lines of `09_structures/fetch_structures.py` were produced in a structure-download task that was carried out by a
sub-agent and resumed part-way. For the few lines that assign `META` and the segment map (around lines 657-660), **it was not
established who wrote them or when**: the main session and another session in the same workspace each said they had not
changed them, and the file's modification time falls inside the period in which the sub-agent resumed. Two things were done to
compensate:

1. **The script was read.** It contacts only `alphafold.ebi.ac.uk`, `data.rcsb.org`, `files.rcsb.org`, `rest.uniprot.org` and `www.ebi.ac.uk`; it uses no `subprocess`, `eval` or `exec`; it reads `pickle` only for its own `api_cache` files; it writes only below `structures/`.
2. **Its output was checked independently** (re-read with gemmi, not using the intermediate results of the run): the C-alpha coordinates of sampled PDB structures equal those of the original mmCIF files; the hashes of the raw files equal those of freshly downloaded copies; the entry metadata agree with the online RCSB and PDBe records; the sequences of the recommended structures agree with the group sequences. See [VALIDATION](VALIDATION.md).

This gap is recorded here as it is. If the structure tables are used for a formal publication, the step should be reviewed
once more as a stand-alone download script that can be re-run.

## Code and instructions not found on the server

These were looked for and not found; they are not in this repository:

- the curation prompts, instructions and model versions used for the target mapping (the curation results are data, in the `work/` directories of `target-uniprot-map-*`);
- the model and the instructions used for the second, independent implementation of the alias expansion (the script itself is `07_structure_repair/alias_expand_independent.py`).

## Code found but incomplete

- The scripts of the evidence audit of the reversed nitriles are included (`07_structure_repair/nitrile_audit/`; they were stored in `structure_repair/nitrile_audit/` of the r6 release). Two parts of the audit are **not scripted**: the step that writes `motif_evidence_geometry.parquet` (the bond geometry of the motifs), which `nitrile_text.py` and `nitrile_finalize.py` read, and the visual review of the image panels, whose result is `visual_review_index.parquet`. The scripts are short one-line-style session scripts with fixed paths; they were not written as reusable tools.
- `10_evaluation/baselines.py` is the r6 script as it ran and has a known defect (the few-shot ECFP4 + random forest predictions are overwritten by the zero-shot ones, see [LIMITATIONS](LIMITATIONS.md)); `10_evaluation/baselines_fixed.py` is the corrected copy.

## What is not in the repository

- Data: the archives, the outputs of the stages and the release directories.
- The curation results of the target mapping (the curation TSVs of v1 and v2): judgements written by the curators after reading the patent text. They are data and live in `target-uniprot-map-*/work/` of the workspace.
- The annotation files of the validation sample (`sar-validation-20261019/annotations/`) and the adjudication results.
- systemd service units, cron jobs and machine-specific launch configuration.
