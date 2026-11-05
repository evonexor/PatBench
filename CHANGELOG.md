# Changelog

All notable changes to this repository are recorded here. The data releases have their own history in
[docs/RELEASES.md](docs/RELEASES.md).

## Unreleased

- MIT License added (`LICENSE`).
- Absolute paths of the build server and a personal user name removed from documentation, comments, the import tool and `provenance/NOT_INCLUDED.tsv`. The four frozen directories (`00_ops`, `02_extract`, `03_recovery`, `04_catalog`) are unchanged: their hashes are part of the record, so they still contain the paths they ran with (see `docs/CODE_PROVENANCE.md`).

Initial standalone repository (snapshot of 2026-11-02), corresponding to PatBench v1.0 revision r6 (frozen 2026-11-01).

- Processing code of all stages imported from the server directories with `tools/import_pipeline.py`: public-data standardisation, extraction, recovery, structure catalogue, target mapping (v1, v2, v2.1), validation, structure repair, build, structures, evaluation, support tables, provenance check, and the legacy first-release builders.
- `provenance/CODE_MANIFEST.tsv` with the original and repository SHA-256 of every pipeline file; `provenance/STAGE_HASHES.json` with the code hashes recorded by the three frozen stages (they reproduce); `provenance/NOT_INCLUDED.tsv` with the server scripts left out and why.
- The evaluator named as recommended in the r6 freeze record (`10_evaluation/evaluate_patent_bootstrap.py`, SHA-256 checked by a test) and the corrected few-shot baseline script (`10_evaluation/baselines_fixed.py`) are included; `patbench release` uses the corrected script.
- Command line `patbench`: `stages`, `steps`, `plan`, `check-workspace`, `run`, `exec`, `release`, `verify-code`, `verify-release`. The runner sets `PATBENCH_SAR_DIR` and `PATBENCH_TARGET_MAP_DIR` to the repaired tables and the v2.1 target map.
- Stage numbering: the structure repair is stage `07`; the build, structures, evaluation, support and verification stages are `08` to `12`.
- Documentation in English only.

### Verification of this snapshot (2026-11-02)

- `patbench verify-code`: all 153 files of `provenance/CODE_MANIFEST.tsv` are consistent; the three recorded stage hashes (`148c103fb46f`, `81fd68912e6f`, `67b554fd4599`) are reproduced from the repository files; no server path is left in the rewritten code.
- `pytest`: 110 passed and 1 skipped on a laptop (Python 3.12, no RDKit); 116 passed and 1 skipped on the original server (Python 3.10, RDKit 2025.09, pyarrow 19), which includes the unit tests shipped with the extraction, recovery and catalogue stages. The skipped test is the one that needs Python 3.12 f-string tokenisation.
- `patbench verify-release` on the frozen r6 release directory: 62 of 62 checks passed (files present, parameters equal to `configs/release_v1.0.json`, one value per (group, molecule), one split per group, no patent or family crossing train / valid / test under `split_time` and `split_random`, no (molecule, target) pair shared by train and test, test rows respect the priority-year cutoff and the grant-year floor, counts equal the manifest).
- Not done: a rebuild of r6 with these scripts compared with the frozen release (see `docs/REPRODUCE.md`).
