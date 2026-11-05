# Release history

Every release directory is kept under the same parent directory. The reference release of this repository is **r6**. Dates and numbers are taken from the README and `release_manifest.json` of each release. The suffix `-r2` to `-r6` marks the revision of the version `patbench-v1.0-20261030`. The date in a directory name need not equal the build date: r5 and r6 carry `20261030` in their names and were built on 2026-10-31.

| Directory | Date | Content |
|---|---|---|
| `affinity-rank-benchmark-20261019/` | 2026-10-19 | First ranking benchmark (`pipeline/legacy_v0/benchmark/`) |
| `affinity-sar-release-20261019/` | 2026-10-19 | First SAR release: all formal SAR with quality and split flags, before curation (`pipeline/legacy_v0/build_release.py`) |
| `patbench-v1.0-20261020/` | 2026-10-20 | PatBench v1.0 original: industrial-grade curation on top of the pre-curation release (13 items in the README of that release) |
| `patbench-v1.0-20261030/`, `-r2/` | 2026-10-30 | Intermediate releases; the priority-year start of `split_time` was 2018, and the test sets were small (`rank_ic50` 84 groups, `rank_ki_kd` 10 groups) |
| `patbench-v1.0-20261030-r3/` | 2026-10-30 | Start changed to 2016, larger test sets (`rank_ic50` 178 groups, `rank_ki_kd` 35 groups); no grant-year restriction; r4 reused its 3D structure directory |
| `patbench-v1.0-20261030-r4/` | 2026-10-30 | Test sets restricted to patents granted in 2020 or later (`rank_ic50` 166 groups, `rank_ki_kd` 33 groups); three first-public-date columns added; relative to r3 only `split_time` and `test_tier` changed in `sar_all.parquet`, and rows already in train or valid are unchanged |
| `patbench-v1.0-20261030-r4/eval_support/` | 2026-10-31 | Add-only evaluation helper tables on r4: domain sequences of the structures, test-group profile, review table for 39 test groups of other tasks |
| `patbench-v1.0-20261030-r5/` | 2026-10-31 | r4 plus three UniProt accession corrections: NDM-1 (`G4WHZ0` to `C7C422`), VIM-1 (`Q8GRG4` to `Q9XAY4`) and Mtb LeuRS (`P9WFV3` to `P9WFV1`). Test sets unchanged. One group (5332, a duplicate series after the correction) removed: `rank_ic50` 62,214 rows / 627 groups, `rank_highcontrast` 74,663 rows / 673 groups; regression task 94,335 to 94,306 rows; `targets.parquet` 1,088 to 1,087 accessions. Baselines were not re-run. Input: `target-uniprot-map-v2.1-20261031/` |
| `patbench-v1.0-20261030-r6/` | built 2026-10-31, frozen 2026-11-01 | r5 plus **structure repair**: the structures of 128,062 rows (21.2%) were re-derived from MOL alias lines; three columns added to `sar_all.parquet` (`structure_repair`, `SMILES_before_repair`, `InChIKey_before_repair`); 1,953 unresolved records keep `structure_alias_suspect`. Test groups: `rank_ic50` 172, `rank_ki_kd` 33, `rank_highcontrast` 182, `rank_complex_mutant` 29, `rank_ec50` 9. Baselines re-run on the four ranking tasks. Input: `uspto-sar-recovery-20261031-structfix/results/SAR_merged.parquet` (upstream copy; the original directory is unchanged) |

## r6, the reference release

- **Frozen on 2026-11-01.** The directory is read-only. `SHA256SUMS` lists the SHA-256 of all 5,788 files (check with `sha256sum -c SHA256SUMS`), the digest of that list is in `SHA256SUMS.sha256`, and the freeze record is `FROZEN.json`. Any later revision has to be a separate release directory. Five input files of the frozen evaluation run (`tasks/rank_ic50.parquet`, `tasks/rank_ki_kd.parquet`, `target_sequences.parquet`, `structures/group_structure_map.parquet`, `structures/target_structures.parquet`) are byte-identical to the frozen release.
- **r6 differs from r5 in its test sets, structures and public-database flags.** Results computed on r5 or r4, including the baseline table of r4, cannot be compared with results on r6.
  - 6,390 rows have a different `family_id` (family membership changed with the stored structures).
  - Test groups: `rank_ic50` 1 group left and 7 entered; `rank_highcontrast` 1 left and 5 entered; the test groups of `rank_ki_kd`, `rank_complex_mutant` and `rank_ec50` are unchanged. Between train and valid, 105 `rank_ic50` groups swapped (the families changed, so the random assignment was redrawn). Of all 13,122 groups, 9 moved from test to train, and none moved from train or valid to test.
  - Counts: parent molecules 214,385 to 211,734; patent families 4,581 to 4,522; rows with `structure_alias_suspect` 15,736 to 3,678; rows with `in_public_data` 252,985 to 321,741 (41.9% to 53.3%); regression task 94,306 to 95,376 rows.
- **Baselines.** The random, heavy-atom, logP and ECFP4 + random forest baselines were re-run on r6 (2026-10-31) for `rank_ic50`, `rank_ki_kd`, `rank_highcontrast` and `rank_complex_mutant` with the scripts of r4 (`scripts/baselines.py`, `scripts/leaderboard.py`, run by `run_baselines_r6.sh`). The eight few-shot ECFP4 + random forest prediction files were invalidated at the freeze (see [LIMITATIONS](LIMITATIONS.md)).
- **Supplement of 2026-10-31.** Review of the nitrile conversion against the original drawings (`structure_repair/nitrile_audit/`, see [VALIDATION](VALIDATION.md)); the main `sar_all.parquet` and the tasks are unchanged, a derived table flags three additional records.
- `README_r5.md` in the r6 directory is the full description of r5. Its sections on construction, splits and technical validation still apply to r6, but its numbers are those of r5. For r6, the r6 README and `release_manifest.json` are authoritative.

## Directories that are not releases of this repository

The directory `patbench-v1.0-20261031-r7` exists next to r6. It has the same rows as r6, but the split of 112 groups differs. It is **not produced by this repository** and is not described in these documents; results computed on it must not be mixed with results on r6.

## Status

No release has been deposited with a DOI. The code is under the MIT License (see `LICENSE`). Splits and task definitions are described in [SPLITS](SPLITS.md).
