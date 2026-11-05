# Splits and task sets

The numbers in this file describe the reference release r6 (built 2026-10-31, frozen 2026-11-01). They are taken from the r6 `release_manifest.json`, the r6 README, the r6 build logs (`build.log`, `chain.log`) and, where stated, `eval_support/README.md`. The split parameters are in `configs/release_v1.0.json`; they are identical to the `config` block of the r6 `release_manifest.json` (which additionally records `BUILT` and `VAL`) and to the constants at the top of `build_patbench.py`.

A statement that was measured on release r4 only is marked "(verified on r4; not repeated on r6)". The r6 test sets differ from those of r5 and r4 (see [RELEASES](RELEASES.md)), so counts of the earlier releases are not comparable with the counts below.

## Split unit

The split unit is the **connected component of patents and patent families**: 6,397 patents and 4,522 patent families form 4,405 components. Under `split_time` and `split_random`, rows of the same patent or the same family never straddle train, valid and test. A group is a patent × assay label × endpoint; all rows of a group lie in the same part under these two splits.

Patent families are a heuristic (same endpoint and molecule overlap of at least 50%, `FAMILY_OVERLAP` = 0.5), not an INPADOC lookup. Because the heuristic compares stored structures, the structure repair of r6 changed the family membership of some rows (see [RELEASES](RELEASES.md)).

## Five splits (one column each in `sar_all.parquet`)

| Column | Rule | Use | Rows train / valid / test (r6) |
|---|---|---|---|
| `split_time` | **Main split.** A component is test if its earliest priority year is 2016 or later **and** every patent in the component was granted in 2020 or later. Of the remaining components, a random 15% is valid and the rest is train. | Formal split for the ranking tasks | 364,238 / 73,879 / 164,648 (1,143 rows not assigned) |
| `split_time_grant` | Control: component whose earliest grant year is 2021 or later is test (rule of an earlier release). | Control | 393,720 / 81,027 / 129,010 (151 rows not assigned) |
| `split_random` | Random 70/15/15 by component. | Control; tasks with a small test set such as `rank_ki_kd` should report it as well | 422,749 / 90,634 / 90,525 |
| `split_scaffold_mol` | By molecular scaffold. | Regression tasks (recommended) | 444,856 / 79,185 / 79,867 |
| `split_scaffold_group` | By the dominant scaffold of the group. | Reference only; patents are **not** guaranteed to be disjoint | 400,463 / 86,001 / 85,302 (32,142 rows not assigned) |

The random seed is `20261019`. For `split_scaffold_group`, 5 patents of `rank_ic50` and 3 patents of `rank_ki_kd` appear in both train and test (verified on r4; not repeated on r6).

## Why priority date and not grant year

Compounds are synthesised around the **priority date**. US applications are published 18 months after filing, and grants often follow 3-5 years later (grant lag: median 5 years, p10 3 years, p90 8 years; verified on r4; not repeated on r6). On the grant-year split of an earlier release, two thirds of the test chemistry had already been filed before 2019 (measured on that earlier release; not repeated on r6). `split_time` therefore judges the age of the chemistry by the priority year and adds a grant-year floor, which restricts the test set to patents granted in 2020-2025.

The priority date is taken, in this order, from the earliest priority claim (5,747 patents), the PCT filing date (85) or the US filing date (142). For 423 patents none of the three exists (mainly old patents granted in 2001-2004; description from the r4 documentation, the count is unchanged in r6); for these the grant year serves as an upper bound. An upper bound before 2016 places the component in train or valid. Units whose date remains unknowable (7 units, 1,143 rows in r6) are excluded from the time split and have an empty `split_time`.

**Grant-year floor.** A component that is test by priority year but contains a patent granted before 2020 is not test; it joins train or valid, with its own 15% valid draw (seed `(SEED, 7)`). The r6 build log reports 1,148 units that are test by priority year alone and 102 units moved to train/valid because a patent in the unit was granted before 2020. Because the unit is a connected component, a moved unit shares no patent and no family with the retained test patents.

**`rank_complex_mutant` exception.** The task has little data: with the common 2016 start the test set is larger than the training set (24 / 8 / 40 groups for train / valid / test; verified on r4; not repeated on r6). The task therefore uses its own start, priority year 2018 (`TASK_TEST_FROM`), together with the same grant-year floor. In r6 this gives 30 / 13 / 29 groups. In `sar_all.parquet`, `split_time` always uses 2016; only `split_time` and `test_tier` in `tasks/rank_complex_mutant.parquet` use 2018. The test set of `rank_complex_mutant` therefore cannot be concatenated with the test sets of the other four tasks.

**First public date** (`first_public_date`, `first_public_year`, `first_public_source`): the earliest of the US application publication date, the PCT publication date and the grant date. It is known for 5,992 of 6,397 patents (sources: US application 3,692; PCT 2,204; grant date 96). For 405 patents none of the three exists; the field stays empty and is not estimated. The field is used for description and stratification only and **not for splitting**: continuations and divisionals are published late although their chemistry was disclosed through the parent.

## Test tiers (`test_tier`)

Defined only for groups with `split_time == test`:

- `novel`: no row of the group is in BindingDB or ChEMBL (the (molecule, target) pairs of the group are all absent from both).
- `reference`: all other groups.
- `unscored`: the group has no target mapped at high or medium confidence.

Among the groups with `split_time == test` in `sar_all.parquet` (common 2016 start), the r6 manifest records 396 `novel`, 920 `reference` and 1,575 `unscored` groups. **A model pretrained on ChEMBL or BindingDB should be reported on the `novel` tier for the main comparison.** BindingDB curates US patents itself, so a model pretrained on it may have seen the same chemistry. In r6, 82-84% of the test rows of `rank_ic50`, `rank_ki_kd`, `rank_highcontrast` and `rank_ec50` (65.3% for `rank_complex_mutant`) are already in BindingDB or ChEMBL, which makes the `novel` tier small (see [LIMITATIONS](LIMITATIONS.md)).

## Task sets

Common conditions for all tasks: exact numeric values, drug-like molecule, not `replicate_disagrees`, and all eight hold-out flags false (`group_is_fragment_suspect`, `pair_under_adjudication`, `pair_validated_error`, `patent_has_validation_issue`, `structure_alias_suspect`, `group_unit_suspect`, `group_is_duplicate`, `group_label_uninformative`). One value per (group, molecule). One group per series. Each patent contributes at least 20 distinct molecules in each task.

The ranking tasks add: biochemical assay, at least 10 molecules in the group, dynamic range of at least 1 log, not graded, and `target_mapped` (`rank_complex_mutant` is exempt from the last condition). `rank_highcontrast` requires a within-group range of at least 2 log.

| Task | Rows | Groups | Patents | Molecules | Targets | Groups train / valid / test | Test groups in the `novel` tier |
|---|---:|---:|---:|---:|---:|---|---:|
| `rank_ic50` | 62,712 | 629 | 493 | 53,534 | 254 | 382 / 75 / 172 | 26 |
| `rank_ki_kd` | 18,358 | 189 | 129 | 14,238 | 87 | 138 / 18 / 33 | 6 |
| `rank_highcontrast` | 76,115 | 682 | 535 | 63,890 | 269 | 421 / 79 / 182 | 30 |
| `rank_complex_mutant` (start 2018) | 5,357 | 72 | 52 | 4,867 | 42 | 30 / 13 / 29 | 4 |
| `rank_ec50` | 2,215 | 23 | 15 | 1,601 | 12 | 10 / 4 / 9 | 1 |

Test sets of the ranking tasks (r6):

| Task | Test rows | Test patents | Test molecules | Test targets | Molecules per group, median (range) | `novel` tier (groups / patents / rows) | Test rows already in BindingDB or ChEMBL | Test rows with a repaired structure |
|---|---:|---:|---:|---:|---|---|---:|---:|
| `rank_ic50` | 16,003 | 136 | 13,903 | 96 | 61 (11-569) | 26 / 24 / 1,959 | 84.3% | 20.1% |
| `rank_ki_kd` | 3,718 | 21 | 2,695 | 26 | 44 (20-673) | 6 / 2 / 203 | 82.1% | 18.6% |
| `rank_highcontrast` | 19,102 | 137 | 15,280 | 103 | 71 (12-673) | 30 / 24 / 2,272 | 83.2% | 17.8% |
| `rank_complex_mutant` | 1,171 | 15 | 896 | 18 | 26 (13-164) | 4 / 4 / 329 | 65.3% | 44.0% |
| `rank_ec50` | 928 | 3 | 351 | 5 | 141 (43-155) | 1 / 1 / 153 | 83.5% | 0.4% |

### Regression task `regression_affinity`

One row is one (parent molecule, target, endpoint). Only biochemical assays are used, and repeated measurements are summarised by their median. The task has 95,376 rows, 79,444 molecules and 335 targets from 967 source patents (IC50 73,279; Ki 17,624; EC50 2,274; Kd 2,199). The recommended split is `split_scaffold_mol`: 67,364 / 11,358 / 16,654 (train / valid / test). The other splits give: `split_time` 60,629 / 10,802 / 23,945; `split_random` 68,344 / 14,991 / 12,041.

A combination that occurs in patents of several splits takes the **earliest** split (train before valid, valid before test), because it was already public at that point.

### Points to note

1. The test sets of `rank_ic50` and `rank_ki_kd` together contain 205 distinct groups (172 + 33). They come from 155 patents and cover 111 UniProt targets (counted on 2026-11-02 from the test rows of `tasks/rank_ic50.parquet` and `tasks/rank_ki_kd.parquet` of the frozen r6 directory; r4 had 199 groups, 151 patents and 107 targets). The `novel` tier of `rank_ki_kd` has only 6 groups from 2 patents, and the test set of `rank_ec50` comes from only 3 patents; both are suitable for exploration only.
2. Several groups of the same patent (for example kinase selectivity panels) share molecules and are strongly rank-correlated; they are not independent samples. Confidence intervals should resample patents. The build-time `evaluate.py` resamples groups (see [LIMITATIONS](LIMITATIONS.md)).
3. The test chemistry was mainly filed in 2016-2018 (priority year). A model pretrained on public databases is more likely to have seen similar chemistry. Results should be stratified by `test_tier` and `priority_year`.

## Models not trained on PatBench

The time-split test is only part of each task. A **model that is not trained on PatBench** (a public pretrained model, zero-shot scoring) does not need this hold-out, and every group can serve as test data. `make_eval_sets.py` writes `eval_all/<task>.parquet` with `split_all` equal to `test` everywhere; `test_tier` is recomputed over all groups (`public_share` equal to 0 gives `novel`, everything else `reference`), `test_tier_time` keeps the tier of the time split, and `priority_year` allows stratification. `novel` only guarantees that the (molecule, target) pairs are not in BindingDB and ChEMBL; a model trained on other sources has to be checked separately.

| Task | Groups | Patents | Targets | Molecules | `novel` groups | `novel` patents | `novel` targets | `novel` molecules |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `rank_ic50` | 629 | 493 | 254 | 53,534 | 189 | 155 | 114 | 12,547 |
| `rank_ki_kd` | 189 | 129 | 87 | 14,238 | 56 | 32 | 39 | 2,711 |
| `rank_ec50` | 23 | 15 | 12 | 1,601 | 8 | 8 | 7 | 958 |
| `rank_highcontrast` | 682 | 535 | 269 | 63,890 | 186 | 150 | 125 | 14,402 |
| `rank_complex_mutant` | 72 | 52 | 42 | 4,867 | 20 | 16 | 18 | 1,337 |
| `rank_ic50` + `rank_ki_kd`, deduplicated | 818 | 618 | 300 | 67,315 | 245 | 185 | 142 | 15,195 |

**A model trained or tuned on PatBench must not use these files**, because they contain the training rows.

## Leakage checks

`patbench verify-release` checks, for the five ranking tasks: the intersection of patents, patent families and groups between train, valid and test is empty (for `split_time` and `split_random`); train and test share no (molecule, target) pair in the task files; and the earliest priority year of every test row is not below the start of the task, and its grant year is not below 2020. Results of these checks for the reference release are recorded in `CHANGELOG.md`; see also [VALIDATION](VALIDATION.md).

In `sar_all.parquet`, among the `target_mapped` rows, 152 of 58,607 (molecule, target) pairs (0.26%) of the test split also occur in train. They come from different patents and families and are **not in the task files** (verified on r4 and r5, whose test sets are identical; not repeated on r6).
