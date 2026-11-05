# Validation: how the data are checked

The aim is that the released data contain no fabricated or altered values. The reference release of this repository is r6 (built 2026-10-31, frozen 2026-11-01). Validation consists of the following parts:

1. **Structure repair (r6):** how the repair of wrongly stored structures was validated.
2. **Data authenticity checks (r4):** row-by-row comparison with the upstream tables and with the original patents. These were run on r4 after its build and were **not repeated on r6**.
3. **Validation sample and unit audit:** an automated structure check and a manual review of 400 records by two annotators, drawn before r4 and not re-run on r6, and a unit audit of the r4 task sets.
4. **Pipeline rebuild (r4, historical):** a rebuild with an earlier snapshot of this repository, compared with r4.
5. **Code:** see the pointer below.

Every statement that rests on r4 or on an earlier stage is labelled as such. This repository does not claim to have rebuilt r6.

## Structure repair (r6)

### Problem

In a MOL file, ChemDraw writes an abbreviated atom (F3C, MeO, NH, ...) as **one placeholder atom whose element is always `C`, plus an alias line** that carries the real atom or group. RDKit discards aliases when it reads a MOL file, so the placeholder became a carbon atom. Releases up to r5 stored structures read this way. The abbreviation check of r5 caught only aliases that introduce several atoms (8,263 patent-molecule pairs) and treated a bare `H` alias as harmless.

| Problem | How it is written | How it was stored | Affected rows (r5) |
|---|---|---|---:|
| NH in a ring | two-line label `H` / `N` on a placeholder with two bonds; the alias text contains a line break (CR) | an extra ring **carbon** instead of **NH** | 98,083 |
| Abbreviated group: `F3^C` (CF3), `Me^O` (OMe), `F3C^O` (OCF3), `D3^C` (CD3), `Et^O`, `Me2^N`, `O2^S`, `Boc`, ... | placeholder atom plus alias | placeholder became CH3, the group was lost | 18,042 |
| Reversed nitrile | charged N(+)#C(-) bonded through N (isonitrile); mostly without any alias | isonitrile | 15,885 (13,205 rows without any alias) |
| Text-only alias: counter-ion (2HCl), annotation (ENT-1, Example 97) | placeholder atom of degree 0 | a free "carbon" fragment | 2,177 (removed; no effect on the parent molecule) |

In total **128,062 rows (21.2% of 603,908) and 46,740 molecules** were affected. The class counts add up to more than this total, so a row can belong to more than one class. The NH class is the largest: 61,174 structure files from 2,154 patents, all of which the r5 check had treated as harmless.

### Method: two implementations

Two implementations were written independently; each starts from the original MOL files.

1. `07_structure_repair/alias_repair.py`: a catalogue of labels, and in-place rewriting of the placeholder atom.
2. `07_structure_repair/alias_expand_independent.py`: written by a separate automated (LLM-based) session that received only the task description and the data and had not seen the code of the first implementation. It reads labels byte by byte, parses labels that contain a caret with a grammar, checks the valence of every group on its own, checks that the whole molecule sanitises, and applies a text-direction check (the label stands left or right of the bond) that does not use any public database. A label that contradicts its direction is not expanded.

r6 uses the **result of the second implementation** (wider coverage, more conservative rules) plus the nitrile step (`07_structure_repair/merge_final.py`); the first implementation serves as a cross-check. Both implementations were written by automated sessions; no review by a domain expert is recorded.

Rules: the placeholder atom is rewritten in place into the first atom of the group, so that wedge stereo information on neighbouring bonds stays valid; a reversed nitrile is turned back into a nitrile. Only labels in the rule table are expanded; nothing is guessed. Labels without a rule (peptide sequences, polymers, antibodies, Markush placeholders such as `X`, `M`, `AA`, codes such as `Pv1`-`Pv6`) stay as drawn and keep the flag `structure_alias_suspect`. Superatoms damaged by the MOL exporter (`HO2^C`, `EtO2^C`, `MeO2^S`: the group is cut in two and part of it is missing; 217 records) are also left unresolved.

The repaired structures are written to an upstream copy (`uspto-sar-recovery-20261031-structfix/results/SAR_merged.parquet`; the original directory is unchanged), which the build reads. `sar_all.parquet` gains three columns: `structure_repair` (rules used: `NH`, `group`, `nitrile`, `text`, joined by semicolons; empty if unchanged), `SMILES_before_repair` and `InChIKey_before_repair`. The task files gain `structure_repair`.

### Check 1: agreement of the two implementations

- For the 69,847 records for which both implementations return a structure, the structures are **identical in 100.0%** (0 disagreements). The comparison excludes the nitrile step, which the second implementation does not perform, and the free carbon fragments that the first implementation keeps.
- The second implementation covers 903 more records (rare labels such as `O\sabs`, `F2^C`, `(i-Pr)^O`); the first covers 487 more (the second refuses to expand them: direction contradicts the label, `#` labels).
- The first reading of the NH class by the first implementation was wrong: the line break (CR) inside an alias text was taken as the end of the line, so only `H` of a two-line label was read. The second implementation read these labels correctly.
- In the report of the second implementation (`structure_repair/second_implementation_README.md`, a run on 9,360 MOL records) the chosen reading of `Me^O` (OMe) matches a public molecule of the same patent for 703 of 759 records (92.6%), the alternative reading (CH2-OH) for 15 of 759 (2.0%); `Et^O`: 26 of 27 against 0 of 27.

### Check 2: recovery of public molecules (independent of the MOL files and of both rule sets)

BindingDB and ChEMBL contain molecules from the same patents. A correct structure should be found among the public molecules of its own patent (comparison of the connectivity block of the InChIKey); a structure with a missing atom is not found.

| Stored structure | Share found |
|---|---:|
| Before repair (35,163 records, covering 1,207 patents that are in the public databases) | **2.6%** |
| After repair | **94.5%** |
| Control: molecules without an alias in the same patents | 92.7% |

By rule:

| Rule | Records | Before | After | Control |
|---|---:|---:|---:|---:|
| NH | 29,629 | 0.4% | 95.4% | 93.4% |
| Abbreviated group | 3,304 | 12.7% | 91.6% | 92.5% |
| Reversed nitrile (molecules that r5 had not flagged) | 1,155 | 15.3% | 63.7% | 71.2% |

For the reversed nitrile, 168 records matched already before the repair, because the public databases carry the same reversed drawing. Only 23 records were found before the repair and not after it. End to end, `in_public_data` among the repaired rows rose from 2.0% to 55.7%, while the unchanged rows stayed at 52.6%; 370 rows lost a public match and 69,126 rows gained one.

### Check 3: review of the nitrile conversion against the original drawings (2026-10-31)

The conversion of isonitrile to nitrile concerns 15,885 SAR records, 7,773 original MOL files (7,938 positions) and 613 patents. All positions were read again from the original MOL files (label `N^C` 7,354 positions; `NC` 571; `CN` 4; `^NC` 4; no corresponding superatom label 5).

- 7,770 files have source evidence for the conversion; for many of them this is the explicit original label and bond geometry, not a visual check of every drawing.
- The original drawings of 583 distinct files (588 positions) were inspected, including all 556 images of the test tasks (561 positions, 35 patents), all rare or unlabelled types and all `NC` orientation anomalies. No wrongly converted nitrile was found in the test-task images. The review covers only the direction of the nitrile bond; it does not confirm other groups, stereochemistry, activity values or targets.
- 3 records are quarantined because drawing and text conflict (two related records of Example 107 in US8263781B2 and US8399681B2, where the drawing shows the C end bonded but the name and starting materials say isocyano; one record of US11285151B2 with an ambiguous orientation). None of them is in a ranking or regression task. They are flagged `structure_alias_suspect` in the derived table `structure_repair/nitrile_audit/sar_all_reviewed.parquet` (3,681 flagged rows against 3,678 in the main `sar_all.parquet`); the main file and all tasks are unchanged by this audit.
- `structure_repair/nitrile_audit/verified_nitrile_guard.py` is an evidence gate for future rebuilds. The global nitrile rule of the conversion script should not be applied to new data without evidence.

### Limits of the repair validation

- The InChIKey connectivity block of CD3 equals that of CH3, so Check 2 cannot distinguish deuterated groups.
- For some patents (US10793535, US11866436 and others), OCF3 molecules match less often than the control; this may only mean that the public databases did not record them.
- About 60% of the patents (half of the records) are not covered by the public databases. There, correctness rests on the rules and on the agreement of the two implementations, both of which were produced by automated sessions.
- `F2^C` read as CHF2 is the least certain rule (6 of 8 public matches, all in one patent). `O\sabs` and `N\sabs` (an atom plus the stereo-flag text "abs") are inferred; 43 of 43 match.
- The public databases can contain the same error as the pipeline (the 168 reversed nitriles above); a shared error is invisible to Check 2.
- The repair covers only MOL files with alias lines and reversed nitriles. A scan of molecules without alias lines against the public databases found no further widespread structure problems; this is not an exhaustive proof.
- 1,953 records could not be expanded and keep the flag (see [LIMITATIONS](LIMITATIONS.md)). About 400 records whose `(S)` stereo text was exported as a free S atom were not processed; this does not affect the parent molecule.
- The target mapping and the three UniProt corrections of r5 were carried over without new verification.
- No review of the repaired structures by a domain expert is recorded in the r6 documentation.

Record-level files are in `structure_repair/`: `repaired_records.tsv` (71,043 MOL files, structure before and after), `unresolved_records.tsv` (1,953 records), `unresolved_labels.tsv` (558 distinct labels), `validation_by_record.tsv`, `second_implementation_README.md`, `second_implementation_match_rates_by_label.tsv`, `second_implementation_hypothesis_tests.tsv` and `nitrile_audit/`.

## Data authenticity checks (r4)

These checks were run on r4 after its build, directly on the released files and the upstream extraction tables (r4 documentation, section on data authenticity). **They were not repeated on r6.** The structure columns of 128,062 rows of r6 differ from the original upstream table by design (see above), so the zero-difference results of checks 1 and 3 do not carry over to r6 unchanged. The r6 build log records one corresponding result: provenance was found for 603,908 of 603,908 rows (100.00%).

1. **Row-by-row comparison with the upstream extraction table (r4):** all 603,908 rows of `sar_all.parquet` were found in the upstream `SAR_merged` (no duplicate record IDs); patent number, target text, endpoint, relation, value, unit, converted nM value, structure SMILES and InChIKey show **zero differences**.
2. **Activity conversion (r4):** for all 526,752 exact values, pAct = 9 - log10(nM) holds with a maximum error of 0.
3. **Row-level provenance of the task files (r4):** every row of the five ranking tasks traces back to `sar_all` through `source_record_id`; `group_id` and `InChIKey` agree; for rows that merge repeated measurements, pAct equals the median of the individual pAct values (maximum difference 0); `n_measurements` equals the number of records.
4. **Difference between r4 and r3 (r4):** `sar_all.parquet` has the same rows, row order and record IDs; only `split_time` and `test_tier` changed, and three first-public-date columns were added. Rows that were already in train or valid are unchanged.
5. **Back-check against the original patent archives (r4):** 900 random measurements (600 test, 300 train) were traced to the original archives. Agreement with the upstream evidence table: 900 of 900. `value_nM` and `pAct` recomputed from value and unit: 900 of 900. The text of the table cell in the patent XML equals `raw_value`: 899 of 900. The structure computed from the patent MOL file equals the stored structure: 900 of 900. The remaining record comes from the patent US6211376B1 in the old SGM format, which the parser cannot read; it was not checked, which is not a failed check. The script is `pipeline/12_verify/provenance_check.py`. Note that the structure comparison uses the same RDKit parse as the pipeline; for alias atoms it could not detect the error that the r6 repair corrects (see the validation sample below).
6. **Independent check of the protein structures (r4):** re-read with gemmi, without using intermediate results of the generation. For 20 PDB structures, the C-alpha coordinates of 8,131 residues in the cleaned files equal the original mmCIF files (maximum distance 0). Five random original mmCIF files have the same SHA-256 as freshly downloaded RCSB files. Ten random PDB entries agree with the RCSB and PDBe records in method, resolution, release date and UniProt accession (10 of 10); five random AlphaFold models agree with the AlphaFold DB records. For 64 groups, the recommended structure sequence was re-aligned to the group sequence: 63 agree with the table within 0.01 identity; the exception is `5DXE_A` of ESR1 Y537S, for which a residue-by-residue comparison of 231 residues shows zero differences. The new structure rows added in r5 and r6 did not receive this independent check; they were checked by `check_structures.py`, which passed all 177 checks for each of three random seeds on r6.
7. **Split invariants (r4):** the minimum grant year of the test set of each task is 2020; the intersections of patents, patent families and assay groups between train, valid and test are empty. These checks are implemented as `patbench verify-release` (see [SPLITS](SPLITS.md)).

**Checks 5 and 6 are random samples, not exhaustive back-checks. The manual review of the target mappings and the adjudication of the 61 annotation disagreements have not been done** (see [LIMITATIONS](LIMITATIONS.md)).

## Validation sample and unit audit (`sar-validation-20261019/VALIDATION_REPORT.md`)

The validation sample was built on 2026-10-19/20, before r4, and was not re-run on r6.

- **Automated structure check:** 20,379 stratified random records from 3,575 patents and eight endpoints; for 100.0% the stored structure was found among the structures parsed from the MOL files of the source patent (95% Wilson interval 99.981-100%). For the 7,293 recovery records that name a specific structure file, the InChIKey of that file equals the claimed one in 100.0%. This shows that the stored structure agrees with a re-parse by the same toolkit. It does not show that the molecule is the right compound for the activity value, and for alias atoms the check is circular: the stored and the recomputed InChIKey come from the same parse, which is why this check did not reveal the problem repaired in r6.
- **Source data:** 134,082 of 1,268,434 distinct MOL files (10.57%) in the sampled archives cannot be parsed by RDKit. No record of the dataset draws on an unparseable file.
- **Unit audit (r4 documentation; r4 task sets, not repeated on r6):** compared with BindingDB and ChEMBL, 97.9% of the shared measurements in the tasks differ by at most 0.5 log; tables that deviate were excluded by the six unit checks.

### Manual review

400 items were drawn from 16 strata (extraction pass × source format × target-mapping confidence; small strata deliberately over-represented). The evidence sheets were built from the patent XML, not from pipeline records, so that an annotator can contradict the pipeline. Two annotators worked in separate sessions without shared context. **Both annotators were automated (LLM-based) sessions; no domain expert took part.** The agreement statistics measure the agreement of two independent automated passes.

| Question | Judgeable items | Population-weighted accuracy (95% CI) | Raw agreement | Cohen's kappa |
|---|---:|---|---:|---:|
| Q1 value (endpoint, number, unit, conversion, compound) | 394 | 97.8% [96.3, 99.4] | 98.5% | 0.72 |
| Q2 structure to compound-label association | 120 | 76.8% [68.2, 85.4] | 89.5% | 0.81 |
| Q3 target protein | 329 | 85.8% [82.4, 89.1] | 95.8% | 0.90 |

- The population-weighted figures should be quoted, because the sample over-represents small strata. The kappa of Q1 is depressed by the very high prevalence of correct items; the raw agreement is the better summary.
- Q2 could not be determined for 238 of the 400 items (59.5%), because the patent gives no usable link between compound label and drawing or the dataset records no label. This share has to be reported together with the accuracy.
- The main Q2 error is a substituent fragment stored as the compound (tables of the layout `Example | R | IC50`). Such groups are now excluded by `group_is_fragment_suspect` (group-median heavy-atom count below 16); the Q2 figure was measured before this exclusion. On the 37 annotated fragment errors the detector has precision 1.00 and recall 0.97. The threshold was chosen on these same items, so the figures are in-sample, and validation on a fresh sample is still open.
- The 61 disagreements are queued for expert adjudication (`adjudication_queue.tsv`) and are not yet adjudicated.

## Pipeline rebuild (r4, earlier snapshot of this repository; not repeated for r6)

This is a historical result. On the original server, a release directory was rebuilt from the same inputs as r4 with the code of an earlier snapshot of this repository (`patbench release`) and compared table by table with r4 using `tools/compare_releases.py` (relative tolerance 1e-9 for floating-point numbers, missing equals missing, row order significant). The files `public_pairs`, `alias_scan` and `novelty_flags` in `work/` were reused from r4 (they depend only on the public databases and the raw archives; r4 had reused them as well); `mol_std.parquet` was regenerated by `standardize`.

| Table | Result |
|---|---|
| `sar_all.parquet` (603,908 rows × 83 columns) | identical values |
| `groups.parquet` (13,122 rows), `targets.parquet` (1,088 rows) | identical values |
| six task files in `tasks/` (`rank_ic50` 62,244 rows, `rank_ki_kd` 17,227, `rank_highcontrast` 74,693, `rank_complex_mutant` 5,235, `rank_ec50` 2,192, `regression_affinity` 94,335) | identical values |
| `label_repaired.tsv` | identical bytes |
| `release_manifest.json` | identical keys except build date, absolute paths and console log |
| `protocols/` (25 files: few-shot support sets, activity-cliff pairs, all seeded) | identical values |
| `eval_all/` (5 files) | identical values |
| `baselines/` (49 prediction files and `leaderboard.csv`, including ECFP4 + random forest, all seeded) | identical values |

`target_audit.parquet`, `target_sequences.parquet`, `sequence_fixes.tsv` and `target_review_test.tsv` are not produced by `release`: r4 reused those of r3 and recomputed only the columns that depend on the split. `carry_over_target_files.py` regenerated the four files from r3, and they were identical to those of r4 by value or by byte.

The first complete rebuild (`patbench release`, from standardize to leaderboard, about 50 minutes) gave 89 identical files; the only difference was one extra key in the manifest, explained below. After the fix, standardize and build were repeated and 11 files (including the manifest) were identical. The steps `eval-support` and `review-page` were regenerated as well: the two parquet files and two TSV files of `eval_support/` were byte-identical to the existing ones, and the 238 rows embedded in the target review page equalled the previously produced page column by column, with an identical template.

Two entry points were run from the command line as a smoke test. `patbench run evaluate` evaluated the heavy-atom baseline of `rank_ki_kd` on r4 and gave a median Spearman of 0.261, as stated in the r4 documentation. `patbench run verify-provenance` back-checked 60 measurements (40 test, 20 train) against the original patent archives: agreement with the upstream evidence table, recomputed value conversion, XML cell text and MOL structure were all 60 of 60, and the pAct of 107 single-measurement rows equals `sar_all` row by row. The 900-measurement sample of check 5 above is the larger one.

Two problems found during the rebuild were in the repository wrapper, not in the data. `build_patbench.py` expects the `tasks/` directory to exist; `patbench release` now creates it (with a test). The variable injected for path rewriting was originally named `WS`; because `build_patbench.py` writes every upper-case global variable into the `config` block of `release_manifest.json`, the first rebuilt manifest had an extra `WS` key. After renaming the variable to `_ws`, the rebuilt manifest equals r4.

## Code

See [CODE_PROVENANCE.md](CODE_PROVENANCE.md); results are recorded in `CHANGELOG.md`.

## What the checks do not cover

- The upstream extraction, recovery and structure-catalogue stages (345,940 scanned patent archives) were not re-run in this repository. The code hashes recorded at run time (see [CODE_PROVENANCE](CODE_PROVENANCE.md)) document which code produced the tables; they do not show that a re-run would give the same tables.
- The curation steps of the target mapping and of the validation annotation depend on curated results that are not part of the repository, so they were not reproduced.
- The back-check against the original patent archives is a sample (900 records), not a full comparison.
- The r4 checks were not repeated on r6, and no comparison of the r6 structure repair with a reference set curated by a domain expert is recorded. Outside the patents covered by BindingDB and ChEMBL, the repaired structures are supported only by the rules and by the agreement of the two implementations.
- A comparison with public databases cannot reveal an error that BindingDB shares with the pipeline: BindingDB curates US patents too, and reading the same table wrongly in the same place (for example US10966980) leaves the two in agreement. Agreement on the overlapping part shows that the extraction is accurate; it is not evidence that the experiments are reproducible.
- The target mappings of the test groups and the 61 annotation disagreements have not been reviewed or adjudicated by experts.
