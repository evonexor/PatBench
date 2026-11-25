# PatBench: an agent-curated patent bioactivity dataset for within-assay affinity model benchmarking

## Overview

PatBench is an open, traceable structure–activity dataset built from 345,940 US patent grants with chemical structure (MOL) files published between 2001 and 2025. Deterministic extraction code parsed 603,908 formal structure–activity records from 6,397 patents (211,734 parent molecules; 14 bioactivity endpoints), after which language-model agents automated UniProt target curation (1,087 targets) and authored deterministic structure-repair code that expanded abbreviated substituent placeholders and reversed nitriles across 128,062 records (raising public-database structure match rates from 2.6% to 94.5%). Quality filters flag unit errors, duplicate tables, substituent fragments and disputed entries without deleting any underlying record, supplemented by expert manual verification of all 205 IC50 and Binding test-group targets and visual inspection of 583 reversed-nitrile patent drawings.

Every evaluation group comprises compounds measured in a single patent under one assay label, endpoint and assigned UniProt target, enabling within-assay molecular ranking benchmarks aligned with lead optimisation. Among 200,101 affinity-record molecules outside substituent tables, 46.4% are absent from BindingDB (release 2026-09) and ChEMBL 37.

## Deposit Contents

All tabular datasets are distributed in Apache Parquet format (total deposit size: ~124.5 MB across 10 files):

| File | Size | Description |
| --- | ---: | --- |
| `sar_all.parquet` | 67.6 MB | Master table of all 603,908 formal records across 86 fields (cell-level provenance, parent and as-drawn structures, pre-repair SMILES/InChIKeys, applied structure-repair rules, pAct values, UniProt target assignments, quality/exclusion flags, public-database novelty flags and five split columns). |
| `tasks.zip` | 26.8 MB | Five within-assay ranking tasks (`rank_ic50`, `rank_ki_kd`, `rank_ec50`, `rank_highcontrast`, `rank_complex_mutant`: 915 groups, 88,774 unique rows) and one pooled regression task (`regression_affinity`: 95,376 rows). |
| `eval_all.zip` | 17.8 MB | Unified evaluation sets packaging all groups across splits for zero-shot benchmarking of external models not trained on PatBench (818 deduplicated IC50 and Binding groups, including 245 novel groups). |
| `baselines.zip` | 8.4 MB | Reference baseline predictions (random, heavy-atom count, Crippen logP and ECFP4 random forest) and leaderboard summary (`leaderboard.csv`). |
| `protocols.zip` | 2.2 MB | Frozen few-shot evaluation splits (`k = 5` and `k = 10`, five random draws each) and fixed activity-cliff pair lists (ECFP4 Tanimoto >= 0.85, potency difference >= 1.0 log). |
| `targets.parquet` | 0.8 MB | Canonical sequences, protein names, gene symbols and organisms for all 1,087 assigned UniProt entries. |
| `groups.parquet` | 0.6 MB | Group-level summary statistics across all 13,122 patent assay groups (compound counts, potency span, assay labels, UniProt accessions, inferred families and group-level flags). |
| `target_sequences.parquet` | 0.3 MB | Exact target sequences for all 915 ranking groups (835 full-length canonical sequences, 47 complex primary subunits, 25 mutation-applied sequences and 8 mature chains). |
| `release_manifest.json` | 22 KB | Machine-readable build parameters, random seeds, quality-flag counts and task statistics. |
| `README.md` | ~6 KB | This documentation file. |

Each `.zip` archive preserves its relative folder name (`tasks/`, `eval_all/`, `baselines/`, `protocols/`) so extracting all archives into a single directory reconstructs the full release tree.

## Using the files with this repository

Download the ten files from figshare into this directory (it is ignored by git), then unpack and verify:

```bash
patbench unpack data --out /tmp/patbench_data
patbench verify-release /tmp/patbench_data
```

## Accompanying Code

Processing, release verification and patent-cluster bootstrap evaluation scripts are openly available under the MIT licence at:

- **Code repository:** <https://github.com/evonexor/PatBench>

## Licence and Third-Party Terms

- **PatBench dataset:** Released under the **Creative Commons Attribution 4.0 International (CC BY 4.0)** licence (<https://creativecommons.org/licenses/by/4.0/>). You are free to share and adapt the dataset for any purpose, provided appropriate credit is given.
- **USPTO source patents:** All bioactivity measurements and chemical structures are derived exclusively from public US Patent and Trademark Office (USPTO) patent grants (2001–2025). PatBench makes no representation regarding patent validity, enforceability or freedom to operate.
- **External reference databases:** BindingDB (release 2026-09) and ChEMBL 37 were used solely for comparison, novelty annotation and unit/structure auditing; no bioactivity records were copied from either database. UniProt sequences and cross-reference identifiers remain governed by the respective terms of UniProt, ChEMBL and BindingDB.

## Citation

If you use PatBench in your research, please cite:

> Wu, W., Xiao, X., Gong, H., Li, A., Deng, S., Dai, B., Yang, G., Ge, H. & Chen, H. PatBench: an agent-curated patent bioactivity dataset for within-assay affinity model benchmarking. *figshare* https://doi.org/10.6084/m9.figshare.34925475 (2026).

```bibtex
@dataset{wu2026patbench,
  title     = {PatBench: an agent-curated patent bioactivity dataset for within-assay affinity model benchmarking},
  author    = {Wu, Wenfan and Xiao, Xianglu and Gong, Hao and Li, Ang and Deng, Shenglong and Dai, Boyuan and Yang, Guang and Ge, Hu and Chen, Hongming},
  year      = {2026},
  publisher = {figshare},
  doi       = {10.6084/m9.figshare.34925475},
  url       = {https://doi.org/10.6084/m9.figshare.34925475}
}
```
