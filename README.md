# PatBench

**Data-processing pipeline of PatBench**, a benchmark for within-series potency ranking built from US patents.

PatBench turns the activity tables of US patents into ranking tasks. Each *group* is one patent, one assay label and one
endpoint, and the question a model has to answer is the one a medicinal chemist faces during lead optimisation: *within
one series, which molecule is more potent?* This repository holds the code that goes from the patent archives to the
released benchmark tables: extraction, recovery of activity links, target mapping, validation sampling, structure repair,
benchmark construction, target structures, evaluation and verification.

The repository contains **code and documentation only**. The patent archives, the intermediate tables and the released
benchmark are not stored here (see [Data](#data)).

## Contents

- [Overview](#overview)
- [Pipeline at a glance](#pipeline-at-a-glance)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Building a release](#building-a-release)
- [Repository layout](#repository-layout)
- [Code provenance and integrity](#code-provenance-and-integrity)
- [Reference release](#reference-release)
- [Data](#data)
- [Documentation](#documentation)
- [Known gaps](#known-gaps)
- [Testing](#testing)
- [Citation](#citation)
- [License](#license)

## Overview

| | |
|---|---|
| Source | 345,940 US patent MOL archives; the 6,397 patents that yield formal SAR were granted 2001-2025 |
| Task | rank the molecules of one assay group by potency |
| Unit | group = patent x assay label x endpoint |
| Splits | time (main), random, scaffold; patents and patent families never cross splits |
| Tasks | `rank_ic50`, `rank_ki_kd`, `rank_ec50`, `rank_highcontrast`, `rank_complex_mutant`, and the regression task `regression_affinity` |
| Provenance | every row carries its patent, table cell, MOL file and the raw value, so it can be traced back to the patent text |
| Nothing deleted | all formal SAR records stay in `sar_all.parquet`; every selection is a flag column, so every filter is reproducible |

## Pipeline at a glance

```
USPTO MOL archives (345,940 ZIPs)                      BindingDB / ChEMBL / HARVEST / SureChEMBL
   |                  |                                               |
 02 extract         04 structure catalogue                       01 standardise
   |                  |                                               |
   +--------+---------+                                               |
            v                                                         |
        03 recovery of activity links  ->  SAR_merged                 |
            |            |                                            |
 05 target mapping    06 validation sample                            |
 (assay label ->      (annotation, agreement,                         |
  UniProt)             failure modes)                                 |
            |            |                                            |
            +-----+------+                                            |
                  v                                                   |
        07 structure repair (MOL alias lines)                         |
                  |                                                   |
                  +------------------------+--------------------------+
                                           v
        08 build: standardize -> novelty -> build -> protocols -> eval-sets -> target audit
                                           v
        release directory: sar_all, groups, targets, tasks/, protocols/, eval_all/, release_manifest.json
                                           v
        09 target structures   10 evaluation (baselines, leaderboard)   11 support tables   12 provenance check
```

Inputs, outputs and the decisions made in each stage are described in [docs/PIPELINE.md](docs/PIPELINE.md).

## Installation

Python 3.9 or newer.

```bash
git clone <repository-url> patbench
cd patbench
pip install -e ".[dev]"
```

The command line and the tests need only `numpy`, `pandas`, `pyarrow` and `scipy`. Running the upstream stages needs more:

| Extra | Needed for | Install |
|---|---|---|
| `chem` | extraction, standardisation, build, protocols, baselines (RDKit) | `pip install -e ".[chem]"` |
| `baselines` | ECFP4 + random-forest baselines (scikit-learn, RDKit) | `pip install -e ".[baselines]"` |
| `structures` | target structure download and checks (gemmi) | `pip install -e ".[structures]"` |

## Quick start

These commands need no data:

```bash
patbench stages          # the processing stages, in order
patbench steps           # the runnable steps, with what each reads and writes
patbench verify-code     # repository files against the recorded hashes
pytest                   # tests that need no data
```

## Building a release

Building a release needs the data workspace, whose layout is described in [docs/REPRODUCE.md](docs/REPRODUCE.md).

```bash
export PATBENCH_WORKSPACE=/path/to/workspace        # layout: docs/REPRODUCE.md
patbench check-workspace                            # are the input tables present?
patbench plan --out /path/to/release                # print the steps of the build chain
patbench release --out /path/to/release --reuse-work /path/to/earlier/release/work
patbench verify-release /path/to/release            # integrity checks on the finished release
```

Individual steps run with `patbench run <step> --out <release> -- <arguments for the script>`. Any script under
`pipeline/` runs with the right environment through `patbench exec <script> -- <arguments>`.

## Repository layout

| Path | Content |
|---|---|
| `pipeline/` | the processing code, one directory per stage ([pipeline/README.md](pipeline/README.md)) |
| `src/patbench/` | command line, workspace layout, step registry, code-manifest check, release check |
| `configs/release_v1.0.json` | build parameters of the reference release (the `config` of its `release_manifest.json`); `verify-release` uses it to detect drift |
| `provenance/` | `CODE_MANIFEST.tsv` (status and hashes of every pipeline file), `STAGE_HASHES.json` (code hashes recorded by the three frozen stages), `NOT_INCLUDED.tsv` (server scripts deliberately left out) |
| `tests/` | tests: code integrity, import tool, command line, metric definitions, release checks, review page; the unit tests shipped with the extraction stages (they need RDKit) |
| `tools/` | `import_pipeline.py` (imports the code from a snapshot of the server directories), `dump_schemas.py`, `make_schema_doc.py`, `compare_releases.py` |
| `docs/` | [PIPELINE](docs/PIPELINE.md), [REPRODUCE](docs/REPRODUCE.md), [CODE_PROVENANCE](docs/CODE_PROVENANCE.md), [SPLITS](docs/SPLITS.md), [VALIDATION](docs/VALIDATION.md), [LIMITATIONS](docs/LIMITATIONS.md), [RELEASES](docs/RELEASES.md), [DATA_SCHEMA](docs/DATA_SCHEMA.md) |
| `data/` | empty; documentation only |

## Code provenance and integrity

The code in `pipeline/` was copied from the server that produced the data, not rewritten.

- The extraction, recovery and structure-catalogue stages (`02_extract`, `03_recovery`, `04_catalog`) are **byte-identical** to the code that ran. Each wrote the SHA-256 of its code into a status file when it ran; `patbench verify-code` recomputes these hashes from the files in this repository and they match.
- In the other scripts the only change is how data paths are written: a string literal that began with the server's workspace prefix now begins with `_ws`, which is read from the environment variable `PATBENCH_WORKSPACE`. Nothing else differs.
- `provenance/CODE_MANIFEST.tsv` records the SHA-256 of the original and of the copy for every file, and `tools/import_pipeline.py` repeats the import.

Details, including the few places where provenance is incomplete, are in [docs/CODE_PROVENANCE.md](docs/CODE_PROVENANCE.md).

## Reference release

The code in this repository corresponds to **PatBench v1.0, revision r6**: built on 2026-10-31 and frozen read-only on 2026-11-01
(a checksum list covers every file of the release directory).

| | r6 |
|---|---:|
| Formal SAR records | 603,908 |
| Source patents | 6,397 |
| Parent molecules | 211,734 |
| Assay groups | 13,122 |
| UniProt targets | 1,087 |
| Endpoints | 14 |

Task sizes, splits and the test tiers are in [docs/SPLITS.md](docs/SPLITS.md); the history of the revisions is in
[docs/RELEASES.md](docs/RELEASES.md). Later directories in the same location, such as `patbench-v1.0-20261031-r7`, are not
produced by this repository.

## Data

No data is distributed with this repository. The released benchmark tables are described column by column in
[docs/DATA_SCHEMA.md](docs/DATA_SCHEMA.md). A public deposit with a persistent identifier has not been made yet; this section
will point to it once it exists.

## Documentation

| Document | Content |
|---|---|
| [docs/PIPELINE.md](docs/PIPELINE.md) | what each stage does: inputs, outputs, decisions |
| [docs/REPRODUCE.md](docs/REPRODUCE.md) | environment, workspace layout, how to build a release, what has not been reproduced |
| [docs/CODE_PROVENANCE.md](docs/CODE_PROVENANCE.md) | how the code got here, what the hashes prove, gaps |
| [docs/SPLITS.md](docs/SPLITS.md) | split rules, task definitions, test tiers, leakage checks |
| [docs/VALIDATION.md](docs/VALIDATION.md) | how the data were checked, including the structure repair |
| [docs/LIMITATIONS.md](docs/LIMITATIONS.md) | known limitations and open items |
| [docs/RELEASES.md](docs/RELEASES.md) | release history |
| [docs/DATA_SCHEMA.md](docs/DATA_SCHEMA.md) | columns of every table |

## Known gaps

- The upstream stages (extraction, recovery, structure catalogue) were **not re-run** from this repository; the code hashes match, which shows the code is unchanged, not that a re-run would give identical tables.
- The curation results of the target mapping and the validation annotations are data and are not in this repository. The prompts and model versions used for the curation were not found and are not included. Two parts of the nitrile evidence audit (the geometry step and the visual review) are not scripted.
- The target mappings of the test groups have not been reviewed group by group, and 61 annotation disagreements of the validation sample have not been adjudicated by an expert ([docs/LIMITATIONS.md](docs/LIMITATIONS.md)).
- A full rebuild of r6 with these scripts, compared table by table with the frozen release, has not been done.
- `pipeline/10_evaluation/baselines.py` is the r6 script as it ran and has a known defect in the few-shot ECFP4 + random-forest baseline; `baselines_fixed.py` is the corrected copy that `patbench` runs ([docs/LIMITATIONS.md](docs/LIMITATIONS.md)).

## Testing

```bash
pytest                      # no data needed; the tests of the extraction stages run when RDKit is installed
patbench verify-code        # code hashes
patbench verify-release /path/to/release    # needs a built release directory
```

The tables used in the tests are small synthetic samples built in memory; they are labelled as not being PatBench data.

## Citation

The data descriptor is in preparation. A `CITATION.cff` will be added when the article and the data deposit have
identifiers. Until then, please cite this repository by its name and commit hash.

## License

The code in this repository is released under the [MIT License](LICENSE). The data release is planned under CC BY 4.0, which is a data license
and does not apply to this code. Third-party databases that the pipeline reads (BindingDB, ChEMBL, HARVEST, SureChEMBL, UniProt, RCSB PDB,
AlphaFold DB, USPTO patents) keep their own terms; they are not redistributed here.
