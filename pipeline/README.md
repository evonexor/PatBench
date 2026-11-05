# pipeline/

The processing code of PatBench, one directory per stage, in the order the data flows.
`patbench stages` prints this table with the command for each stage; [docs/PIPELINE.md](../docs/PIPELINE.md)
describes inputs, outputs and decisions in detail.

| Directory | What it does | State of the files |
|---|---|---|
| `00_ops/` | copy the patent archives to the server's SSD; launch the recovery service | verbatim, site specific |
| `01_public_data/` | standardise BindingDB / ChEMBL / HARVEST / SureChEMBL and merge them | path literals rewritten |
| `02_extract/` | SAR tables from 345,940 USPTO MOL archives | **verbatim**, hash checked |
| `03_recovery/` | recover activity links the first pass skipped; `SAR_merged` | **verbatim**, hash checked |
| `04_catalog/` | catalogue of every MOL and the patent metadata, all domains | **verbatim**, hash checked |
| `05_target_map/{v1,v2,v2_1}/` | assay label to UniProt accession for every (patent, label) pair | path literals rewritten |
| `06_validation/` | validation sample, annotation pack, agreement, failure analysis | path literals rewritten |
| `07_structure_repair/` | repair of structures stored wrongly from MOL alias lines; evidence audit of the nitrile conversion (`nitrile_audit/`) | path literals rewritten |
| `08_build/` | molecule standardisation, public novelty, the PatBench build, protocols, evaluation sets, target audit | path literals rewritten |
| `09_structures/` | protein 3D structures and their independent check | path literals rewritten |
| `10_evaluation/` | evaluators (build-time and the recommended patent-bootstrap one), reference baselines (as run and corrected), leaderboard | verbatim; `evaluate_patent_bootstrap.py` and `baselines_fixed.py` come from the evaluation runs |
| `11_support/` | evaluation-support tables, target-name table, target-mapping review page, release profile | imported scripts plus `review_page/` written for this repository |
| `12_verify/` | provenance check against the original patent archives | written for this repository |
| `legacy_v0/` | the first release and benchmark builders PatBench grew out of | history only |

"Path literals rewritten" means exactly one kind of edit, made only in files that named the server: a string
literal that began with the server's workspace prefix now begins with `_ws`, defined once near the top of the file
from the environment variable `PATBENCH_WORKSPACE`. Nothing else differs from the original.
`provenance/CODE_MANIFEST.tsv` records, for every file, its status and the SHA-256 of the original and of the copy
here; `patbench verify-code` checks them and recomputes the code hashes that the three frozen stages recorded when they
ran. Scripts that were deliberately left out are listed in `provenance/NOT_INCLUDED.tsv`.

The scripts are research code and run as scripts (several do their work at import time). Use them through
`patbench run <step>` or `patbench exec <script>`, which set the environment they expect.
