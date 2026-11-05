# Reproduction

## Environment

- Python 3.9 or newer (the original runs used 3.10; the tests of this repository were also run on 3.12). `pip install -e ".[dev]"` installs the command line and the test tools. The extraction, standardisation, build and baseline stages also need **RDKit**; the baselines need scikit-learn; the structure stage needs **gemmi**.
- On the original server `pyarrow 19` lived in a separate dependency directory that was put on `PYTHONPATH`; any other way of installing pyarrow works as well.
- The scripts are written to be run as scripts (some start working at import time). Start them through `patbench run` or `patbench exec`: these set `PATBENCH_WORKSPACE` and `PYTHONPATH` (with `pipeline/02_extract` on it, because the later stages import the extraction modules by flat module name).
- Two scripts of the structure stage (`fetch_structures.py`, `check_structures.py`) insert optional dependency directories (`struct-deps-20261030`, `outputs/uspto-sar-demo-20261003/runtime/deps`) into `sys.path` if you provide them in the workspace; a normal environment does not need them.

## Workspace

Every data path in the scripts is written as `_ws + "<relative path>"`, and `_ws` is taken from the environment variable `PATBENCH_WORKSPACE`. The workspace must have the layout below (the dates in the directory names identify the original runs and are kept as they were):

```
$PATBENCH_WORKSPACE/
  raw/patent-data/                                     original archives and the public databases
    standardized/{BindingDB,ChEMBL}.parquet            output of 01_public_data
  outputs/
    uspto-sar-full-20261004/results/                   output of stage 02 (SAR_with_evidence, patents ...)
    uspto-all-structures-20261004/results/             output of stage 04
    uspto-sar-recovery-20261006/results/               output of stage 03 (SAR_merged, new_SAR, patents ...), before the repair
    uspto-sar-recovery-20261031-structfix/results/     output of stage 07: the same tables with repaired structures
    target-uniprot-map-20261016/work/                  stage 05 v1 working files
    target-uniprot-map-v2-20261019/                    stage 05 v2: target_map.parquet, work/
    target-uniprot-map-v2.1-20261031/                  stage 05 v2.1: target_map.parquet (read by the build)
    sar-validation-20261019/                           stage 06: merged_annotations.tsv, sample.parquet ...
    nitrile-evidence-audit-20261031/                   working files of the nitrile evidence audit (stage 07)
    affinity-sar-release-20261019/                     the first release (input of the first alias scan)
    <release directory>/                               PatBench itself, any name, given with --out
  struct-deps-20261030/                                optional: gemmi installed with pip --target
```

`patbench check-workspace` lists whether the input tables of the build chain are present. Input tables may be placed in this layout as symbolic links.

The build scripts take two input directories from environment variables, with a trailing slash: `PATBENCH_SAR_DIR` (default in the scripts: the **unrepaired** recovery results) and `PATBENCH_TARGET_MAP_DIR` (default in the scripts: target map v2). `patbench` sets them to the repaired tables and to v2.1; set them yourself first to override.

## How the reference release (r6) was built

The order below is the one that was run on the original server; each line is one registered step (`patbench steps`).

1. Stages 02, 04, 03, 05, 06 produced `SAR_merged`, the target map and the validation annotations (see [PIPELINE](PIPELINE.md)).
2. A release was built from them without the structure repair (r5). It is the input of stage 07.
3. **Structure repair** (stage 07), on the r5 release directory: `alias-extract`, `alias-scan-all`, `alias-pubkeys`, `alias-repair` (implementation A), `alias-expand` (implementation B, used), `alias-validate`, `alias-merge`, `patch-sar`, `alias-scan-final`, `list-unresolved`, and the nitrile evidence audit `nitrile-inventory`, `nitrile-extract`, `nitrile-text`, `nitrile-panels`, `nitrile-finalize`, `nitrile-verify` (two parts of the audit are not scripted, see [CODE_PROVENANCE](CODE_PROVENANCE.md)).
4. **Release chain** on the repaired tables: `standardize`, `novelty`, `build`, `protocols`, `eval-sets`.
5. **Target audit and structures**: `join-target-map`, `target-audit`, `audit-targets`, then `fetch-structures` for each stage `meta entries gql select download process groups write`, then `check-structures` with several seeds.
6. **Baselines**: `baselines`, `leaderboard`. (On the original server the as-run script, now `baselines-r6-asrun`, was used; `baselines` runs the corrected copy, see [LIMITATIONS](LIMITATIONS.md).)
7. **Support tables**: `eval-support`, `target-names`, `review-table`, `review-page`, `profile`.
8. The directory was frozen: a checksum list over every file, then write permission removed (see [RELEASES](RELEASES.md)).

## Build a release directory

```bash
export PATBENCH_WORKSPACE=/path/to/workspace
patbench check-workspace --steps standardize,novelty,build
patbench plan --out /path/to/release                  # inspect the chain first
patbench release --out /path/to/release \
    --reuse-work /path/to/earlier/release/work        # copies public_pairs.parquet and alias_scan.parquet
patbench verify-release /path/to/release
```

`release` runs `standardize -> novelty -> build -> protocols -> eval-sets -> baselines -> leaderboard`, as `run_r6.sh` and `run_baselines_r6.sh` did on the original server (with the corrected baselines script); the log of each step goes to `<release>/<step>.log`. `build_patbench.py` expects `<release>/tasks/` to exist; `release` creates it.

`<release>/work/` must contain two files that no step of the chain produces:

- `public_pairs.parquet` depends only on the public databases. Produce it with `patbench run public-pairs --out /path/to/release`.
- `alias_scan.parquet` depends on the structure repair (step `alias-scan-final`). **Reuse it only from a release built on the same repaired tables.**

Environment variables you can tune: `PATBENCH_JOBS` (parallel processes of the baselines, default 4), `OMP_NUM_THREADS` (2 for `build`, 1 for the others). Add `--nice` on a shared server.

## Other steps

```bash
# Target audit (network: three UniProt requests)
patbench run join-target-map --out <release> -- <target_map.parquet> <join_output.parquet>
patbench run target-audit    --out <release> -- <join_output.parquet> <release>/target_audit.parquet
patbench run audit-targets   -- <release>/target_audit.parquet <release>/targets.parquet <release>

# Target structures (network: public files; --max-gb caps the raw files at 1.8 GB by default; stages in order)
for s in meta entries gql select download process groups write; do
  patbench run fetch-structures --out <release> -- $s
done
patbench run check-structures --out <release>

# Support tables and the target-mapping review page
patbench run eval-support --out <release>
patbench run review-page  --out <release>

# Evaluate a model on a ranking task
patbench run evaluate -- <release>/tasks/rank_ic50.parquet predictions.parquet --split split_time --part test

# Trace a random sample back to the original patent archives
patbench run verify-provenance --out <release> -- 600 300
```

### Structure repair (stage 07)

The stage works on a directory of its own (`<alias_dir>`) and takes the release built **before** the repair as input (`<pre-repair release>`). Each script documents its arguments and outputs in its docstring; `patbench steps` lists them.

| Step | Arguments after `--` | Writes |
|---|---|---|
| `alias-extract` | `<pre-repair release> <alias_dir> [workers]` | `mols/`, `records.tsv`, `alias_atoms.tsv` |
| `alias-scan-all` | `<pre-repair release> <alias_dir> [workers]` | `mols_all/`, `all_records.tsv`, `all_alias_atoms.tsv` |
| `alias-pubkeys` | `<alias_dir> [records.tsv] [out.tsv]` | `public_keys_by_patent.tsv` |
| `alias-repair` | `<alias_dir> <records.tsv> <mols_dir> <out.tsv> [workers]` | repair table of implementation A |
| `alias-expand` | none; copy `alias_expand_independent.py` into `<alias_dir>` as `expand.py` and run it there (`--eval` for the self-check) | `expansions.tsv` |
| `alias-validate` | `<alias_dir> <records.tsv> <repairs.tsv> <public_keys.tsv> <report_prefix> [pre-repair release]` | validation tables |
| `alias-merge` | `<alias_dir> <expansions.tsv> <records.tsv> <out.tsv>` | the final repair table (adds the nitrile step) |
| `patch-sar` | `<pre-repair release> <repairs.tsv> <out_results_dir>` | `SAR_merged.parquet` with repaired structures (the original directory is untouched) |
| `alias-scan-final` | `<pre-repair release> <new release> <alias_dir> <all_records.tsv> <repairs.tsv>` | `work/alias_scan.parquet` of the new release |
| `list-unresolved` | `<alias_dir> <unresolved_records.tsv> <out.tsv>` | table of unresolved labels |

## Upstream stages

Extraction, recovery and the structure catalogue are long jobs on 345,940 archives, each with its own command line (`patbench stages` lists them). In short:

```bash
python pipeline/02_extract/run_full.py   --root <archive root> --manifest <manifest.json> --output <output dir> --workers 20 --full
python pipeline/04_catalog/run_structure_catalog.py --root ... --manifest ... --output ... --prior <results dir of the earlier extraction> --workers 20 --full
python pipeline/03_recovery/run_recovery.py --manifest ... --output ... --workers 20
```

The code of these three directories **must not be edited**: when they ran, they wrote the SHA-256 of their code into a status file, and `patbench verify-code` recomputes and compares it (see [CODE_PROVENANCE](CODE_PROVENANCE.md)).

### Path constants of the recovery stage

`03_recovery/run_recovery.py` and `03_recovery/activity_recovery.py` contain four hard-coded paths. They were not changed so that the code hash stays the same. To run the stage elsewhere, create the same directory layout (or symbolic links) at these locations:

| Constant | File | Value |
|---|---|---|
| `BASELINE` | `run_recovery.py` | `<workspace>/outputs/uspto-sar-full-20261004/results/SAR.parquet` |
| `CAT` | `activity_recovery.py` | `<workspace>/outputs/uspto-all-structures-20261004/results/checkpoint_cache` |
| `OLD` | `activity_recovery.py` | `<workspace>/outputs/uspto-sar-full-20261004/results/checkpoint_cache` |
| `RAW` | `activity_recovery.py` | `<workspace>/raw/patent-data` |

## What has not been reproduced

- The upstream extraction, recovery and structure-catalogue stages were not re-run from this repository (345,940 archives). The code hashes match, which shows the code is unchanged; it does not show that a re-run would give the same tables.
- The curation steps of the target mapping and the validation annotations depend on curation results that are data, not code. They are not in this repository and are not regenerated.
- A full rebuild of r6 with this repository's scripts, compared table by table with the frozen release, has not been done. See [VALIDATION](VALIDATION.md) for the checks that have been run.
