"""The workspace the pipeline reads from and writes to.

The scripts under pipeline/ were written for one directory tree. Rather than rewrite every path, the
repository keeps that tree as its contract: `PATBENCH_WORKSPACE` points at a directory with the same
relative layout, and every path literal in the scripts is `_ws + "<relative path>"`.

    $PATBENCH_WORKSPACE/
      raw/patent-data/                              USPTO MOL archives and the public databases
        standardized/BindingDB.parquet              public data, standardised by pipeline/01_public_data
        standardized/ChEMBL.parquet
      outputs/
        uspto-sar-full-20261004/results/            stage 02 output
        uspto-all-structures-20261004/results/      stage 04 output
        uspto-sar-recovery-20261006/results/        stage 03 output (structures as extracted, before the repair)
        uspto-sar-recovery-20261031-structfix/results/   stage 07 output: the same tables with repaired structures
        target-uniprot-map-20261016/work/           stage 05 (v1) working files
        target-uniprot-map-v2-20261019/             stage 05 (v2) output
        target-uniprot-map-v2.1-20261031/           stage 05 (v2.1) output: target_map.parquet, read by the build
        sar-validation-20261019/                    stage 06 output
        nitrile-evidence-audit-20261031/            working files of the nitrile evidence audit (stage 07)
        affinity-sar-release-20261019/              the first release (input of the first alias scan)
        <release directory>/                        PatBench itself (any name; passed as --out)
      struct-deps-20261030/                         optional: gemmi installed with `pip --target`

The build scripts read two directories from environment variables (the second default of the
original run is the pre-repair one): `PATBENCH_SAR_DIR` and `PATBENCH_TARGET_MAP_DIR`. The runner sets
them to the repaired SAR tables and to the v2.1 target map unless they are already set.
"""
import os
from dataclasses import dataclass
from pathlib import Path

REL = "outputs/"
SAR_PRE = REL + "uspto-sar-recovery-20261006/results/"
SAR = REL + "uspto-sar-recovery-20261031-structfix/results/"
FULL = REL + "uspto-sar-full-20261004/results/"
MAP = REL + "target-uniprot-map-v2.1-20261031/"
VAL = REL + "sar-validation-20261019/"
PREV = REL + "affinity-sar-release-20261019/"
STD = "raw/patent-data/standardized/"


@dataclass(frozen=True)
class Requirement:
    path: str          # relative to the workspace
    what: str
    needed_by: tuple   # step names


REQUIREMENTS = [
    Requirement(SAR + "SAR_merged.parquet", "all formal SAR rows, structures repaired (stage 07 output)",
                ("standardize", "novelty", "build")),
    Requirement(SAR + "new_SAR.parquet", "rows added by the recovery stage, with provenance", ("build",)),
    Requirement(SAR + "patents.parquet", "per-archive status of the recovery run", ("build",)),
    Requirement(SAR_PRE + "SAR_merged.parquet", "the same table before the structure repair (input of the repair)",
                ("patch-sar",)),
    Requirement(FULL + "SAR_with_evidence.parquet", "first-run SAR rows with the source cell", ("build",)),
    Requirement(FULL + "patents.parquet", "patent metadata: priority and publication dates", ("build",)),
    Requirement(MAP + "target_map.parquet", "(patent, assay label) -> UniProt mapping",
                ("novelty", "build", "join-target-map")),
    Requirement(VAL + "merged_annotations.tsv", "annotations of the validation sample", ("build",)),
    Requirement(VAL + "sample.parquet", "the validation sample", ("build",)),
    Requirement(STD + "BindingDB.parquet", "BindingDB, standardised", ("novelty", "public-pairs")),
    Requirement(STD + "ChEMBL.parquet", "ChEMBL, standardised", ("novelty", "public-pairs")),
    Requirement(PREV + "sar_all.parquet", "first-release table the first alias scan was run on", ("alias-scan",)),
    Requirement("raw/patent-data", "original patent archives (the alias scans and the provenance check open them)",
                ("alias-scan", "alias-extract", "alias-scan-all", "verify-provenance")),
]

# files build_patbench.py expects in <release>/work/ that no step of the release chain produces
# (the chain itself writes mol_std.parquet and novelty_flags.parquet)
WORK_INPUTS = ("public_pairs.parquet", "alias_scan.parquet")


def workspace(arg=None):
    ws = arg or os.environ.get("PATBENCH_WORKSPACE")
    if not ws:
        raise SystemExit("set PATBENCH_WORKSPACE or pass --workspace (see docs/REPRODUCE.md)")
    return Path(ws).expanduser().resolve()


def check_workspace(ws, steps=None):
    """[(Requirement, present)] for the requirements of `steps` (all when None)."""
    ws = Path(ws)
    out = []
    for r in REQUIREMENTS:
        if steps is None or set(steps) & set(r.needed_by):
            out.append((r, (ws / r.path).exists()))
    return out
