"""High-level evaluation and release-unpacking helpers for PatBench."""
import importlib.util
import shutil
import zipfile
from pathlib import Path

from patbench import PIPELINE

_EVAL_SCRIPT = PIPELINE / "10_evaluation" / "evaluate_patent_bootstrap.py"

RELEASE_TOP_FILES = (
    "sar_all.parquet",
    "groups.parquet",
    "targets.parquet",
    "target_sequences.parquet",
    "release_manifest.json",
)

RELEASE_ZIPS = (
    "tasks.zip",
    "protocols.zip",
    "eval_all.zip",
    "baselines.zip",
)


def _load_evaluator_module():
    """Import the frozen patent-cluster bootstrap evaluator without modifying it."""
    spec = importlib.util.spec_from_file_location("_patbench_eval_bootstrap", _EVAL_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def evaluate_predictions(
    task_path,
    pred_path,
    *,
    split="split_time",
    part="test",
    tier=None,
    reference_path=None,
    profile_path=None,
    collapse_stereo=False,
    cliffs_path=None,
    fewshot_path=None,
    bootstrap=None,
):
    """Score a prediction table on a PatBench ranking task with patent-cluster bootstrap CIs.

    Returns the same JSON-serialisable report dictionary produced by
    ``pipeline/10_evaluation/evaluate_patent_bootstrap.py``.
    """
    mod = _load_evaluator_module()
    if bootstrap is not None:
        mod.BOOTSTRAP = int(bootstrap)

    full = mod.read(str(task_path))
    task = full[full[split] == part].copy() if (split in full.columns and part) else full.copy()
    scope_rows = len(task)
    scope_groups = int(task.group_id.nunique())

    if tier and "test_tier" in task.columns:
        task = task[task.test_tier == tier].copy()

    pred = mod.read(str(pred_path))
    protocol = mod.read(str(fewshot_path)) if fewshot_path else None
    gm = mod.per_group_fewshot(task, pred, protocol, collapse_stereo) if fewshot_path else mod.per_group(task, pred, collapse_stereo)
    if not len(gm):
        raise ValueError("No scorable groups in this slice")

    report = {
        "summary": mod.summary(gm),
        "coverage": {
            "task_rows": len(task),
            "task_groups": int(task.group_id.nunique()),
            "scored_groups": len(gm),
            "score_semantics": "higher means more potent",
            "scope_rows_before_filter": scope_rows,
            "scope_groups_before_filter": scope_groups,
        },
    }

    if reference_path:
        ref = mod.read(str(reference_path))
        rg = mod.per_group_fewshot(task, ref, protocol, collapse_stereo) if fewshot_path else mod.per_group(task, ref, collapse_stereo)
        report["paired_vs_reference"] = mod.paired_comparison(gm, rg)

    if cliffs_path:
        if fewshot_path:
            raise ValueError("Cliff evaluation requires zero-shot predictions; run separately")
        report["cliffs"] = mod.cliff_accuracy(task, pred, mod.read(str(cliffs_path)), set(task.group_id))

    profile = mod.read(str(profile_path)) if profile_path else None
    if profile is not None and "task" in profile.columns:
        profile = profile[profile.task == Path(task_path).stem]
    report["strata"] = mod.stratified_summaries(gm, profile)
    return mod.json_safe(report)


def unpack_figshare_release(src_dir, out_dir):
    """Unpack a downloaded figshare PatBench deposit into a standard release directory tree.

    Returns the resolved ``out_dir`` ``Path``.
    """
    src = Path(src_dir).expanduser().resolve()
    if not src.is_dir():
        raise FileNotFoundError(f"Source directory does not exist: {src}")

    found_top = [src / name for name in RELEASE_TOP_FILES if (src / name).is_file()]
    found_zips = [src / name for name in RELEASE_ZIPS if (src / name).is_file()]
    if not found_top and not found_zips:
        raise FileNotFoundError(
            f"No PatBench deposit files found in {src}. "
            "Download the dataset from https://doi.org/10.6084/m9.figshare.34925475 first."
        )

    out = Path(out_dir).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    for p in found_top:
        dst = out / p.name
        if not dst.exists() or dst.stat().st_size != p.stat().st_size:
            shutil.copy2(p, dst)

    for zp in found_zips:
        with zipfile.ZipFile(zp, "r") as zf:
            zf.extractall(out)

    return out
