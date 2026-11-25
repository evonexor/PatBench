"""Self-contained quickstart demonstrating PatBench within-assay evaluation.

Run with no arguments for a fast synthetic demo:
    python examples/quickstart_eval.py

Or evaluate real PatBench reference baselines on an unpacked figshare release:
    patbench unpack data --out /tmp/patbench_data
    python examples/quickstart_eval.py --release /tmp/patbench_data
"""
import argparse
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from patbench import evaluate_predictions, verify_release


def _make_synthetic_demo(root: Path):
    """Create a minimal synthetic multi-patent ranking task and two prediction files."""
    rng = np.random.default_rng(42)
    rows, model_preds, rand_preds = [], [], []

    for p_idx in range(12):
        patent = f"US99000{p_idx:02d}"
        gid = 1000 + p_idx
        n_mols = 25
        true_pact = rng.normal(7.0, 0.9, size=n_mols)
        noisy_pred = true_pact + rng.normal(0.0, 0.7, size=n_mols)
        rand_pred = rng.normal(0.0, 1.0, size=n_mols)
        for m_idx in range(n_mols):
            ikey = f"SYNTH{p_idx:02d}MOL{m_idx:04d}-UHFFFAOYSA-N"
            rows.append({
                "group_id": gid,
                "patent_core": patent,
                "uniprot_accession": f"P00{p_idx:03d}",
                "Activity Type": "IC50",
                "InChIKey": ikey,
                "pAct": float(true_pact[m_idx]),
                "heavy_atoms": int(rng.integers(20, 40)),
                "split_time": "test",
                "test_tier": "novel" if p_idx < 4 else "reference",
            })
            model_preds.append({"group_id": gid, "InChIKey": ikey, "score": float(noisy_pred[m_idx])})
            rand_preds.append({"group_id": gid, "InChIKey": ikey, "score": float(rand_pred[m_idx])})

    task_path = root / "synthetic_ic50.parquet"
    pred_path = root / "model_preds.parquet"
    ref_path = root / "random_preds.parquet"
    pd.DataFrame(rows).to_parquet(task_path, index=False)
    pd.DataFrame(model_preds).to_parquet(pred_path, index=False)
    pd.DataFrame(rand_preds).to_parquet(ref_path, index=False)
    return task_path, pred_path, ref_path


def _print_report(label: str, report: dict):
    cov = report["coverage"]
    s = report["summary"]
    print(f"\n=== {label} ===")
    print(
        f"Coverage: {cov['scored_groups']}/{cov['task_groups']} groups "
        f"({cov['task_rows']:,} rows, {s['patents']} patents, {s['bootstrap_replicates']} bootstrap resamples)"
    )
    for metric in ("spearman", "kendall", "pair_accuracy", "ndcg10", "top10pct_hit"):
        med = s.get(f"{metric}_median")
        lo, hi = s.get(f"{metric}_ci_lo"), s.get(f"{metric}_ci_hi")
        if med is not None:
            print(f"  {metric:<16} {med:+.3f}  95% CI [{lo:+.3f}, {hi:+.3f}]")
    paired = report.get("paired_vs_reference")
    if paired and paired.get("delta_median") is not None:
        print(
            f"  delta Spearman   {paired['delta_median']:+.3f}  "
            f"95% CI [{paired['ci_lo']:+.3f}, {paired['ci_hi']:+.3f}] "
            f"(vs reference, {paired['groups']} shared groups)"
        )


def main():
    ap = argparse.ArgumentParser(description="PatBench within-assay evaluation quickstart.")
    ap.add_argument("--release", help="Path to an unpacked PatBench release directory (optional).")
    ap.add_argument("--task", default="rank_ic50", help="Task name when --release is given (default: rank_ic50).")
    ap.add_argument("--bootstrap", type=int, default=500, help="Bootstrap resamples for the demo (default: 500).")
    args = ap.parse_args()

    if args.release:
        rel = Path(args.release).expanduser().resolve()
        checks = verify_release(rel)
        passed = sum(c.ok for c in checks)
        print(f"Release verification: {passed}/{len(checks)} checks passed ({rel})")
        task_path = rel / "tasks" / f"{args.task}.parquet"
        pred_path = rel / "baselines" / f"{args.task}_ecfp_rf.parquet"
        ref_path = rel / "baselines" / f"{args.task}_random.parquet"
        rep = evaluate_predictions(task_path, pred_path, reference_path=ref_path, bootstrap=args.bootstrap)
        _print_report(f"Real PatBench ({args.task}: ECFP4 + RF vs random)", rep)
    else:
        with tempfile.TemporaryDirectory() as td:
            task_path, pred_path, ref_path = _make_synthetic_demo(Path(td))
            rep = evaluate_predictions(task_path, pred_path, reference_path=ref_path, bootstrap=args.bootstrap)
            _print_report("Synthetic Demo (noisy oracle vs random baseline)", rep)


if __name__ == "__main__":
    main()
