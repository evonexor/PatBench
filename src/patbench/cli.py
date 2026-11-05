"""patbench command line."""
import argparse
import sys

from patbench import PIPELINE, __version__
from patbench import layout, manifest, runner, verify
from patbench.steps import RELEASE_CHAIN, STAGES, STEPS, Ctx


def _split_extra(argv):
    if "--" in argv:
        i = argv.index("--")
        return argv[:i], argv[i + 1:]
    return argv, []


def cmd_stages(a):
    for num, d, what, how in STAGES:
        print(f"{num}  {d:<26} {what}\n    {how}")
    return 0


def cmd_steps(a):
    w = max(map(len, STEPS))
    for s in STEPS.values():
        flags = ",".join(f for f, on in (("chain", s.name in RELEASE_CHAIN), ("heavy", s.heavy), ("network", s.network)) if on)
        print(f"{s.name:<{w}}  {s.script:<34} {s.summary}" + (f"  [{flags}]" if flags else ""))
        if s.notes:
            print(f"{'':<{w}}  note: {s.notes}")
    return 0


def cmd_plan(a):
    ws = layout.workspace(a.workspace) if (a.workspace or a.check) else "<workspace>"
    out = a.out or "<release>"
    print(f"workspace {ws}\nrelease   {out}\n")
    for i, name in enumerate(RELEASE_CHAIN, 1):
        s = STEPS[name]
        cmd = runner.command(s.script, s.args(Ctx(out=out, workers=a.workers)), "python")
        print(f"{i}. {name:<12} {s.summary}\n   reads : {', '.join(s.reads)}\n   writes: {', '.join(s.writes)}\n   $ "
              + " ".join(cmd).replace(str(PIPELINE), "pipeline"))
    print("\nbefore step 1, <release>/work/ needs public_pairs.parquet and alias_scan.parquet "
          "(patbench run public-pairs / alias-scan-final, or --reuse-work)")
    if a.check:
        return cmd_check(a)
    return 0


def cmd_check(a):
    ws = layout.workspace(a.workspace)
    res = layout.check_workspace(ws, None if not getattr(a, "steps", None) else a.steps.split(","))
    for r, ok in res:
        print(f"{'OK     ' if ok else 'MISSING'} {r.path:<62} {r.what}")
    miss = [r for r, ok in res if not ok]
    print(f"\n{len(res) - len(miss)} of {len(res)} inputs present under {ws}")
    return 1 if miss else 0


def cmd_run(a, rest):
    if a.step not in STEPS:
        raise SystemExit(f"unknown step '{a.step}'; see `patbench steps`")
    ws = layout.workspace(a.workspace)
    return runner.run_step(a.step, ws, a.out or "", rest, python=a.python, workers=a.workers, nice=a.nice,
                           dry_run=a.dry_run)


def cmd_exec(a, rest):
    ws = layout.workspace(a.workspace)
    cmd = runner.command(a.script, rest, a.python)
    return runner.run(cmd, runner.environment(ws), nice=a.nice, dry_run=a.dry_run)


def cmd_release(a):
    ws = layout.workspace(a.workspace)
    return runner.release(ws, a.out, a.reuse_work, a.python, a.workers, a.nice, a.dry_run, a.skip.split(",") if a.skip else ())


def cmd_verify_code(a):
    bad = 0
    files = manifest.check_files()
    for p, why in files:
        print(f"FAIL {p}: {why}")
    bad += len(files)
    for d, rec, got in manifest.check_stage_hashes():
        ok = rec == got
        print(f"{'OK  ' if ok else 'FAIL'} {d}: recorded {rec[:12]}  recomputed {got[:12]}")
        bad += not ok
    left = manifest.server_paths_left()
    for rel, ln, s in left:
        print(f"FAIL {rel}:{ln}: server path left in code: {s}")
    bad += len(left)
    n = len(manifest.read_manifest())
    print(f"\n{n} files in provenance/CODE_MANIFEST.tsv; {'all consistent' if not bad else f'{bad} problems'}")
    return 1 if bad else 0


def cmd_verify_release(a):
    ok = verify.report(verify.verify_release(a.release, a.config))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    argv, rest = _split_extra(argv)
    ap = argparse.ArgumentParser(prog="patbench", description=__doc__)
    ap.add_argument("--version", action="version", version=f"patbench {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p, out=False):
        p.add_argument("--workspace", help="the data workspace (default: $PATBENCH_WORKSPACE)")
        p.add_argument("--python", help="interpreter for the scripts (default: the current one)")
        p.add_argument("--workers", type=int, default=20)
        p.add_argument("--nice", action="store_true", help="run with nice -n 19")
        p.add_argument("--dry-run", action="store_true", help="print the commands, run nothing")
        if out:
            p.add_argument("--out", help="the release directory")

    sub.add_parser("stages", help="the processing stages, in order").set_defaults(fn=cmd_stages)
    sub.add_parser("steps", help="the runnable steps").set_defaults(fn=cmd_steps)
    p = sub.add_parser("plan", help="print the release chain")
    common(p, True)
    p.add_argument("--check", action="store_true", help="also check the workspace inputs")
    p.set_defaults(fn=cmd_plan)
    p = sub.add_parser("check-workspace", help="are the input tables where the scripts expect them?")
    p.add_argument("--workspace")
    p.add_argument("--steps", help="comma separated step names (default: all)")
    p.set_defaults(fn=cmd_check)
    p = sub.add_parser("run", help="run one step; arguments for the script go after --")
    p.add_argument("step")
    common(p, True)
    p.set_defaults(fn=cmd_run, takes_rest=True)
    p = sub.add_parser("exec", help="run any script under pipeline/ with the right environment; arguments after --")
    p.add_argument("script", help="path relative to pipeline/, e.g. 05_target_map/v2/v2_universe.py")
    common(p)
    p.set_defaults(fn=cmd_exec, takes_rest=True)
    p = sub.add_parser("release", help="build a release directory: standardize, novelty, build, protocols, eval-sets, baselines, leaderboard")
    common(p)
    p.add_argument("--out", required=True)
    p.add_argument("--reuse-work", help="directory holding the public_pairs and alias_scan parquet files to copy")
    p.add_argument("--skip", help="comma separated steps to skip")
    p.set_defaults(fn=cmd_release)
    sub.add_parser("verify-code", help="repository files against provenance/CODE_MANIFEST.tsv and the recorded stage hashes"
                   ).set_defaults(fn=cmd_verify_code)
    p = sub.add_parser("verify-release", help="integrity checks on a built release directory")
    p.add_argument("release")
    p.add_argument("--config", help="parameters to compare with (default: configs/release_v1.0.json)")
    p.set_defaults(fn=cmd_verify_release)

    a = ap.parse_args(argv)
    if getattr(a, "takes_rest", False):
        return a.fn(a, rest)
    if rest:
        ap.error("arguments after -- are only for `run` and `exec`")
    return a.fn(a)


if __name__ == "__main__":
    raise SystemExit(main())
