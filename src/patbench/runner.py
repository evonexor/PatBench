"""Run pipeline scripts with the environment they expect."""
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

from patbench import PIPELINE
from patbench.layout import MAP, SAR, WORK_INPUTS
from patbench.steps import RELEASE_CHAIN, STEPS, Ctx


def environment(ws, omp_threads=1, extra_pythonpath=()):
    """PATBENCH_WORKSPACE, the two input directories of the build, and a PYTHONPATH for the stage scripts.

    The extraction modules (uspto_sar, xml_source, ...) are imported by flat name from the stages
    that build on them, exactly as on the original server, so 02_extract goes on the path."""
    env = dict(os.environ)
    env["PATBENCH_WORKSPACE"] = str(ws)
    env.setdefault("PATBENCH_JOBS", "4")
    # the build scripts take these two directories from the environment (trailing slash included)
    env.setdefault("PATBENCH_SAR_DIR", str(ws).rstrip("/") + "/" + SAR)
    env.setdefault("PATBENCH_TARGET_MAP_DIR", str(ws).rstrip("/") + "/" + MAP)
    env["OMP_NUM_THREADS"] = str(omp_threads)
    parts = [str(PIPELINE / "02_extract"), *map(str, extra_pythonpath)]
    if env.get("PYTHONPATH"):
        parts.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = os.pathsep.join(parts)
    return env


def command(script, args, python=None):
    path = PIPELINE / script
    if not path.is_file():
        raise SystemExit(f"no such pipeline script: {script}")
    return [python or sys.executable, str(path), *map(str, args)]


def run(cmd, env, nice=False, dry_run=False, log=None):
    shown = " ".join(shlex.quote(c) for c in cmd)
    print("$", shown, flush=True)
    if dry_run:
        return 0
    if nice and shutil.which("nice"):
        cmd = ["nice", "-n", "19", *cmd]
    if log:
        Path(log).parent.mkdir(parents=True, exist_ok=True)
        with open(log, "ab") as fh:
            fh.write(f"\n$ {shown}\n".encode())
            return subprocess.call(cmd, env=env, stdout=fh, stderr=subprocess.STDOUT)
    return subprocess.call(cmd, env=env)


def run_step(name, ws, out="", extra=(), python=None, workers=20, nice=False, dry_run=False, log=None):
    step = STEPS[name]
    if step.needs_out and not out:
        raise SystemExit(f"step '{name}' needs --out <release directory>")
    ctx = Ctx(out=str(Path(out).resolve()) if out else "", workers=workers)
    cmd = command(step.script, step.args(ctx) + list(extra), python)
    return run(cmd, environment(ws, step.omp_threads), nice=nice, dry_run=dry_run, log=log)


def ensure_work_inputs(out, reuse_from=None, dry_run=False):
    """build_patbench.py reads three files from <release>/work/ that other steps produce.

    Returns the list still missing. With `reuse_from`, copies them from that directory. public_pairs.parquet
    depends only on the public databases; alias_scan.parquet is written by `alias-scan-final` after the
    structure repair (docs/REPRODUCE.md)."""
    work = Path(out) / "work"
    missing = [f for f in WORK_INPUTS if not (work / f).exists()]
    if missing and reuse_from:
        for f in list(missing):
            src = Path(reuse_from) / f
            if src.exists():
                print(f"reuse {src} -> {work / f}")
                if not dry_run:
                    work.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(src, work / f)
                missing.remove(f)
    return missing


def prepare_release_dir(out):
    """build_patbench.py writes tasks/*.parquet without creating tasks/, and reads work/ for its inputs."""
    for d in ("work", "tasks"):
        (Path(out) / d).mkdir(parents=True, exist_ok=True)


def release(ws, out, reuse_work=None, python=None, workers=20, nice=False, dry_run=False, skip=()):
    out = Path(out).resolve()
    if not dry_run:
        prepare_release_dir(out)
    missing = ensure_work_inputs(out, reuse_work, dry_run)
    if missing and not dry_run:
        raise SystemExit(
            f"{out}/work lacks {', '.join(missing)}.\n"
            "Produce them with `patbench run public-pairs` and `patbench run alias-scan-final`, or copy them "
            "from an earlier release with --reuse-work <dir> (docs/REPRODUCE.md).")
    log = out / "chain.log"
    for name in RELEASE_CHAIN:
        if name in skip:
            continue
        print(f"=== {name}", flush=True)
        rc = run_step(name, ws, str(out), python=python, workers=workers, nice=nice, dry_run=dry_run,
                      log=None if dry_run else out / f"{name}.log")
        if rc:
            raise SystemExit(f"step '{name}' failed (exit {rc}); see {out}/{name}.log")
    return 0
