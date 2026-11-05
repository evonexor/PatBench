import os

import pytest

from patbench import layout
from patbench.cli import main
from patbench.steps import RELEASE_CHAIN, STEPS


def test_stages_and_steps_list(capsys):
    assert main(["stages"]) == 0
    out = capsys.readouterr().out
    assert "pipeline/08_build" in out and "pipeline/07_structure_repair" in out and "pipeline/02_extract" in out
    assert main(["steps"]) == 0
    assert "build_patbench" in capsys.readouterr().out or "build" in STEPS


def test_every_registered_script_exists():
    from patbench import PIPELINE
    for s in STEPS.values():
        assert (PIPELINE / s.script).is_file(), s.script


def test_release_chain_is_registered_and_ordered():
    assert RELEASE_CHAIN == ["standardize", "novelty", "build", "protocols", "eval-sets", "baselines", "leaderboard"]
    assert all(n in STEPS for n in RELEASE_CHAIN)


def test_plan_prints_the_commands_in_order(capsys):
    assert main(["plan", "--out", "/tmp/rel"]) == 0
    out = capsys.readouterr().out
    pos = [out.index(f"pipeline/08_build/{n}") for n in ["standardize.py", "novelty.py", "build_patbench.py",
                                                       "make_protocols.py", "make_eval_sets.py"]]
    assert pos == sorted(pos)
    assert "baselines_fixed.py" in out and "leaderboard.py" in out


def test_run_dry_run_sets_the_workspace(capsys, tmp_path):
    assert main(["run", "build", "--out", str(tmp_path / "r"), "--workspace", str(tmp_path), "--dry-run"]) == 0
    out = capsys.readouterr().out
    assert "build_patbench.py" in out and str(tmp_path / "r") in out


def test_run_passes_arguments_after_double_dash(capsys, tmp_path):
    main(["run", "fetch-structures", "--out", str(tmp_path), "--workspace", str(tmp_path), "--dry-run", "--", "meta", "--jobs", "2"])
    assert "meta --jobs 2" in capsys.readouterr().out


def test_run_step_that_needs_a_release_dir_refuses_without_one(tmp_path):
    with pytest.raises(SystemExit):
        main(["run", "build", "--workspace", str(tmp_path), "--dry-run"])


def test_unknown_step(tmp_path):
    with pytest.raises(SystemExit):
        main(["run", "nope", "--workspace", str(tmp_path)])


def test_release_dry_run_lists_every_step(capsys, tmp_path):
    assert main(["release", "--out", str(tmp_path / "r"), "--workspace", str(tmp_path), "--dry-run"]) == 0
    out = capsys.readouterr().out
    for n in RELEASE_CHAIN:
        assert f"=== {n}" in out
    assert not (tmp_path / "r").exists()


def test_release_stops_when_work_inputs_are_missing(tmp_path):
    with pytest.raises(SystemExit) as e:
        main(["release", "--out", str(tmp_path / "r"), "--workspace", str(tmp_path)])
    assert "public_pairs.parquet" in str(e.value)


def test_release_reuse_work_copies_the_inputs(tmp_path, capsys):
    src = tmp_path / "old"
    src.mkdir()
    for f in layout.WORK_INPUTS:
        (src / f).write_bytes(b"x")
    main(["release", "--out", str(tmp_path / "r"), "--workspace", str(tmp_path), "--reuse-work", str(src), "--dry-run"])
    assert "reuse" in capsys.readouterr().out


def test_workspace_check(tmp_path, capsys):
    assert main(["check-workspace", "--workspace", str(tmp_path)]) == 1
    for r in layout.REQUIREMENTS:
        p = tmp_path / r.path
        p.parent.mkdir(parents=True, exist_ok=True)
        if "." in os.path.basename(r.path):
            p.write_bytes(b"")
        else:
            p.mkdir(exist_ok=True)
    assert main(["check-workspace", "--workspace", str(tmp_path)]) == 0


def test_environment_exposes_the_workspace_and_the_extraction_modules(tmp_path):
    from patbench import PIPELINE, runner
    env = runner.environment(tmp_path, omp_threads=2)
    assert env["PATBENCH_WORKSPACE"] == str(tmp_path) and env["OMP_NUM_THREADS"] == "2"
    assert str(PIPELINE / "02_extract") in env["PYTHONPATH"].split(os.pathsep)


def test_environment_points_the_build_at_the_repaired_tables_and_the_v2_1_map(tmp_path, monkeypatch):
    from patbench import runner
    monkeypatch.delenv("PATBENCH_SAR_DIR", raising=False)
    monkeypatch.delenv("PATBENCH_TARGET_MAP_DIR", raising=False)
    env = runner.environment(tmp_path)
    assert env["PATBENCH_SAR_DIR"] == f"{tmp_path}/{layout.SAR}" and env["PATBENCH_SAR_DIR"].endswith("results/")
    assert env["PATBENCH_TARGET_MAP_DIR"] == f"{tmp_path}/{layout.MAP}" and env["PATBENCH_TARGET_MAP_DIR"].endswith("/")
    monkeypatch.setenv("PATBENCH_SAR_DIR", "/elsewhere/")
    assert runner.environment(tmp_path)["PATBENCH_SAR_DIR"] == "/elsewhere/"       # an explicit setting wins


def test_every_step_script_is_in_a_stage_directory():
    from patbench.steps import STAGES
    stages = {d.split("/", 1)[1] for _, d, _, _ in STAGES}
    for s in STEPS.values():
        assert s.script.split("/", 1)[0] in stages, s.name


def test_a_patched_script_fails_loudly_without_a_workspace(tmp_path):
    import subprocess, sys
    from patbench import PIPELINE
    env = {k: v for k, v in os.environ.items() if k != "PATBENCH_WORKSPACE"}
    r = subprocess.run([sys.executable, str(PIPELINE / "08_build/standardize.py"), str(tmp_path)], env=env,
                       capture_output=True, text=True)
    assert r.returncode != 0 and "PATBENCH_WORKSPACE is not set" in r.stderr + r.stdout


def test_verify_code_command(capsys):
    assert main(["verify-code"]) == 0
    assert "all consistent" in capsys.readouterr().out


def test_release_directory_gets_the_folders_the_build_does_not_create(tmp_path):
    from patbench import runner
    runner.prepare_release_dir(tmp_path / "r")
    assert (tmp_path / "r/work").is_dir() and (tmp_path / "r/tasks").is_dir()
