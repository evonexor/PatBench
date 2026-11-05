"""The imported pipeline code is what the manifest says it is."""
import py_compile
from pathlib import Path

import pytest

from patbench import PIPELINE, REPO, manifest


def test_files_match_manifest():
    assert manifest.check_files() == []


@pytest.mark.parametrize("directory,recorded,recomputed", manifest.check_stage_hashes())
def test_frozen_stage_reproduces_the_hash_its_run_recorded(directory, recorded, recomputed):
    assert recomputed == recorded, directory


def test_no_server_path_left_in_rewritten_code():
    assert manifest.server_paths_left() == []


def test_every_python_file_compiles(tmp_path):
    for i, p in enumerate(sorted(q for q in PIPELINE.rglob("*.py") if not q.name.startswith("._"))):
        py_compile.compile(str(p), cfile=str(tmp_path / f"{i}.pyc"), doraise=True)


def test_recovery_vendor_copies_are_the_extraction_modules():
    for f in (PIPELINE / "03_recovery/vendor").glob("*.py"):
        assert f.read_bytes() == (PIPELINE / "02_extract" / f.name).read_bytes(), f.name


def test_catalog_stage_shares_the_extraction_modules():
    for name in ["uspto_sar.py", "xml_source.py", "patent_documents.py", "patent_metadata.py", "select_samples.py"]:
        assert (PIPELINE / "04_catalog" / name).read_bytes() == (PIPELINE / "02_extract" / name).read_bytes(), name


def test_release_config_equals_the_constants_of_the_build_script():
    """build_patbench.py writes every upper-case global into release_manifest.json; the config file must match them."""
    import ast
    import json
    ref = json.loads((REPO / "configs/release_v1.0.json").read_text())
    tree = ast.parse((PIPELINE / "08_build/build_patbench.py").read_text(encoding="utf-8"))
    consts = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name) \
                and node.targets[0].id.isupper():
            try:
                consts[node.targets[0].id] = ast.literal_eval(node.value)
            except (ValueError, SyntaxError):
                pass
    shared = {k for k in ref if k in consts}
    assert len(shared) >= 20, sorted(shared)
    assert {k: ref[k] for k in shared} == {k: json.loads(json.dumps(consts[k])) for k in shared}     # tuples read back as lists


def test_a_changed_imported_file_is_caught_and_a_changed_authored_file_is_not(tmp_path, monkeypatch):
    import shutil
    from patbench import manifest as M
    repo = tmp_path / "repo"
    shutil.copytree(REPO / "pipeline", repo / "pipeline")
    shutil.copytree(REPO / "provenance", repo / "provenance")
    monkeypatch.setattr(M, "REPO", repo)
    monkeypatch.setattr(M, "PIPELINE", repo / "pipeline")
    monkeypatch.setattr(M, "PROVENANCE", repo / "provenance")
    assert M.check_files() == []
    (repo / "pipeline/README.md").write_text("edited\n")                       # authored: allowed
    assert M.check_files() == []
    with open(repo / "pipeline/10_evaluation/evaluate.py", "a") as fh:         # imported: caught
        fh.write("# edit\n")
    assert [p for p, _ in M.check_files()] == ["pipeline/10_evaluation/evaluate.py"]
    (repo / "pipeline/10_evaluation/extra.py").write_text("x = 1\n")           # unregistered file: caught
    assert "pipeline/10_evaluation/extra.py" in [p for p, _ in M.check_files()]
