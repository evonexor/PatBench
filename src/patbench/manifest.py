"""Code integrity: the repository copy of every pipeline file against the recorded hashes."""
import hashlib
import io
import json
import tokenize
from pathlib import Path

from patbench import PIPELINE, PROVENANCE, REPO

SERVER_PREFIXES = ("/home/",)
# directories that are byte-for-byte what ran on the server; literals there are part of the record
FROZEN = ("pipeline/00_ops", "pipeline/02_extract", "pipeline/03_recovery", "pipeline/04_catalog")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_manifest():
    rows = []
    with open(PROVENANCE / "CODE_MANIFEST.tsv", encoding="utf-8") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            rows.append(dict(zip(head, line.rstrip("\n").split("\t"))))
    return rows


def hash_stage(path, formula):
    """The two formulas the original runs used for `code_sha256`."""
    path = Path(path)
    # `._name` files are macOS AppleDouble metadata that tar leaves behind on other systems; the
    # original directories had none, so ignoring them changes nothing there and keeps a copied
    # checkout verifiable
    real = lambda ps: [p for p in sorted(ps) if not p.name.startswith("._")]
    if formula == "glob_bytes":              # run_full.py, run_structure_catalog.py
        return hashlib.sha256(b"".join(p.read_bytes() for p in real(path.glob("*.py")))).hexdigest()
    if formula == "rglob_named":             # launch_service.py of the recovery run
        return hashlib.sha256(b"".join(
            str(p.relative_to(path)).encode() + b"\0" + p.read_bytes() for p in real(path.rglob("*.py")))).hexdigest()
    raise ValueError(formula)


def check_files():
    """[(path, problem)] for files that differ from, or are missing in, the manifest."""
    problems, seen = [], set()
    for r in read_manifest():
        p = REPO / r["repo_path"]
        seen.add(r["repo_path"])
        if not p.is_file():
            problems.append((r["repo_path"], "missing"))
        elif r["status"] != "authored" and sha256(p) != r["sha256_repo"]:
            # files written in this repository are expected to change; imported code is not
            problems.append((r["repo_path"], "content differs from the manifest"))
        elif r["status"] == "verbatim" and r["sha256_repo"] != r["sha256_original"]:
            problems.append((r["repo_path"], "marked verbatim but the two hashes differ"))
    for p in sorted(PIPELINE.rglob("*")):
        rel = str(p.relative_to(REPO))
        if (p.is_file() and rel not in seen and "__pycache__" not in p.parts
                and p.suffix != ".pyc" and not p.name.startswith("._") and p.name != ".DS_Store"):
            problems.append((rel, "not in the manifest"))
    return problems


def check_stage_hashes():
    """[(directory, recorded, recomputed)] for the stage directories whose original run recorded a hash."""
    rec = json.loads((PROVENANCE / "STAGE_HASHES.json").read_text(encoding="utf-8"))
    return [(d, v["sha256"], hash_stage(REPO / d, v["formula"])) for d, v in rec.items()]


def _docstring_lines(tree):
    import ast
    lines = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.ClassDef, ast.AsyncFunctionDef)):
            body = getattr(node, "body", [])
            if body and isinstance(body[0], ast.Expr) and isinstance(getattr(body[0], "value", None), ast.Constant) \
                    and isinstance(body[0].value.value, str):
                lines.update(range(body[0].lineno, body[0].end_lineno + 1))
    return lines


def server_paths_left():
    """String literals (outside docstrings) in the patched directories that still name the server."""
    import ast
    hits = []
    for p in sorted(PIPELINE.rglob("*.py")):
        if p.name.startswith("._"):
            continue
        rel = str(p.relative_to(REPO))
        if rel.startswith(FROZEN):
            continue
        text = p.read_text(encoding="utf-8")
        doc = _docstring_lines(ast.parse(text))
        for t in tokenize.generate_tokens(io.StringIO(text).readline):
            if t.type in (tokenize.STRING, getattr(tokenize, "FSTRING_MIDDLE", -1)) and any(s in t.string for s in SERVER_PREFIXES) \
                    and t.start[0] not in doc:
                hits.append((rel, t.start[0], t.string[:80]))
    return hits
