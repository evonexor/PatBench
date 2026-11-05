"""Import the PatBench processing code from a snapshot of the original server layout.

The code in `pipeline/` was written and run on one server, in several directories, one per stage.
This tool is how those files became this repository, so the step can be repeated and audited:

  * `verbatim` files are copied byte for byte. Three of them are frozen stage directories whose
    SHA-256 was recorded by the original runs; `patbench verify-code` recomputes those hashes.
  * `patched` files differ from the original only in how a path is written: a string literal that
    began with the server's workspace prefix now begins with `_ws`, which is read from the
    environment variable PATBENCH_WORKSPACE (a line defining `_ws` is added after the module
    docstring; lower case on purpose: build_patbench.py writes every upper-case global into release_manifest.json). Nothing else is edited. Every changed file is listed in provenance/CODE_MANIFEST.tsv
    with the SHA-256 of the original and of the repository copy.

    PATBENCH_IMPORT_SERVER_PREFIX=<absolute workspace path of the server> PATBENCH_IMPORT_PUBLIC_PREFIX=<absolute path of the public-data directory> \
    python tools/import_pipeline.py <snapshot_dir>

The snapshot directory mirrors the server: `outputs/<stage directory>/...`, `raw_std/` (the public
data standardisation scripts), `ops/` (one-off operations scripts). See SNAPSHOT_LAYOUT below.
"""
import argparse
import ast
import hashlib
import io
import os
import re
import shutil
import sys
import tokenize
from pathlib import Path

# the two absolute prefixes of the original server are not recorded in this repository: give them through the environment (the defaults are placeholders)
SERVER_PREFIX = os.environ.get("PATBENCH_IMPORT_SERVER_PREFIX", "/path/to/server-workspace")
PUBLIC_PREFIX = os.environ.get("PATBENCH_IMPORT_PUBLIC_PREFIX", "/path/to/server-public-data")      # where 01_public_data originally read and wrote
V1_SCRIPTS = SERVER_PREFIX + "/outputs/target-uniprot-map-20261016/scripts"
V2_SCRIPTS = SERVER_PREFIX + "/outputs/target-uniprot-map-v2-20261019/scripts"
HERE = "__import__('os').path.dirname(__import__('os').path.abspath(__file__))"
HERE_V1 = "__import__('os').path.join(" + HERE + ", '..', 'v1')"
WS_DEF = ('_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit('
          '"PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")')

R6 = "outputs/patbench-v1.0-20261030-r6/scripts/"
ALIAS = "outputs/alias-expansion-20261031/scripts/"
RANK = "outputs/affinity-rank-benchmark-20261019/scripts/"
NITRILE = "outputs/patbench-v1.0-20261030-r6/structure_repair/nitrile_audit/"


def _r6(names, dst):
    return [(R6 + n, dst + n, "patched") for n in names]


# (source relative to the snapshot, destination relative to the repository, mode)
# mode: verbatim = byte for byte; patched = path literals rewritten; dir_* = every file of a directory
SNAPSHOT_LAYOUT = (
    [
        ("ops/copy_patent_data_20261001.py", "pipeline/00_ops/copy_patent_data.py", "verbatim"),
        ("outputs/uspto-sar-recovery-20261006/launch_service.py", "pipeline/00_ops/launch_recovery_service.py", "verbatim"),
        ("raw_std", "pipeline/01_public_data", "dir_patched"),
        ("outputs/uspto-sar-full-20261004/code", "pipeline/02_extract", "dir_verbatim"),
        ("outputs/uspto-sar-recovery-20261006/code", "pipeline/03_recovery", "dir_verbatim"),
        ("outputs/uspto-all-structures-20261004/code", "pipeline/04_catalog", "dir_verbatim"),
        ("outputs/target-uniprot-map-20261016/scripts", "pipeline/05_target_map/v1", "dir_patched"),
        ("outputs/target-uniprot-map-v2-20261019/scripts", "pipeline/05_target_map/v2", "dir_patched"),
        ("outputs/target-uniprot-map-v2.1-20261031/scripts", "pipeline/05_target_map/v2_1", "dir_patched"),
        ("outputs/sar-validation-20261019/scripts", "pipeline/06_validation", "dir_patched"),
    ]
    # 07: structure repair
    + _r6(["alias_extract.py", "alias_scan_all.py", "alias_scan.py", "alias_pubkeys.py", "alias_repair.py",
           "alias_expand_independent.py", "alias_validate2.py", "merge_final.py", "make_patched_sar.py",
           "make_alias_scan_r6.py"], "pipeline/07_structure_repair/")
    + [(ALIAS + "list_unresolved_labels.py", "pipeline/07_structure_repair/list_unresolved_labels.py", "patched")]
    # evidence audit of the reversed-nitrile conversion (scripts stored inside the r6 release directory)
    + [(NITRILE + "scripts", "pipeline/07_structure_repair/nitrile_audit", "dir_patched"),
       (NITRILE + "verified_nitrile_guard.py", "pipeline/07_structure_repair/nitrile_audit/verified_nitrile_guard.py", "verbatim")]
    # 08: build, protocols, target audit
    + _r6(["standardize.py", "novelty.py", "build_patbench.py", "join_target_map.py", "make_target_audit.py",
           "audit_targets.py", "make_protocols.py", "make_eval_sets.py"], "pipeline/08_build/")
    # public_pairs.parquet, which build_patbench.py reads from work/, is written by this script
    + [(RANK + "public_overlap.py", "pipeline/08_build/public_overlap.py", "patched")]
    + _r6(["fetch_structures.py", "check_structures.py"], "pipeline/09_structures/")
    # 10: evaluation. evaluate_patent_bootstrap.py is the evaluator named as recommended in the freeze record of r6
    # (patent-cluster bootstrap, paired comparison); baselines_fixed.py is baselines.py with the few-shot overwrite removed
    + _r6(["evaluate.py", "baselines.py", "leaderboard.py"], "pipeline/10_evaluation/")
    + [("frozen_run/evaluate.py", "pipeline/10_evaluation/evaluate_patent_bootstrap.py", "verbatim"),
       ("fewshot/baselines_fixed.py", "pipeline/10_evaluation/baselines_fixed.py", "patched")]
    + _r6(["build_eval_support.py", "build_target_names.py", "make_review_table.py", "profile_r6.py", "profile2_r6.py"],
          "pipeline/11_support/")
    # the first (2026-10-19) release and benchmark builders; PatBench v1.0 grew out of them
    + [("outputs/affinity-sar-release-20261019/scripts/build_release.py", "pipeline/legacy_v0/build_release.py", "patched")]
    + [(RANK + n, "pipeline/legacy_v0/benchmark/" + n, "patched")
       for n in ["baselines.py", "build_benchmark.py", "evaluate.py", "novelty.py", "rebuild_slices.py"]]
)

# server scripts that are deliberately not in the repository: (server path, path in the snapshot or None, reason);
# written to provenance/NOT_INCLUDED.tsv
_P = SERVER_PREFIX + "/"
NOT_INCLUDED = [
    (_P + R6 + "freeze_r6.py", R6 + "freeze_r6.py", "release housekeeping: writes the freeze block and notes into the data release README (Chinese text)"),
    (_P + R6 + "freeze_r6.sh", R6 + "freeze_r6.sh", "release housekeeping: SHA256SUMS and the read-only switch of the r6 directory (see docs/RELEASES.md)"),
    (_P + R6 + "note_unresolved_r6.py", R6 + "note_unresolved_r6.py", "release housekeeping: appends the extension-trial note to the data release README"),
    (_P + R6 + "stamp_baselines_r6.py", R6 + "stamp_baselines_r6.py", "release housekeeping: writes baseline results into the data release README and manifest"),
    (_P + R6 + "stamp_revision.py", R6 + "stamp_revision.py", "release housekeeping: adds a revision note to release_manifest.json"),
    (_P + ALIAS + "alias_expand.py", ALIAS + "alias_expand.py", "superseded: first rule-based expander, replaced by alias_repair.py and alias_expand_independent.py"),
    (_P + ALIAS + "alias_validate.py", ALIAS + "alias_validate.py", "superseded by alias_validate2.py"),
    (_P + ALIAS + "make_expand_v7.py", ALIAS + "make_expand_v7.py", "trial rule extension that was not merged into r6 (see docs/RELEASES.md)"),
    (_P + ALIAS + "draw_labels.py", ALIAS + "draw_labels.py", "visual inspection helper for alias labels"),
    (_P + "outputs/patbench-v1.0-20261030-r4/scripts/carry_over_target_files.py", None, "r4 only: reused the r3 target files; not used by r6"),
    ("runs/r6_frozen_20261101/scripts/ (all files except evaluate.py)", None,
     "model-run harness of the frozen evaluation run (adapters, docking and DeepPurpose workers, scoring); not part of the data pipeline"),
    (_P + "outputs/patbench-r6-fewshot-baselines-20261101/scripts/ (score.py, verify.py, run.sh, evaluate_patent.py)", None,
     "scoring and verification helpers of the corrected few-shot baseline run; evaluate_patent.py is byte-identical to evaluate_patent_bootstrap.py"),
]

# recorded by the original runs (job_config.json / completion.json of each stage directory)
STAGE_HASHES = {
    "pipeline/02_extract": ("148c103fb46fd807b6d7f64613ac36b3ece45528112a02d4097a6a1184b407f7", "glob_bytes",
                            "uspto-sar-full-20261004, version 0.8.0-full"),
    "pipeline/03_recovery": ("81fd68912e6f5bb92bb3629ba3fbc8f1ba68e40ebd9431f91598d1b8a3b89fd2", "rglob_named",
                             "uspto-sar-recovery-20261006"),
    "pipeline/04_catalog": ("67b554fd4599c2cde756f0689263214bb22fdb97ef46e7631c49970d4f7e0214", "glob_bytes",
                            "uspto-all-structures-20261004, version 0.9.1-sharded-aggregation"),
}


def sha(b):
    return hashlib.sha256(b).hexdigest()


def patch_source(text, rel):
    """Rewrite string literals that start with the server prefix. Returns (new_text, n_changes)."""
    toks = list(tokenize.generate_tokens(io.StringIO(text).readline))
    lines = text.splitlines(keepends=True)
    offs = [0]
    for ln in lines:
        offs.append(offs[-1] + len(ln))
    pos = lambda rc: offs[rc[0] - 1] + rc[1]
    edits = []
    in_v2 = "/05_target_map/v2/" in rel or rel.startswith("pipeline/05_target_map/v2/")
    for i, t in enumerate(toks):
        if t.type == tokenize.STRING:
            m = re.fullmatch(r"([rRbBuU]*)(['\"])(.*)\2", t.string, flags=re.S)
            if not m or t.string.lstrip("rRbBuU").startswith(('"""', "'''")):
                continue
            pre, q, body = m.groups()
            public = body.startswith(PUBLIC_PREFIX)
            if not (body.startswith(SERVER_PREFIX) or public):
                continue
            if public:
                # the public-database scripts were written on another machine, where the raw data lived here
                new = f"_ws + {pre}{q}/raw/patent-data{body[len(PUBLIC_PREFIX):]}{q}"
            elif body == V1_SCRIPTS:
                new = HERE_V1 if in_v2 else HERE
            elif body == V2_SCRIPTS:
                new = HERE
            else:
                new = f"_ws + {pre}{q}{body[len(SERVER_PREFIX):]}{q}"
            edits.append((pos(t.start), pos(t.end), new))
        elif t.type == getattr(tokenize, "FSTRING_MIDDLE", -1) and t.string.startswith(SERVER_PREFIX):
            # a path at the start of an f-string: 'f"<server prefix>..."' -> 'f"{WS}..."'
            if toks[i - 1].type == getattr(tokenize, "FSTRING_START", -1) and len(toks[i - 1].string) <= 3:
                edits.append((pos(t.start), pos(t.start) + len(SERVER_PREFIX), "{_ws}"))
    if not edits:
        return text, 0
    for a, b, new in sorted(edits, reverse=True):
        text = text[:a] + new + text[b:]
    if in_v2:
        # v2_assemble put the v1 *data* directory plus 'scripts' on sys.path; the v1 scripts live
        # beside the v2 ones in this repository
        text, k = re.subn(r"\bV1\s*\+\s*(['\"])scripts\1", lambda _: HERE_V1, text)
        edits += [None] * k
    # define WS after the module docstring and any __future__ imports
    tree = ast.parse(text)
    at = 0
    for node in tree.body:
        is_doc = isinstance(node, ast.Expr) and isinstance(getattr(node, "value", None), ast.Constant) \
            and isinstance(node.value.value, str)
        if is_doc or (isinstance(node, ast.ImportFrom) and node.module == "__future__"):
            at = node.end_lineno
            continue
        break
    ls = text.splitlines(keepends=True)
    ls.insert(at, WS_DEF + "\n")
    return "".join(ls), len(edits)


def hash_stage(path, how):
    files = sorted(path.glob("*.py")) if how == "glob_bytes" else sorted(path.rglob("*.py"))
    if how == "glob_bytes":
        return sha(b"".join(p.read_bytes() for p in files))
    return sha(b"".join(str(p.relative_to(path)).encode() + b"\0" + p.read_bytes() for p in files))


def main():
    if SERVER_PREFIX == "/path/to/server-workspace" or PUBLIC_PREFIX == "/path/to/server-public-data":
        sys.exit("set PATBENCH_IMPORT_SERVER_PREFIX and PATBENCH_IMPORT_PUBLIC_PREFIX (see the docstring of this file)")
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("snapshot")
    ap.add_argument("--repo", default=str(Path(__file__).resolve().parent.parent))
    a = ap.parse_args()
    snap, repo = Path(a.snapshot).resolve(), Path(a.repo).resolve()
    skip = {"__pycache__", ".DS_Store"}
    rows = []
    for src, dst, mode in SNAPSHOT_LAYOUT:
        s = snap / src
        pairs = []
        if mode.startswith("dir_"):
            for f in sorted(p for p in s.rglob("*") if p.is_file() and not (set(p.parts) & skip)
                            and not p.name.startswith("._") and p.suffix != ".pyc"):
                pairs.append((f, repo / dst / f.relative_to(s)))
        else:
            pairs.append((s, repo / dst))
        for f, d in pairs:
            raw = f.read_bytes()
            body, n = raw, 0
            is_py = f.suffix == ".py"
            if mode.endswith("patched") and is_py:
                text, n = patch_source(raw.decode("utf-8"), str(d.relative_to(repo)))
                body = text.encode("utf-8")
            d.parent.mkdir(parents=True, exist_ok=True)
            d.write_bytes(body)
            rows.append((str(d.relative_to(repo)), src if not mode.startswith("dir_") else str(f.relative_to(snap)),
                         "patched-paths" if n else "verbatim", n, sha(raw), sha(body)))
    # files written in this repository (not copied from the server) are listed too, so that the manifest
    # covers every file under pipeline/
    have = {r[0] for r in rows}
    for f in sorted((repo / "pipeline").rglob("*")):
        rel = str(f.relative_to(repo))
        if f.is_file() and rel not in have and "__pycache__" not in f.parts and f.suffix != ".pyc" \
                and not f.name.startswith("._") and f.name != ".DS_Store":
            rows.append((rel, "-", "authored", 0, "-", sha(f.read_bytes())))
    # frozen stage directories: the recorded hash must come out of the copied files
    for d, (rec, how, what) in STAGE_HASHES.items():
        got = hash_stage(repo / d, how)
        print(f"{d}: recorded {rec[:12]} recomputed {got[:12]} {'OK' if got == rec else 'MISMATCH'}")
        if got != rec:
            sys.exit(f"{d} does not reproduce its recorded code hash")
    out = repo / "provenance" / "CODE_MANIFEST.tsv"
    out.parent.mkdir(exist_ok=True)
    import json
    (repo / "provenance" / "STAGE_HASHES.json").write_text(json.dumps(
        {d: {"sha256": h, "formula": how, "run": what} for d, (h, how, what) in STAGE_HASHES.items()},
        indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    with out.open("w") as fh:
        fh.write("repo_path\tsnapshot_path\tstatus\tpath_edits\tsha256_original\tsha256_repo\n")
        for r in sorted(rows):
            fh.write("\t".join(map(str, r)) + "\n")
    with (repo / "provenance" / "NOT_INCLUDED.tsv").open("w") as fh:
        fh.write("path_in_workspace\tsha256\treason\n")
        for server, rel, why in NOT_INCLUDED:
            f = snap / rel if rel else None
            server = server[len(SERVER_PREFIX) + 1:] if server.startswith(SERVER_PREFIX + "/") else server      # relative to the workspace
            fh.write(f"{server}\t{sha(f.read_bytes()) if f is not None and f.is_file() else '-'}\t{why}\n")
    patched = sum(1 for r in rows if r[2] == "patched-paths")
    authored = sum(1 for r in rows if r[2] == "authored")
    print(f"{len(rows)} files in the manifest: {patched} with rewritten path literals ({sum(r[3] for r in rows)} literals), "
          f"{authored} authored here, {len(rows) - patched - authored} verbatim")


if __name__ == "__main__":
    main()
