"""The review-page generator on a small synthetic release (fixture tables, not PatBench data)."""
import importlib.util
import json
import re

import pandas as pd
import pytest

from patbench import PIPELINE

spec = importlib.util.spec_from_file_location("make_review_page", PIPELINE / "11_support/review_page/make_review_page.py")
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)


def _release(tmp_path):
    root = tmp_path / "rel"
    (root / "eval_support").mkdir(parents=True)
    ids = [1, 2, 3, 4]       # 1-2 are in the review table; 3-4 are test groups of another task only
    rv = pd.DataFrame({c: "" for c in M.ORIG_COLS}, index=range(2))
    rv["group_id"] = ["1", "2"]
    rv["task"] = "rank_ic50"
    rv["uniprot"] = ["P00001", "P00002"]
    rv["patent_title"] = ["<b>bold</b> title", "plain"]
    rv.to_csv(root / "target_review_test.tsv", sep="\t", index=False, encoding="utf-8-sig")
    au = pd.DataFrame({"group_id": ids, "patent_core": [f"US{i}" for i in ids], "Target": [f"T{i}" for i in ids],
                       "activity_header": "", "mapped_identity": [f"G{i}" for i in ids], "protein_name": "Protein",
                       "organism": "Homo sapiens", "uniprot_accession": [f"P0000{i}" for i in ids],
                       "curator_note": "", "patent_title": "t", **{c: "" for c in M.AUDIT_EXTRA}})
    au.to_parquet(root / "target_audit.parquet", index=False)
    pd.DataFrame({"group_id": ids, "sequence_id": [f"P0000{i}" for i in ids], "sequence_kind": "full_length", "seq_len": 100,
                  "confidence": "medium", "label_check": "label_names_target", "xref_status": "agrees"}
                 ).to_parquet(root / "target_sequences.parquet", index=False)
    pd.DataFrame({"group_id": ids, **{c: 1 for c in M.GROUP_EXTRA}}).to_parquet(root / "groups.parquet", index=False)
    pd.DataFrame({"task": ["rank_ic50", "rank_ic50", "rank_complex_mutant", "rank_ec50", "rank_highcontrast"],
                  "group_id": [1, 2, 3, 4, 4], "test_tier": ["novel", "reference", "reference", "novel", "novel"]}
                 ).to_parquet(root / "eval_support/test_group_profile.parquet", index=False)
    return root


def test_records_split_into_reviewed_table_and_other_tasks(tmp_path):
    main, other = M.build_records(_release(tmp_path))
    assert list(main.group_id) == [1, 2] or list(main.group_id.astype(int)) == [1, 2]
    assert list(other.group_id.astype(int)) == [3, 4]
    assert other.set_index(other.group_id.astype(int)).loc[4, "task"] == "rank_ec50;rank_highcontrast"
    assert (other.reviewer_verdict == "").all() and (other.reviewer_note == "").all()


def test_page_embeds_every_group_and_keeps_markup_inert(tmp_path, monkeypatch):
    root = _release(tmp_path)
    out = tmp_path / "page.html"
    monkeypatch.setattr("sys.argv", ["make_review_page.py", str(root), "--out", str(out)])
    M.main()
    html = out.read_text(encoding="utf-8")
    m = re.search(r'<script id="data" type="application/json">(.*?)</script>', html, flags=re.S)
    rows = json.loads(m.group(1))
    assert [int(r["group_id"]) for r in rows] == [1, 2, 3, 4]
    assert [r.get("extra") for r in rows] == [None, None, "1", "1"]
    assert "<b>bold</b>" not in m.group(1) and "\\u003cb>bold" in m.group(1)       # '<' is escaped inside the script block
    assert rows[0]["patent_title"] == "<b>bold</b> title"                          # and decodes back to the original text
    assert "__DATA__" not in html
    other = pd.read_csv(root / "eval_support/target_review_test_other_tasks.tsv", sep="\t", encoding="utf-8-sig", dtype=str)
    assert list(other.columns) == M.ORIG_COLS and list(other.group_id) == ["3", "4"]


def test_the_template_is_self_contained():
    html = (PIPELINE / "11_support/review_page/template.html").read_text(encoding="utf-8")
    assert not re.search(r'(src|href)="https?://(?!patents\.google|www\.uniprot)', html)
    assert "<script src=" not in html
