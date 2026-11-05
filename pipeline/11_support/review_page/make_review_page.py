"""Build the target-mapping review page: one self-contained HTML file a person can open in a browser.

    python make_review_page.py <release_dir> [--out target_review.html]

The page lists every test group of the five ranking tasks. For each group a reviewer sees the label and
column header printed in the patent next to the protein the pipeline mapped it to, and records one
verdict (correct / wrong / unsure), an optional corrected UniProt accession and a note. Nothing is
pre-filled. The reviewer exports a TSV with the original 20 columns of target_review_test.tsv plus
`reviewer_corrected_uniprot`.

Inputs (all in the release directory):
  target_review_test.tsv              the test groups of rank_ic50 and rank_ki_kd, reviewer columns blank
  target_audit.parquet, groups.parquet, target_sequences.parquet
  eval_support/test_group_profile.parquet   (patbench run eval-support)

Test groups of the other tasks (rank_highcontrast, rank_complex_mutant, rank_ec50) that are not in
target_review_test.tsv get the same fields from the audit tables and are written to
eval_support/target_review_test_other_tasks.tsv, then added to the page.
"""
import argparse
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ORIG_COLS = ["group_id", "task", "tier", "patent_core", "patent_title", "assay_label_as_printed", "activity_header",
             "mapped_gene_or_identity", "protein_name", "organism", "uniprot", "confidence", "curator_note",
             "sequence_id", "sequence_kind", "seq_len", "label_check", "xref_status", "reviewer_verdict",
             "reviewer_note"]
AUDIT_EXTRA = ["mapping_kind", "curator_confidence", "sequence_caveat", "organism_stated", "organism_detected",
               "organism_resolution", "mutation", "xref_db", "xref_accession", "xref_symbol", "xref_value_match",
               "target_symbols", "target_class"]
GROUP_EXTRA = ["endpoint", "n_molecules", "pAct_range", "year", "first_public_year", "public_share"]


def clean(df):
    for c in df.columns:
        df[c] = df[c].astype(str).replace({"<NA>": "", "nan": "", "None": ""})
    return df


def build_records(root):
    root = Path(root)
    rv = pd.read_csv(root / "target_review_test.tsv", sep="\t", encoding="utf-8-sig", dtype=str, keep_default_na=False)
    au = pd.read_parquet(root / "target_audit.parquet")
    sq = pd.read_parquet(root / "target_sequences.parquet")
    gr = pd.read_parquet(root / "groups.parquet")
    prof = pd.read_parquet(root / "eval_support/test_group_profile.parquet")

    rv["group_id"] = rv.group_id.astype(int)
    main = rv.merge(au[["group_id"] + AUDIT_EXTRA], on="group_id", how="left").merge(
        gr[["group_id"] + GROUP_EXTRA], on="group_id", how="left")
    assert len(main) == len(rv) and main.group_id.is_unique

    extra = prof[~prof.group_id.isin(set(rv.group_id))]
    ids = sorted(extra.group_id.unique())
    tasks = extra.groupby("group_id").task.apply(lambda s: ";".join(sorted(s)))
    tier = extra.sort_values("task").groupby("group_id").test_tier.first()
    a, s, g = au.set_index("group_id").loc[ids], sq.set_index("group_id").loc[ids], gr.set_index("group_id").loc[ids]
    oth = pd.DataFrame({
        "group_id": ids, "task": tasks.loc[ids].values, "tier": tier.loc[ids].values,
        "patent_core": a.patent_core.values, "patent_title": a.patent_title.values,
        "assay_label_as_printed": a.Target.values, "activity_header": a.activity_header.values,
        "mapped_gene_or_identity": a.mapped_identity.values, "protein_name": a.protein_name.values,
        "organism": a.organism.values, "uniprot": a.uniprot_accession.values, "confidence": s.confidence.values,
        "curator_note": a.curator_note.values, "sequence_id": s.sequence_id.values,
        "sequence_kind": s.sequence_kind.values, "seq_len": s.seq_len.values, "label_check": s.label_check.values,
        "xref_status": s.xref_status.values, "reviewer_verdict": "", "reviewer_note": "",
        **{c: a[c].values for c in AUDIT_EXTRA}, **{c: g[c].values for c in GROUP_EXTRA}})
    return clean(main.copy()), clean(oth)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("release")
    ap.add_argument("--out", default=None, help="HTML file to write (default: <release>/review/target_review.html)")
    a = ap.parse_args()
    root = Path(a.release)
    main_df, oth = build_records(root)
    assert not set(main_df.group_id.astype(str)) & set(oth.group_id.astype(str))
    oth.iloc[:, :20].to_csv(root / "eval_support/target_review_test_other_tasks.tsv", sep="\t", index=False,
                            encoding="utf-8-sig")
    rows = main_df.to_dict("records") + [dict(r, extra="1") for r in oth.to_dict("records")]
    assert all(r["reviewer_verdict"] == "" and r["reviewer_note"] == "" for r in rows), "reviewer columns must start blank"
    data = json.dumps(rows, ensure_ascii=False, separators=(",", ":")).replace("<", "\\u003c")
    tpl = (HERE / "template.html").read_text(encoding="utf-8")
    assert tpl.count("__DATA__") == 1
    out = Path(a.out) if a.out else root / "review/target_review.html"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(tpl.replace("__DATA__", data), encoding="utf-8")
    print(f"{len(rows)} groups ({len(main_df)} from target_review_test.tsv, {len(oth)} from other tasks) -> {out}")


if __name__ == "__main__":
    main()
