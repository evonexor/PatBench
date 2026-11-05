"""Manual-review tables for the time-split test groups, rebuilt from target_audit.parquet / target_sequences.parquet / tasks.

    python make_review_table.py <release_dir> <out_dir>

target_review_test.tsv             test groups of rank_ic50 and rank_ki_kd (one row per group)
target_review_test_other_tasks.tsv test groups of rank_highcontrast, rank_complex_mutant, rank_ec50 that are not in the first table
reviewer_verdict and reviewer_note are left blank.
"""
import sys
import pandas as pd
R, OUT = sys.argv[1].rstrip("/"), sys.argv[2].rstrip("/")
au = pd.read_parquet(f"{R}/target_audit.parquet").set_index("group_id")
sq = pd.read_parquet(f"{R}/target_sequences.parquet").set_index("group_id")
COLS = ["group_id", "task", "tier", "patent_core", "patent_title", "assay_label_as_printed", "activity_header", "mapped_gene_or_identity", "protein_name",
        "organism", "uniprot", "confidence", "curator_note", "sequence_id", "sequence_kind", "seq_len", "label_check", "xref_status", "reviewer_verdict", "reviewer_note"]


def rows(tasks, skip=()):
    out, seen = [], set(skip)
    for t in tasks:
        d = pd.read_parquet(f"{R}/tasks/{t}.parquet", columns=["group_id", "split_time", "test_tier"])
        d = d[d.split_time == "test"].drop_duplicates("group_id")
        for g, tier in zip(d.group_id, d.test_tier):
            if g in seen:
                continue
            seen.add(g)
            a, s = au.loc[g], sq.loc[g]
            out.append(dict(group_id=g, task=t, tier=tier, patent_core=a.patent_core, patent_title=a.patent_title, assay_label_as_printed=a.Target,
                            activity_header=a.activity_header, mapped_gene_or_identity=a.mapped_identity, protein_name=a.protein_name, organism=a.organism,
                            uniprot=a.uniprot_accession, confidence=a.confidence, curator_note=a.curator_note, sequence_id=s.sequence_id, sequence_kind=s.sequence_kind,
                            seq_len=s.seq_len, label_check=s.label_check, xref_status=s.xref_status, reviewer_verdict="", reviewer_note=""))
    return pd.DataFrame(out, columns=COLS).sort_values(["task", "group_id"]).reset_index(drop=True)


main = rows(["rank_ic50", "rank_ki_kd"])
other = rows(["rank_highcontrast", "rank_complex_mutant", "rank_ec50"], skip=set(main.group_id))
# a group that is a test group of several tasks lists them all, alphabetically
member = {}
for t in ["rank_highcontrast", "rank_complex_mutant", "rank_ec50"]:
    d = pd.read_parquet(f"{R}/tasks/{t}.parquet", columns=["group_id", "split_time"])
    for g in set(d[d.split_time == "test"].group_id):
        member.setdefault(g, []).append(t)
other["task"] = other.group_id.map(lambda g: ";".join(sorted(member[g])))
for name, d in (("target_review_test.tsv", main), ("target_review_test_other_tasks.tsv", other)):
    d.astype(str).to_csv(f"{OUT}/{name}", sep="\t", index=False, encoding="utf-8-sig")
    print(name, len(d), d.task.value_counts().to_dict())
