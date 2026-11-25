"""Target-name table for PatBench r4. Additive: reads the release, writes only eval_support/.

  group_target_names.{parquet,tsv}   one row per group_id used by the five rank_* tasks : what the
        target IS (UniProt accession, official gene symbol, protein name, organism, mutation) next to what the
        patent PRINTED (assay label / column header), so `target_label` ("GR binding") is not mistaken for the
        target name. Regression-only targets are covered by structures/target_name_map.tsv.

Nothing here changes a split, a label or a task file.

    python build_target_names.py <release_dir>
"""
import re
import sys
from pathlib import Path

import pandas as pd

R = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
OUT = R / "eval_support"
OUT.mkdir(exist_ok=True)
TASKS = ["rank_ic50", "rank_ki_kd", "rank_highcontrast", "rank_complex_mutant", "rank_ec50"]
log = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    log.append(s)


grp = pd.read_parquet(R / "groups.parquet", columns=["group_id", "target_label", "endpoint"])
ta = pd.read_parquet(R / "target_audit.parquet")
sq = pd.read_parquet(R / "target_sequences.parquet",
                     columns=["group_id", "sequence_id", "sequence_kind", "confidence", "label_check", "xref_status"])
nm = pd.read_csv(R / "structures/target_name_map.tsv", sep="\t", encoding="utf-8-sig", dtype=str, keep_default_na=False)
assert nm.uniprot_accession.is_unique

in_tasks = {}
for t in TASKS:
    for g in pd.read_parquet(R / f"tasks/{t}.parquet", columns=["group_id"]).group_id.unique():
        in_tasks.setdefault(int(g), []).append(t)
ids = sorted(in_tasks)
assert ids == sorted(ta.group_id.astype(int)), "target_audit must cover exactly the groups of the five rank tasks"

d = ta.merge(grp.rename(columns={"target_label": "printed_label"}), on="group_id", how="left") \
      .merge(sq, on="group_id", how="left", suffixes=("", "_seq"))
# protein name and organism come from the release's own target table (same source as the sequences);
# the gene symbol only exists in structures/target_name_map.tsv. For 8 non-human accessions the name map spells the
# protein name out in full (EC numbers, aliases) and the organism with a different strain wording: same record,
# same length, so the two are not reconciled here.
d = d.merge(nm[["uniprot_accession", "gene_symbol", "taxid", "reviewed", "uniprot_length"]],
            on="uniprot_accession", how="left")
assert len(d) == len(ids) and d.group_id.is_unique
assert d.gene_symbol.notna().all() and d.taxid.notna().all(), "accession missing from target_name_map"
assert (d.protein_name.fillna("") != "").all()
full = d.sequence_kind == "full_length"
assert (d.loc[full, "seq_len"].astype(int) == d.loc[full, "uniprot_length"].astype(int)).all(), \
    "full-length sequence length differs from the UniProt length in target_name_map"
say("groups whose protein-name text differs from target_name_map (wording only, same length):",
    int((d.protein_name != d.merge(nm[["uniprot_accession", "protein_name"]], on="uniprot_accession",
                                   suffixes=("", "_nm")).protein_name_nm).sum()))


def short(n):
    return re.split(r" \(| \[", n, maxsplit=1)[0].strip()


def tokens(s):
    return [re.sub(r"[^A-Z0-9]", "", t.upper()) for t in re.split(r"[ ;,/|]+", s) if t]


# curated names that are free text ("ORG:organism|name") carry no gene symbol, so they cannot be compared by
# symbol. These were read against the UniProt record name by hand (2026-10-31): the three accessions below are
# named for one protein but the accession is another one.
CONFLICT = set()   # r4 had G4WHZ0 (NDM-1), Q8GRG4 (VIM-1), P9WFV3 (Mtb LeuRS) here; corrected in r5 (C7C422, Q9XAY4, P9WFV1)
check, note = [], []
for r in d.itertuples(index=False):
    cur = r.mapped_identity
    if cur.startswith("ORG:"):
        if r.uniprot_accession in CONFLICT:
            check.append("free_text_conflicts_uniprot")
            alt = f"; ChEMBL 编号 {r.xref_accession}" if r.xref_status == "differs" else ""
            note.append(f"标注为「{cur.split('|', 1)[1]}」,该编号在 UniProt 里是「{short(r.protein_name)}」{alt}")
        else:
            check.append("free_text_consistent")
            note.append("")
    else:
        ok = any(t in tokens(r.gene_symbol) for t in tokens(cur))
        check.append("symbol_matches_uniprot" if ok else "symbol_not_in_uniprot_gene_names")
        note.append("")
d["name_check"], d["name_note"] = check, note

out = pd.DataFrame({
    "group_id": d.group_id.astype(int),
    "patent_core": d.patent_core,
    "printed_label": d.printed_label,
    "printed_label_raw": d.Target,
    "activity_header": d.activity_header,
    "endpoint": d.endpoint,
    "uniprot_accession": d.uniprot_accession,
    "uniprot_entry": d.uniprot_entry,
    "gene_symbol": d.gene_symbol,
    "curated_symbol": d.mapped_identity,
    "protein_name": d.protein_name,
    "protein_name_short": d.protein_name.map(short),
    "organism": d.organism,
    "taxid": d.taxid,
    "uniprot_reviewed": d.reviewed,
    "target_class": d.target_class,
    "mutation": d.mutation,
    "sequence_kind": d.sequence_kind,
    "sequence_id": d.sequence_id,
    "confidence": d.confidence,
    "xref_status": d.xref_status,
    "name_check": d.name_check,
    "name_note": d.name_note,
    "in_tasks": d.group_id.astype(int).map(lambda g: ";".join(in_tasks[g])),
}).sort_values("group_id").reset_index(drop=True)
out.to_parquet(OUT / "group_target_names.parquet", index=False)
out.to_csv(OUT / "group_target_names.tsv", sep="\t", index=False, encoding="utf-8-sig")

say("group_target_names rows:", len(out), " distinct accessions:", out.uniprot_accession.nunique())
say("name_check:", out.name_check.value_counts().to_dict())
say("empty gene_symbol:", int((out.gene_symbol == "").sum()), " empty protein_name:", int((out.protein_name == "").sum()))
bad = out[out.name_check == "free_text_conflicts_uniprot"]
say("conflict groups:", bad.group_id.tolist(), " accessions:", sorted(bad.uniprot_accession.unique()))
(OUT / "build_target_names.log").write_text("\n".join(log) + "\n", encoding="utf-8")
