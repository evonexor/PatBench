"""Join each ranking-task group to the target-map pair it was mapped from (input of audit_targets.py).

    python join_target_map.py <release_dir> <target_map.parquet> <out.parquet>
"""
import sys

import pandas as pd

R, MAP, OUT = sys.argv[1].rstrip("/"), sys.argv[2], sys.argv[3]
tasks = ["rank_ic50", "rank_ki_kd", "rank_ec50", "rank_highcontrast", "rank_complex_mutant"]
t = pd.read_parquet(MAP)
ids = sorted(set().union(*[set(pd.read_parquet(f"{R}/tasks/{n}.parquet", columns=["group_id"]).group_id)
                           for n in tasks]))
s = pd.read_parquet(R + "/sar_all.parquet", columns=["group_id", "patent_core", "Target"])
s = s[s.group_id.isin(ids)]
# a group can join several raw labels; the pair with the most rows is the one it was mapped from
pair = (s.groupby(["group_id", "patent_core", "Target"]).size().reset_index(name="n")
        .sort_values(["group_id", "n"], ascending=[True, False]).drop_duplicates("group_id"))
cols = ["patent_core", "Target", "uniprot_accession", "mapped_identity", "protein_name", "uniprot_entry",
        "organism", "organism_stated", "organism_detected", "organism_resolution", "mutation", "target_class",
        "mapping_kind", "confidence", "curator_confidence", "mappability", "sequence_caveat", "target_symbols",
        "xref_db", "xref_accession", "xref_symbol", "xref_value_match", "patent_title", "activity_header",
        "curator_note", "seq_len"]
a = pair.merge(t[cols].drop_duplicates(["patent_core", "Target"]), on=["patent_core", "Target"], how="left")
a.to_parquet(OUT, index=False)
print(len(a), "groups,", int(a.uniprot_accession.notna().sum()), "with a mapping")
