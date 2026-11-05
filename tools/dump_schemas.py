"""Dump row counts and column names / types of the main tables of a release as JSON.

    python tools/dump_schemas.py <release_dir> > schemas.json

docs/DATA_SCHEMA.md is generated from this by tools/make_schema_doc.py."""
import json
import os
import sys

import pyarrow.parquet as pq

R = sys.argv[1].rstrip("/") + "/"
out = {}
files = ["sar_all.parquet", "groups.parquet", "targets.parquet", "tasks/rank_ic50.parquet", "tasks/rank_ki_kd.parquet", "tasks/rank_ec50.parquet",
         "tasks/rank_highcontrast.parquet", "tasks/rank_complex_mutant.parquet", "tasks/regression_affinity.parquet", "target_sequences.parquet",
         "target_audit.parquet", "structures/target_structures.parquet", "structures/group_structure_map.parquet", "eval_support/test_group_profile.parquet",
         "eval_support/target_domain_sequences.parquet", "eval_all/rank_ic50.parquet", "baselines/rank_ic50_heavy_atoms.parquet",
         "protocols/rank_ic50_time_k5.parquet", "protocols/rank_ic50_cliffs.parquet", "structures/structure_variant_metrics.parquet"]
for f in files:
    p = R + f
    if not os.path.exists(p): continue
    pf = pq.ParquetFile(p)
    out[f] = {"rows": pf.metadata.num_rows, "cols": [(n, str(t)) for n, t in zip(pf.schema_arrow.names, pf.schema_arrow.types)]}
print(json.dumps(out))
