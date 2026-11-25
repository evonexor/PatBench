# Standardised public activity tables

`build_standardized.py` turns four public sources into one shared schema, `merge_sar.py` merges three of them into one
SAR table, and `filter_invalid_smiles.py` writes `*.valid.parquet` copies that keep only the rows whose SMILES RDKit can parse. PatBench uses the output only to
ask how much of the extracted patent data is already public (`08_build/novelty.py`, `08_build/public_overlap.py`).

The standard columns always come first:

`Patent`, `Target`, `Activity Type`, `Relation`, `Value`, `Unit`, `SMILES`, `UniProt`

The remaining source columns follow. There is no merging across databases and no de-duplication in the standardised files.

## Output (`$PATBENCH_WORKSPACE/raw/patent-data/standardized/`)

| File | Content |
|---|---|
| `BindingDB.parquet` | the full BindingDB patent TSV; Ki / IC50 / Kd / EC50 / kon / koff split into one row per activity |
| `ChEMBL.parquet` | ChEMBL 37 patent activity sources SRC_ID 37 and 38; empty kon / koff placeholder rows without a value removed |
| `HARVEST.parquet` | the public H-bench CSV files (48 targets); Ki / IC50 / Kd / EC50 split into one row per activity |
| `SureChEMBL.parquet` | patent-structure pairs from claims for molecules of weight 200 or more; the snapshot has no quantitative activity and the full compound map (1.54 billion rows) is not expanded |
| `demo/` | ten rows per source as CSV and JSON |
| `stats.json` | row counts |

PatBench's public comparison set is BindingDB and ChEMBL only. HARVEST is left out because its public snapshot covers just
the H-bench subset (48 targets).

## Field conventions

- **Relation**: for BindingDB the relation is split from prefixes such as `<100` and `>10000`; HARVEST keeps its own `relation` (including `range`); ChEMBL uses `standard_relation`.
- **Value / Unit**: value and unit are separate columns. For BindingDB and HARVEST the Ki, IC50, Kd and EC50 units are `nM`. HARVEST's original range text is kept in `original_range`.
- **UniProt**: BindingDB uses the SwissProt accession of chain 1 and falls back to TrEMBL; ChEMBL uses the accession of the target component and joins complex components with `|`; HARVEST uses `uniprot_acc`.
- The original activity columns (such as `Ki (nM)`, `IC50 (nM)`, `Kd (nM)`, `EC50 (nM)`, and ChEMBL's `original_type`, `original_relation`, `original_value`, `original_units`) and the original ligand structure columns (`Ligand SMILES`, `SMILES`, `InChI`, `InChIKey`) are dropped after conversion. BindingDB's empty target chain 2-50 columns are not kept; the chain 1 sequence and UniProt accession are.
- SureChEMBL has no target and no activity, so `Target`, `Activity Type`, `Relation`, `Value`, `Unit` and `UniProt` are empty.

## Running

The scripts read the raw database files from `$PATBENCH_WORKSPACE/raw/patent-data/` (the layout they were written for) and
write below it:

```bash
export PATBENCH_WORKSPACE=/path/to/workspace
patbench exec 01_public_data/build_standardized.py
patbench exec 01_public_data/merge_sar.py
patbench exec 01_public_data/filter_invalid_smiles.py
```

`build_standardized.py --source {all,bindingdb,chembl,harvest,surechembl}` processes one source or all of them (default `all`).
The raw database files themselves are not distributed with this repository.
