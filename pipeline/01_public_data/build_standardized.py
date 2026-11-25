#!/usr/bin/env python3
"""Standardize BindingDB, ChEMBL, HARVEST, and SureChEMBL into a shared SAR schema.

Standard columns (always first):
  Patent, Target, Activity Type, Relation, Value, Unit, SMILES, UniProt
Remaining source columns are kept after those eight.
"""

from __future__ import annotations
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")

import argparse
import csv
import io
import json
import re
import sqlite3
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(_ws + "/raw/patent-data")
OUT_DIR = ROOT / "standardized"
DEMO_DIR = OUT_DIR / "demo"

STANDARD_COLS = [
    "Patent",
    "Target",
    "Activity Type",
    "Relation",
    "Value",
    "Unit",
    "SMILES",
    "UniProt",
]

BINDINGDB_ACTIVITY = [
    ("Ki (nM)", "Ki", "nM"),
    ("IC50 (nM)", "IC50", "nM"),
    ("Kd (nM)", "Kd", "nM"),
    ("EC50 (nM)", "EC50", "nM"),
    ("kon (M-1-s-1)", "kon", "M-1-s-1"),
    ("koff (s-1)", "koff", "s-1"),
]

HARVEST_ACTIVITY = [
    ("Ki (nM)", "Ki", "nM"),
    ("IC50 (nM)", "IC50", "nM"),
    ("Kd (nM)", "Kd", "nM"),
    ("EC50 (nM)", "EC50", "nM"),
]

HARVEST_DROP_COLS = {
    "Ki (nM)",
    "IC50 (nM)",
    "Kd (nM)",
    "EC50 (nM)",
    "SMILES",
    "InChIKey",
}

CHAIN_GE2_RE = re.compile(
    r"(?:Target Chain Sequence|of Target Chain) (?:[2-9]|[1-4][0-9]|50)\b"
)
LEADING_REL_RE = re.compile(r"^(<=|>=|<|>|~|=)\s*(.*)$")


def log(msg: str) -> None:
    print(msg, flush=True)


def as_text(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if value != value:  # NaN
            return ""
        text = f"{value:.12g}"
        return text
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace").strip()
    return str(value).strip()


def parse_bindingdb_cell(raw: str) -> tuple[str, str] | None:
    text = as_text(raw)
    if not text:
        return None
    match = LEADING_REL_RE.match(text)
    if match:
        relation, rest = match.group(1), match.group(2).strip()
        if not rest:
            return None
        return relation, rest
    return "=", text


BINDINGDB_DROP_COLS = {
    "Ki (nM)",
    "IC50 (nM)",
    "Kd (nM)",
    "EC50 (nM)",
    "kon (M-1-s-1)",
    "koff (s-1)",
    "Ligand SMILES",
    "Ligand InChI",
    "Ligand InChI Key",
}


def keep_bindingdb_column(name: str) -> bool:
    if name in BINDINGDB_DROP_COLS or name in STANDARD_COLS:
        return False
    return not bool(CHAIN_GE2_RE.search(name))


def make_schema(columns: list[str]) -> pa.Schema:
    return pa.schema([pa.field(col, pa.string()) for col in columns])


class ParquetBatchWriter:
    def __init__(self, path: Path, columns: list[str], batch_size: int = 20000):
        self.path = path
        self.columns = columns
        self.batch_size = batch_size
        self.schema = make_schema(columns)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.writer = pq.ParquetWriter(
            str(path),
            self.schema,
            compression="zstd",
            use_dictionary=True,
        )
        self.buffer: list[dict[str, str]] = []
        self.rows = 0

    def add(self, row: dict[str, str]) -> None:
        self.buffer.append(row)
        if len(self.buffer) >= self.batch_size:
            self.flush()

    def flush(self) -> None:
        if not self.buffer:
            return
        arrays = {
            col: pa.array([row.get(col, "") for row in self.buffer], type=pa.string())
            for col in self.columns
        }
        table = pa.Table.from_pydict(arrays, schema=self.schema)
        self.writer.write_table(table)
        self.rows += len(self.buffer)
        self.buffer.clear()

    def close(self) -> int:
        self.flush()
        self.writer.close()
        return self.rows


def write_demo(rows: list[dict[str, str]], dataset: str, source: str, selection: str) -> None:
    if not rows:
        raise RuntimeError(f"No demo rows for {dataset}")
    columns = list(rows[0].keys())
    csv_path = DEMO_DIR / f"{dataset}.csv"
    json_path = DEMO_DIR / f"{dataset}.json"
    DEMO_DIR.mkdir(parents=True, exist_ok=True)

    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in columns})

    records = []
    for row in rows:
        record = {}
        for col in columns:
            value = row.get(col, "")
            record[col] = None if value == "" else value
        records.append(record)

    payload = {
        "dataset": dataset,
        "source": source,
        "selection": selection,
        "row_count": len(records),
        "standard_columns": STANDARD_COLS,
        "records": records,
    }
    json_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    log(f"  demo -> {csv_path.name}, {json_path.name} ({len(records)} rows)")


def iter_parquet_rows(path: Path, limit: int | None = None):
    parquet_file = pq.ParquetFile(path)
    yielded = 0
    for batch in parquet_file.iter_batches(batch_size=2048):
        columns = batch.column_names
        arrays = [batch.column(i).to_pylist() for i in range(batch.num_columns)]
        for i in range(batch.num_rows):
            yield {columns[j]: as_text(arrays[j][i]) for j in range(len(columns))}
            yielded += 1
            if limit is not None and yielded >= limit:
                return


def process_bindingdb() -> dict:
    zip_path = ROOT / "bindingdb/202609/BindingDB_Patents_202609_tsv.zip"
    out_path = OUT_DIR / "BindingDB.parquet"
    log(f"[BindingDB] reading {zip_path}")
    demo_rows: list[dict[str, str]] = []
    source_rows = 0
    skipped_no_activity = 0
    writer = None
    extra_cols: list[str] = []

    with zipfile.ZipFile(zip_path) as zf:
        with zf.open(zf.namelist()[0]) as raw:
            text = io.TextIOWrapper(raw, encoding="utf-8", errors="replace", newline="")
            reader = csv.DictReader(text, delimiter="\t")
            if reader.fieldnames is None:
                raise RuntimeError("BindingDB TSV has no header")
            extra_cols = [name for name in reader.fieldnames if keep_bindingdb_column(name)]
            columns = STANDARD_COLS + extra_cols
            writer = ParquetBatchWriter(out_path, columns)
            for src in reader:
                source_rows += 1
                patent = as_text(src.get("Patent Number"))
                target = as_text(src.get("Target Name"))
                smiles = as_text(src.get("Ligand SMILES"))
                uniprot = as_text(src.get("UniProt (SwissProt) Primary ID of Target Chain 1"))
                if not uniprot:
                    uniprot = as_text(src.get("UniProt (TrEMBL) Primary ID of Target Chain 1"))
                extras = {name: as_text(src.get(name)) for name in extra_cols}
                emitted = 0
                for src_col, act_type, unit in BINDINGDB_ACTIVITY:
                    parsed = parse_bindingdb_cell(src.get(src_col, ""))
                    if parsed is None:
                        continue
                    relation, value = parsed
                    row = {
                        "Patent": patent,
                        "Target": target,
                        "Activity Type": act_type,
                        "Relation": relation,
                        "Value": value,
                        "Unit": unit,
                        "SMILES": smiles,
                        "UniProt": uniprot,
                    }
                    row.update(extras)
                    writer.add(row)
                    emitted += 1
                    if len(demo_rows) < 10 and patent and smiles and value:
                        demo_rows.append(row)
                if emitted == 0:
                    skipped_no_activity += 1
                if source_rows % 200000 == 0:
                    log(f"  scanned {source_rows:,} source rows, written {writer.rows + len(writer.buffer):,}")

    assert writer is not None
    n = writer.close()
    write_demo(
        demo_rows,
        "BindingDB",
        "BindingDB 202609 / BindingDB_Patents.tsv",
        "First 10 standardized activity rows with Patent, SMILES and Value. Ki/IC50/Kd/EC50/kon/koff unpivoted into Activity Type/Relation/Value/Unit; Ligand SMILES/InChI/InChI Key and original activity columns dropped.",
    )
    stats = {
        "dataset": "BindingDB",
        "source_rows": source_rows,
        "output_rows": n,
        "skipped_no_activity": skipped_no_activity,
        "output": str(out_path),
        "extra_columns": len(extra_cols),
    }
    log(f"[BindingDB] done source={source_rows:,} output={n:,}")
    return stats


def process_harvest() -> dict:
    csv_dir = ROOT / "harvest/repository/data/h_bench"
    out_path = OUT_DIR / "HARVEST.parquet"
    files = sorted(csv_dir.glob("*.csv"))
    log(f"[HARVEST] {len(files)} target CSVs")
    extra_cols: list[str] | None = None
    writer = None
    demo_rows: list[dict[str, str]] = []
    source_rows = 0
    skipped = 0

    for csv_path in files:
        with csv_path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                continue
            if extra_cols is None:
                extra_cols = [
                    name
                    for name in reader.fieldnames
                    if name not in STANDARD_COLS and name not in HARVEST_DROP_COLS
                ] + ["source_csv"]
                writer = ParquetBatchWriter(out_path, STANDARD_COLS + extra_cols)
            for src in reader:
                source_rows += 1
                extras = {name: as_text(src.get(name)) for name in extra_cols if name != "source_csv"}
                extras["source_csv"] = csv_path.name
                patent = as_text(src.get("patent_number"))
                target = as_text(src.get("protein_target_name"))
                smiles = as_text(src.get("SMILES"))
                uniprot = as_text(src.get("uniprot_acc")) or as_text(src.get("UniProt ID"))
                relation = as_text(src.get("relation")) or "="
                emitted = 0
                for src_col, act_type, unit in HARVEST_ACTIVITY:
                    value = as_text(src.get(src_col))
                    if not value:
                        continue
                    row = {
                        "Patent": patent,
                        "Target": target,
                        "Activity Type": act_type,
                        "Relation": relation,
                        "Value": value,
                        "Unit": unit,
                        "SMILES": smiles,
                        "UniProt": uniprot,
                    }
                    row.update(extras)
                    writer.add(row)
                    emitted += 1
                    if len(demo_rows) < 10 and patent and smiles and value:
                        demo_rows.append(row)
                if emitted == 0:
                    skipped += 1

    assert writer is not None and extra_cols is not None
    n = writer.close()
    write_demo(
        demo_rows,
        "HARVEST",
        "HARVEST commit 61fd705edbdfcd212f0c9ac23786412b8fdcede2 / public H-bench",
        "First 10 standardized activity rows. Ki/IC50/Kd/EC50 unpivoted into Activity Type/Relation/Value/Unit; original activity columns, SMILES and InChIKey dropped. Unit is nM.",
    )
    stats = {
        "dataset": "HARVEST",
        "source_rows": source_rows,
        "output_rows": n,
        "skipped_no_activity": skipped,
        "target_files": len(files),
        "output": str(out_path),
    }
    log(f"[HARVEST] done source={source_rows:,} output={n:,}")
    return stats


def load_chembl_lookups(cur: sqlite3.Cursor) -> dict:
    log("[ChEMBL] loading lookups")
    src_name = {
        int(src_id): as_text(desc)
        for src_id, desc in cur.execute(
            "SELECT src_id, src_description FROM source WHERE src_id IN (37, 38)"
        )
    }

    docs = {}
    for row in cur.execute(
        """
        SELECT doc_id, patent_id, chembl_id, title, year, authors, src_id
        FROM docs
        WHERE src_id IN (37, 38)
        """
    ):
        docs[int(row[0])] = {
            "patent_id": as_text(row[1]),
            "document_chembl_id": as_text(row[2]),
            "document_title": as_text(row[3]),
            "year": as_text(row[4]),
            "authors": as_text(row[5]),
            "doc_src_id": as_text(row[6]),
        }

    assays = {}
    for row in cur.execute(
        """
        SELECT assay_id, tid, description, chembl_id, assay_type, src_id, confidence_score
        FROM assays
        WHERE src_id IN (37, 38)
        """
    ):
        assays[int(row[0])] = {
            "tid": int(row[1]) if row[1] is not None else None,
            "assay_description": as_text(row[2]),
            "assay_chembl_id": as_text(row[3]),
            "assay_type": as_text(row[4]),
            "assay_src_id": as_text(row[5]),
            "confidence_score": as_text(row[6]),
        }

    targets = {}
    for row in cur.execute(
        """
        SELECT tid, pref_name, target_type, organism, chembl_id
        FROM target_dictionary
        """
    ):
        targets[int(row[0])] = {
            "pref_name": as_text(row[1]),
            "target_type": as_text(row[2]),
            "organism": as_text(row[3]),
            "target_chembl_id": as_text(row[4]),
        }

    uniprot_lists: dict[int, list[str]] = defaultdict(list)
    for tid, accession in cur.execute(
        """
        SELECT tc.tid, seq.accession
        FROM target_components tc
        JOIN component_sequences seq ON seq.component_id = tc.component_id
        WHERE seq.accession IS NOT NULL AND trim(seq.accession) != ''
        """
    ):
        acc = as_text(accession)
        if acc and acc not in uniprot_lists[int(tid)]:
            uniprot_lists[int(tid)].append(acc)
    uniprot = {tid: "|".join(vals) for tid, vals in uniprot_lists.items()}

    molecules = {}
    for molregno, chembl_id in cur.execute(
        "SELECT molregno, chembl_id FROM molecule_dictionary"
    ):
        molecules[int(molregno)] = as_text(chembl_id)

    structures = {}
    for molregno, smiles, inchikey in cur.execute(
        "SELECT molregno, canonical_smiles, standard_inchi_key FROM compound_structures"
    ):
        structures[int(molregno)] = {
            "smiles": as_text(smiles),
            "inchikey": as_text(inchikey),
        }

    log(
        f"[ChEMBL] lookups docs={len(docs):,} assays={len(assays):,} "
        f"targets={len(targets):,} molecules={len(molecules):,} structures={len(structures):,}"
    )
    return {
        "src_name": src_name,
        "docs": docs,
        "assays": assays,
        "targets": targets,
        "uniprot": uniprot,
        "molecules": molecules,
        "structures": structures,
    }


def process_chembl() -> dict:
    db_path = ROOT / "chembl/chembl_37/extracted/chembl_37.db"
    out_path = OUT_DIR / "ChEMBL.parquet"
    extra_cols = [
        "source_dataset",
        "activity_id",
        "molecule_chembl_id",
        "assay_chembl_id",
        "target_chembl_id",
        "target_type",
        "organism",
        "pchembl_value",
        "data_validity_comment",
        "potential_duplicate",
        "document_chembl_id",
        "document_title",
        "year",
        "authors",
        "assay_description",
        "assay_type",
        "confidence_score",
        "src_id",
        "src_description",
        "doc_id",
        "assay_id",
        "molregno",
        "activity_comment",
        "standard_upper_value",
        "standard_text_value",
    ]
    columns = STANDARD_COLS + extra_cols
    log(f"[ChEMBL] opening {db_path}")
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    cur = con.cursor()
    lookups = load_chembl_lookups(cur)
    writer = ParquetBatchWriter(out_path, columns, batch_size=25000)
    demo_rows: list[dict[str, str]] = []
    source_rows = 0
    skipped_no_value = 0
    skipped_no_patent = 0

    cur.execute(
        """
        SELECT
            activity_id, assay_id, doc_id, molregno,
            standard_relation, standard_value, standard_units, standard_type,
            relation, value, units, type,
            pchembl_value, data_validity_comment, potential_duplicate,
            activity_comment, src_id, standard_text_value, text_value,
            standard_upper_value
        FROM activities
        WHERE src_id IN (37, 38)
        ORDER BY activity_id
        """
    )
    while True:
        batch = cur.fetchmany(25000)
        if not batch:
            break
        for rec in batch:
            source_rows += 1
            doc = lookups["docs"].get(int(rec["doc_id"])) if rec["doc_id"] is not None else None
            patent = doc["patent_id"] if doc else ""
            if not patent:
                skipped_no_patent += 1
                continue
            act_type = as_text(rec["standard_type"]) or as_text(rec["type"])
            relation = as_text(rec["standard_relation"]) or as_text(rec["relation"])
            value = as_text(rec["standard_value"]) or as_text(rec["value"]) or as_text(rec["standard_text_value"]) or as_text(rec["text_value"])
            unit = as_text(rec["standard_units"]) or as_text(rec["units"])
            if not value:
                skipped_no_value += 1
                continue
            if not relation:
                relation = "="
            assay = lookups["assays"].get(int(rec["assay_id"])) if rec["assay_id"] is not None else None
            tid = assay["tid"] if assay else None
            target_info = lookups["targets"].get(tid) if tid is not None else None
            molregno = int(rec["molregno"]) if rec["molregno"] is not None else None
            struct = lookups["structures"].get(molregno) if molregno is not None else None
            src_id = rec["src_id"]
            row = {
                "Patent": patent,
                "Target": target_info["pref_name"] if target_info else "",
                "Activity Type": act_type,
                "Relation": relation,
                "Value": value,
                "Unit": unit,
                "SMILES": struct["smiles"] if struct else "",
                "UniProt": lookups["uniprot"].get(tid, "") if tid is not None else "",
                "source_dataset": "ChEMBL",
                "activity_id": as_text(rec["activity_id"]),
                "molecule_chembl_id": lookups["molecules"].get(molregno, "") if molregno is not None else "",
                "assay_chembl_id": assay["assay_chembl_id"] if assay else "",
                "target_chembl_id": target_info["target_chembl_id"] if target_info else "",
                "target_type": target_info["target_type"] if target_info else "",
                "organism": target_info["organism"] if target_info else "",
                "pchembl_value": as_text(rec["pchembl_value"]),
                "data_validity_comment": as_text(rec["data_validity_comment"]),
                "potential_duplicate": as_text(rec["potential_duplicate"]),
                "document_chembl_id": doc["document_chembl_id"] if doc else "",
                "document_title": doc["document_title"] if doc else "",
                "year": doc["year"] if doc else "",
                "authors": doc["authors"] if doc else "",
                "assay_description": assay["assay_description"] if assay else "",
                "assay_type": assay["assay_type"] if assay else "",
                "confidence_score": assay["confidence_score"] if assay else "",
                "src_id": as_text(src_id),
                "src_description": lookups["src_name"].get(int(src_id), "") if src_id is not None else "",
                "doc_id": as_text(rec["doc_id"]),
                "assay_id": as_text(rec["assay_id"]),
                "molregno": as_text(rec["molregno"]),
                "activity_comment": as_text(rec["activity_comment"]),
                "standard_upper_value": as_text(rec["standard_upper_value"]),
                "standard_text_value": as_text(rec["standard_text_value"]),
            }
            writer.add(row)
            if len(demo_rows) < 10 and row["SMILES"] and row["Value"] and row["Activity Type"]:
                demo_rows.append(row)
        if source_rows % 250000 == 0:
            log(f"  scanned {source_rows:,} activities, written {writer.rows + len(writer.buffer):,}")

    n = writer.close()
    con.close()
    write_demo(
        demo_rows,
        "ChEMBL",
        "ChEMBL 37 / chembl_37.db / SRC_ID 37+38 patent bioactivity",
        "First 10 standardized patent activities with Patent, SMILES, Activity Type and Value. Original type/relation/value/units and InChIKey dropped; empty kon/k_off placeholders without a value were omitted.",
    )
    stats = {
        "dataset": "ChEMBL",
        "source_rows": source_rows,
        "output_rows": n,
        "skipped_no_value": skipped_no_value,
        "skipped_no_patent": skipped_no_patent,
        "output": str(out_path),
    }
    log(f"[ChEMBL] done source={source_rows:,} output={n:,} skipped_no_value={skipped_no_value:,}")
    return stats


def _arrow_to_text(value) -> str:
    if isinstance(value, list):
        return " | ".join(as_text(x) for x in value if x not in (None, ""))
    return as_text(value)


def _table_by_id(table, id_col: str, keep_cols: list[str]) -> dict[int, dict[str, str]]:
    ids = table.column(id_col).to_pylist()
    cols = {name: table.column(name).to_pylist() for name in keep_cols}
    out: dict[int, dict[str, str]] = {}
    for i, raw_id in enumerate(ids):
        out[int(raw_id)] = {name: _arrow_to_text(cols[name][i]) for name in keep_cols}
    return out


def process_surechembl() -> dict:
    """SureChEMBL has no paired activity. Write a standardized extract plus demo.

    The full patent_compound_map is 1.54 billion rows and is not fully joined.
    This extract takes unique claims (field_id=2) patent-compound pairs with
    mol_weight >= 200, in mapping order, capped at 100,000 rows.
    """
    base = ROOT / "surechembl/2026-09-15"
    out_path = OUT_DIR / "SureChEMBL.parquet"
    EXTRACT_LIMIT = 100_000
    CANDIDATE_LIMIT = 140_000
    log("[SureChEMBL] no quantitative activity; extracting unique claims pairs")

    fields_tbl = pq.read_table(base / "fields.parquet")
    field_name = {
        int(i): as_text(n)
        for i, n in zip(fields_tbl.column("id").to_pylist(), fields_tbl.column("field_name").to_pylist())
    }

    extra_cols = [
        "source_dataset",
        "patent_title",
        "compound_id",
        "mol_weight",
        "field_name",
        "patent_id",
        "field_id",
        "country",
        "publication_date",
        "family_id",
        "assignee",
    ]
    columns = STANDARD_COLS + extra_cols

    map_file = pq.ParquetFile(base / "patent_compound_map.parquet")
    pairs: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    source_rows = 0
    log("[SureChEMBL] scanning patent_compound_map for unique field_id=2 pairs")
    for rg in range(map_file.num_row_groups):
        table = map_file.read_row_group(rg, columns=["patent_id", "compound_id", "field_id"])
        pids = table.column("patent_id").to_pylist()
        cids = table.column("compound_id").to_pylist()
        fids = table.column("field_id").to_pylist()
        for patent_id, compound_id, field_id in zip(pids, cids, fids):
            source_rows += 1
            if field_id != 2:
                continue
            pair = (int(patent_id), int(compound_id))
            if pair in seen:
                continue
            seen.add(pair)
            pairs.append(pair)
            if len(pairs) >= CANDIDATE_LIMIT:
                break
        log(f"  map row_group {rg}/{map_file.num_row_groups} scanned={source_rows:,} unique_claims={len(pairs):,}")
        if len(pairs) >= CANDIDATE_LIMIT:
            break

    patent_ids = sorted({p for p, _ in pairs})
    compound_ids = sorted({c for _, c in pairs})
    log(f"[SureChEMBL] lookup patents={len(patent_ids):,} compounds={len(compound_ids):,}")

    patents = pq.read_table(
        base / "patents.parquet",
        columns=["id", "patent_number", "country", "publication_date", "family_id", "assignee", "title"],
        filters=[("id", "in", patent_ids)],
    )
    compounds = pq.read_table(
        base / "compounds.parquet",
        columns=["id", "smiles", "inchi", "inchi_key", "mol_weight"],
        filters=[("id", "in", compound_ids)],
    )
    patent_by_id = _table_by_id(
        patents, "id", ["patent_number", "country", "publication_date", "family_id", "assignee", "title"]
    )
    compound_by_id = _table_by_id(compounds, "id", ["smiles", "inchi", "inchi_key", "mol_weight"])
    del patents, compounds

    writer = ParquetBatchWriter(out_path, columns, batch_size=20000)
    demo_rows: list[dict[str, str]] = []
    skipped_mw = 0
    skipped_missing = 0
    for patent_id, compound_id in pairs:
        compound = compound_by_id.get(compound_id)
        patent = patent_by_id.get(patent_id)
        if not compound or not patent:
            skipped_missing += 1
            continue
        try:
            mw = float(compound["mol_weight"]) if compound["mol_weight"] else 0.0
        except ValueError:
            mw = 0.0
        if mw < 200:
            skipped_mw += 1
            continue
        row = {
            "Patent": patent["patent_number"],
            "Target": "",
            "Activity Type": "",
            "Relation": "",
            "Value": "",
            "Unit": "",
            "SMILES": compound["smiles"],
            "UniProt": "",
            "source_dataset": "SureChEMBL",
            "patent_title": patent["title"],
            "compound_id": as_text(compound_id),
            "mol_weight": compound["mol_weight"],
            "field_name": field_name.get(2, "clms"),
            "patent_id": as_text(patent_id),
            "field_id": "2",
            "country": patent["country"],
            "publication_date": patent["publication_date"],
            "family_id": patent["family_id"],
            "assignee": patent["assignee"],
        }
        writer.add(row)
        if len(demo_rows) < 10 and row["Patent"] and row["SMILES"]:
            demo_rows.append(row)
        if writer.rows + len(writer.buffer) >= EXTRACT_LIMIT:
            break

    n = writer.close()
    write_demo(
        demo_rows,
        "SureChEMBL",
        "SureChEMBL full snapshot 2026-09-15",
        "First 10 unique patent-compound pairs from claims (field_id=2) with mol_weight>=200. InChI/InChIKey dropped after mapping SMILES. Target/Activity Type/Relation/Value/Unit/UniProt are empty because this snapshot has no paired quantitative activity.",
    )
    stats = {
        "dataset": "SureChEMBL",
        "note": "No quantitative activity in this snapshot. Full patent_compound_map has 1,542,719,188 rows and was not fully materialized. Parquet extract is unique claims pairs with MW>=200, capped at 100000.",
        "map_rows_scanned": source_rows,
        "candidate_pairs": len(pairs),
        "output_rows": n,
        "skipped_mw_lt_200": skipped_mw,
        "skipped_missing_join": skipped_missing,
        "extract_limit": EXTRACT_LIMIT,
        "output": str(out_path),
    }
    log(f"[SureChEMBL] done extract={n:,} scanned_map_rows={source_rows:,}")
    return stats


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source",
        choices=["all", "bindingdb", "chembl", "harvest", "surechembl"],
        default="all",
    )
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DEMO_DIR.mkdir(parents=True, exist_ok=True)

    jobs = {
        "harvest": process_harvest,
        "bindingdb": process_bindingdb,
        "chembl": process_chembl,
        "surechembl": process_surechembl,
    }
    selected = list(jobs) if args.source == "all" else [args.source]
    report_path = OUT_DIR / "stats.json"
    existing: list[dict] = []
    if report_path.exists() and args.source != "all":
        try:
            loaded = json.loads(report_path.read_text(encoding="utf-8"))
            if isinstance(loaded, list):
                existing = loaded
        except json.JSONDecodeError:
            existing = []
    by_name = {item.get("dataset", "").lower(): item for item in existing if isinstance(item, dict)}
    for name in selected:
        item = jobs[name]()
        by_name[item.get("dataset", name).lower()] = item
    order = ["harvest", "bindingdb", "chembl", "surechembl"]
    stats = [by_name[k] for k in order if k in by_name]
    stats.extend(item for key, item in by_name.items() if key not in order)
    report_path.write_text(json.dumps(stats, ensure_ascii=False, indent=2), encoding="utf-8")
    log(f"stats -> {report_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
