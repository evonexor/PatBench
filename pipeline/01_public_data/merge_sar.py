#!/usr/bin/env python3
"""Merge BindingDB + ChEMBL src38 + HARVEST into one SAR table."""

from __future__ import annotations
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")

import csv
import json
import re
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = Path(_ws + "/raw/patent-data")
STD = ROOT / "standardized"
OUT = STD / "merged"
DEMO = STD / "demo"

KIND = re.compile(r"(?:A\d?|B\d?|C\d?|U\d?|S\d?|E\d?)$")
REL_MAP = {
    "=": "=",
    "==": "=",
    "<": "<",
    ">": ">",
    "<=": "<=",
    ">=": ">=",
    "~": "~",
    ">>": ">",
    "<<": "<",
    "range": "range",
}
ACT_MAP = {
    "IC50": "IC50",
    "KI": "Ki",
    "KD": "Kd",
    "EC50": "EC50",
    "AC50": "AC50",
    "KON": "kon",
    "K_OFF": "koff",
    "KOFF": "koff",
    "PIC50": "pIC50",
    "PKI": "pKi",
}
UNIT_NANO = {
    "NM": ("nM", 1.0),
    "NANOMOLAR": ("nM", 1.0),
    "UM": ("uM", 1000.0),
    "µM": ("uM", 1000.0),
    "ΜM": ("uM", 1000.0),
    "MICROMOLAR": ("uM", 1000.0),
    "MM": ("mM", 1_000_000.0),
    "MMOLAR": ("mM", 1_000_000.0),
    "M": ("M", 1_000_000_000.0),
    "PM": ("pM", 0.001),
}

OUT_COLS = [
    "source",
    "Patent",
    "patent_core",
    "Target",
    "Activity Type",
    "Relation",
    "Value",
    "value_nM",
    "Unit",
    "SMILES",
    "SMILES_original",
    "InChIKey",
    "UniProt",
    "organism",
    "gene",
    "final_label",
    "source_record_id",
]


def log(msg: str) -> None:
    print(msg, flush=True)


def as_text(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if value != value:
            return ""
        return f"{value:.12g}"
    return str(value).strip()


def patent_core(patent: str) -> str:
    s = re.sub(r"[^A-Z0-9]", "", patent.upper())
    s2 = KIND.sub("", s)
    return s2 if len(s2) > 4 else s


def norm_relation(raw: str) -> str:
    text = as_text(raw)
    if not text:
        return "="
    return REL_MAP.get(text, text)


def norm_activity(raw: str) -> str:
    text = as_text(raw)
    if not text:
        return ""
    key = re.sub(r"[^A-Z0-9]", "", text.upper())
    return ACT_MAP.get(key, text)


def parse_number(raw: str) -> float | None:
    text = as_text(raw).replace(",", "")
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def to_nM(value: str, unit: str) -> str:
    num = parse_number(value)
    if num is None:
        return ""
    key = as_text(unit).upper().replace(" ", "")
    key = key.replace("Μ", "U").replace("µ", "U")
    mapped = UNIT_NANO.get(key)
    if mapped is None:
        if as_text(unit) == "nM":
            return as_text(num)
        return ""
    return as_text(num * mapped[1])


def canonical_smiles(smiles: str, cache: dict[str, tuple[str, str] | None]) -> tuple[str, str] | None:
    if smiles in cache:
        return cache[smiles]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        cache[smiles] = None
        return None
    can = Chem.MolToSmiles(mol, canonical=True)
    key = Chem.MolToInchiKey(mol)
    cache[smiles] = (can, key)
    return cache[smiles]


def iter_parquet(path: Path, columns: list[str] | None = None):
    pf = pq.ParquetFile(path)
    names = columns or pf.schema_arrow.names
    names = [c for c in names if c in set(pf.schema_arrow.names)]
    for batch in pf.iter_batches(batch_size=16384, columns=names):
        arrays = {c: batch.column(c).to_pylist() for c in batch.column_names}
        for i in range(batch.num_rows):
            yield {c: as_text(arrays[c][i]) for c in batch.column_names}


def write_parquet(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        raise RuntimeError(f"no rows for {path}")
    arrays = {c: pa.array([r.get(c, "") for r in rows], type=pa.string()) for c in OUT_COLS}
    table = pa.Table.from_pydict(arrays, schema=pa.schema([pa.field(c, pa.string()) for c in OUT_COLS]))
    pq.write_table(table, path, compression="zstd")


def write_demo(rows: list[dict[str, str]]) -> None:
    DEMO.mkdir(parents=True, exist_ok=True)
    csv_path = DEMO / "merged_SAR.csv"
    json_path = DEMO / "merged_SAR.json"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUT_COLS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "") for c in OUT_COLS})
    payload = {
        "dataset": "merged_SAR",
        "sources": ["BindingDB", "ChEMBL_src38", "HARVEST"],
        "selection": "First 8 complete rows from each source after merge processing.",
        "row_count": len(rows),
        "records": [{k: (v if v else None) for k, v in row.items()} for row in rows],
    }
    json_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def process() -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    cache: dict[str, tuple[str, str] | None] = {}
    seen: set[tuple[str, ...]] = set()
    demo_by_source: dict[str, list[dict[str, str]]] = {"BindingDB": [], "ChEMBL": [], "HARVEST": []}
    stats = {
        "input": {},
        "kept": {},
        "dropped_incomplete": {},
        "dropped_invalid_smiles": {},
        "dropped_duplicate": {},
    }
    out_path = OUT / "SAR.parquet"
    schema = pa.schema([pa.field(c, pa.string()) for c in OUT_COLS])
    writer = pq.ParquetWriter(str(out_path), schema, compression="zstd", use_dictionary=True)
    buffer: list[dict[str, str]] = []
    n_out = 0
    patents: set[str] = set()
    inchis: set[str] = set()

    def flush() -> None:
        nonlocal buffer, n_out
        if not buffer:
            return
        arrays = {c: pa.array([r.get(c, "") for r in buffer], type=pa.string()) for c in OUT_COLS}
        writer.write_table(pa.Table.from_pydict(arrays, schema=schema))
        n_out += len(buffer)
        buffer = []

    streams = [
        (
            "BindingDB",
            STD / "BindingDB.parquet",
            None,
            {
                "organism": "Target Source Organism According to Curator or DataSource",
                "source_record_id": "BindingDB Reactant_set_id",
            },
        ),
        (
            "ChEMBL",
            STD / "ChEMBL.parquet",
            "38",
            {
                "organism": "organism",
                "source_record_id": "activity_id",
            },
        ),
        (
            "HARVEST",
            STD / "HARVEST.parquet",
            None,
            {
                "organism": "organism_scientific",
                "gene": "gene",
                "final_label": "final_label",
                "source_record_id": "source_csv",
            },
        ),
    ]

    for source, path, src_filter, fieldmap in streams:
        log(f"[{source}] reading {path.name}")
        n_in = n_incomplete = n_bad_smi = n_dup = n_keep = 0
        for src in iter_parquet(path):
            if src_filter is not None and src.get("src_id") != src_filter:
                continue
            n_in += 1
            patent = src.get("Patent", "")
            smiles0 = src.get("SMILES", "")
            act = norm_activity(src.get("Activity Type", ""))
            value = as_text(src.get("Value", ""))
            if not (patent and smiles0 and act and value):
                n_incomplete += 1
                continue
            can = canonical_smiles(smiles0, cache)
            if can is None:
                n_bad_smi += 1
                continue
            smiles, inchikey = can
            uniprot = as_text(src.get("UniProt", "")).split("|")[0]
            relation = norm_relation(src.get("Relation", ""))
            unit = as_text(src.get("Unit", ""))
            core = patent_core(patent)
            dedup = (core, inchikey, uniprot, act, relation, value, unit)
            if dedup in seen:
                n_dup += 1
                continue
            seen.add(dedup)
            row = {
                "source": source if source != "ChEMBL" else "ChEMBL_src38",
                "Patent": patent,
                "patent_core": core,
                "Target": src.get("Target", ""),
                "Activity Type": act,
                "Relation": relation,
                "Value": value,
                "value_nM": to_nM(value, unit),
                "Unit": unit,
                "SMILES": smiles,
                "SMILES_original": smiles0,
                "InChIKey": inchikey,
                "UniProt": uniprot,
                "organism": src.get(fieldmap.get("organism", ""), "") if fieldmap.get("organism") else "",
                "gene": src.get(fieldmap.get("gene", ""), "") if fieldmap.get("gene") else "",
                "final_label": src.get(fieldmap.get("final_label", ""), "") if fieldmap.get("final_label") else "",
                "source_record_id": src.get(fieldmap.get("source_record_id", ""), ""),
            }
            if source == "HARVEST" and fieldmap.get("source_record_id") == "source_csv":
                rec = src.get("source_csv", "")
                chem = src.get("chemical_id", "") or src.get("original_alias", "")
                row["source_record_id"] = f"{rec}:{chem}" if chem else rec
            buffer.append(row)
            patents.add(core)
            if inchikey:
                inchis.add(inchikey)
            n_keep += 1
            if len(buffer) >= 20000:
                flush()
            demo_key = "ChEMBL" if source == "ChEMBL" else source
            if len(demo_by_source[demo_key]) < 8:
                demo_by_source[demo_key].append(row)
            if n_in % 200000 == 0:
                log(f"  {source} scanned={n_in:,} kept={n_keep:,}")
        stats["input"][source] = n_in
        stats["kept"][source] = n_keep
        stats["dropped_incomplete"][source] = n_incomplete
        stats["dropped_invalid_smiles"][source] = n_bad_smi
        stats["dropped_duplicate"][source] = n_dup
        log(f"[{source}] in={n_in:,} keep={n_keep:,} incomplete={n_incomplete:,} invalid_smi={n_bad_smi:,} dup={n_dup:,}")

    flush()
    writer.close()
    demo_rows = demo_by_source["BindingDB"] + demo_by_source["ChEMBL"] + demo_by_source["HARVEST"]
    write_demo(demo_rows)
    stats["output_rows"] = n_out
    stats["output_patents"] = len(patents)
    stats["output_inchikeys"] = len(inchis)
    stats["output"] = str(out_path)
    stats["priority"] = "BindingDB > ChEMBL_src38 > HARVEST"
    stats["dedup_key"] = "patent_core + InChIKey + UniProt + Activity Type + Relation + Value + Unit"
    (OUT / "stats.json").write_text(json.dumps(stats, ensure_ascii=False, indent=2), encoding="utf-8")
    log(f"done rows={n_out:,} patents={len(patents):,} inchikeys={len(inchis):,}")
    return stats


if __name__ == "__main__":
    process()
