_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
#!/usr/bin/env python3
from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
STD = Path(_ws + "/raw/patent-data/standardized")
cache = {}


def smiles_ok(smi) -> bool:
    if smi is None:
        return False
    text = str(smi).strip()
    if not text:
        return False
    if text not in cache:
        mol = Chem.MolFromSmiles(text)
        cache[text] = mol is not None and mol.GetNumHeavyAtoms() >= 1
    return cache[text]


def filter_file(path: Path, dest: Path) -> None:
    pf = pq.ParquetFile(path)
    writer = pq.ParquetWriter(str(dest), pf.schema_arrow, compression="zstd", use_dictionary=True)
    kept = dropped = 0
    for batch in pf.iter_batches(batch_size=8192):
        smis = batch.column("SMILES").to_pylist()
        mask = [smiles_ok(s) for s in smis]
        n_keep = sum(1 for x in mask if x)
        dropped += len(mask) - n_keep
        if n_keep == 0:
            continue
        if n_keep == len(mask):
            writer.write_batch(batch)
        else:
            writer.write_table(pa.Table.from_batches([batch]).filter(pa.array(mask)))
        kept += n_keep
    writer.close()
    print(f"{path.name} -> {dest.name}: kept={kept:,} dropped={dropped:,}", flush=True)


if __name__ == "__main__":
    jobs = [
        (STD / "BindingDB.parquet", STD / "BindingDB.valid.parquet"),
        (STD / "ChEMBL.parquet", STD / "ChEMBL.valid.parquet"),
        (STD / "HARVEST.parquet", STD / "HARVEST.valid.parquet"),
        (STD / "SureChEMBL.parquet", STD / "SureChEMBL.valid.parquet"),
        (STD / "merged/SAR.parquet", STD / "merged/SAR.valid.parquet"),
    ]
    for src, dest in jobs:
        if not src.exists():
            print(f"missing {src}", flush=True)
            continue
        filter_file(src, dest)
    print("unique_checked", len(cache), flush=True)
