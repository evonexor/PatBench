"""Alias inventory of EVERY structure file behind PatBench r5, not only the ones flagged structure_alias_suspect.

    python alias_scan_all.py <release_dir> <out_dir> [workers]

Why: the first scan (scripts/alias_scan.py) only looked for labels that "add atoms". A bare "H" label on a two-bond placeholder was counted as
harmless, but the placeholder is always written as element C, and the label is the real atom (NH). So the damage is wider than the flag says.

out_dir/all_alias_atoms.tsv   one row per alias line of any MOL file used by a record: patent_core, structure_file, atom, element, degree, label
out_dir/all_records.tsv       one row per (patent_core, structure_file) that has >= 1 alias line: stored SMILES_as_drawn / InChIKey_as_drawn, n_rows
out_dir/mols_all/<structure_file>   the MOL text of those files (unchanged)
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import multiprocessing as mp
import re
import sys
import zipfile
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

R = Path(sys.argv[1]).resolve()
OUT = Path(sys.argv[2]).resolve()
WORKERS = int(sys.argv[3]) if len(sys.argv) > 3 else 16
RAW = Path(_ws + "/raw/patent-data/")
RECOVERY = _ws + "/outputs/uspto-sar-recovery-20261006/results/"
ALIAS = re.compile(r"A\s+(\d+)\s*$")


def job(args):
    patent, archive, files = args
    atoms, texts = [], {}
    try:
        with zipfile.ZipFile(RAW / archive) as z:
            names = set(z.namelist())
            for sf in files:
                if sf not in names:
                    continue
                block = z.read(sf).decode("latin1")
                lines = block.splitlines()
                hits = [(k, ALIAS.match(l)) for k, l in enumerate(lines) if ALIAS.match(l) and k + 1 < len(lines)]
                if not hits:
                    continue
                texts[sf] = block
                try:
                    na = int(lines[3][:3]); nb = int(lines[3][3:6])
                except Exception:
                    continue
                elem, deg = {}, {}
                for i in range(na):
                    elem[i + 1] = lines[4 + i][31:34].strip() if len(lines[4 + i]) >= 34 else ""
                for b in lines[4 + na:4 + na + nb]:
                    try:
                        a1, a2 = int(b[:3]), int(b[3:6])
                    except ValueError:
                        continue
                    deg[a1] = deg.get(a1, 0) + 1
                    deg[a2] = deg.get(a2, 0) + 1
                for k, m in hits:
                    idx = int(m.group(1))
                    atoms.append(dict(structure_file=sf, atom=idx, element=elem.get(idx, "?"), degree=deg.get(idx, 0), label=lines[k + 1].rstrip("\r\n")))
    except Exception as exc:
        atoms.append(dict(structure_file="", atom=0, element="", degree=0, label=f"ERROR {patent}: {exc!s:.60}"))
    return texts, atoms


def main():
    s = pq.read_table(R / "sar_all.parquet", columns=["Patent", "patent_core", "structure_file", "compound_alias", "InChIKey_as_drawn", "SMILES_as_drawn", "InChIKey",
                                                       "structure_alias_suspect", "group_id"]).to_pandas()
    s = s[s.structure_file.ne("")]
    rec = (s.groupby(["Patent", "patent_core", "structure_file"], as_index=False)
             .agg(compound_alias=("compound_alias", "first"), InChIKey_as_drawn=("InChIKey_as_drawn", "first"), SMILES_as_drawn=("SMILES_as_drawn", "first"),
                  n_rows=("group_id", "size"), flagged=("structure_alias_suspect", "max")))
    print(f"records with a structure file: {len(rec):,}; patents {rec.Patent.nunique():,}", flush=True)
    arch = dict(zip(*[pq.read_table(RECOVERY + "patents.parquet", columns=["Patent", "archive"]).to_pandas()[c] for c in ("Patent", "archive")]))
    jobs = [(p, arch[p], list(g.structure_file)) for p, g in rec.groupby("Patent") if p in arch]
    (OUT / "mols_all").mkdir(parents=True, exist_ok=True)
    atoms, n = [], 0
    with mp.Pool(WORKERS) as pool:
        for texts, at in pool.imap_unordered(job, jobs, chunksize=4):
            for sf, t in texts.items():
                p = OUT / "mols_all" / sf
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_bytes(t.encode("latin1"))
            atoms.extend(at)
            n += 1
            if n % 500 == 0:
                print(f"  {n}/{len(jobs)} patents", flush=True)
    a = pd.DataFrame(atoms)
    files = set(a.structure_file) - {""}
    rec = rec[rec.structure_file.isin(files)]
    rec.to_csv(OUT / "all_records.tsv", sep="\t", index=False)
    a.to_csv(OUT / "all_alias_atoms.tsv", sep="\t", index=False)
    print("records with >=1 alias:", len(rec), " alias lines:", len(a), " of which flagged records:", int(rec.flagged.sum()))


if __name__ == "__main__":
    main()
