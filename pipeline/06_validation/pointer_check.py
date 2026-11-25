"""Does the structure file the pipeline named actually hold the claimed structure?

The first pass compared file names, which is wrong: a patent often draws one molecule in several
MOL files, so the name recorded can differ from the first file that happens to produce that
InChIKey while still being correct. The right test parses the named file itself.

This is the automated counterpart of the manual Q2 (compound label to structure), so it is worth
getting right: it tests the pointer from a record to a drawing, at a scale no annotator can reach.
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import json
import multiprocessing as mp
import sys
import warnings
import zipfile
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")
RAW = Path(_ws + "/raw/patent-data/")
W = sys.argv[1].rstrip("/")
WORKERS = int(sys.argv[2]) if len(sys.argv) > 2 else 18
Z = 1.959963985


def wilson(k, n, z=Z):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def job(args):
    archive, items = args
    out = []
    try:
        with zipfile.ZipFile(RAW / archive) as z:
            names = {n.split("/")[-1].upper(): n for n in z.namelist()
                     if n.upper().endswith(".MOL")}
            n_mol = len(names)
            cache = {}
            for i in items:
                want = i["structure_file"].split("/")[-1].upper()
                rec = dict(i, n_mol_in_archive=n_mol)
                if want not in names:
                    rec.update(verdict="named_file_absent", named_inchikey="")
                    out.append(rec)
                    continue
                if want not in cache:
                    try:
                        mol = Chem.MolFromMolBlock(z.read(names[want]).decode("utf-8", "replace"))
                        cache[want] = Chem.MolToInchiKey(mol) if mol is not None else ""
                    except Exception:
                        cache[want] = ""
                got = cache[want]
                claim = i["InChIKey"]
                rec["named_inchikey"] = got
                if not got:
                    rec["verdict"] = "named_file_unparseable"
                elif got == claim:
                    rec["verdict"] = "exact"
                elif got.split("-")[0] == claim.split("-")[0]:
                    rec["verdict"] = "skeleton"
                else:
                    rec["verdict"] = "different_molecule"
                out.append(rec)
    except Exception as exc:
        return [dict(i, verdict="archive_error", named_inchikey=str(exc)[:60],
                     n_mol_in_archive=0) for i in items]
    return out


def main():
    d = pq.read_table(W + "/structure_check_records.parquet").to_pandas()
    d = d[d.structure_file.fillna("").ne("") & d.InChIKey.fillna("").ne("")]
    print(f"records with a named structure file: {len(d):,} "
          f"across {d.archive.nunique():,} archives", flush=True)
    cols = ["source_record_id", "Patent", "InChIKey", "structure_file", "stratum",
            "source_format", "Activity Type", "archive"]
    jobs = [(a, g[cols].to_dict("records")) for a, g in d.groupby("archive", sort=False)]
    res = []
    with mp.Pool(WORKERS) as pool:
        for n, batch in enumerate(pool.imap_unordered(job, jobs, chunksize=4), start=1):
            res.extend(batch)
            if n % 500 == 0:
                print(f"  {n:,}/{len(jobs):,} archives", flush=True)
    r = pd.DataFrame(res)
    r.to_parquet(W + "/pointer_check_records.parquet", index=False)

    def rates(frame):
        n = len(frame)
        k = int(frame.verdict.eq("exact").sum())
        lo, hi = wilson(k, n)
        return {"records": n, "exact": k, "exact_rate": round(k / n, 5) if n else None,
                "wilson_95": [round(lo, 5), round(hi, 5)],
                **{v: int((frame.verdict == v).sum()) for v in
                   ["skeleton", "different_molecule", "named_file_absent",
                    "named_file_unparseable", "archive_error"]}}

    # distinct MOL files per archive, so the integrity figure is not counted once per record
    per_archive = r.drop_duplicates("archive")
    summary = {"overall": rates(r),
               "by_stratum": {s: rates(g) for s, g in r.groupby("stratum")},
               "by_source_format": {s: rates(g) for s, g in r.groupby("source_format")},
               "distinct_mol_files_in_sampled_archives": int(per_archive.n_mol_in_archive.sum()),
               "sampled_archives": int(r.archive.nunique())}
    with open(W + "/pointer_check.json", "w") as fh:
        json.dump(summary, fh, indent=1, default=float)
    o = summary["overall"]
    print(f"\n=== the named structure file holds the claimed structure")
    print(f"  exact {o['exact']:,}/{o['records']:,} = {o['exact_rate']:.4%}  "
          f"95% CI [{o['wilson_95'][0]:.4%}, {o['wilson_95'][1]:.4%}]")
    for k in ["skeleton", "different_molecule", "named_file_absent",
              "named_file_unparseable", "archive_error"]:
        if o[k]:
            print(f"  {k}: {o[k]:,}")
    print("\n=== by source format")
    print(pd.DataFrame([dict(fmt=s, **{k: v for k, v in x.items()
                                       if k in ("records", "exact_rate", "different_molecule",
                                                "skeleton", "named_file_absent")})
                        for s, x in summary["by_source_format"].items()]).to_string(index=False))
    if (r.verdict == "different_molecule").any():
        print("\n=== records whose named file holds a different molecule")
        print(r[r.verdict == "different_molecule"].head(10)[
            ["Patent", "structure_file", "InChIKey", "named_inchikey"]].to_string(index=False))
    print(f"\ndistinct MOL files in the {summary['sampled_archives']:,} sampled archives: "
          f"{summary['distinct_mol_files_in_sampled_archives']:,}")
    print(f"wrote {W}/pointer_check.json")


if __name__ == "__main__":
    main()
