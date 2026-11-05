"""Automated structure verification at scale.

For a stratified random sample of formal SAR records, recompute the InChIKey of every structure
drawn in the patent archive and ask whether the record's claimed structure is among them. This
needs no human judgement, so it can run over tens of thousands of records and stands on its own
as a technical-validation measure.

Reported three ways per record:
  exact       the claimed InChIKey is among the structures the patent draws
  skeleton    only the connectivity layer (first 14 characters) matches, so the patent's drawing
              and the stored structure differ in stereochemistry, isotopes or protonation
  none        neither matches
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import json
import multiprocessing as mp
import sys
import warnings
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
from rdkit import Chem, RDLogger

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")

R = _ws + "/outputs/uspto-sar-recovery-20261006/results/"
RAW = Path(_ws + "/raw/patent-data/")
OUT = sys.argv[1].rstrip("/")
N_TARGET = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
WORKERS = int(sys.argv[3]) if len(sys.argv) > 3 else 20
SEED = 20261019
Z = 1.959963985


def wilson(k, n, z=Z):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def archive_job(args):
    """Parse every MOL in one archive and judge each sampled record from that patent."""
    archive, items = args
    keys, skeletons, n_mol, n_bad = {}, {}, 0, 0
    try:
        with zipfile.ZipFile(RAW / archive) as z:
            for name in z.namelist():
                if not name.upper().endswith(".MOL"):
                    continue
                n_mol += 1
                try:
                    mol = Chem.MolFromMolBlock(z.read(name).decode("utf-8", "replace"))
                    if mol is None:
                        n_bad += 1
                        continue
                    ik = Chem.MolToInchiKey(mol)
                    if not ik:
                        n_bad += 1
                        continue
                    keys.setdefault(ik, name.split("/")[-1].upper())
                    skeletons.setdefault(ik.split("-")[0], name.split("/")[-1].upper())
                except Exception:
                    n_bad += 1
    except Exception as exc:
        return [dict(i, verdict="archive_error", detail=str(exc)[:80],
                     n_mol_in_archive=0, n_mol_unparseable=0) for i in items]
    out = []
    for i in items:
        claim = i["InChIKey"] if isinstance(i["InChIKey"], str) else ""
        rec = dict(i, n_mol_in_archive=n_mol, n_mol_unparseable=n_bad)
        if not claim:
            rec.update(verdict="no_claim", detail="")
        elif claim in keys:
            rec.update(verdict="exact", detail=keys[claim])
        elif claim.split("-")[0] in skeletons:
            rec.update(verdict="skeleton", detail=skeletons[claim.split("-")[0]])
        else:
            rec.update(verdict="none", detail="")
        # when the pipeline named a structure file, did that file hold the claimed structure?
        # a missing value arrives as NaN, which is truthy, so `x or ""` will not do
        sf = i.get("structure_file")
        sf = sf if isinstance(sf, str) else ""
        rec["named_file_matches"] = bool(sf) and keys.get(claim, "") == sf.split("/")[-1].upper()
        out.append(rec)
    return out


def main():
    sar = ds.dataset(R + "SAR_merged.parquet").to_table(
        columns=["source", "source_record_id", "Patent", "patent_core", "Activity Type",
                 "InChIKey"]).to_pandas()
    pt = ds.dataset(R + "patents.parquet").to_table(
        columns=["Patent", "archive", "source_format"]).to_pandas()
    sar = sar.merge(pt, on="Patent", how="left")
    ev = ds.dataset(R + "new_SAR.parquet").to_table(
        columns=["source_record_id", "structure_file"]).to_pandas().set_index("source_record_id")
    sar = sar.join(ev, on="source_record_id")
    sar["pass_name"] = np.where(sar.source.eq("USPTO_MOL_RECOVERY"), "recovery", "baseline")
    sar["stratum"] = sar.pass_name + " | " + sar.source_format.fillna("?")
    print(f"population {len(sar):,} rows", flush=True)
    pop = sar.stratum.value_counts()
    print(pop.to_string(), flush=True)

    rng = np.random.default_rng(SEED)
    frames = []
    floor = 500
    for s, n_pop in pop.items():
        want = max(floor, int(round(N_TARGET * n_pop / len(sar))))
        frames.append(sar[sar.stratum == s].sample(n=min(want, n_pop), random_state=SEED))
    sample = pd.concat(frames).reset_index(drop=True)
    print(f"\nsampled {len(sample):,} records from {sample.archive.nunique():,} archives", flush=True)
    print(sample.stratum.value_counts().to_string(), flush=True)

    cols = ["source_record_id", "Patent", "patent_core", "Activity Type", "InChIKey",
            "structure_file", "pass_name", "source_format", "stratum", "archive"]
    jobs = [(arch, grp[cols].to_dict("records"))
            for arch, grp in sample.groupby("archive", sort=False)]
    total_mols = 0
    results = []
    with mp.Pool(WORKERS) as pool:
        for n, batch in enumerate(pool.imap_unordered(archive_job, jobs, chunksize=4), start=1):
            results.extend(batch)
            if n % 500 == 0:
                print(f"  {n:,}/{len(jobs):,} archives, {len(results):,} records", flush=True)
    d = pd.DataFrame(results)
    d.to_parquet(OUT + "/structure_check_records.parquet", index=False)

    def rates(frame):
        n = len(frame)
        k_exact = int(frame.verdict.eq("exact").sum())
        k_any = int(frame.verdict.isin(["exact", "skeleton"]).sum())
        lo, hi = wilson(k_exact, n)
        return {"records": n, "exact": k_exact, "exact_rate": round(k_exact / n, 5) if n else None,
                "exact_wilson_95": [round(lo, 5), round(hi, 5)],
                "skeleton_only": int(frame.verdict.eq("skeleton").sum()),
                "exact_or_skeleton_rate": round(k_any / n, 5) if n else None,
                "no_match": int(frame.verdict.eq("none").sum()),
                "archive_error": int(frame.verdict.eq("archive_error").sum()),
                "no_claim": int(frame.verdict.eq("no_claim").sum())}

    summary = {"sampled_records": len(d), "archives": int(d.archive.nunique()),
               "target_n": N_TARGET, "seed": SEED, "floor_per_stratum": floor,
               "overall": rates(d),
               "by_stratum": {s: rates(g) for s, g in d.groupby("stratum")},
               "by_pass": {s: rates(g) for s, g in d.groupby("pass_name")},
               "by_source_format": {s: rates(g) for s, g in d.groupby("source_format")},
               "by_endpoint": {s: rates(g) for s, g in
                               d.groupby("Activity Type") if len(g) >= 50},
               "mol_file_integrity": {
                   "mol_files_parsed": int(d.n_mol_in_archive.sum()),
                   "distinct_archives": int(d.archive.nunique()),
                   "unparseable_mol_files": int(d.n_mol_unparseable.sum()),
                   "records_whose_named_structure_file_matched":
                       int(d.named_file_matches.sum()),
                   "records_with_a_named_structure_file":
                       int(d.structure_file.fillna("").ne("").sum())}}
    with open(OUT + "/structure_check.json", "w") as fh:
        json.dump(summary, fh, indent=1, default=float)

    print("\n=== overall")
    o = summary["overall"]
    print(f"  exact match {o['exact']:,}/{o['records']:,} = {o['exact_rate']:.4%}  "
          f"95% CI [{o['exact_wilson_95'][0]:.4%}, {o['exact_wilson_95'][1]:.4%}]")
    print(f"  skeleton-only {o['skeleton_only']:,}   no match {o['no_match']:,}   "
          f"archive errors {o['archive_error']:,}")
    print("\n=== by extraction pass and source format")
    rows = [dict(stratum=s, **{k: v for k, v in r.items() if k in
                               ("records", "exact", "exact_rate", "skeleton_only", "no_match")})
            for s, r in summary["by_stratum"].items()]
    print(pd.DataFrame(rows).to_string(index=False))
    if (d.verdict == "none").any():
        print("\n=== sample of records with no matching structure in the patent")
        print(d[d.verdict == "none"].head(12)[
            ["Patent", "stratum", "Activity Type", "InChIKey", "n_mol_in_archive"]
        ].to_string(index=False))
    print(f"\nwrote {OUT}/structure_check.json and structure_check_records.parquet")


if __name__ == "__main__":
    main()
