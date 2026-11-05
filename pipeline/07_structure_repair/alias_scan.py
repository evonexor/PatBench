"""Find stored structures that rest on a MOL abbreviation RDKit could not expand.

ChemDraw writes an abbreviated group ("Ph", "OMe", "HN", "CF3") as one atom plus an alias line.
RDKit keeps the atom's element and drops the group, so the stored structure is missing atoms.
A record is suspect when every MOL in its patent that yields its InChIKey carries such an alias:
if some clean drawing yields the same key, the stored structure is right whatever the others say.
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import json
import multiprocessing as mp
import re
import sys
import warnings
import zipfile
from collections import Counter
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")
D = _ws + "/outputs/affinity-sar-release-20261019/"
RAW = Path(_ws + "/raw/patent-data/")
WORKERS = int(sys.argv[1]) if len(sys.argv) > 1 else 20

ELEMENTS = set("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Ti V Cr Mn Fe Co Ni Cu Zn Ga "
               "Ge As Se Br Kr Rb Sr Zr Mo Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba Pt Au Hg Tl Pb Bi".split())
# a label that only names the atom's own element plus hydrogens, or a methyl on a carbon,
# describes the atom exactly; anything else implies atoms the stored structure does not have
METHYL = {"Me", "CH3", "H3C", "C", "CH", "CH2"}


# stereo and drawing annotations ChemDraw attaches as aliases; they describe, they add no atoms
ANNOTATIONS = {"chiral", "abs", "absolute", "rel", "relative", "rac", "racemic", "racemate", "cis",
               "trans", "z", "e", "r", "s", "or", "and", "enantiomer", "diastereomer", "isomer",
               "mixture", "endo", "exo", "syn", "anti", "meso"}


def clean_label(label):
    # ChemDraw text codes: \s..\n wraps a subscript, ^ marks the bonding side
    lab = re.sub(r"\\[sn]", "", label).replace("^", "")
    return re.sub(r"[\s()\[\]+\-−0-9·•.,;:#$&*'\"]", "", lab)


def alias_changes_atom(label, symbol):
    lab = clean_label(label)
    if not lab or lab == symbol:
        return False
    if lab.lower() in ANNOTATIONS or lab.lower().rstrip("l") in ANNOTATIONS:
        return False
    if symbol == "C" and lab in METHYL:
        return False
    heavy = 0
    for tok in re.findall(r"[A-Z][a-z]?|[a-z]+", lab):
        if tok == "H":
            continue
        if tok in ELEMENTS:
            if tok != symbol:
                return True
            heavy += 1
        else:
            return True          # Ph, Et, Boc, Ts ... : a group, never one atom
    return heavy > 1


def job(args):
    patent, archive, wanted = args
    out = {}
    try:
        with zipfile.ZipFile(RAW / archive) as z:
            for name in z.namelist():
                if not name.upper().endswith(".MOL"):
                    continue
                block = z.read(name).decode("utf-8", "replace")
                mol = Chem.MolFromMolBlock(block)
                if mol is None:
                    continue
                try:
                    ik = Chem.MolToInchiKey(mol)
                except Exception:
                    continue
                if ik not in wanted:
                    continue
                raw = Chem.MolFromMolBlock(block, sanitize=False)
                bad, labels = False, []
                lines = block.splitlines()
                for k, line in enumerate(lines):
                    if line.startswith("A  ") and k + 1 < len(lines):
                        try:
                            idx = int(line[3:].strip()) - 1
                        except ValueError:
                            continue
                        lab = lines[k + 1].strip()
                        sym = raw.GetAtomWithIdx(idx).GetSymbol() if raw is not None and 0 <= idx < raw.GetNumAtoms() else "?"
                        labels.append(f"{sym}:{lab}")
                        if alias_changes_atom(lab, sym):
                            bad = True
                rec = out.setdefault(ik, {"n_mol": 0, "n_bad": 0, "labels": Counter()})
                rec["n_mol"] += 1
                rec["n_bad"] += int(bad)
                rec["labels"].update(labels)
    except Exception as exc:
        return [{"Patent": patent, "InChIKey": "", "error": str(exc)[:80]}]
    return [{"Patent": patent, "InChIKey": ik, "n_mol": v["n_mol"], "n_bad": v["n_bad"],
             "labels": "|".join(f"{l}" for l, _ in v["labels"].most_common(6)), "error": ""}
            for ik, v in out.items()]


def main():
    sar = pq.read_table(D + "sar_all.parquet", columns=[
        "Patent", "patent_core", "InChIKey", "is_affinity", "is_exact", "is_druglike"]).to_pandas()
    sar = sar[sar.is_affinity & sar.is_exact & sar.is_druglike]
    pt = pq.read_table(_ws + "/outputs/uspto-sar-recovery-20261006/results/patents.parquet",
                       columns=["Patent", "archive"]).to_pandas()
    arch = dict(zip(pt.Patent, pt.archive))
    wanted = sar.groupby("Patent").InChIKey.agg(lambda s: frozenset(s.dropna()))
    jobs = [(p, arch[p], w) for p, w in wanted.items() if p in arch]
    print(f"scanning {len(jobs):,} patents for {sar.InChIKey.nunique():,} molecules", flush=True)
    res = []
    with mp.Pool(WORKERS) as pool:
        for n, b in enumerate(pool.imap_unordered(job, jobs, chunksize=4), start=1):
            res.extend(b)
            if n % 1000 == 0:
                print(f"  {n:,}/{len(jobs):,}", flush=True)
    r = pd.DataFrame(res)
    r.to_parquet(D + "work/alias_scan.parquet", index=False)
    ok = r[r.error.eq("")]
    ok = ok.assign(suspect=ok.n_bad.eq(ok.n_mol) & ok.n_bad.gt(0))
    print(f"\n(patent, molecule) pairs found: {len(ok):,}; errors: {int(r.error.ne('').sum())}")
    print(f"drawn with any changing alias: {int(ok.n_bad.gt(0).sum()):,}")
    print(f"SUSPECT (every drawing of it carries a changing alias): {int(ok.suspect.sum()):,} "
          f"= {ok.suspect.mean():.2%}")
    lab = Counter()
    for s in ok.loc[ok.suspect, "labels"]:
        lab.update(x for x in s.split("|") if x)
    print("most common labels on suspect structures:", lab.most_common(25))
    keep = Counter()
    for s in ok.loc[~ok.suspect & ok.labels.ne(""), "labels"]:
        keep.update(x for x in s.split("|") if x)
    print("labels judged harmless (sample):", [k for k, _ in keep.most_common(40)
                                               if not alias_changes_atom(k.split(':', 1)[1], k.split(':', 1)[0])][:20])


if __name__ == "__main__":
    main()
