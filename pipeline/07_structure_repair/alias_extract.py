"""Collect the MOL files of every record flagged structure_alias_suspect in PatBench r5, plus an inventory of their alias lines.

    python alias_extract.py <release_dir> <out_dir> [workers]

out_dir/mols/<structure_file>      the original MOL text, unchanged
out_dir/records.tsv                one row per (patent_core, structure_file): stored SMILES / InChIKey, compound label, group ids
out_dir/alias_atoms.tsv            one row per alias line: file, atom number (1-based), element in the atom block, degree, bond orders,
                                   raw label, and whether alias_scan.py's rule counts the label as changing the atom
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import json
import multiprocessing as mp
import re
import sys
import zipfile
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parent))
RDLogger.DisableLog("rdApp.*")
R = Path(sys.argv[1]).resolve()
OUT = Path(sys.argv[2]).resolve()
WORKERS = int(sys.argv[3]) if len(sys.argv) > 3 else 16
RAW = Path(_ws + "/raw/patent-data/")
RECOVERY = _ws + "/outputs/uspto-sar-recovery-20261006/results/"
ELEMENTS = set("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Zr Mo Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba Pt Au Hg Tl Pb Bi".split())
METHYL = {"Me", "CH3", "H3C", "C", "CH", "CH2"}
ANNOTATIONS = {"chiral", "abs", "absolute", "rel", "relative", "rac", "racemic", "racemate", "cis", "trans", "z", "e", "r", "s", "or",
               "and", "enantiomer", "diastereomer", "isomer", "mixture", "endo", "exo", "syn", "anti", "meso"}


def clean_label(label):          # same cleaning as scripts/alias_scan.py
    lab = re.sub(r"\\[sn]", "", label).replace("^", "")
    return re.sub(r"[\s()\[\]+\-−0-9·•.,;:#$&*'\"]", "", lab)


def changes_atom(label, symbol):  # same rule as scripts/alias_scan.py
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
            return True
    return heavy > 1


def job(args):
    patent, archive, files = args
    recs, atoms, texts = [], [], {}
    try:
        with zipfile.ZipFile(RAW / archive) as z:
            for sf in files:
                try:
                    block = z.read(sf).decode("latin1")
                except KeyError:
                    atoms.append(dict(structure_file=sf, atom=0, label="", note="file not in archive"))
                    continue
                texts[sf] = block
                raw = Chem.MolFromMolBlock(block, sanitize=False)
                lines = block.splitlines()
                for k, line in enumerate(lines):
                    if line.startswith("A  ") and k + 1 < len(lines):
                        try:
                            idx = int(line[3:].strip())
                        except ValueError:
                            continue
                        lab = lines[k + 1].rstrip("\r\n")
                        if raw is None or not (1 <= idx <= raw.GetNumAtoms()):
                            atoms.append(dict(structure_file=sf, atom=idx, label=lab, note="atom index out of range"))
                            continue
                        a = raw.GetAtomWithIdx(idx - 1)
                        atoms.append(dict(structure_file=sf, atom=idx, element=a.GetSymbol(), degree=a.GetDegree(),
                                          bond_orders="".join(str(int(b.GetBondTypeAsDouble())) for b in a.GetBonds()),
                                          neighbours=",".join(f"{n.GetIdx()+1}{n.GetSymbol()}" for n in a.GetNeighbors()),
                                          charge=a.GetFormalCharge(), label=lab, changes_atom=bool(changes_atom(lab, a.GetSymbol())), note=""))
    except Exception as exc:
        atoms.append(dict(structure_file="", atom=0, label="", note=f"{patent}: {exc!s:.80}"))
    return texts, atoms


def main():
    s = pq.read_table(R / "sar_all.parquet", columns=["Patent", "patent_core", "structure_file", "compound_alias", "InChIKey", "InChIKey_as_drawn",
                                                       "SMILES", "SMILES_as_drawn", "group_id", "structure_alias_suspect", "is_affinity", "is_exact"]).to_pandas()
    a = s[s.structure_alias_suspect == True]
    rec = (a.groupby(["Patent", "patent_core", "structure_file"], as_index=False)
             .agg(compound_alias=("compound_alias", "first"), InChIKey=("InChIKey", "first"), InChIKey_as_drawn=("InChIKey_as_drawn", "first"),
                  SMILES=("SMILES", "first"), SMILES_as_drawn=("SMILES_as_drawn", "first"), n_rows=("group_id", "size"),
                  groups=("group_id", lambda v: ";".join(map(str, sorted(set(v)))))))
    print(f"flagged rows {len(a):,}; distinct (patent, MOL file) {len(rec):,}; patents {rec.Patent.nunique():,}", flush=True)
    arch = dict(pq.read_table(RECOVERY + "patents.parquet", columns=["Patent", "archive"]).to_pandas().pipe(lambda d: zip(d.Patent, d.archive)))
    jobs = [(p, arch[p], list(g.structure_file)) for p, g in rec.groupby("Patent") if p in arch]
    print("patents with an archive:", len(jobs), "of", rec.Patent.nunique(), flush=True)
    (OUT / "mols").mkdir(parents=True, exist_ok=True)
    atoms, n = [], 0
    with mp.Pool(WORKERS) as pool:
        for texts, at in pool.imap_unordered(job, jobs, chunksize=4):
            for sf, t in texts.items():
                p = OUT / "mols" / sf
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_bytes(t.encode("latin1"))
            atoms.extend(at)
            n += 1
            if n % 200 == 0:
                print(f"  {n}/{len(jobs)} patents", flush=True)
    rec.to_csv(OUT / "records.tsv", sep="\t", index=False)
    pd.DataFrame(atoms).to_csv(OUT / "alias_atoms.tsv", sep="\t", index=False)
    print("written", len(rec), "records,", len(atoms), "alias lines")


if __name__ == "__main__":
    main()
