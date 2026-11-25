"""Parent form and chemistry flags for every molecule in the formal SAR.

The parent is the largest organic fragment, neutralised: a hydrochloride and its free base are one
compound, and counting them twice in a series double-weights it. Tautomers are left as drawn
(canonicalising them rewrites structures and is not reliable enough to apply blind).

Writes work/mol_std.parquet keyed by the InChIKey as drawn.
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import multiprocessing as mp
import os
import sys
import warnings

import pandas as pd
import pyarrow.dataset as ds
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")
R = os.environ.get("PATBENCH_SAR_DIR", _ws + "/outputs/uspto-sar-recovery-20261006/results/")
OUT = sys.argv[1].rstrip("/")
WORKERS = int(sys.argv[2]) if len(sys.argv) > 2 else 20
ORGANIC = set("C H N O S P F Cl Br I B Si Se".split())
# covalent warheads whose IC50 depends on incubation time; boronic acids bind reversibly-covalently
WARHEADS = {
    "acrylamide": "[CH2]=[CH]C(=O)[NX3]",
    "substituted_acrylamide": "[CX3;!R;!$([CH2])]=[CX3;!R]C(=O)[NX3]",
    "butynamide": "C#CC(=O)[NX3]",
    "chloroacetamide": "Cl[CH2]C(=O)[NX3]",
    "vinyl_sulfone": "[CX3;!R]=[CX3;!R]S(=O)(=O)",
    "boronic_acid": "[BX3]([OX2])[OX2]",
}
_PAT = None


def init():
    global _PAT, _CHOOSE, _UNCHARGE
    _PAT = {k: Chem.MolFromSmarts(v) for k, v in WARHEADS.items()}
    _CHOOSE = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
    _UNCHARGE = rdMolStandardize.Uncharger()


def job(item):
    key, smi = item
    mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
    if mol is None:
        return (key, smi, key, -1, True, "", 0, False, False, "", 1)
    n_comp = len(Chem.GetMolFrags(mol))
    try:
        par = _UNCHARGE.uncharge(_CHOOSE.choose(mol))
        psmi = Chem.MolToSmiles(par)
        pkey = Chem.MolToInchiKey(par) or key
    except Exception:
        par, psmi, pkey = mol, smi, key
    try:
        core = MurckoScaffold.GetScaffoldForMol(par)
        sc = Chem.MolToSmiles(core) if core.GetNumAtoms() else ""
    except Exception:
        sc = ""
    centres = Chem.FindMolChiralCenters(par, includeUnassigned=True, useLegacyImplementation=False)
    undefined = any(lbl == "?" for _, lbl in centres)
    isotope = any(a.GetIsotope() for a in par.GetAtoms())
    hits = [k for k, p in _PAT.items() if par.HasSubstructMatch(p)]
    return (key, psmi, pkey, par.GetNumHeavyAtoms(),
            any(a.GetSymbol() not in ORGANIC for a in par.GetAtoms()), sc, len(centres),
            undefined, isotope, ";".join(hits), n_comp)


def main():
    sar = ds.dataset(R + "SAR_merged.parquet").to_table(columns=["InChIKey", "SMILES"]).to_pandas()
    uniq = sar.dropna(subset=["InChIKey"]).drop_duplicates("InChIKey")
    print(f"standardising {len(uniq):,} molecules on {WORKERS} workers", flush=True)
    with mp.Pool(WORKERS, initializer=init) as pool:
        recs = pool.map(job, list(zip(uniq.InChIKey, uniq.SMILES)), chunksize=500)
    out = pd.DataFrame(recs, columns=[
        "InChIKey", "parent_smiles", "parent_inchikey", "heavy_atoms", "unusual_elements", "scaffold",
        "n_stereocentres", "has_undefined_stereo", "is_isotope_labelled", "covalent_warhead",
        "n_components"])
    out.to_parquet(OUT + "/work/mol_std.parquet", index=False)
    changed = out.parent_inchikey.ne(out.InChIKey)
    print(f"multi-component as drawn: {int((out.n_components > 1).sum()):,}; parent key differs from "
          f"drawn key: {int(changed.sum()):,}; distinct parents {out.parent_inchikey.nunique():,}")
    print(f"undefined stereo {int(out.has_undefined_stereo.sum()):,}; isotope {int(out.is_isotope_labelled.sum()):,}; "
          f"warheads {out.covalent_warhead.ne('').sum():,}: "
          f"{out.covalent_warhead[out.covalent_warhead.ne('')].str.split(';').explode().value_counts().to_dict()}")


if __name__ == "__main__":
    main()
