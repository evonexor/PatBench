"""Expand ChemDraw abbreviations (MOL alias lines) that RDKit drops, for the records flagged structure_alias_suspect.

    python alias_repair.py <alias_dir> <records.tsv> <mols_dir> <out.tsv> [workers]

Reads  <records.tsv> and <mols_dir>/<structure_file>  (written by alias_extract.py / alias_scan_all.py)
Writes <out.tsv>, one row per record:
    patent_core, structure_file, status, smiles, inchikey, rules, labels_expanded, labels_unresolved, note
    status: expanded (a group or NH was restored) | counterion_only | as_drawn (nothing to restore) | unresolved (record stays flagged)

Repairs applied to a drawing
  1. abbreviation groups (F3C, MeO, CD3, Boc ...) from the catalogue below
  2. a bare "H" label on a placeholder with two bonds is an N-H (the placeholder is always written as element C)
  3. a nitrile that ChemDraw wrote reversed, N(+)#C(-) bonded through N, is turned into C#N bonded through C

Reading a MOL block
  ChemDraw writes an abbreviated group as ONE placeholder atom (always element C in these files) plus an alias entry
  "A  n" / label. RDKit keeps the placeholder and drops the label. Here the placeholder atom is changed in place into the
  first atom of the group (so stereo perceived from wedge bonds at its neighbours stays valid) and the other atoms of the
  group are added to it.

Label codes (checked against the data):  "\\s" ... "\\n" wrap a subscript run; "^" marks the atom that bonds to the molecule
  (the group is written with the bonding atom last: F3^C = CF3, Me^O = OMe, F3C^O = OCF3, Me2^N = NMe2, O2^S = SO2).
  Labels without a caret are matched whole (MeO, OMe, CF3, F3C, O2N, NO2 ...). Anything not in the catalogue below is left
  unresolved and the record stays flagged: nothing is guessed.
"""
import multiprocessing as mp
import re
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
O = Path(sys.argv[1]).resolve()
RECORDS = Path(sys.argv[2]).resolve()
MOLS = Path(sys.argv[3]).resolve()
OUTFILE = Path(sys.argv[4]).resolve()
WORKERS = int(sys.argv[5]) if len(sys.argv) > 5 else 16

# Each entry: the whole group as SMILES; [*] marks a bond to the rest of the molecule. The placeholder atom becomes the atom
# bonded to the first [*]; the number of [*] must equal the degree of the placeholder in the drawing.
G1 = {   # monovalent groups
    # fluoroalkyl / deuterated / alkoxy: label spellings seen in the files, all meaning the same group
    "CF3": "[*]C(F)(F)F", "F3^C": "[*]C(F)(F)F", "F3C": "[*]C(F)(F)F",
    "OCF3": "[*]OC(F)(F)F", "F3C^O": "[*]OC(F)(F)F", "F3CO": "[*]OC(F)(F)F",
    "CHF2": "[*]C(F)F", "F2H^C": "[*]C(F)F", "F2HC": "[*]C(F)F", "HF2C": "[*]C(F)F",
    "OCHF2": "[*]OC(F)F", "F2HC^O": "[*]OC(F)F", "F2HCO": "[*]OC(F)F",
    "CD3": "[*]C([2H])([2H])[2H]", "D3^C": "[*]C([2H])([2H])[2H]", "D3C": "[*]C([2H])([2H])[2H]",
    "OCD3": "[*]OC([2H])([2H])[2H]", "D3C^O": "[*]OC([2H])([2H])[2H]", "D3CO": "[*]OC([2H])([2H])[2H]",
    "C2F5": "[*]C(F)(F)C(F)(F)F",
    "OMe": "[*]OC", "Me^O": "[*]OC", "MeO": "[*]OC", "OCH3": "[*]OC", "H3C^O": "[*]OC", "CH3O": "[*]OC",
    "OEt": "[*]OCC", "Et^O": "[*]OCC", "EtO": "[*]OCC",
    "OiPr": "[*]OC(C)C", "OtBu": "[*]OC(C)(C)C", "OBn": "[*]OCc1ccccc1", "Bn^O": "[*]OCc1ccccc1",
    "OAc": "[*]OC(C)=O", "Ac^O": "[*]OC(C)=O", "AcO": "[*]OC(C)=O",
    "Me^S": "[*]SC", "SMe": "[*]SC", "F5^S": "[*]S(F)(F)(F)(F)F", "SF5": "[*]S(F)(F)(F)(F)F", "F3C^S": "[*]SC(F)(F)F", "SCF3": "[*]SC(F)(F)F",
    "NMe2": "[*]N(C)C", "Me2^N": "[*]N(C)C", "Me2N": "[*]N(C)C",
    "MeH^N": "[*]NC", "NHMe": "[*]NC", "AcH^N": "[*]NC(C)=O", "NHAc": "[*]NC(C)=O",
    "BocH^N": "[*]NC(=O)OC(C)(C)C", "NHBoc": "[*]NC(=O)OC(C)(C)C",
    "(PMB)2^N": "[*]N(Cc1ccc(OC)cc1)Cc1ccc(OC)cc1",
    "NO2": "[*][N+](=O)[O-]", "O2N": "[*][N+](=O)[O-]", "O2^N": "[*][N+](=O)[O-]",
    "Ph": "[*]c1ccccc1", "tBu": "[*]C(C)(C)C", "nC4H9": "[*]CCCC", "CO2Et": "[*]C(=O)OCC",
    "+N(CH3)3": "[*][N+](C)(C)C", "+N(CH2CH3)3": "[*][N+](CC)(CC)CC",
    "O-": "[*][O-]", "O^-": "[*][O-]", "O_": "[*][O-]", "H2N": "[*]N", "NH2": "[*]N",
    "H^O|": "[*]O", "O|": "[*]O",     # hydroxyl; "|" is a ChemDraw line-end code (tested against BindingDB/ChEMBL, see the validation)
}
G2 = {   # divalent groups (placeholder with two bonds)
    "O2^S": "[*]S(=O)(=O)[*]", "SO2": "[*]S(=O)(=O)[*]", "N": "[*]N[*]",
    "Boc^N": "[*]N(C(=O)OC(C)(C)C)[*]", "Me^N": "[*]N(C)[*]", "NMe": "[*]N(C)[*]",
    "O|": "[*]O[*]",
    "HN": "[*]N[*]", "NH": "[*]N[*]",    # stacked "H" over "N": NH (the placeholder atom is written as C); checked against BindingDB/ChEMBL
    "MeN": "[*]N(C)[*]", "O2S": "[*]S(=O)(=O)[*]", "BocN": "[*]N(C(=O)OC(C)(C)C)[*]",
}
SALT = re.compile(r"^[\s·•\x95.,;+\d]*((H?Cl|HBr|HI|TFA|TEA|HCOOH|HOAc|H2O|THF|DMF|DMSO|MeOH|EtOH|Et2O|HNO3|H2SO4|H3PO4|iPr2NEt|NEt3|Et3N)[\s·•.,;^]*)+$")


ELEMENTS = set("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Zr Mo Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba Pt Au Hg Tl Pb Bi".split())
METHYL = {"Me", "CH3", "H3C", "C", "CH", "CH2"}
ANNOTATIONS = {"chiral", "abs", "absolute", "rel", "relative", "rac", "racemic", "racemate", "cis", "trans", "z", "e", "r", "s", "or", "and",
               "enantiomer", "diastereomer", "isomer", "mixture", "endo", "exo", "syn", "anti", "meso"}


def changes_atom(label, symbol="C"):
    """The rule scripts/alias_scan.py uses: does the label imply atoms the placeholder atom does not have?
    A bare 'H', a stereo note ('abs', 'rac'), a methyl on a carbon etc. describe the atom as drawn and need no expansion."""
    lab = re.sub(r"\\[sn]", "", label).replace("^", "")
    lab = re.sub(r"[\s()\[\]+\-\u2212 0-9\u00b7\u2022.,;:#$&*'\"]", "", lab)
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


def norm(label):
    lab = re.sub(r"\\[sn]", "", label.strip())
    lab = re.sub(r"\s+", "", lab)            # also drops the CR of a stacked label: "H\rN" -> "HN"
    return lab


def lookup(lab, degree):
    if degree == 1 and lab in G1:
        return G1[lab]
    if degree == 2 and lab in G2:
        return G2[lab]
    return None


def apply_group(rw, idx, spec):
    """Turn the placeholder atom idx into the group `spec` in place."""
    frag = Chem.MolFromSmiles(spec)
    dummies = [a.GetIdx() for a in frag.GetAtoms() if a.GetAtomicNum() == 0]
    core = frag.GetAtomWithIdx(dummies[0]).GetNeighbors()[0].GetIdx()
    if any(frag.GetAtomWithIdx(d).GetNeighbors()[0].GetIdx() != core for d in dummies):
        raise ValueError("all attachment points must sit on the core atom")
    if len(dummies) != rw.GetAtomWithIdx(idx).GetDegree():
        raise ValueError("attachment points do not match the drawing")
    mapping = {core: idx}
    for a in frag.GetAtoms():
        if a.GetIdx() in dummies or a.GetIdx() == core:
            continue
        na = Chem.Atom(a.GetAtomicNum())
        na.SetFormalCharge(a.GetFormalCharge()); na.SetIsotope(a.GetIsotope()); na.SetIsAromatic(a.GetIsAromatic())
        na.SetNumExplicitHs(a.GetNumExplicitHs()); na.SetNoImplicit(a.GetNoImplicit())
        mapping[a.GetIdx()] = rw.AddAtom(na)
    ca, fa = rw.GetAtomWithIdx(idx), frag.GetAtomWithIdx(core)
    ca.SetAtomicNum(fa.GetAtomicNum()); ca.SetFormalCharge(fa.GetFormalCharge()); ca.SetIsotope(fa.GetIsotope())
    ca.SetIsAromatic(fa.GetIsAromatic()); ca.SetNumExplicitHs(fa.GetNumExplicitHs()); ca.SetNoImplicit(fa.GetNoImplicit())
    for b in frag.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in dummies or j in dummies:
            continue
        rw.AddBond(mapping[i], mapping[j], b.GetBondType())


NITRILE = Chem.MolFromSmarts("[#6][N+]#[C-]")


def fix_nitrile(mol):
    """ChemDraw writes a reversed cyano label as N(+)#C(-) bonded through N: turn it into C#N bonded through C. Returns True if changed."""
    changed = False
    for _ in range(10):
        m = mol.GetSubstructMatch(NITRILE)
        if not m:
            break
        _, n, c = m
        mol.GetAtomWithIdx(n).SetAtomicNum(6); mol.GetAtomWithIdx(n).SetFormalCharge(0)
        mol.GetAtomWithIdx(c).SetAtomicNum(7); mol.GetAtomWithIdx(c).SetFormalCharge(0)
        Chem.SanitizeMol(mol)
        changed = True
    return changed


def process(rec):
    sf = rec["structure_file"]
    base = dict(patent_core=rec["patent_core"], structure_file=sf, status="unresolved", smiles="", inchikey="", rules="", labels_expanded="",
                labels_unresolved="", note="")
    try:
        block = (MOLS / sf).read_bytes().decode("latin1")
    except Exception as exc:
        return {**base, "note": f"cannot read file: {exc!s:.60}"}
    raw = Chem.MolFromMolBlock(block, sanitize=False)
    if raw is None:
        return {**base, "note": "RDKit cannot parse the MOL block"}
    lines = block.split("\n")           # LF ends a file line; a CR inside the label text is a stacked-label line break ("H" over "N")
    aliases = []
    for k, line in enumerate(lines):
        m = re.match(r"A\s+(\d+)\s*\r?$", line)
        if m and k + 1 < len(lines):
            aliases.append((int(m.group(1)) - 1, lines[k + 1].rstrip("\r")))
    rw = Chem.RWMol(raw)
    done, left, ignored, rules = [], [], [], set()
    for idx, label in aliases:
        if not 0 <= idx < raw.GetNumAtoms():
            left.append(f"{idx+1}:{label}")
            continue
        at = raw.GetAtomWithIdx(idx)
        lab = norm(label)
        if at.GetDegree() == 0 and SALT.match(lab):
            ignored.append(f"{idx+1}:{label}")
            continue
        spec = lookup(lab, at.GetDegree())
        if spec is not None:
            try:
                apply_group(rw, idx, spec)
                done.append(f"{idx+1}:{label}->{spec}")
                rules.add("NH" if lab in ("HN", "NH") else "group")
            except Exception as exc:
                left.append(f"{idx+1}:{label} ({exc!s:.40})")
            continue
        if lab == "H" or changes_atom(lab, at.GetSymbol()):
            left.append(f"{idx+1}:{label}")      # a bare H on anything but a two-bond atom, or a label we cannot expand
        # else: the label describes the atom as drawn (stereo note, methyl on carbon ...) -- nothing to do
    if left:
        return {**base, "labels_expanded": ";".join(done), "labels_unresolved": ";".join(left), "rules": ";".join(sorted(rules)),
                "note": "some alias not in the catalogue" if not done else "partly expanded; record stays flagged"}
    try:
        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
        if fix_nitrile(mol):
            rules.add("nitrile")
        smi = Chem.MolToSmiles(mol)
        ik = Chem.MolToInchiKey(mol)
    except Exception as exc:
        return {**base, "labels_expanded": ";".join(done), "note": f"sanitize/InChI failed: {exc!s:.60}"}
    status = "expanded" if (done or "nitrile" in rules) else ("counterion_only" if ignored else "as_drawn")
    return {**base, "status": status, "smiles": smi, "inchikey": ik, "rules": ";".join(sorted(rules)), "labels_expanded": ";".join(done),
            "note": ("counter-ion/solvent labels ignored: " + ";".join(ignored)) if ignored else ""}


def main():
    rec = pd.read_csv(RECORDS, sep="\t", dtype=str, keep_default_na=False)
    rows = rec.to_dict("records")
    with mp.Pool(WORKERS) as pool:
        out = pool.map(process, rows, chunksize=100)
    d = pd.DataFrame(out)
    for c in ("labels_expanded", "labels_unresolved", "note"):
        d[c] = d[c].str.replace("\r", "<CR>", regex=False)       # stacked labels contain CR; keep the table one line per record
    d.to_csv(OUTFILE, sep="\t", index=False)
    print(len(d), "records;", d.status.value_counts().to_dict(), " rules:", d.rules.value_counts().head(6).to_dict())


if __name__ == "__main__":
    main()
