#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Expand ChemDraw abbreviations ("aliases") in the MOL files of this directory.

Run:  PYTHONPATH=<directory with the dependencies> \
      python expand.py [--eval]

Input  (all relative to this script):  records.tsv, mols/<structure_file>
Output: expansions.tsv   (one row per record of records.tsv, same order)
        --eval additionally writes match_rates_by_label.tsv (self-check against
        public_keys_by_patent.tsv; the public keys are NEVER used to choose a structure).

Design (see README.md):
  * Alias labels are read from the raw MOL bytes.  A label may contain CR characters
    (multi-line ChemDraw labels such as "H<CR>N" = NH); alias_atoms.tsv truncates at the CR.
  * The placeholder carbon is converted IN PLACE to the root atom of the group and the rest
    of the group is appended, so the bonds (and thus wedge-derived chirality / E-Z flags) of
    the drawing are untouched.
  * Every group is defined by a closed set of rules (explicit tokens + the "caret grammar").
    Anything not covered by a rule is unresolved.  Valence is verified with RDKit for the
    group itself (dummy atoms stand for the neighbours, with the real bond orders) and for the
    whole molecule (sanitisation + no new radicals).
"""
import os
import re
import sys
import argparse
from collections import defaultdict, Counter
from multiprocessing import Pool

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.*')

HERE = os.path.dirname(os.path.abspath(__file__))
N_WORKERS = 16

# --------------------------------------------------------------------------------------
# 1. reading the MOL files
# --------------------------------------------------------------------------------------
ALIAS_RE = re.compile(r'^A\s+(\d+)\s*$')


def read_lines(path):
    """Lines of a MOL file.  Lines end with LF or CRLF; a lone CR inside a line is data
    (multi-line ChemDraw label) and is kept.  latin-1 keeps every byte."""
    raw = open(path, 'rb').read().decode('latin-1')
    lines = raw.split('\n')
    return [l[:-1] if l.endswith('\r') else l for l in lines]


def parse_molfile(path):
    """returns (mol, aliases{atom_no: raw_label}, lines); mol is unsanitised, Hs kept,
    alias lines removed, substance groups cleared."""
    lines = read_lines(path)
    aliases = {}
    clean = []
    skip = False
    for i, l in enumerate(lines):
        if skip:
            skip = False
            continue
        mt = ALIAS_RE.match(l)
        if mt and i + 1 < len(lines):
            aliases[int(mt.group(1))] = lines[i + 1]
            skip = True
            continue
        clean.append(l)
    mol = Chem.MolFromMolBlock('\n'.join(clean), sanitize=False, removeHs=False)
    if mol is not None:
        Chem.ClearMolSubstanceGroups(mol)
    return mol, aliases, lines


# --------------------------------------------------------------------------------------
# 2. label normalisation
# --------------------------------------------------------------------------------------
def decode_line(s):
    """Decode the ChemDraw text codes we understand.
       \\s<digits>\\n  -> digits (subscript run);  \\s- / \\s^- -> - / ^- (subscript charge)
       \\S<letters>\\n -> letters (superscript run, e.g. the t of tBu)"""
    s = re.sub(r'\\s([0-9]+)\\n', r'\1', s)
    s = re.sub(r'\\s([0-9]+)$', r'\1', s)
    s = re.sub(r'\\s(\^?[+-])(?:\\n)?', r'\1', s)
    s = re.sub(r'\\S([A-Za-z])\\n', r'\1', s)
    return s


def norm_label(raw):
    """Normalised key: codes decoded, blanks removed, lines (CR) kept as '\r'."""
    lines = [decode_line(x).replace(' ', '').replace('\t', '') for x in raw.split('\r')]
    lines = [x for x in lines if x != '']
    return '\r'.join(lines)


def show_label(raw):
    """label as written to the tsv (no tab / newline / semicolon)."""
    return raw.replace('\r', '<CR>').replace(';', '<SC>').replace('\t', ' ')


# --------------------------------------------------------------------------------------
# 3. the rules
# --------------------------------------------------------------------------------------
# Monovalent substituent tokens; SMILES start with the atom that carries the bond.
SUBS = {
    'Me': 'C', 'Et': 'CC', 'nPr': 'CCC', 'n-Pr': 'CCC', 'iPr': 'C(C)C', 'i-Pr': 'C(C)C',
    'Bu': 'CCCC', 'nBu': 'CCCC', 'n-Bu': 'CCCC', 'tBu': 'C(C)(C)C', 't-Bu': 'C(C)(C)C',
    'iBu': 'CC(C)C', 'i-Bu': 'CC(C)C', 'sBu': 'C(C)CC', 'Hex': 'CCCCCC',
    'Ph': 'c1ccccc1', 'Bn': 'Cc1ccccc1', 'PMB': 'Cc1ccc(OC)cc1',
    'Ac': 'C(C)=O', 'Boc': 'C(=O)OC(C)(C)C', 'Cbz': 'C(=O)OCc1ccccc1', 'Ms': 'S(C)(=O)=O',
    'F': 'F', 'Cl': 'Cl', 'Br': 'Br', 'I': 'I', 'D': '[2H]',
    # written "backwards" in front of a caret
    'F3C': 'C(F)(F)F', 'F2HC': 'C(F)F', 'D3C': 'C([2H])([2H])[2H]', 'H3C': 'C',
    'C2D5': 'C([2H])([2H])C([2H])([2H])[2H]',
    'MeO2C': 'C(=O)OC', 'BnO2C': 'C(=O)OCc1ccccc1', 'MeO2S': 'S(=O)(=O)C',
}
# H and O are handled specially (explicit hydrogens on the root / oxo on S,P)
TOKEN_NAMES = sorted(list(SUBS) + ['H', 'O'], key=len, reverse=True)
CORE_ELEMENTS = ('C', 'N', 'O', 'S', 'P')


def tokenize_prefix(p):
    """'Me2' -> [('Me',2)] ; 'BocH' -> [('Boc',1),('H',1)] ; '(PMB)2' ; '3F' ; None if not understood"""
    out = []
    i = 0
    while i < len(p):
        lead = re.match(r'\d+', p[i:])
        if lead:
            cnt_lead = int(lead.group(0))
            i += lead.end()
        else:
            cnt_lead = None
        if i < len(p) and p[i] == '(':
            j = p.find(')', i)
            if j < 0:
                return None
            name = p[i + 1:j]
            i = j + 1
            if name not in SUBS:
                return None
        else:
            name = None
            for t in TOKEN_NAMES:
                if p.startswith(t, i):
                    name = t
                    break
            if name is None:
                return None
            i += len(name)
        trail = re.match(r'\d+', p[i:])
        cnt_trail = None
        if trail:
            cnt_trail = int(trail.group(0))
            i += trail.end()
        if cnt_lead is not None and cnt_trail is not None:
            return None
        out.append((name, cnt_lead or cnt_trail or 1))
    return out or None


def caret_template(label):
    """<substituents>^<root element>: the atom behind the caret carries the bond(s) to the
    drawing, everything in front of it is attached to that atom.  Returns template or None."""
    if label.count('^') != 1:
        return None
    pre, core = label.split('^')
    if core not in CORE_ELEMENTS or not pre:
        return None
    toks = tokenize_prefix(pre)
    if toks is None:
        return None
    nH = 0
    branches = ''
    for name, cnt in toks:
        if name == 'H':
            nH += cnt
        elif name == 'O':
            if core not in ('S', 'P'):
                return None
            branches += '(=O)' * cnt
        else:
            branches += ('(' + SUBS[name] + ')') * cnt
    root = '[' + core + ('H' + (str(nH) if nH > 1 else '') if nH else '') + ']'
    return root + branches + '{D}'


# Equivalent spellings -> (canonical caret spelling, orientation).
# orientation: text order as drawn.  'rev' = written backwards (bond at the LAST atom, label sits
# to the left of its bond), 'fwd' = bond at the FIRST atom (label to the right of its bond),
# None = vertically stacked / not order-sensitive.
EQUIV = {
    'H\rN': ('H^N', None), 'N\rH': ('H^N', None),
    'O2\rS': ('O2^S', None), 'SO2': ('O2^S', None),
    'Boc\rN': ('Boc^N', None), 'Me\rN': ('Me^N', None), 'N\rMe': ('Me^N', None), 'NMe': ('Me^N', None),
    'OMe': ('Me^O', 'fwd'), 'MeO': ('Me^O', 'rev'), 'OCH3': ('Me^O', 'fwd'),
    'OEt': ('Et^O', 'fwd'), 'EtO': ('Et^O', 'rev'), 'H5C2O': ('Et^O', 'rev'),
    'OiPr': ('iPr^O', 'fwd'), 'OtBu': ('tBu^O', 'fwd'), 'tBuO': ('tBu^O', 'rev'),
    'OBn': ('Bn^O', 'fwd'), 'OCH2C6H5': ('Bn^O', 'fwd'), 'PhO': ('Ph^O', 'rev'), 'OPh': ('Ph^O', 'fwd'),
    'SMe': ('Me^S', 'fwd'), 'MeS': ('Me^S', 'rev'),
    'NMe2': ('Me2^N', 'fwd'), 'Me2N': ('Me2^N', 'rev'), 'HNMe': ('MeH^N', 'fwd'),
    'NMs2': ('Ms2^N', 'fwd'), '(MeO2S)2N': ('(MeO2S)2^N', 'rev'),
    'NHHex': ('HexH^N', 'fwd'),
}

# Groups that are a single formula; template with {D}; orientation as above
DIRECT = {
    'Me': '[CH3]{D}', 'Et': '[CH2]{D}C', 'Ph': '[c]{D}1ccccc1', 'tBu': '[C]{D}(C)(C)C',
    'iPr': '[CH]{D}(C)C', 'i-Pr': '[CH]{D}(C)C', 'nPr': '[CH2]{D}CC', 'n-Pr': '[CH2]{D}CC',
    'Bu': '[CH2]{D}CCC', 'nBu': '[CH2]{D}CCC', 'n-Bu': '[CH2]{D}CCC', 'n-butyl': '[CH2]{D}CCC',
    'nC4H9': '[CH2]{D}CCC', 'n-C4H9': '[CH2]{D}CCC', 'sBu': '[CH]{D}(C)CC', 'octyl': '[CH2]{D}CCCCCCC',
    'CF3': '[C]{D}(F)(F)F', 'F3C': '[C]{D}(F)(F)F', 'C2F5': '[C]{D}(F)(F)C(F)(F)F',
    'CH2CH2C6H5': '[CH2]{D}[CH2]c1ccccc1',
    'CO2Et': '[C]{D}(=O)OCC', 'COOEt': '[C]{D}(=O)OCC', 'CO2tBu': '[C]{D}(=O)OC(C)(C)C',
    'COOiPr': '[C]{D}(=O)OC(C)C', 'HO2C': '[C]{D}(=O)O',
    'NO2': '[N+]{D}(=O)[O-]', 'O2N': '[N+]{D}(=O)[O-]',
    'F': '[F]{D}', 'Bpin': '[B]{D}1OC(C)(C)C(C)(C)O1',
    'Sn(nBu)3': '[Sn]{D}(CCCC)(CCCC)CCCC', 'Sn(n-Bu)3': '[Sn]{D}(CCCC)(CCCC)CCCC',
}
DIRECT_ORIENT = {'CO2Et': 'fwd', 'COOEt': 'fwd', 'CO2tBu': 'fwd', 'COOiPr': 'fwd', 'HO2C': 'rev',
                 'NO2': 'fwd', 'O2N': 'rev'}
# groups whose special shape is not covered by the grammar
SPECIAL = {
    'O-': '[O-]{D}', 'O^-': '[O-]{D}',            # O(-) on an N(+): N-oxide
    'H^O|': '[OH]{D}',                              # '|' = ChemDraw artefact after O (see README)
    'O|': 'O{D}',                                   # fill valence; only O= or ring -O-
    'O\\sabs': 'O{D}', 'N\\sabs': 'N{D}',            # atom + stereo flag text "abs" (see README)
    'NH': '[NH]{D}',                                # plain NH, two single bonds
    'F2^C': '[CH](F)(F){D}',                        # CHF2 (the H is not written in the label)
}
SPECIAL_ALLOWED = {'NH': {(1, 1)}, 'O|': {(2,), (1, 1)}, 'O\\sabs': {(1, 1)}, 'N\\sabs': {(1, 2)}}
DIRECT_ALLOWED = {(1,)}                            # bond order pattern of the placeholder (one single bond)

# quaternary onium groups: explicit charge in the label
_ONIUM = {
    'N(CH3)3': '[N+]{D}(C)(C)C', 'N(CH2CH3)3': '[N+]{D}(CC)(CC)CC',
    'N((CH2)3CH3)3': '[N+]{D}(CCCC)(CCCC)CCCC', 'N(C6H5)3': '[N+]{D}(c1ccccc1)(c1ccccc1)c1ccccc1',
    'P(CH2CH3)3': '[P+]{D}(CC)(CC)CC', 'P(C6H5)3': '[P+]{D}(c1ccccc1)(c1ccccc1)c1ccccc1',
    'S(CH2CH3)2': '[S+]{D}(CC)CC',
}
for _k, _v in _ONIUM.items():
    for _form in ('+' + _k, _k + '+', _k + '\r+', '+\r' + _k):
        SPECIAL[_form] = _v

_NALK = re.compile(r'^n-?C(\d+)H(\d+)$')

OVERRIDE = {}          # used only by the hypothesis tests of --eval


def nalkyl_template(label):
    m = _NALK.match(label)
    if not m:
        return None
    n, h = int(m.group(1)), int(m.group(2))
    if h != 2 * n + 1 or n < 3 or n > 30:
        return None
    return '[CH2]{D}' + 'C' * (n - 1)


def resolve_template(label):
    """normalised label -> (template, allowed_bond_patterns or None, orientation) or None"""
    if label in OVERRIDE:
        return OVERRIDE[label], None, None
    if label in SPECIAL:
        return SPECIAL[label], SPECIAL_ALLOWED.get(label), None
    if label in EQUIV:
        canon, orient = EQUIV[label]
        t = caret_template(canon)
        if t:
            return t, None, orient
    t = caret_template(label)
    if t:
        return t, None, 'rev'
    if label in DIRECT:
        return DIRECT[label], DIRECT_ALLOWED, DIRECT_ORIENT.get(label)
    t = nalkyl_template(label)
    if t:
        return t, DIRECT_ALLOWED, None
    return None


# --------------------------------------------------------------------------------------
# 4. building the group and merging it into the drawing
# --------------------------------------------------------------------------------------
BOND_SYM = {1: '', 2: '=', 3: '#'}


def build_group(template, bond_orders):
    """RDKit mol of the group with one dummy atom per bond of the placeholder (real bond
    orders).  Atom 0 is the root.  Returns (mol, reason)."""
    dummies = ''.join('(%s*)' % BOND_SYM[b] for b in bond_orders)
    smi = template.replace('{D}', dummies)
    g = Chem.MolFromSmiles(smi)
    if g is None:
        return None, 'valence: group cannot bind with these bonds'
    for a in g.GetAtoms():
        if a.GetAtomicNum() != 0 and a.GetNumRadicalElectrons() > 0:
            return None, 'valence: group would be a radical with these bonds'
    Chem.Kekulize(g, clearAromaticFlags=True)
    return g, ''


def merge_group(rw, idx, g):
    """Convert placeholder atom idx of rw in place into the root of g and append the rest."""
    mapping = {0: idx}
    root = rw.GetAtomWithIdx(idx)
    ga = g.GetAtomWithIdx(0)
    root.SetAtomicNum(ga.GetAtomicNum())
    root.SetFormalCharge(ga.GetFormalCharge())
    root.SetIsotope(ga.GetIsotope())
    root.SetNoImplicit(True)
    root.SetNumExplicitHs(ga.GetTotalNumHs())
    root.SetNumRadicalElectrons(0)
    root.SetIsAromatic(False)
    root.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    for a in g.GetAtoms():
        i = a.GetIdx()
        if i == 0 or a.GetAtomicNum() == 0:
            continue
        na = Chem.Atom(a.GetAtomicNum())
        na.SetFormalCharge(a.GetFormalCharge())
        na.SetIsotope(a.GetIsotope())
        na.SetNoImplicit(True)
        na.SetNumExplicitHs(a.GetTotalNumHs())
        mapping[i] = rw.AddAtom(na)
    for b in g.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in mapping and j in mapping:
            rw.AddBond(mapping[i], mapping[j], b.GetBondType())


def group_smiles(g):
    return Chem.MolToSmiles(g)


def sup_defects(lines, mol):
    """Defects of superatom (SUP) groups that the exporter already expanded.  Not alias atoms,
    so they are NOT repaired -- only reported in the note."""
    sg = {}
    for l in lines:
        if l.startswith('M  STY'):
            p = l.split()
            for i in range(int(p[2])):
                sg.setdefault(int(p[3 + 2 * i]), {})['type'] = p[4 + 2 * i]
        elif l.startswith('M  SAL'):
            p = l.split()
            sg.setdefault(int(p[2]), {})['atoms'] = [int(x) - 1 for x in p[4:4 + int(p[3])]]
        elif l.startswith('M  SMT'):
            p = l.split(None, 3)
            sg.setdefault(int(p[2]), {})['label'] = p[3] if len(p) > 3 else ''
    out = []
    n = mol.GetNumAtoms()
    for g in sg.values():
        if g.get('type') != 'SUP' or not g.get('atoms') or any(a >= n for a in g['atoms']):
            continue
        label = norm_label(g.get('label', ''))
        aset = set(g['atoms'])
        seen, stack = set(), [g['atoms'][0]]
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            stack.extend(nb.GetIdx() for nb in mol.GetAtomWithIdx(x).GetNeighbors()
                         if nb.GetIdx() in aset and nb.GetIdx() not in seen)
        cross = [b for b in mol.GetBonds() if (b.GetBeginAtomIdx() in aset) != (b.GetEndAtomIdx() in aset)]
        if seen != aset and cross and '^' in label:
            out.append('SUP "%s" is broken in the file (atoms not connected: spurious fragment, group missing)' % label)
        elif label == 'N^C' and len(aset) == 2:
            for b in cross:
                inner = b.GetBeginAtom() if b.GetBeginAtomIdx() in aset else b.GetEndAtom()
                if inner.GetSymbol() == 'N':
                    out.append('SUP "N^C" is attached through N in the file (isocyanide R-N#C instead of nitrile)')
        elif len(aset) == 1 and not cross and re.match(r'^\\[sS]?\(?[RS]\^?\)?$', g.get('label', '')):
            out.append('SUP "%s" stereo text became a stray atom' % label)
    return sorted(set(out))


# --------------------------------------------------------------------------------------
# 5. deg-0 text labels (counter-ions, solvents, annotations)
# --------------------------------------------------------------------------------------
_SALT = re.compile(r'^[\x95·]*(\d+(\.\d+)?|½|1/2)?[\x95·]*(HCl|HBr|TFA|MsOH|HCOOH|NH4|NH3|H2O|citricacid|HO2C|CF3CO2-|'
                   r'CO2|AcOH|oxalicacid|fumaricacid|maleicacid|H2SO4|HOAc|Na|K|Li)', re.I)


def classify_text(label):
    n = label.replace(' ', '')
    return 'counterion/solvent' if _SALT.match(n) else 'annotation text'


# --------------------------------------------------------------------------------------
# 6. one record
# --------------------------------------------------------------------------------------
def heavy_atoms(m):
    return sum(1 for a in m.GetAtoms() if a.GetAtomicNum() > 1)


def inchikey_of(m):
    try:
        return Chem.MolToInchiKey(m)
    except Exception:
        return ''


def process(args):
    """never raises: an internal error makes the record unresolved"""
    try:
        return _process(args)
    except Exception as e:                                   # pragma: no cover
        return dict(patent_core=args[0], structure_file=args[1], status='unresolved', smiles='', inchikey='',
                    labels_expanded='', labels_unresolved='', note='internal error: %s' % str(e)[:80])


def _process(args):
    patent_core, sfile = args
    out = dict(patent_core=patent_core, structure_file=sfile, status='unresolved', smiles='',
               inchikey='', labels_expanded='', labels_unresolved='', note='')
    notes = []
    try:
        mol, aliases, lines = parse_molfile(os.path.join(HERE, 'mols', sfile))
    except Exception as e:                                   # pragma: no cover
        out['note'] = 'read error: %s' % e
        return out
    if mol is None:
        out['note'] = 'RDKit cannot read the MOL file'
        return out

    expanded, unresolved, plan, texts = [], [], [], []
    for atom_no in sorted(aliases):
        raw = aliases[atom_no]
        tag = '%d:%s' % (atom_no, show_label(raw))
        if atom_no < 1 or atom_no > mol.GetNumAtoms():
            unresolved.append(tag)
            notes.append('%s: alias refers to a missing atom' % tag)
            continue
        at = mol.GetAtomWithIdx(atom_no - 1)
        if at.GetAtomicNum() != 6:
            unresolved.append(tag)
            notes.append('%s: alias atom is not carbon' % tag)
            continue
        deg = at.GetDegree()
        if deg == 0:
            texts.append((atom_no, raw))
            continue
        if any(b.GetBondType() not in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE)
               for b in at.GetBonds()):
            unresolved.append(tag)
            notes.append('%s: placeholder has a non-single/double/triple bond' % tag)
            continue
        bos = tuple(sorted(int(round(b.GetBondTypeAsDouble())) for b in at.GetBonds()))
        label = norm_label(raw)
        res = resolve_template(label)
        if res is None:
            unresolved.append(tag)
            continue
        template, allowed, orient = res
        if allowed is not None and bos not in allowed:
            unresolved.append(tag)
            notes.append('%s: bond pattern %s not allowed for this group' % (tag, bos))
            continue
        if orient is not None and deg == 1:
            # the text of an order-sensitive label must run away from its bond
            nb = at.GetBonds()[0].GetOtherAtom(at).GetIdx()
            conf = mol.GetConformer()
            dx = conf.GetAtomPosition(nb).x - conf.GetAtomPosition(atom_no - 1).x
            dy = conf.GetAtomPosition(nb).y - conf.GetAtomPosition(atom_no - 1).y
            vertical = abs(dx) < 0.3 * abs(dy)
            if not vertical and ((orient == 'rev' and dx < 0) or (orient == 'fwd' and dx > 0)):
                unresolved.append(tag)
                notes.append('%s: label written %s but its bond runs the other way' % (tag, orient))
                continue
        g, why = build_group(template, list(bos))
        if g is None:
            unresolved.append(tag)
            notes.append('%s: %s' % (tag, why))
            continue
        # trial: this group alone must sanitise and add no radical
        trial = Chem.RWMol(mol)
        merge_group(trial, atom_no - 1, g)
        try:
            Chem.SanitizeMol(trial)
        except Exception as e:
            unresolved.append(tag)
            notes.append('%s: sanitisation fails (%s)' % (tag, str(e).split('\n')[0][:60]))
            continue
        plan.append((atom_no, g))
        expanded.append('%s->%s' % (tag, group_smiles(g)))

    out['labels_expanded'] = ';'.join(expanded)
    out['labels_unresolved'] = ';'.join(unresolved)
    if texts:
        notes.append('deg-0 text dropped: ' + ', '.join(
            '%d:%s[%s]' % (a, show_label(r), classify_text(r)) for a, r in texts))
    if unresolved:
        out['note'] = '; '.join(notes)
        return out

    rw = Chem.RWMol(mol)
    for atom_no, g in plan:
        merge_group(rw, atom_no - 1, g)
    drop = sorted((a - 1 for a, _ in texts), reverse=True)
    rw.BeginBatchEdit()
    for i in drop:
        rw.RemoveAtom(i)
    rw.CommitBatchEdit()
    try:
        Chem.SanitizeMol(rw)
        m2 = Chem.RemoveHs(rw)
        Chem.AssignStereochemistry(m2, cleanIt=True, force=True)
        if any(a.GetNumRadicalElectrons() > 0 for a in m2.GetAtoms()):
            raise ValueError('radical atom after expansion')
    except Exception as e:
        out['note'] = '; '.join(notes + ['sanitisation of the expanded drawing failed: %s' % str(e).split('\n')[0][:80]])
        out['labels_unresolved'] = ';'.join('%d:%s' % (a, show_label(aliases[a])) for a, _ in plan)
        out['labels_expanded'] = ''
        return out
    if heavy_atoms(m2) == 0:
        out['note'] = '; '.join(notes + ['no structure left after removing text labels'])
        return out
    m2.RemoveAllConformers()
    out['smiles'] = Chem.MolToSmiles(m2)
    out['inchikey'] = inchikey_of(m2)
    out['status'] = 'expanded' if plan else 'counterion_only'
    try:
        notes.extend('FILE-DEFECT (not repaired): ' + d for d in sup_defects(lines, mol))
    except Exception:
        pass
    out['note'] = '; '.join(notes)
    return out


# --------------------------------------------------------------------------------------
# 7. self-check against public molecules (never used to choose structures)
# --------------------------------------------------------------------------------------
def keys14_of_smiles(smi):
    """set of inchikey14 of the whole drawing and of each fragment with >= 6 heavy atoms, plus largest"""
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return set()
    ks = set()
    ks.add(inchikey_of(m)[:14])
    frags = Chem.GetMolFrags(m, asMols=True)
    frags = sorted(frags, key=lambda f: -heavy_atoms(f))
    if frags:
        ks.add(inchikey_of(frags[0])[:14])
    for f in frags:
        if heavy_atoms(f) >= 6:
            ks.add(inchikey_of(f)[:14])
    return ks


def _eval_one(args):
    idx, smi = args
    return idx, keys14_of_smiles(smi) if smi else set()


def evaluate(df, rec):
    pk = pd.read_csv(os.path.join(HERE, 'public_keys_by_patent.tsv'), sep='\t', dtype=str)
    pub = defaultdict(set)
    for d, k in zip(pk.patent_digits, pk.inchikey14):
        pub[d].add(k)
    pub = dict(pub)                      # plain dict: 'd in pub' must mean "patent has public molecules"
    digits = df.patent_core.str.replace('US', '', regex=False).str.lstrip('0')
    df = df.assign(digits=digits)
    df['in_db'] = df.digits.isin(pub)
    with Pool(N_WORKERS) as p:
        res = p.map(_eval_one, list(zip(df.index, df.smiles)), chunksize=50)
    keys = dict(res)
    df['match'] = [bool(keys[i] & pub.get(d, set())) if s else False for i, d, s in zip(df.index, df.digits, df.smiles)]
    # as drawn (placeholder structure) for comparison
    drawn = rec.set_index('structure_file').SMILES_as_drawn
    with Pool(N_WORKERS) as p:
        res2 = p.map(_eval_one, list(zip(df.index, [drawn[s] for s in df.structure_file])), chunksize=50)
    keys2 = dict(res2)
    df['match_drawn'] = [bool(keys2[i] & pub.get(d, set())) for i, d in zip(df.index, df.digits)]
    rows = []
    ok = df[df.status != 'unresolved']
    # label sets per record
    lab_sets = []
    for s in ok.labels_expanded:
        labs = set()
        for e in s.split(';'):
            if e:
                labs.add(e.split('->')[0].split(':', 1)[1])
        lab_sets.append(labs)
    ok = ok.assign(labs=lab_sets)
    stat = defaultdict(lambda: Counter())
    for _, r in ok.iterrows():
        for lab in r.labs:
            c = stat[lab]
            c['records'] += 1
            c['in_db'] += int(r.in_db)
            c['match'] += int(r.match)
            c['match_drawn'] += int(r.match_drawn)
            if len(r.labs) == 1:
                c['solo_records'] += 1
                c['solo_in_db'] += int(r.in_db)
                c['solo_match'] += int(r.match)
                c['solo_match_drawn'] += int(r.match_drawn)
    for lab, c in sorted(stat.items(), key=lambda x: -x[1]['records']):
        rows.append(dict(label=lab, records=c['records'], in_db=c['in_db'], match=c['match'],
                         match_rate=round(c['match'] / c['in_db'], 3) if c['in_db'] else '',
                         drawn_match=c['match_drawn'],
                         solo_records=c['solo_records'], solo_in_db=c['solo_in_db'], solo_match=c['solo_match'],
                         solo_rate=round(c['solo_match'] / c['solo_in_db'], 3) if c['solo_in_db'] else '',
                         solo_drawn_match=c['solo_match_drawn']))
    pd.DataFrame(rows).to_csv(os.path.join(HERE, 'match_rates_by_label.tsv'), sep='\t', index=False)
    tot = ok[ok.in_db]
    print('records expanded/counterion_only:', len(ok), ' in patents with public molecules:', len(tot),
          ' matching:', int(tot.match.sum()), ' as-drawn matching:', int(tot.match_drawn.sum()))
    hypothesis_tests(df, pub)
    geometry_table(rec)


def _label_records(df, key):
    """(patent_core, structure_file) of the resolved records that carry the normalised label key"""
    recs = []
    for pc, sf, st, le in zip(df.patent_core, df.structure_file, df.status, df.labels_expanded):
        if st == 'unresolved' or not le:
            continue
        for e in le.split(';'):
            if e and norm_label(e.split('->')[0].split(':', 1)[1].replace('<CR>', '\r').replace('<SC>', ';')) == key:
                recs.append((pc, sf))
                break
    return recs


def hypothesis_tests(df, pub):
    """Alternative readings of a few labels, run on exactly the same records; match against the
    public molecules of the patent.  (General rules only: one reading per label for all records.)"""
    rows = []
    byfile = df.set_index('structure_file')

    def digits(pc):
        return pc.replace('US', '').lstrip('0')

    def score(smiles_by_file, recs):
        n = m = 0
        for pc, sf in recs:
            d = digits(pc)
            s = smiles_by_file.get(sf)
            if d in pub and s:
                n += 1
                m += bool(keys14_of_smiles(s) & pub[d])
        return n, m

    tests = [
        ('Me^O', 'chosen: O bonds, Me on O (OMe)', None),
        ('Me^O', 'alt: Me carbon bonds (CH2-OH)', '[CH2]{D}O'),
        ('Et^O', 'chosen: O bonds (OEt)', None),
        ('Et^O', 'alt: Et carbon bonds (CH2CH2-OH)', '[CH2]{D}CO'),
        ('H\rN', 'chosen: NH', None),
        ('H\rN', 'alt: O', 'O{D}'),
        ('H\rN', 'alt: S', 'S{D}'),
        ('H\rN', 'alt: CH2 (= placeholder kept)', '[CH2]{D}'),
        ('O\\sabs', 'chosen: O', None),
        ('O\\sabs', 'alt: S', 'S{D}'),
        ('O\\sabs', 'alt: CH2 (= placeholder kept)', '[CH2]{D}'),
        ('O|', 'chosen: O', None),
        ('F2^C', 'chosen: CHF2', None),
    ]
    for key, desc, tmpl in tests:
        recs = _label_records(df, key)
        if tmpl is None:
            res = {sf: byfile.loc[sf, 'smiles'] for _, sf in recs}
        else:
            OVERRIDE.clear()
            OVERRIDE[key] = tmpl
            with Pool(N_WORKERS) as p:
                out = p.map(process, recs, chunksize=20)
            OVERRIDE.clear()
            res = {o['structure_file']: o['smiles'] for o in out}
        n, m = score(res, recs)
        rows.append(dict(label=key.replace('\r', '<CR>'), reading=desc, records=len(recs), in_db=n, match=m,
                         rate=round(m / n, 3) if n else ''))
    # natural experiment: superatom N^C exported attached through N (R-N#C) vs the caret rule (R-C#N)
    iso = df[df.note.str.contains('isocyanide R-N#C') & (df.smiles != '')]
    from rdkit.Chem import AllChem
    rxn = AllChem.ReactionFromSmarts('[*:1][N+:2]#[C-:3]>>[*:1][C+0:3]#[N+0:2]')
    n = m_exp = m_fix = 0
    for pc, sf, smi in zip(iso.patent_core, iso.structure_file, iso.smiles):
        d = digits(pc)
        if d not in pub:
            continue
        mol = Chem.MolFromSmiles(smi)
        for _ in range(10):
            ps = rxn.RunReactants((mol,))
            if not ps:
                break
            mol = ps[0][0]
            Chem.SanitizeMol(mol)
        n += 1
        m_exp += bool(keys14_of_smiles(smi) & pub[d])
        m_fix += bool(keys14_of_smiles(Chem.MolToSmiles(mol)) & pub[d])
    rows.append(dict(label='SUP N^C (already expanded in file)', reading='as exported: R-N#C (bond at first atom)',
                     records=len(iso), in_db=n, match=m_exp, rate=round(m_exp / n, 3) if n else ''))
    rows.append(dict(label='SUP N^C (already expanded in file)', reading='caret rule: R-C#N (bond at atom after caret)',
                     records=len(iso), in_db=n, match=m_fix, rate=round(m_fix / n, 3) if n else ''))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(HERE, 'hypothesis_tests.tsv'), sep='\t', index=False)
    print(t.to_string(index=False))


def geometry_table(rec):
    """Which side of its placeholder does the bonded neighbour lie?  (label text runs away from its bond)"""
    want = ['F3^C', 'Me^O', 'F3C^O', 'MeO', 'OMe', 'Me2^N', 'Me2N', 'NMe2', 'O2N', 'NO2', 'EtO', 'OEt', 'CO2Et', 'OtBu',
            'CF3', 'Ph', 'tBu']
    res = defaultdict(lambda: [0, 0, 0])
    for sf in rec.structure_file:
        mol, aliases, _ = parse_molfile(os.path.join(HERE, 'mols', sf))
        if mol is None:
            continue
        conf = mol.GetConformer()
        for a, raw in aliases.items():
            lab = norm_label(raw)
            at = mol.GetAtomWithIdx(a - 1)
            if lab not in want or at.GetDegree() != 1:
                continue
            nb = at.GetNeighbors()[0].GetIdx()
            dx = conf.GetAtomPosition(nb).x - conf.GetAtomPosition(a - 1).x
            dy = conf.GetAtomPosition(nb).y - conf.GetAtomPosition(a - 1).y
            res[lab][2 if abs(dx) < 0.3 * abs(dy) else (0 if dx > 0 else 1)] += 1
    print('label: bond partner to the RIGHT of the label / LEFT / (roughly vertical)')
    for lab in want:
        if lab in res:
            print('  %-8s %5d %5d %5d' % (lab, *res[lab]))


# --------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--eval', action='store_true')
    ap.add_argument('--limit', type=int, default=0)
    args = ap.parse_args()
    rec = pd.read_csv(os.path.join(HERE, 'records.tsv'), sep='\t', dtype=str, keep_default_na=False)
    if args.limit:
        rec = rec.head(args.limit)
    jobs = list(zip(rec.patent_core, rec.structure_file))
    with Pool(N_WORKERS) as p:
        res = p.map(process, jobs, chunksize=40)
    cols = ['patent_core', 'structure_file', 'status', 'smiles', 'inchikey', 'labels_expanded',
            'labels_unresolved', 'note']
    df = pd.DataFrame(res, columns=cols)
    if not args.limit:
        df.to_csv(os.path.join(HERE, 'expansions.tsv'), sep='\t', index=False)
    print(df.status.value_counts().to_string())
    if args.eval:
        evaluate(df, rec)


if __name__ == '__main__':
    main()
