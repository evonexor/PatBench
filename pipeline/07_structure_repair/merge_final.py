"""Turn the output of the alias expander (expand.py) into the repair table that make_patched_sar.py consumes.

    python merge_final.py <alias_dir> <expansions.tsv> <records.tsv> <out.tsv>

Per record of <records.tsv> (every MOL file that has an alias line):
  status  expanded    at least one alias was restored (group, NH), or text atoms (counter-ions, annotations) were removed, and nothing is left open
          unresolved  an alias could not be expanded, or the file itself is broken (an exporter superatom that cut a group in two)
  rules   group | NH | text | nitrile  (joined with ';')
The reversed nitrile (N(+)#C(-) bonded through N) is turned into C#N here, so one table holds the final structure of the record.
Records the expander left unresolved stay unresolved.
"""
import re
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
O = Path(sys.argv[1]).resolve()
EXP, REC, OUT = O / sys.argv[2], O / sys.argv[3], O / sys.argv[4]
NITRILE = Chem.MolFromSmarts("[#6][N+]#[C-]")
e = pd.read_csv(EXP, sep="\t", dtype=str, keep_default_na=False)
rec = pd.read_csv(REC, sep="\t", dtype=str, keep_default_na=False)
assert len(e) == len(rec) and (e.structure_file.values == rec.structure_file.values).all(), "expansions must be in the order of the records"
BROKEN = re.compile(r'FILE-DEFECT \(not repaired\): SUP "[^"]*" is broken')
rows = []
for r in e.itertuples(index=False):
    d = dict(patent_core=r.patent_core, structure_file=r.structure_file, status="unresolved", smiles="", inchikey="", rules="",
             labels_expanded=r.labels_expanded, labels_unresolved=r.labels_unresolved, note=r.note)
    if r.status == "unresolved":
        rows.append(d); continue
    if BROKEN.search(r.note):
        d["note"] = "file defect: an exporter superatom was cut in two, the group is missing; " + r.note[:200]
        rows.append(d); continue
    rules = set()
    for item in filter(None, r.labels_expanded.split(";")):
        label, frag = item.split("->", 1)
        label = label.split(":", 1)[1].replace("<CR>", "").replace(" ", "")
        rules.add("NH" if frag == "*N*" and label in ("HN", "NH") else "group")
    if re.search(r"counterion/solvent|annotation text", r.note):
        rules.add("text")
    m = Chem.MolFromSmiles(r.smiles)
    if m is None:
        d["note"] = "expanded SMILES does not parse"; rows.append(d); continue
    for _ in range(10):
        mt = m.GetSubstructMatch(NITRILE)
        if not mt:
            break
        _, n, c = mt
        m.GetAtomWithIdx(n).SetAtomicNum(6); m.GetAtomWithIdx(n).SetFormalCharge(0)
        m.GetAtomWithIdx(c).SetAtomicNum(7); m.GetAtomWithIdx(c).SetFormalCharge(0)
        Chem.SanitizeMol(m)
        rules.add("nitrile")
    d.update(status="expanded", smiles=Chem.MolToSmiles(m), inchikey=Chem.MolToInchiKey(m), rules=";".join(sorted(rules)))
    rows.append(d)
out = pd.DataFrame(rows)
for c in ("labels_expanded", "labels_unresolved", "note"):
    out[c] = out[c].str.replace("\r", "<CR>", regex=False)
out.to_csv(OUT, sep="\t", index=False)
print(len(out), "records;", out.status.value_counts().to_dict())
print("rules:", out[out.status == "expanded"].rules.value_counts().to_dict())
print("unresolved because of a file defect:", int(out.note.str.startswith("file defect").sum()))
