"""List every unresolved alias label of the r6 unresolved records with its context.

    python list_unresolved_labels.py <alias_dir> <unresolved_records.tsv> <out.tsv>

One row per distinct normalised label: records, patents, bond patterns of the placeholder atom, the engine's reason, one example record.
"""
import os
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

O = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(O / "indep_all"))
import expand as X  # noqa: E402

u = pd.read_csv(sys.argv[2], sep="\t", dtype=str, keep_default_na=False)
agg = defaultdict(lambda: dict(records=set(), patents=set(), bonds=defaultdict(int), reasons=defaultdict(int), example="", raw=""))
n_missing = 0
for r in u.itertuples(index=False):
    p = O / "indep_all" / "mols" / r.structure_file
    if not p.exists():
        n_missing += 1
        continue
    mol, aliases, _ = X.parse_molfile(str(p))
    if mol is None:
        continue
    for atom_no, raw in aliases.items():
        if atom_no < 1 or atom_no > mol.GetNumAtoms():
            continue
        at = mol.GetAtomWithIdx(atom_no - 1)
        if at.GetAtomicNum() != 6 or at.GetDegree() == 0:
            continue
        label = X.norm_label(raw)
        res = X.resolve_template(label)
        bos = tuple(sorted(int(round(b.GetBondTypeAsDouble())) for b in at.GetBonds()))
        if res is not None:
            template, allowed, orient = res
            if allowed is not None and bos not in allowed:
                why = "bond pattern not allowed"
            else:
                g, w = X.build_group(template, list(bos))
                if g is None:
                    why = w
                else:
                    continue            # this label resolves fine; the record is open for another label
        else:
            why = "no rule"
        a = agg[label]
        a["records"].add((r.patent_core, r.structure_file))
        a["patents"].add(r.patent_core)
        a["bonds"][bos] += 1
        a["reasons"][why] += 1
        a["raw"] = X.show_label(raw)
        if not a["example"]:
            a["example"] = r.structure_file
rows = []
for label, a in agg.items():
    rows.append(dict(label=X.show_label(label), records=len(a["records"]), patents=len(a["patents"]),
                     bonds=";".join(f"{k}x{v}" for k, v in sorted(a["bonds"].items(), key=lambda t: -t[1])),
                     reason=";".join(f"{k}x{v}" for k, v in a["reasons"].items()), example=a["example"]))
out = pd.DataFrame(rows).sort_values("records", ascending=False)
out.to_csv(sys.argv[3], sep="\t", index=False)
print("records without mol file:", n_missing, " distinct labels:", len(out), " label occurrences:", int(out.records.sum()))
