"""Assemble the adjudication queue for the browser page.

Each item carries what the two annotators saw (the claim and the patent evidence from its sheet),
both annotators' verdicts and reasons, and a drawing of the stored structure. Pipeline-side hints
the annotators were blind to (mapping confidence, stratum, the fragment flag) are left out, so the
adjudicator is blind to them too.
"""
import json
import re
import sys

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.Draw import rdMolDraw2D

RDLogger.DisableLog("rdApp.*")
W = sys.argv[1].rstrip("/")
OUT = sys.argv[2]

q = pd.read_csv(W + "/adjudication_queue.tsv", sep="\t", dtype=str).fillna("")


def sections(md):
    parts = {}
    for chunk in ("\n" + md).split("\n## ")[1:]:
        head, _, body = chunk.partition("\n")
        parts[head.strip()] = body.strip()
    return parts


def bullet_pairs(block):
    pairs = []
    for line in block.splitlines():
        m = re.match(r"\s*-\s*([^:]+):\s*(.*)$", line)
        if m:
            pairs.append([m.group(1).strip(), m.group(2).strip()])
        elif line.strip() and pairs:
            pairs[-1][1] += " " + line.strip()
    return [p for p in pairs if p[1] and p[1] not in ("-",)]


def draw(smiles):
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return ""
    d = rdMolDraw2D.MolDraw2DSVG(360, 260)
    o = d.drawOptions()
    o.clearBackground = False
    o.bondLineWidth = 1.6
    o.fixedBondLength = 26
    d.DrawMolecule(mol)
    d.FinishDrawing()
    svg = d.GetDrawingText()
    svg = re.sub(r"<\?xml[^>]*\?>\s*", "", svg)
    svg = re.sub(r"<!--.*?-->", "", svg, flags=re.S)
    return svg.strip()


def gp_url(patent):
    # Google Patents drops the zero padding after RE (USRE037556E1 -> USRE37556E1)
    p = re.sub(r"^USRE0+", "USRE", patent)
    return f"https://patents.google.com/patent/{p}/en"


items = []
for _, r in q.iterrows():
    md = open(f"{W}/annotation_pack/{r.item_id}.md", encoding="utf-8").read()
    s = sections(md)
    patent_block = bullet_pairs(s.get("Patent", ""))
    title = dict(patent_block).get("Title", "")
    claim = bullet_pairs(s.get("Claim made by the dataset (this is what you are checking)", ""))
    tables = s.get("Evidence from the patent: activity table(s)", "")
    tables = re.sub(r"^```\s*|\s*```$", "", tables.strip())
    struct = bullet_pairs(s.get("Evidence from the patent: structure file", ""))
    ctx = s.get("Evidence from the patent: assay context", "")
    ctx = [c.lstrip("- ").strip() for c in ctx.splitlines() if c.strip().startswith("-")]
    smiles = dict(claim).get("Structure (SMILES)", "")
    open_qs = [k for k in ("q1", "q2", "q3") if r[f"{k}_consensus"] == "DISAGREE"]
    items.append({
        "id": r.item_id,
        "patent": r.Patent,
        "patent_url": gp_url(r.Patent),
        "title": "" if title.startswith("(not found") else title,
        "archive": dict(patent_block).get("Archive", ""),
        "endpoint": r["Activity Type"],
        "source_format": r.source_format,
        "open": open_qs,
        "verdicts": {k: {"A": r[f"{k}_A"], "B": r[f"{k}_B"],
                         "rA": r[f"{k}_reason_A"], "rB": r[f"{k}_reason_B"],
                         "agreed": r[f"{k}_consensus"] != "DISAGREE"}
                     for k in ("q1", "q2", "q3")},
        "claim": claim,
        "smiles": smiles,
        "svg": draw(smiles),
        "tables": tables,
        "structure": struct,
        "context": ctx[:4],
    })

with open(OUT, "w", encoding="utf-8") as fh:
    json.dump({"built": "2026-10-20", "n_items": len(items),
               "n_open": sum(len(i["open"]) for i in items), "items": items},
              fh, ensure_ascii=False, separators=(",", ":"))
print(f"{len(items)} items, {sum(len(i['open']) for i in items)} open questions")
print("open by question:", pd.Series([k for i in items for k in i["open"]]).value_counts().to_dict())
print("items with a structure drawing:", sum(1 for i in items if i["svg"]))
print("items with no table evidence:", sum(1 for i in items if i["tables"].startswith("<<")))
import os
print(f"wrote {OUT} ({os.path.getsize(OUT) / 1024:.0f} KB)")
