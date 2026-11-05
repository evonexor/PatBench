"""Draw a stratified validation sample and build blinded evidence sheets from the patent XML.

Evidence comes from the patent archive itself, never from the extraction's own evidence record,
so an annotator can contradict the pipeline. Sheets carry no pipeline confidence, no mapping run
and no verdict hint, so the two annotation passes are blind to each other and to the pipeline.
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import html
import json
import re
import sys
import unicodedata
import warnings
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")

R = _ws + "/outputs/uspto-sar-recovery-20261006/results/"
V = _ws + "/outputs/target-uniprot-map-v2-20261019/"
RAW = Path(_ws + "/raw/patent-data/")
OUT = Path(sys.argv[1])
N_TARGET = 400
SEED = 20261019
MIN_PER_CELL = 8
rng = np.random.default_rng(SEED)

# ----------------------------------------------------------------- population
sar = ds.dataset(R + "SAR_merged.parquet").to_table(
    columns=["source", "source_record_id", "Patent", "patent_core", "Target", "Activity Type",
             "Relation", "Value", "Unit", "value_nM", "SMILES", "InChIKey"]).to_pandas()
sar["Target"] = sar.Target.fillna("").str.strip()
pt = ds.dataset(R + "patents.parquet").to_table(
    columns=["Patent", "archive", "source_format"]).to_pandas()
sar = sar.merge(pt, on="Patent", how="left")
tm = pq.read_table(V + "target_map.parquet", columns=[
    "patent_core", "Target", "uniprot_accession", "mapped_identity", "protein_name", "organism",
    "confidence", "target_class", "mappability", "source_run", "measured_in_cells"]).to_pandas()
tm = tm.drop_duplicates(["patent_core", "Target"])
sar = sar.merge(tm, on=["patent_core", "Target"], how="left")
sar["conf"] = sar.confidence.fillna("none")
print(f"population: {len(sar):,} formal SAR rows, {sar.archive.notna().sum():,} with an archive path")
sar = sar[sar.archive.notna()]

# evidence the recovery pass recorded; used only to locate the cell, never shown as truth
ev = ds.dataset(R + "new_SAR.parquet").to_table(columns=[
    "source_record_id", "activity_locator", "activity_header", "raw_value", "compound_alias",
    "structure_file", "source_format"]).to_pandas().set_index("source_record_id")
sar = sar.join(ev.drop(columns=["source_format"]), on="source_record_id")

# ----------------------------------------------------------------- strata
sar["endpoint_grp"] = np.where(sar["Activity Type"].isin(["IC50", "Ki", "Kd", "EC50"]),
                               sar["Activity Type"], "other")
sar["stratum"] = (sar.source.str.replace("USPTO_MOL_RECOVERY", "recovery").str.replace("USPTO_MOL", "baseline")
                  + " | " + sar.source_format.fillna("?")
                  + " | " + sar.conf)
counts = sar.stratum.value_counts()
print(f"\n{len(counts)} non-empty strata")
alloc = {}
remaining = N_TARGET
for s in counts.index:                                   # floor every stratum first
    take = min(MIN_PER_CELL, counts[s])
    alloc[s] = take
    remaining -= take
if remaining > 0:                                        # spread the rest proportionally
    weight = counts / counts.sum()
    for s in counts.index:
        extra = int(round(remaining * weight[s]))
        alloc[s] = min(alloc[s] + extra, counts[s])
rows = []
for s, n in alloc.items():
    pool = sar[sar.stratum == s]
    rows.append(pool.sample(n=min(n, len(pool)), random_state=SEED))
sample = pd.concat(rows).sample(frac=1.0, random_state=SEED).reset_index(drop=True)
sample["item_id"] = [f"V{i + 1:04d}" for i in range(len(sample))]
print(f"sampled {len(sample)} records")
print(sample.stratum.value_counts().to_string())


# ----------------------------------------------------------------- xml helpers
SGML_ENTITIES = {
    "mgr": "\u00b5", "agr": "\u03b1", "bgr": "\u03b2", "ggr": "\u03b3", "dgr": "\u03b4",
    "egr": "\u03b5", "zgr": "\u03b6", "eegr": "\u03b7", "thgr": "\u03b8", "igr": "\u03b9",
    "kgr": "\u03ba", "lgr": "\u03bb", "ngr": "\u03bd", "pgr": "\u03c0", "rgr": "\u03c1",
    "sgr": "\u03c3", "tgr": "\u03c4", "phgr": "\u03c6", "khgr": "\u03c7", "psgr": "\u03c8",
    "ohgr": "\u03c9", "OHgr": "\u03a9", "Dgr": "\u0394", "Sgr": "\u03a3",
    "plusmn": "\u00b1", "ap": "\u2248", "ges": "\u2265", "les": "\u2264",
    "deg": "\u00b0", "times": "\u00d7", "minus": "\u2212", "prime": "\u2032",
    "ndash": "\u2013", "mdash": "\u2014", "cir": "\u25cb", "squf": "\u25aa",
}


def sgml_entities(text):
    return re.sub(r"&([A-Za-z][A-Za-z0-9]*);",
                  lambda m: SGML_ENTITIES.get(m.group(1), m.group(0)), text)


def clean(text):
    text = re.sub(r"</?(?:PTEXT|PDAT|HIL|SB|SP|ITALIC|BOLD)\b[^>]*>", "", text, flags=re.I)
    text = re.sub(r"<sub>(.*?)</sub>", r"\1", text, flags=re.S)
    text = re.sub(r"<sup>(.*?)</sup>", r"^\1", text, flags=re.S)
    text = re.sub(r"<[^>]+>", " ", text)
    text = sgml_entities(html.unescape(text))
    text = unicodedata.normalize("NFKC", text)
    return re.sub(r"\s+", " ", text).strip()


def render_table(block, alias="", forms=(), max_rows=40):
    """Render a <tables> block as a text grid.

    Long activity tables run to hundreds of rows, so showing the first N would often hide the
    very row under validation. Header rows are always kept, and the rest is centred on the rows
    that mention the compound label or the claimed value, with the hit marked.
    """
    out = []
    for gi, tg in enumerate(re.findall(r"<tgroup\b.*?</tgroup>", block, re.S | re.I), start=1):
        is_head = bool(re.search(r"<thead\b", tg, re.I))
        rows = []
        for ri, row in enumerate(re.findall(r"<row\b.*?</row>", tg, re.S | re.I), start=1):
            cells = [clean(m.group(1) or "") for m in re.finditer(
                r"<entry\b[^>]*>(.*?)</entry>|<entry\b[^>]*/>", row, re.S | re.I)]
            if any(cells):
                rows.append((ri, [c if c else "-" for c in cells]))
        if not rows:
            continue
        # Marking has to be strict. A loose match flags "Ex. no. 6" as the row for a claimed
        # value of 6, which would point an annotator at the wrong line.
        alias_s = str(alias).strip() if alias and str(alias) != "nan" else ""
        strong, medium = set(), set()
        for idx, (ri, cells) in enumerate(rows):
            stripped = [c.strip() for c in cells]
            has_alias = bool(alias_s) and alias_s in stripped
            has_value = any(f and f in stripped for f in forms)
            if has_alias and has_value:
                strong.add(idx)
            elif has_alias:
                medium.add(idx)
            elif has_value and not alias_s:
                medium.add(idx)
        hits = strong or medium
        if is_head or len(rows) <= max_rows:
            keep = set(range(len(rows)))
        elif hits:
            keep = set()
            for h in hits:
                keep |= set(range(max(0, h - 4), min(len(rows), h + 5)))
            keep |= set(range(min(6, len(rows))))          # keep the top of the table for context
        else:
            keep = set(range(min(max_rows, len(rows))))
        out.append(f"  tgroup-{gi}:" + ("  (header)" if is_head else ""))
        last = -2
        for idx, (ri, cells) in enumerate(rows):
            if idx not in keep:
                continue
            if idx != last + 1 and last >= 0:
                out.append("    ...")
            mark = ("  <<< claimed compound and value" if idx in strong else
                    ("  <<< claimed compound" if idx in medium and alias_s else
                     ("  <<< claimed value" if idx in medium else "")))
            out.append(f"    row-{ri}: " + " | ".join(cells) + mark)
            last = idx
        if max(keep) < len(rows) - 1:
            out.append(f"    ... ({len(rows) - 1 - max(keep)} further rows not shown)")
    return "\n".join(out)


def value_strings(value, raw):
    """Plausible printed forms of the stored value, for locating it in the patent."""
    forms = set()
    if isinstance(raw, str) and raw.strip():
        forms.add(raw.strip())
    if pd.notna(value):
        for f in (f"{value:g}", f"{value:.1f}", f"{value:.2f}", f"{value:.3f}",
                  f"{value:.0f}", str(value)):
            forms.add(f.rstrip("0").rstrip(".") if "." in f else f)
            forms.add(f)
    return {f for f in forms if f}


ARCHIVE_CACHE = {}


def read_patent(archive):
    if archive in ARCHIVE_CACHE:
        return ARCHIVE_CACHE[archive]
    path = RAW / archive
    xml, mols = None, {}
    try:
        with zipfile.ZipFile(path) as z:
            names = z.namelist()
            xn = [n for n in names if n.upper().endswith((".XML", ".SGM", ".SGML", ".TXT"))]
            if xn:
                xml = z.read(xn[0]).decode("utf-8", "replace")
            for n in names:
                if n.upper().endswith(".MOL"):
                    mols[n.split("/")[-1].upper()] = z.read(n).decode("utf-8", "replace")
    except Exception as exc:
        xml = f"<<could not open archive: {exc}>>"
    ARCHIVE_CACHE.clear()
    ARCHIVE_CACHE[archive] = (xml, mols)
    return xml, mols


def mol_inchikey(block):
    try:
        m = Chem.MolFromMolBlock(block)
        return (Chem.MolToInchiKey(m), Chem.MolToSmiles(m)) if m is not None else (None, None)
    except Exception:
        return (None, None)


# ----------------------------------------------------------------- build sheets
pack = OUT / "annotation_pack"
pack.mkdir(parents=True, exist_ok=True)
records = []
for _, r in sample.sort_values("archive").iterrows():
    xml, mols = read_patent(r.archive)
    title = ""
    if xml:
        for pat in (r"<invention-title[^>]*>(.*?)</invention-title>",
                    r"<B540>.*?<STEXT>(.*?)</STEXT>", r"<TTL>(.*?)</TTL>",
                    r"<B542>(.*?)</B542>", r"<title>(.*?)</title>"):
            t = re.search(pat, xml, re.S | re.I)
            if t and clean(t.group(1)):
                title = clean(t.group(1))
                break
    # locate candidate tables: the recorded locator first, then any table printing the value
    blocks = {}
    if xml:
        for pattern in (r'<tables\b[^>]*\bid="([^"]+)"[^>]*>.*?</tables>',
                        r'<TABLE-US\b[^>]*\bID="([^"]+)"[^>]*>.*?</TABLE-US>'):
            for m in re.finditer(pattern, xml, re.S | re.I):
                blocks.setdefault(m.group(1), m.group(0))
    want, forms = None, value_strings(r.Value, r.raw_value)
    if isinstance(r.activity_locator, str) and r.activity_locator:
        want = r.activity_locator.split("/")[0]

    def table_score(block):
        """Rank candidate tables. Substring matching alone picks absurd tables: '56.5' occurs
        inside '256.58' in an LCMS table, so a cell has to equal the value, and naming the
        endpoint or the assay counts for as much."""
        cells = {c.strip() for c in re.findall(
            r"<entry\b[^>]*>(.*?)</entry>", block, re.S | re.I)}
        cells = {clean(c) for c in cells}
        text = clean(block).lower()
        score = 0
        if any(f and f in cells for f in forms):
            score += 3
        endpoint = str(r["Activity Type"]).lower()
        if endpoint and endpoint in text.replace(" ", ""):
            score += 2
        tgt = str(r.Target).strip().lower()
        if len(tgt) >= 3 and tgt[:24] in text:
            score += 2
        if any(w in text for w in ("ic50", "ic 50", "ki", "kd", "ec50", "activity", "inhibition")):
            score += 1
        return score

    if want in blocks:
        chosen = [want]
    else:
        ranked = sorted(((table_score(b), tid) for tid, b in blocks.items()), reverse=True)
        chosen = [tid for sc, tid in ranked[:2] if sc >= 3]
    tables_text = "\n\n".join(f"TABLE {t}\n{render_table(blocks[t], r.compound_alias, forms)}" for t in chosen)
    if not tables_text:
        # no table block printed the value: show where the number appears in the running text,
        # otherwise the item is unjudgeable and the annotator can only say "uncertain"
        plain_all = clean(xml) if xml else ""
        windows = []
        for f in sorted(forms, key=len, reverse=True):
            for m in re.finditer(re.escape(f), plain_all):
                windows.append("..." + plain_all[max(0, m.start() - 320):m.start() + 320] + "...")
                if len(windows) >= 3:
                    break
            if windows:
                break
        tables_text = ("<<no table block printed this value. Passages where the number appears "
                       "in the patent text:>>\n\n" + "\n\n".join(windows)) if windows else \
                      "<<neither a table nor the running text prints this value>>"
    # the structure as the patent drew it
    claim_ik = r.InChIKey if isinstance(r.InChIKey, str) else ""
    mol_key = (r.structure_file.split("/")[-1].upper()
               if isinstance(r.structure_file, str) and r.structure_file else "")
    mol_block = mols.get(mol_key, "")
    mol_ik, mol_smiles = mol_inchikey(mol_block) if mol_block else (None, None)
    # baseline rows carry no structure file. Searching every MOL in the archive for the claimed
    # structure is a stronger check anyway: it asks whether the patent drew this molecule at all.
    found_by_search = ""
    if not mol_ik and claim_ik:
        for name, block in mols.items():
            ik, sm = mol_inchikey(block)
            if ik == claim_ik:
                found_by_search, mol_ik, mol_smiles = name, ik, sm
                break
    # assay context: sentences naming the endpoint or the target text
    ctx = []
    if xml:
        plain = clean(xml)
        for kw in {str(r["Activity Type"]), str(r.Target)[:28]} - {"", "nan"}:
            for m in re.finditer(re.escape(kw), plain):
                ctx.append(plain[max(0, m.start() - 260):m.start() + 260])
                if len(ctx) >= 4:
                    break
    sheet = f"""# Validation item {r.item_id}

## Patent
- Number: {r.Patent}
- Title: {title or "(not found)"}
- Archive: {r.archive}

## Claim made by the dataset (this is what you are checking)
- Compound label in the patent: {r.compound_alias if isinstance(r.compound_alias, str) and r.compound_alias else "(not recorded)"}
- Structure (SMILES): {r.SMILES}
- Assay / target text: {r.Target!r}
- Endpoint: {r["Activity Type"]}
- Relation: {r.Relation}
- Value as printed: {r.raw_value if isinstance(r.raw_value, str) and r.raw_value else "(not recorded)"}
- Value stored: {r.Value} {r.Unit}  (= {r.value_nM:,.4g} nM)
- Protein assigned: {r.uniprot_accession or "(none)"} {("- " + str(r.mapped_identity)) if pd.notna(r.mapped_identity) else ""}
  {("Protein name: " + str(r.protein_name)) if pd.notna(r.protein_name) and str(r.protein_name) else ""}
  {("Organism: " + str(r.organism)) if pd.notna(r.organism) and str(r.organism) else ""}

## Evidence from the patent: activity table(s)
```
{tables_text}
```

## Evidence from the patent: structure file
- MOL file in the archive: {mol_key or found_by_search or "(none recorded)"}{" (located by searching every MOL in the archive for the claimed structure)" if found_by_search else ""}
- InChIKey recomputed from that MOL: {mol_ik or "(could not parse)"}
- InChIKey claimed by the dataset: {claim_ik}
- Match: {"YES" if mol_ik and claim_ik and mol_ik == claim_ik else ("NO" if mol_ik and claim_ik else "cannot tell")}
- SMILES recomputed from that MOL: {mol_smiles or "(n/a)"}

## Evidence from the patent: assay context
{chr(10).join("- ..." + c + "..." for c in ctx[:4]) or "- (no matching passage found)"}

## Your judgements

Answer each with one of: correct / incorrect / uncertain, then one short reason.

- **Q1 value**: does the patent report this endpoint and this value, with this unit, for this
  compound label? `Value` and `Unit` are meant to be the pair as printed in the patent, and the
  nM figure in brackets is the conversion of that pair. So check the printed pair against the
  table, and check that the conversion is right. A thousands separator (`12,000`) or scientific
  notation (`1.44E-10`) in the patent is fine and is not an error.
- **Q2 structure**: does the structure above correspond to that compound label in this patent?
- **Q3 target**: is the assigned protein the entity this assay measured? Answer `not applicable`
  if no protein is assigned, and `incorrect` if a protein is assigned but the assay is a cell,
  tissue or organism readout with no single protein target.
"""
    (pack / f"{r.item_id}.md").write_text(sheet)
    records.append({
        "item_id": r.item_id, "stratum": r.stratum, "source": r.source,
        "source_format": r.source_format, "confidence": r.conf, "target_class": r.target_class,
        "mappability": r.mappability, "mapping_run": r.source_run, "Patent": r.Patent,
        "patent_core": r.patent_core, "Target": r.Target, "endpoint": r["Activity Type"],
        "Relation": r.Relation, "Value": r.Value, "Unit": r.Unit, "value_nM": r.value_nM,
        "InChIKey": claim_ik, "SMILES": r.SMILES, "uniprot_accession": r.uniprot_accession,
        "mapped_identity": r.mapped_identity, "mol_file": mol_key or found_by_search, "mol_found_by_search": bool(found_by_search),
        "mol_inchikey": mol_ik,
        "mol_inchikey_matches": bool(mol_ik and claim_ik and mol_ik == claim_ik),
        "table_found": not tables_text.startswith("<<"),
        "has_recorded_locator": bool(isinstance(r.activity_locator, str) and r.activity_locator),
    })

frame = pd.DataFrame(records)
frame.to_parquet(OUT / "sample.parquet", index=False)
summary = {
    "built": "2026-10-19", "seed": SEED, "n_items": len(frame),
    "target_n": N_TARGET, "min_per_stratum": MIN_PER_CELL,
    "strata": {k: int(v) for k, v in frame.stratum.value_counts().items()},
    "by_source": {k: int(v) for k, v in frame.source.value_counts().items()},
    "by_source_format": {k: int(v) for k, v in frame.source_format.value_counts().items()},
    "by_confidence": {k: int(v) for k, v in frame.confidence.value_counts().items()},
    "by_endpoint": {k: int(v) for k, v in frame.endpoint.value_counts().items()},
    "automated_checks": {
        "activity_table_located": int(frame.table_found.sum()),
        "mol_file_present": int(frame.mol_file.ne("").sum()),
        "mol_located_by_structure_search": int(frame.mol_found_by_search.sum()),
        "mol_inchikey_recomputed": int(frame.mol_inchikey.notna().sum()),
        "mol_inchikey_matches_claim": int(frame.mol_inchikey_matches.sum()),
    },
}
(OUT / "sample_manifest.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False, default=str))
print("\n=== automated checks on the sample")
for k, v in summary["automated_checks"].items():
    print(f"  {k:34s} {v:>4,d} / {len(frame):,}")
print(f"\nwrote {OUT}/sample.parquet, sample_manifest.json and {len(frame)} sheets in annotation_pack/")
