"""Accession corrections for target-uniprot-map-v2 -> v2.1 (2026-10-31). The v2 directory is left as it is.

    python apply_target_overrides.py <v2_dir> <out_dir> <target_overrides.tsv>

Why: fetch_unresolved.py resolved free-text identities ("ORG:species|name") by a full-text UniProt search and
kept the longest hit, without checking that the returned protein is the one named. Three identities came back
as the wrong protein (NDM-1 -> GroEL, VIM-1 -> an aminoglycoside phosphotransferase, Mtb LeuRS -> IleRS).

target_overrides.tsv has one line per (old accession, mapped identity):
    old_accession  mapped_identity  new_accession  expected_rows  must_contain  evidence
`must_contain` is a word that has to occur in the new record's protein name or gene names, so a mistyped
accession cannot slip through. Every other row of target_map.parquet is copied unchanged.
"""
import io
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

V2, OUT, OVR = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)
DATE = "20261031"
FIELDS = "accession,id,protein_name,gene_names,organism_name,organism_id,length,reviewed"  # the search endpoint answers 500 when `sequence` is requested (2026-10-31); the FASTA endpoint is used for it


def get(url, tries=8):
    """The UniProt REST service answers 500 now and then; the same request succeeds a few seconds later."""
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "PatBench-target-fix/1.0 (research use)"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode("utf8")
        except urllib.error.HTTPError as e:
            if e.code < 500 or k == tries - 1:
                raise
            time.sleep(3 * (k + 1))


def uniprot(acc):
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
        {"query": f"accession:{acc}", "fields": FIELDS, "format": "tsv", "size": 2})
    df = pd.read_csv(io.StringIO(get(url)), sep="\t", dtype=str, keep_default_na=False)
    assert len(df) == 1 and df.Entry.iloc[0] == acc, (acc, len(df))
    row = df.iloc[0].copy()
    furl = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    lines = get(furl).strip().split("\n")
    assert lines[0].startswith(">") and f"|{acc}|" in lines[0], lines[0]
    row["Sequence"] = "".join(lines[1:])
    return row


t = pd.read_parquet(V2 / "target_map.parquet")
ovr = pd.read_csv(OVR, sep="\t", dtype=str, keep_default_na=False)
t["accession_override"] = ""
t["override_evidence"] = ""
done = []
for o in ovr.itertuples(index=False):
    mask = (t.uniprot_accession == o.old_accession) & (t.mapped_identity == o.mapped_identity)
    assert int(mask.sum()) == int(o.expected_rows), (o.old_accession, int(mask.sum()), o.expected_rows)
    # nothing else may still carry the old accession, or the correction would be half done
    assert not ((t.uniprot_accession == o.old_accession) & ~mask).any(), o.old_accession
    r = uniprot(o.new_accession)
    seq = r["Sequence"]
    assert len(seq) == int(r["Length"]), "sequence length differs from the UniProt length"
    text = (r["Protein names"] + " " + r["Gene Names"]).lower()
    assert o.must_contain.lower() in text, f"{o.new_accession}: '{o.must_contain}' not in '{text[:120]}'"
    # an accession already in the map keeps the values of the rows it already has (the map spells some protein
    # names short and some long, and leaves the entry name empty for some); the sequence must be identical to UniProt today
    have = t[(t.uniprot_accession == o.new_accession)]
    if len(have):
        h = have.iloc[0]
        assert h.sequence == seq, f"{o.new_accession}: sequence in the map differs from UniProt"
        assert h.uniprot_entry in ("", r["Entry Name"]), f"{o.new_accession}: entry name in the map is {h.uniprot_entry!r}, UniProt says {r['Entry Name']!r}"
        assert r["Protein names"].startswith(h.protein_name), (h.protein_name, r["Protein names"])
        new = {c: h[c] for c in ("uniprot_accession", "uniprot_entry", "protein_name", "organism", "tax_id", "seq_len", "sequence")}
    else:
        new = {"uniprot_accession": o.new_accession, "uniprot_entry": r["Entry Name"], "protein_name": r["Protein names"],
               "organism": r["Organism"], "tax_id": r["Organism (ID)"], "seq_len": len(seq), "sequence": seq}
    for c, v in new.items():
        t.loc[mask, c] = type(t[c].iloc[0])(v) if c in ("seq_len",) else v
    t.loc[mask, "organism_resolution"] = "accession_override"
    t.loc[mask, "method"] = t.loc[mask, "method"].astype(str) + "+accession_override_" + DATE
    t.loc[mask, "accession_override"] = f"{o.old_accession}->{o.new_accession}"
    t.loc[mask, "override_evidence"] = o.evidence
    for p in t[mask].patent_core.unique():
        done.append(dict(patent_core=p, Target=";".join(sorted(t[mask & (t.patent_core == p)].Target.unique())),
                         old_accession=o.old_accession, new_accession=o.new_accession,
                         new_entry=r["Entry Name"], new_protein_name=r["Protein names"], new_gene_names=r["Gene Names"],
                         new_organism=r["Organism"], new_length=len(seq), reviewed=r["Reviewed"]))
    print(f"{o.old_accession} -> {o.new_accession}: {int(mask.sum())} rows; {r['Entry Name']} | {r['Protein names'][:70]} | "
          f"{r['Gene Names'][:40]} | {r['Organism']} | {len(seq)} aa | {r['Reviewed']}")

t.to_parquet(OUT / "target_map.parquet", index=False)
pd.DataFrame(done).to_csv(OUT / "overrides_applied.tsv", sep="\t", index=False)
print("rows:", len(t), " changed:", int((t.accession_override != "").sum()))
