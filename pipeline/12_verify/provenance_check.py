"""Trace a random sample of benchmark rows back to the original patent archive.

    python provenance_check.py <release_dir> <n_test_rows> <n_train_rows>

For every sampled measurement: (a) the stored row equals the upstream evidence table, (b) value_nM follows
from Value and Unit, (c) the cell of the original patent XML holds raw_value, (d) the original MOL file
gives the stored structure."""
import os, sys, zipfile, re, math, random
import numpy as np, pandas as pd, pyarrow.dataset as ds
from collections import defaultdict
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
import xml.etree.ElementTree as ET
def parse_xml(raw):
    raw = re.sub(rb"<!DOCTYPE[^>\[]*(\[.*?\])?\s*>", b"", raw, flags=re.S)
    raw = re.sub(rb"&(?!#|amp;|lt;|gt;|quot;|apos;)[A-Za-z0-9._-]+;", b" ", raw)
    return ET.fromstring(raw)
WS = (os.environ.get("PATBENCH_WORKSPACE") or sys.exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
BASE = WS + "/raw/patent-data/"; O = WS + "/outputs/"
D, N_TEST, N_OTHER = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
os.chdir(D)   # the checks read sar_all.parquet and tasks/ from the release directory
sar = pd.read_parquet("sar_all.parquet", columns=["source", "source_record_id", "Patent", "Activity Type", "raw_value", "Unit", "Value", "Relation", "value_nM", "pAct", "structure_file",
                                                  "compound_alias", "activity_locator", "InChIKey_as_drawn", "split_time"]).set_index("source_record_id")
tk = pd.concat([pd.read_parquet(f"tasks/{n}.parquet", columns=["group_id", "InChIKey", "pAct", "n_measurements", "source_record_id", "split_time"]).assign(task=n) for n in ["rank_ic50", "rank_ki_kd", "rank_highcontrast"]])
tk["rid"] = tk.source_record_id.str.split(" \\| ")
ex = tk.explode("rid")
first_ids = lambda d: d.rid.drop_duplicates()
test_ids = ex[ex.split_time == "test"].rid.drop_duplicates()
train_ids = ex[ex.split_time == "train"].rid.drop_duplicates()
rs = np.random.default_rng(11)
pick = list(rs.choice(test_ids.to_numpy(), N_TEST, replace=False)) + list(rs.choice(train_ids.to_numpy(), N_OTHER, replace=False))
S = sar.loc[pick].copy(); S["rid"] = S.index; S["part"] = ["test"] * N_TEST + ["train"] * N_OTHER
# upstream evidence
need = set(pick)
cols = ["source_record_id", "archive", "xml_file", "structure_file", "activity_locator", "raw_value", "Unit", "Value", "value_nM", "compound_alias", "Activity Type", "Patent"]
ev = []
for p in ["uspto-sar-full-20261004/results/SAR_with_evidence.parquet", "uspto-sar-recovery-20261006/results/new_SAR.parquet"]:
    ev.append(ds.dataset(O + p).to_table(columns=cols, filter=ds.field("source_record_id").isin(list(need))).to_pandas())
ev = pd.concat(ev).drop_duplicates("source_record_id").set_index("source_record_id")
print("sampled", len(S), "| found in upstream evidence:", int(S.rid.isin(ev.index).sum()))
UNIT = {"nM": 1.0, "uM": 1e3, "μM": 1e3, "µM": 1e3, "mM": 1e6, "M": 1e9, "pM": 1e-3, "fM": 1e-6}
def norm(s): return re.sub(r"\s+", " ", str(s).replace("−", "-").replace("−", "-")).strip()
def cell_text(root, loc):
    tab, tg, rw, col = loc.split("/")
    if re.fullmatch(r"table-\d+", tab):          # recovery source: the n-th table element of the document
        alltabs = [e for e in root.iter() if e.tag == "table"]
        n = int(tab.split("-")[1]); tabs = [alltabs[n - 1]] if 0 < n <= len(alltabs) else []
    else:
        tabs = [e for e in root.iter() if e.tag in ("tables", "table") and e.get("id") == tab]
    if not tabs: return None
    t = tabs[0]; tgs = list(t.iter("tgroup"))
    k, r, c = (int(x.split("-")[1]) for x in (tg, rw, col))
    if k > len(tgs): return None
    rows = list(tgs[k - 1].iter("row"))
    if r > len(rows): return None
    ents = rows[r - 1].findall("entry")
    if c > len(ents): return None
    return "".join(ents[c - 1].itertext())
res = defaultdict(list)
for arch, g in S.groupby(S.rid.map(ev.archive)):
    try:
        z = zipfile.ZipFile(BASE + arch)
    except Exception as e:
        for rid in g.rid: res["archive_missing"].append(rid)
        continue
    xml_cache = {}
    for rid, r in g.iterrows():
        e = ev.loc[rid]
        same = all(str(r[c]) == str(e[c]) or (pd.isna(r[c]) and pd.isna(e[c])) for c in ["Patent", "Activity Type", "raw_value", "Unit", "activity_locator", "compound_alias", "structure_file"]) \
            and abs(float(r.value_nM) - float(e.value_nM)) <= 1e-9 * max(1, abs(float(e.value_nM)))
        res["a_matches_upstream"].append(same)
        f = UNIT.get(r.Unit); v = float(r.Value)
        res["b_nM_from_value_unit"].append(f is not None and abs(v * f - float(r.value_nM)) <= 1e-6 * max(1.0, abs(v * f)) and abs(9 - math.log10(float(r.value_nM)) - float(r.pAct)) < 1e-9)
        xf = e.xml_file
        if xf not in xml_cache:
            raw = z.read(xf)
            try:
                xml_cache[xf] = parse_xml(raw)
            except Exception as err:
                xml_cache[xf] = None; res["xml_parse_error"].append((xf, str(err)[:60]))
        txt = cell_text(xml_cache[xf], e.activity_locator) if xml_cache[xf] is not None else None
        ok = txt is not None and (norm(txt) == norm(e.raw_value) or norm(e.raw_value) in norm(txt))
        res["c_cell_holds_raw_value"].append(ok)
        if not ok: res["c_fail_examples"].append((rid, e.activity_locator, e.raw_value, None if txt is None else norm(txt)[:60]))
        try:
            mol = Chem.MolFromMolBlock(z.read(e.structure_file).decode("utf-8", "replace"), sanitize=True, removeHs=True)
            k = Chem.MolToInchiKey(mol).split("-")[0] if mol is not None else None
        except Exception:
            k = None
        okm = k is not None and k == str(r.InChIKey_as_drawn).split("-")[0]
        res["d_mol_gives_structure"].append(okm)
        if not okm: res["d_fail_examples"].append((rid, e.structure_file, k, r.InChIKey_as_drawn))
for k in ["a_matches_upstream", "b_nM_from_value_unit", "c_cell_holds_raw_value", "d_mol_gives_structure"]:
    v = np.array(res[k]); print(f"{k:28s} {int(v.sum())}/{len(v)} = {v.mean()*100:.2f}%")
print("archive missing:", len(res["archive_missing"]), "| xml parse errors:", len(res["xml_parse_error"]), res["xml_parse_error"][:3])
print("cell fail examples:", res["c_fail_examples"][:8]); print("mol fail examples:", res["d_fail_examples"][:8])
# task level: pAct of single-measurement rows reproduces the evidence value
one = ex.merge(tk[["group_id", "InChIKey"]].assign(), on=["group_id", "InChIKey"], how="left")
t1 = tk[tk.n_measurements == 1].copy(); t1["rid"] = t1.source_record_id
t1 = t1[t1.rid.isin(S.rid)]
d = (t1.set_index("rid").pAct - S.loc[t1.rid, "pAct"]).abs()
print("task pAct vs sar_all pAct, single-measurement rows in sample:", len(t1), "max abs diff", float(d.max()) if len(d) else None)
