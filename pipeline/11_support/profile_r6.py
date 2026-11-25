"""Profile of a PatBench release: size, tasks, splits, tiers, molecular weight, and checks for remaining problems.

    python profile_r6.py <release dir> <out dir>

Writes <out dir>/profile.json (all numbers) and <out dir>/mw_hist.tsv (histogram bins).  Everything is computed from the files of the release.
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")
D = Path(sys.argv[1]).resolve()
OUT = Path(sys.argv[2]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
TASKS = ["rank_ic50", "rank_ki_kd", "rank_highcontrast", "rank_complex_mutant", "rank_ec50"]
R = {}


def q(x, ps=(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)):
    x = pd.Series(x).dropna()
    d = {"n": int(len(x)), "mean": round(float(x.mean()), 2), "min": round(float(x.min()), 2), "max": round(float(x.max()), 2)}
    for p in ps:
        d[f"p{int(p * 100)}"] = round(float(x.quantile(p)), 2)
    return d


# ---------------------------------------------------------------- whole release
sar = pq.read_table(D / "sar_all.parquet", columns=[
    "patent_core", "group_id", "InChIKey", "target_label", "uniprot_accession", "year", "pAct", "heavy_atoms", "in_public_data",
    "structure_alias_suspect", "structure_repair", "split_time", "test_tier", "Activity Type", "source", "SMILES"]).to_pandas()
R["sar_all"] = dict(rows=len(sar), patents=int(sar.patent_core.nunique()), molecules=int(sar.InChIKey.nunique()),
                    groups=int(sar.group_id.nunique()), accessions=int(sar.uniprot_accession.nunique()),
                    year_min=int(sar.year.min()), year_max=int(sar.year.max()),
                    rows_by_year=sar.year.value_counts().sort_index().to_dict(),
                    rows_by_activity=sar["Activity Type"].value_counts().head(10).to_dict(),
                    rows_by_source=sar.source.value_counts().to_dict(),
                    alias_suspect_rows=int(sar.structure_alias_suspect.sum()), repaired_rows=int(sar.structure_repair.fillna("").ne("").sum()))
grp = pd.read_parquet(D / "groups.parquet")
R["groups"] = dict(groups=len(grp), usable=int(grp.has_usable_rows.sum()), n_mol_quantiles=q(grp.n_molecules),
                   range_quantiles=q(grp.pAct_range))

# ---------------------------------------------------------------- per task
tasks = {t: pd.read_parquet(D / "tasks" / f"{t}.parquet") for t in TASKS}
reg = pd.read_parquet(D / "tasks" / "regression_affinity.parquet")
allsmi = set(reg.SMILES)
for t in tasks.values():
    allsmi |= set(t.SMILES)
cache = {}
for s in allsmi:
    m = Chem.MolFromSmiles(s)
    if m is None:
        cache[s] = None
        continue
    frags = Chem.GetMolFrags(m, asMols=True)
    big = max(frags, key=lambda f: f.GetNumHeavyAtoms())
    cache[s] = dict(mw=Descriptors.MolWt(big), mw_full=Descriptors.MolWt(m), nfrag=len(frags), ha=big.GetNumHeavyAtoms(),
                    logp=Descriptors.MolLogP(big), hbd=rdMolDescriptors.CalcNumHBD(big), rot=rdMolDescriptors.CalcNumRotatableBonds(big),
                    charge=sum(a.GetFormalCharge() for a in big.GetAtoms()),
                    iso=int(any(a.GetIsotope() for a in big.GetAtoms())),
                    odd=sorted({a.GetSymbol() for a in big.GetAtoms()} - {"C", "N", "O", "S", "F", "Cl", "Br", "I", "P", "H"}))
R["unparsable_smiles"] = int(sum(v is None for v in cache.values()))


def add(df):
    df = df.copy()
    for k in ("mw", "mw_full", "nfrag", "logp", "hbd", "rot", "charge", "iso"):
        df[k] = [cache[s][k] if cache[s] else np.nan for s in df.SMILES]
    df["odd"] = [",".join(cache[s]["odd"]) if cache[s] else "" for s in df.SMILES]
    return df


tabs = {}
rows = []
hist_rows = []
bins = list(range(0, 1500, 25))
for name, df in tasks.items():
    df = add(df)
    tabs[name] = df
    te = df[df.split_time == "test"]
    r = dict(task=name, rows=len(df), groups=int(df.group_id.nunique()), patents=int(df.patent_core.nunique()),
             molecules=int(df.InChIKey.nunique()), targets=int(df.uniprot_accession.nunique()),
             rows_by_split=df.split_time.value_counts().to_dict(),
             groups_by_split=df.drop_duplicates("group_id").split_time.value_counts().to_dict(),
             test_groups_by_tier=te.drop_duplicates("group_id").test_tier.value_counts().to_dict(),
             test_patents=int(te.patent_core.nunique()), test_targets=int(te.uniprot_accession.nunique()),
             test_molecules=int(te.InChIKey.nunique()),
             group_size=q(df.groupby("group_id").size()), pAct=q(df.pAct), pAct_group_range=q(df.groupby("group_id").pAct.agg(lambda x: x.max() - x.min())),
             mw=q(df.mw), mw_test=q(te.mw), mw_train=q(df[df.split_time == "train"].mw), mw_valid=q(df[df.split_time == "valid"].mw),
             ha=q(df.heavy_atoms), logp=q(df.logp), hbd=q(df.hbd), rot=q(df.rot),
             mw_gt_500=round(float((df.mw > 500).mean()), 3), mw_gt_700=round(float((df.mw > 700).mean()), 3), mw_lt_250=round(float((df.mw < 250).mean()), 3),
             rows_multifragment=int((df.nfrag > 1).sum()), rows_charged=int((df.charge != 0).sum()), rows_isotope=int(df.iso.sum()),
             rows_unusual_elements=int((df.odd != "").sum()), unusual_elements=Counter(x for s in df.odd for x in s.split(",") if x).most_common(8),
             in_public_data_rows=round(float(df.in_public_data.mean()), 3), in_public_test=round(float(te.in_public_data.mean()), 3),
             novel_test_rows=int((te.test_tier == "novel").sum()),
             dup_group_molecule=int(df.duplicated(["group_id", "InChIKey"]).sum()),
             top5_target_row_share=round(float(df.target_label.value_counts().head(5).sum() / len(df)), 3),
             top_targets=df.target_label.value_counts().head(8).to_dict(),
             top10_patent_row_share=round(float(df.patent_core.value_counts().head(10).sum() / len(df)), 3),
             max_patent_share_test=round(float(te.patent_core.value_counts().iloc[0] / len(te)), 3),
             groups_gt_300=int((df.groupby("group_id").size() > 300).sum()),
             groups_lt_10=int((df.groupby("group_id").size() < 10).sum()),
             at_floor=int(df.at_assay_floor.sum()) if "at_assay_floor" in df else None,
             at_ceiling=int(df.at_assay_ceiling.sum()) if "at_assay_ceiling" in df else None)
    # within-group MW structure
    gm = df.groupby("group_id").mw.agg(lambda x: x.max() - x.min())
    r["within_group_mw_range"] = q(gm)
    gs = df.groupby("group_id").mw.std()
    r["within_group_mw_std"] = q(gs)
    # molecule overlap between splits
    tr_set = set(df[df.split_time == "train"].InChIKey)
    r["test_molecules_also_in_train"] = int(te.InChIKey.isin(tr_set).sum())
    r["test_scaffold_in_train_rows"] = round(float(te.scaffold.isin(set(df[df.split_time == "train"].scaffold)).mean()), 3)
    # test molecules: how many groups of the task are the same target as a train group
    trt = set(df[df.split_time == "train"].uniprot_accession)
    r["test_groups_target_seen_in_train"] = round(float(te.drop_duplicates("group_id").uniprot_accession.isin(trt).mean()), 3)
    rows.append(r)
    h = np.histogram(df.mw.dropna(), bins=bins)[0]
    ht = np.histogram(te.mw.dropna(), bins=bins)[0]
    for b, a, c in zip(bins[:-1], h, ht):
        hist_rows.append((name, b, int(a), int(c)))
R["tasks"] = rows

# regression
reg = add(reg)
R["regression"] = dict(rows=len(reg), molecules=int(reg.InChIKey.nunique()), accessions=int(reg.uniprot_accession.nunique()),
                       mw=q(reg.mw), by_split=reg.split_time.value_counts().to_dict(), pAct=q(reg.pAct))

# unique molecules over all ranking tasks, and the union MW curve
um = pd.concat([tabs[t][["InChIKey", "mw", "split_time"]] for t in tabs]).drop_duplicates("InChIKey")
R["union_ranking_molecules"] = dict(n=len(um), mw=q(um.mw))
for b, a in zip(bins[:-1], np.histogram(um.mw.dropna(), bins=bins)[0]):
    hist_rows.append(("all_ranking_unique_molecules", b, int(a), 0))

# eval_all
ea = {}
for t in ("rank_ic50", "rank_ki_kd"):
    e = pd.read_parquet(D / "eval_all" / f"{t}.parquet")
    ea[t] = dict(rows=len(e), groups=int(e.group_id.nunique()), novel_groups=int(e[e.test_tier == "novel"].group_id.nunique()),
                 novel_rows=int((e.test_tier == "novel").sum()))
R["eval_all"] = ea

# reference: ChEMBL / BindingDB molecular weight (sample of rows)
ref = {}
S = _ws + "/raw/patent-data/standardized/"
for src in ("ChEMBL", "BindingDB"):
    try:
        t = pq.read_table(S + src + ".parquet", columns=["SMILES"]).to_pandas().dropna().drop_duplicates()
        t = t.sample(min(len(t), 60000), random_state=1)
        mws = []
        for s in t.SMILES:
            m = Chem.MolFromSmiles(s)
            if m is None:
                continue
            f = max(Chem.GetMolFrags(m, asMols=True), key=lambda x: x.GetNumHeavyAtoms())
            mws.append(Descriptors.MolWt(f))
        ref[src] = q(mws)
        for b, a in zip(bins[:-1], np.histogram(mws, bins=bins)[0]):
            hist_rows.append((f"ref_{src}_sample", b, int(a), 0))
    except Exception as e:
        ref[src] = str(e)[:80]
R["reference_mw"] = ref

pd.DataFrame(hist_rows, columns=["set", "bin_lo", "count_all", "count_test"]).to_csv(OUT / "mw_hist.tsv", sep="\t", index=False)
json.dump(R, open(OUT / "profile.json", "w"), indent=1, ensure_ascii=False, default=str)
print("done", OUT)
