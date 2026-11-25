"""Evaluation-support tables for PatBench r4. Additive: reads the release, writes only eval_support/.

  target_domain_sequences.{parquet,tsv}  one row per group_id : the stretch of the group's
        target sequence that the recommended 3D structure actually models (first to last aligned
        residue, gaps inside kept), so sequence models and structure models can see the same region.
  test_group_profile.{parquet,tsv}       one row per (task, test group): target seen in train or not,
        mapping-check status, sequence-length / truncation flags, structure eligibility, domain span.

Nothing here changes a split, a label or a task file.

    python build_eval_support.py <release_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

R = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
OUT = R / "eval_support"
OUT.mkdir(exist_ok=True)
TASKS = ["rank_ic50", "rank_ki_kd", "rank_highcontrast", "rank_complex_mutant", "rank_ec50"]
log = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    log.append(s)


seq = pd.read_parquet(R / "target_sequences.parquet")
sm = pd.read_parquet(R / "structures/group_structure_map.parquet")
ts = pd.read_parquet(R / "structures/target_structures.parquet")
grp = pd.read_parquet(R / "groups.parquet")
assert seq.group_id.is_unique and sm.group_id.is_unique and len(seq) == len(sm)

# ------------------------------------------------------------------ domain sequences
seq_by_gid = seq.set_index("group_id")
rows = []
for r in sm.itertuples(index=False):
    s = seq_by_gid.loc[r.group_id]
    full = s.sequence
    base = dict(group_id=int(r.group_id), uniprot_accession=r.uniprot_accession, sequence_id=s.sequence_id,
                sequence_kind=s.sequence_kind, full_len=len(full), structure_id=r.recommended_structure_id or "",
                structure_source=r.recommended_source or "", is_holo=r.recommended_is_holo)
    if not r.recommended_structure_id:
        rows.append({**base, "status": "no_structure"})
        continue
    p = Path(r.recommended_clean_path)
    p = p if p.is_absolute() else R / "structures" / p
    rm = p.with_suffix("").as_posix() + ".residue_map.tsv"
    if not Path(rm).exists():
        rows.append({**base, "status": "residue_map_missing"})
        continue
    m = pd.read_csv(rm, sep="\t", dtype=str, keep_default_na=False)
    m = m[m.group_pos_main != ""].copy()
    m["pos"] = m.group_pos_main.astype(int)
    m = m[m.in_crop.str.lower() == "true"] if "in_crop" in m else m
    if m.empty:
        rows.append({**base, "status": "no_aligned_residues"})
        continue
    # each group position is counted once even when the file holds several chains / insertion codes
    m = m.drop_duplicates("pos")
    lo, hi = int(m.pos.min()), int(m.pos.max())
    assert 1 <= lo <= hi <= len(full), (r.group_id, lo, hi, len(full))
    mism = int(sum(1 for q, a in zip(m.pos, m.uniprot_aa) if a and full[q - 1] != a))
    st = "alphafold_full_length" if r.recommended_source == "AlphaFold" else "ok"
    rows.append({**base, "status": st, "domain_start": lo, "domain_end": hi, "domain_len": hi - lo + 1,
                 "domain_fraction": (hi - lo + 1) / len(full), "n_modelled": int(len(m)),
                 "n_gap_inside": int(hi - lo + 1 - len(m)), "n_aa_mismatch_vs_group_sequence": mism,
                 "domain_sequence": full[lo - 1:hi]})
dom = pd.DataFrame(rows)
for c in ["domain_start", "domain_end", "domain_len", "n_modelled", "n_gap_inside", "n_aa_mismatch_vs_group_sequence"]:
    dom[c] = dom[c].astype("Int64")
say("domain table:", len(dom), dom.status.value_counts().to_dict())
ok = dom[dom.domain_len.notna()]
say("domain_len median/q90/max:", ok.domain_len.median(), ok.domain_len.quantile(.9), ok.domain_len.max(),
    " fraction of full length median:", round(ok.domain_fraction.median(), 3))
say("groups with aa mismatch vs group sequence:", int((ok.n_aa_mismatch_vs_group_sequence > 0).sum()),
    " max mismatches:", int(ok.n_aa_mismatch_vs_group_sequence.max()))
dom.to_parquet(OUT / "target_domain_sequences.parquet", index=False)
dom.to_csv(OUT / "target_domain_sequences.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ test group profile
frames, train_targets = [], {}
for t in TASKS:
    d = pd.read_parquet(R / f"tasks/{t}.parquet", columns=[
        "group_id", "patent_core", "uniprot_accession", "target_label", "Activity Type", "split_time", "test_tier",
        "group_n_molecules", "group_pAct_range", "InChIKey"])
    tr = d[d.split_time == "train"].drop_duplicates("group_id")
    train_targets[t] = tr.groupby("uniprot_accession").size()
    te = d[d.split_time == "test"]
    g = te.groupby("group_id").agg(patent_core=("patent_core", "first"), uniprot_accession=("uniprot_accession", "first"),
                                   target_label=("target_label", "first"), endpoint=("Activity Type", "first"),
                                   test_tier=("test_tier", "first"), n_test_molecules=("InChIKey", "nunique"),
                                   pAct_range=("group_pAct_range", "first")).reset_index()
    g.insert(0, "task", t)
    frames.append(g)
prof = pd.concat(frames, ignore_index=True)
prof["n_train_groups_same_target_same_task"] = [
    int(train_targets[t].get(u, 0)) for t, u in zip(prof.task, prof.uniprot_accession)]
prof["target_seen_in_train_same_task"] = prof.n_train_groups_same_target_same_task > 0
gtrain = grp[grp.split_time == "train"].groupby("uniprot_accession").size()
prof["n_train_groups_same_target_global"] = [int(gtrain.get(u, 0)) for u in prof.uniprot_accession]

q = seq[["group_id", "confidence", "label_check", "xref_status", "seq_len"]]
prof = prof.merge(q, on="group_id", how="left")
prof["mapping_strict"] = (prof.confidence == "high") & (prof.label_check == "label_names_target") & (prof.xref_status != "differs")
prof["review_priority"] = np.select(
    [prof.xref_status == "differs", prof.label_check == "label_does_not_name_target", prof.confidence == "medium"],
    [1, 2, 3], default=4)
rev = pd.read_csv(R / "target_review_test.tsv", sep="\t", encoding="utf-8-sig", dtype=str, keep_default_na=False)
prof["in_review_table"] = prof.group_id.astype(str).isin(set(rev.group_id))
prof["seq_gt_545"] = prof.seq_len > 545
prof["seq_gt_1000"] = prof.seq_len > 1000
prof["seq_gt_1200"] = prof.seq_len > 1200

s2 = sm[["group_id", "recommended_structure_id", "recommended_source", "recommended_is_holo", "structure_matches_sequence",
         "structure_is_homolog_standin", "seq_coverage", "recommended_ligand_in_test_same_target",
         "recommended_ligand_in_ic50_kikd_test_same_target"]]
prof = prof.merge(s2, on="group_id", how="left")
prof["has_structure"] = prof.recommended_structure_id.fillna("") != ""
prof["structure_ligand_leak"] = prof.recommended_ligand_in_test_same_target.fillna(False).astype(bool)
prof["docking_eligible"] = (prof.has_structure & (prof.recommended_is_holo == True)
                            & (prof.structure_matches_sequence != False) & ~prof.structure_ligand_leak)
prof = prof.merge(dom[["group_id", "status", "domain_start", "domain_end", "domain_len", "domain_fraction"]]
                  .rename(columns={"status": "domain_status"}), on="group_id", how="left")
prof = prof.sort_values(["task", "group_id"]).reset_index(drop=True)
assert not prof.duplicated(["task", "group_id"]).any()
prof.to_parquet(OUT / "test_group_profile.parquet", index=False)
prof.to_csv(OUT / "test_group_profile.tsv", sep="\t", index=False)

say("\ntest_group_profile rows:", len(prof), prof.task.value_counts().to_dict())
for t in TASKS:
    p = prof[prof.task == t]
    say(f"  {t}: groups={len(p)} seen_in_train(same task)={int(p.target_seen_in_train_same_task.sum())} "
        f"strict_mapping={int(p.mapping_strict.sum())} docking_eligible={int(p.docking_eligible.sum())} "
        f"seq>1000={int(p.seq_gt_1000.sum())} not_in_review_table={int((~p.in_review_table).sum())}")
two = prof[prof.task.isin(["rank_ic50", "rank_ki_kd"])].drop_duplicates("group_id")
say("ic50+ki_kd distinct test groups:", len(two), " strict:", int(two.mapping_strict.sum()),
    " docking_eligible:", int(two.docking_eligible.sum()), " in review table:", int(two.in_review_table.sum()))
allg = prof.drop_duplicates("group_id")
say("distinct test groups over all five tasks:", len(allg), " not in review table:", int((~allg.in_review_table).sum()))
(OUT / "build.log").write_text("\n".join(log) + "\n", encoding="utf-8")
