"""Agreement and accuracy from the two independent annotation passes.

The sample is stratified with a floor of 8 per cell, so small strata (legacy SGML, low-confidence
mappings) are deliberately over-represented. A plain proportion over the 400 items would therefore
be biased against the dataset. Both estimates are reported: the unweighted sample proportion, and
the population-weighted estimate with its stratified-sampling interval, which is the one to quote.
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import json
import sys

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import pyarrow.parquet as pq

pd.set_option("display.width", 220)
W = sys.argv[1].rstrip("/")
R = _ws + "/outputs/uspto-sar-recovery-20261006/results/"
V = _ws + "/outputs/target-uniprot-map-v2-20261019/"
QS = ["q1", "q2", "q3"]
VALID = {"correct", "incorrect", "uncertain", "not_applicable"}
Z = 1.959963985


def wilson(k, n, z=Z):
    """Wilson score interval; behaves sensibly at k = 0 and k = n, unlike the normal interval."""
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def cohen_kappa(a, b):
    cats = sorted(set(a) | set(b))
    idx = {c: i for i, c in enumerate(cats)}
    n = len(a)
    m = np.zeros((len(cats), len(cats)))
    for x, y in zip(a, b):
        m[idx[x], idx[y]] += 1
    po = np.trace(m) / n
    pe = float((m.sum(axis=0) / n) @ (m.sum(axis=1) / n))
    return (po - pe) / (1 - pe) if pe < 1 else np.nan, po


# --------------------------------------------------------------- load passes
sample = pq.read_table(W + "/sample.parquet").to_pandas()
passes = {}
for tag in ["A", "B"]:
    p = pd.read_csv(f"{W}/annotations/annotator_{tag}.tsv", sep="\t", dtype=str).fillna("")
    p.columns = [c.strip() for c in p.columns]
    bad = {q: sorted(set(p[q]) - VALID) for q in QS if q in p}
    print(f"annotator {tag}: {len(p)} rows; unexpected verdicts: "
          f"{ {k: v for k, v in bad.items() if v} or 'none'}")
    missing = set(sample.item_id) - set(p.item_id)
    if missing:
        print(f"  WARNING: {len(missing)} items missing, e.g. {sorted(missing)[:5]}")
    passes[tag] = p.set_index("item_id")

d = sample.set_index("item_id")
for tag, p in passes.items():
    for q in QS:
        d[f"{q}_{tag}"] = p[q].reindex(d.index)
        d[f"{q}_{tag}_reason"] = p[f"{q}_reason"].reindex(d.index)

# --------------------------------------------------------- stratum weights
sar = ds.dataset(R + "SAR_merged.parquet").to_table(
    columns=["source", "patent_core", "Target"]).to_pandas()
sar["Target"] = sar.Target.fillna("").str.strip()
pt = ds.dataset(R + "patents.parquet").to_table(columns=["Patent", "source_format"]).to_pandas()
sar2 = ds.dataset(R + "SAR_merged.parquet").to_table(columns=["Patent"]).to_pandas()
sar = pd.concat([sar, sar2], axis=1).merge(pt, on="Patent", how="left")
tm = pq.read_table(V + "target_map.parquet",
                   columns=["patent_core", "Target", "confidence",
                            "uniprot_accession"]).to_pandas().drop_duplicates(
    ["patent_core", "Target"])
sar = sar.merge(tm, on=["patent_core", "Target"], how="left")
sar["stratum"] = (sar.source.str.replace("USPTO_MOL_RECOVERY", "recovery")
                  .str.replace("USPTO_MOL", "baseline") + " | "
                  + sar.source_format.fillna("?") + " | " + sar.confidence.fillna("none"))
pop = sar.stratum.value_counts()
N = int(pop.sum())
# Q3 asks whether an ASSIGNED protein is right, so its population is the rows that have one
pop_q3 = sar[sar.uniprot_accession.fillna("").ne("")].stratum.value_counts()
print(f"\npopulation: {N:,} rows across {len(pop)} strata; sample covers {d.stratum.nunique()}")

report = {"population_rows": N, "sample_items": len(d), "questions": {}}

# ------------------------------------------------------------- per question
for q in QS:
    a, b = d[f"{q}_A"], d[f"{q}_B"]
    both = a.notna() & b.notna() & a.ne("") & b.ne("")
    kappa, po = cohen_kappa(a[both].tolist(), b[both].tolist())
    # the accuracy claim rests on "correct" versus anything else
    abin = a[both].eq("correct")
    bbin = b[both].eq("correct")
    kappa_bin, po_bin = cohen_kappa(abin.astype(str).tolist(), bbin.astype(str).tolist())

    agreed = both & a.eq(b)
    scope = both if q != "q3" else both & a.ne("not_applicable") & b.ne("not_applicable")
    # Accuracy is only defined where a verdict was possible. Counting agreed "uncertain" items as
    # failures understates accuracy by the undeterminable share, so the two are kept apart.
    judgeable = agreed & scope & a.isin(["correct", "incorrect"])
    cons = d[judgeable]
    k_cons = int(cons[f"{q}_A"].eq("correct").sum())
    n_cons = len(cons)
    lo, hi = wilson(k_cons, n_cons)
    n_undet = int((agreed & scope & a.eq("uncertain")).sum())

    full = d[scope]
    n_full = len(full)
    k_opt = int((full[f"{q}_A"].eq("correct") | full[f"{q}_B"].eq("correct")).sum())
    k_pes = int((full[f"{q}_A"].eq("correct") & full[f"{q}_B"].eq("correct")).sum())

    # Population weighting over the strata that actually have judgeable items, renormalised.
    # Without renormalising, a stratum that drops out (e.g. every "none" stratum for Q3, where no
    # protein is assigned) keeps its weight in the denominator and drags the estimate down.
    weights = pop_q3 if q == "q3" else pop
    present = [st for st in weights.index if (cons.stratum == st).any()]
    w_tot = float(sum(weights[st] for st in present))
    num, var = 0.0, 0.0
    per_stratum = {}
    for st in present:
        cell = cons[cons.stratum == st]
        p_h = cell[f"{q}_A"].eq("correct").mean()
        w = weights[st] / w_tot
        num += w * p_h
        var += (w ** 2) * p_h * (1 - p_h) / len(cell)
        per_stratum[st] = {"n": len(cell), "accuracy": round(float(p_h), 4),
                           "population_rows": int(weights[st])}
    se = float(np.sqrt(var))
    report["questions"][q] = {
        "scope_items": n_full, "items_both_annotated": int(both.sum()),
        "raw_agreement": round(float(po), 4), "cohen_kappa": round(float(kappa), 4),
        "raw_agreement_binary": round(float(po_bin), 4),
        "cohen_kappa_binary": round(float(kappa_bin), 4),
        "undeterminable_agreed": n_undet,
        "consensus": {"n": n_cons, "correct": k_cons,
                      "accuracy": round(k_cons / n_cons, 4) if n_cons else None,
                      "wilson_95": [round(lo, 4), round(hi, 4)]},
        "bounds_over_full_scope": {
            "optimistic": round(k_opt / n_full, 4) if n_full else None,
            "pessimistic": round(k_pes / n_full, 4) if n_full else None},
        "population_weighted": {"accuracy": round(num, 4),
                                "se": round(se, 4),
                                "ci_95": [round(max(0, num - Z * se), 4),
                                          round(min(1, num + Z * se), 4)]},
        "per_stratum": per_stratum,
        "verdict_counts": {"A": a.value_counts().to_dict(), "B": b.value_counts().to_dict()},
    }
    print(f"\n=== {q}")
    print(f"  raw agreement {po:.3f}   Cohen's kappa {kappa:.3f} "
          f"(binary correct-vs-rest: {po_bin:.3f} / {kappa_bin:.3f})")
    print(f"  undeterminable (both said uncertain): {n_undet}")
    print(f"  consensus accuracy {k_cons}/{n_cons} = "
          f"{(k_cons / n_cons if n_cons else float('nan')):.3f}  95% CI [{lo:.3f}, {hi:.3f}]")
    print(f"  bounds over all {n_full} in scope: pessimistic {k_pes / n_full:.3f}, "
          f"optimistic {k_opt / n_full:.3f}")
    print(f"  population-weighted {num:.3f} +/- {Z * se:.3f} "
          f"-> [{max(0, num - Z * se):.3f}, {min(1, num + Z * se):.3f}]")
    print(f"  A: {a.value_counts().to_dict()}")
    print(f"  B: {b.value_counts().to_dict()}")

# ------------------------------------------------------- disagreements
rows = []
for q in QS:
    m = d[f"{q}_A"].ne(d[f"{q}_B"]) & d[f"{q}_A"].notna() & d[f"{q}_B"].notna()
    for item, r in d[m].iterrows():
        rows.append({"item_id": item, "question": q, "stratum": r.stratum,
                     "Patent": r.Patent, "endpoint": r.endpoint,
                     "A": r[f"{q}_A"], "A_reason": r[f"{q}_A_reason"],
                     "B": r[f"{q}_B"], "B_reason": r[f"{q}_B_reason"]})
dis = pd.DataFrame(rows)
dis.to_csv(W + "/disagreements_for_adjudication.tsv", sep="\t", index=False)
print(f"\n{len(dis)} disagreements written for human adjudication "
      f"({dis.item_id.nunique() if len(dis) else 0} distinct items, "
      f"{(dis.item_id.nunique() / len(d) * 100) if len(dis) else 0:.1f}% of the sample)")
if len(dis):
    print(dis.question.value_counts().to_string())

# ------------------------------------------- items either annotator called wrong
wrong = []
for q in QS:
    m = d[f"{q}_A"].eq("incorrect") | d[f"{q}_B"].eq("incorrect")
    for item, r in d[m].iterrows():
        wrong.append({"item_id": item, "question": q, "Patent": r.Patent,
                      "stratum": r.stratum, "endpoint": r.endpoint,
                      "both_agree": r[f"{q}_A"] == r[f"{q}_B"],
                      "A": r[f"{q}_A"], "B": r[f"{q}_B"],
                      "A_reason": r[f"{q}_A_reason"], "B_reason": r[f"{q}_B_reason"]})
wr = pd.DataFrame(wrong)
wr.to_csv(W + "/flagged_incorrect.tsv", sep="\t", index=False)
print(f"{len(wr)} (item, question) pairs flagged incorrect by at least one annotator; "
      f"{int(wr.both_agree.sum()) if len(wr) else 0} with both agreeing")

report["automated_structure_check"] = {
    "items": len(sample),
    "inchikey_recomputed_from_patent_mol": int(sample.mol_inchikey.notna().sum()),
    "matches_claim": int(sample.mol_inchikey_matches.sum()),
    "located_by_searching_archive": int(sample.mol_found_by_search.sum()),
}
report["disagreements"] = {"pairs": len(dis), "distinct_items": int(dis.item_id.nunique()) if len(dis) else 0}
with open(W + "/validation_report.json", "w") as fh:
    json.dump(report, fh, indent=1, ensure_ascii=False, default=float)
print(f"\nwrote {W}/validation_report.json, disagreements_for_adjudication.tsv, flagged_incorrect.tsv")
