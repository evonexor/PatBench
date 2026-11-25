_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import sys,pathlib,json,hashlib
sys.path.insert(0,_ws + '/outputs/uspto-sar-demo-20261003/runtime/deps')
r=pathlib.Path(_ws + '/outputs/patbench-v1.0-20261030-r6');o=r/'structure_repair/nitrile_audit';sys.path.insert(0,str(o))
import pandas as pd
from verified_nitrile_guard import NitrileEvidenceGuard,RepairNotApproved
rec=pd.read_parquet(o/'record_decisions.parquet');guard=NitrileEvidenceGuard();accepted=rejected=0
for t in rec.itertuples():
 try:guard.require_approved(t.patent_core,t.structure_file,t.smiles_before,t.smiles_after);accepted+=1;assert t.decision=='supported_nitrile'
 except RepairNotApproved:rejected+=1;assert t.decision=='holdout_unresolved'
neg=[]
t=rec[rec.decision=='supported_nitrile'].iloc[0]
for args in [('US_UNKNOWN','unknown.MOL','C[N+]#[C-]','CC#N'),(t.patent_core,t.structure_file,'CC',t.smiles_after),(t.patent_core,t.structure_file,t.smiles_before,'CC#N')]:
 try:guard.require_approved(*args);neg.append(False)
 except RepairNotApproved:neg.append(True)
assert all(neg)
s=pd.read_parquet(r/'sar_all.parquet');d=pd.read_parquet(o/'sar_all_reviewed.parquet');assert len(s)==len(d)
for col in s:
 if col!='structure_alias_suspect':assert s[col].equals(d[col]),col
assert (s.structure_alias_suspect!=d.structure_alias_suspect).sum()==3
assert d.nitrile_audit_holdout.sum()==3
m=pd.read_parquet(o/'motif_decisions.parquet');assert m[m.test_tasks.ne('')].image_visually_reviewed.all();assert (m[m.test_tasks.ne('')].decision=='supported_nitrile').all()
result={'approved_candidates_passed':accepted,'unresolved_candidates_rejected':rejected,'negative_checks_passed':len(neg),'derived_sar_rows':len(d),'molecular_or_activity_columns_changed':0,'additional_structure_holdouts':3,'source_images_missing':sum(not (o/'original_images'/p).exists() for p in m[m.image_visually_reviewed].image_file.unique())}
assert result['source_images_missing']==0
(o/'verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2));print(result)
