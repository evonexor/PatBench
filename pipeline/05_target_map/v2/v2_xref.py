"""Cross-reference every (patent_core, Target) pair against the curated patent reference built in v1
(BindingDB 2026-09 patents + ChEMBL 37 src 37/38). Now over ALL endpoint types: molecules are matched by
InChIKey connectivity block within the patent family; value agreement is only claimed when the endpoint
type matches (so GI50/MIC/etc. pairs get molecule overlap but no value corroboration)."""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pyarrow.dataset as ds, pandas as pd, numpy as np
V1W=_ws + '/outputs/target-uniprot-map-20261016/work/'
W=_ws + '/outputs/target-uniprot-map-v2-20261019/work/'
R=_ws + '/outputs/uspto-sar-recovery-20261006/results/'
REF=pd.read_parquet(V1W+'ref_patent_activities.parquet')
m=ds.dataset(R+'SAR_merged.parquet').to_table(columns=['patent_core','Target','Activity Type','InChIKey','value_nM']).to_pandas()
m=m[m.InChIKey!=''].copy(); m['ik14']=m.InChIKey.str[:14]
tot=m.groupby(['patent_core','Target']).agg(n_mol=('ik14','nunique')).reset_index()
mm=m[m.patent_core.isin(set(REF.patent_core))]
J=mm.merge(REF,on=['patent_core','ik14'])
same=J['Activity Type']==J.act
J['vmatch']=same&(np.abs(np.log10(J.value_nM.clip(lower=1e-9))-np.log10(J.ref_val.clip(lower=1e-9)))<=0.02)
J['ik_v']=np.where(J.vmatch,J.ik14,None)
key=['patent_core','Target','ref_db','ref_target','ref_organism','ref_acc','ref_type']
C=J.groupby(key).agg(n_overlap=('ik14','nunique'),n_vmatch=('ik_v','nunique'),ref_desc=('ref_desc','first'),
    acts=('act',lambda s:','.join(sorted(set(s))))).reset_index()
C=C.merge(tot,on=['patent_core','Target'])
C['vfrac']=C.n_vmatch/C.n_mol
C.to_parquet(W+'xref_candidates_all.parquet',index=False)
P=pd.read_parquet(W+'all_pairs.parquet')
top=C.sort_values(['n_vmatch','n_overlap'],ascending=False).groupby(['patent_core','Target']).first().reset_index()
x=P.merge(top,on=['patent_core','Target'],how='left')
v1=set(zip(*[pd.read_parquet(_ws + '/outputs/target-uniprot-map-20261016/target_map.parquet',columns=['patent_core','Target'])[c] for c in ['patent_core','Target']]))
x['new']=[(p,t) not in v1 for p,t in zip(x.patent_core,x.Target)]
nw=x[x.new]
print('new pairs',len(nw),'with xref overlap',nw.ref_db.notna().sum(),'strong (vm>=5,vfrac>=0.5)',((nw.n_vmatch>=5)&(nw.vfrac>=0.5)).sum(),
      'any value match',(nw.n_vmatch.fillna(0)>0).sum())
print(C.shape)
