"""Enumerate every (patent_core, Target) pair in SAR_merged across all endpoint types, with group sizes,
and check which are already covered by the v1 map."""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pyarrow.dataset as ds, pandas as pd, os
R=_ws + '/outputs/uspto-sar-recovery-20261006/results/'
V1=_ws + '/outputs/target-uniprot-map-20261016/'
W=_ws + '/outputs/target-uniprot-map-v2-20261019/work/'
os.makedirs(W,exist_ok=True)
m=ds.dataset(R+'SAR_merged.parquet').to_table(columns=['patent_core','Target','Activity Type','InChIKey','value_nM','Unit']).to_pandas()
print('rows',len(m),'patents',m.patent_core.nunique(),'pairs',m.groupby(['patent_core','Target']).ngroups)
AFF=['IC50','Ki','Kd','EC50']
grp=m.groupby(['patent_core','Target','Activity Type']).agg(n_mol=('InChIKey','nunique'),n_rows=('InChIKey','size')).reset_index()
grp.to_parquet(W+'all_groups.parquet',index=False)
print('groups(all endpoints,all sizes)',len(grp),'rows',grp.n_rows.sum())
pairs=m.groupby(['patent_core','Target']).agg(n_mol=('InChIKey','nunique'),n_rows=('InChIKey','size'),
    endpoints=('Activity Type',lambda s:','.join(sorted(set(s))))).reset_index()
pairs.to_parquet(W+'all_pairs.parquet',index=False)
print('distinct pairs',len(pairs))
v1=pd.read_parquet(V1+'target_map.parquet',columns=['patent_core','Target','confidence','target_class'])
key=set(zip(v1.patent_core,v1.Target))
pairs['in_v1']=[ (p,t) in key for p,t in zip(pairs.patent_core,pairs.Target)]
print('pairs in v1',pairs.in_v1.sum(),'new',(~pairs.in_v1).sum(), 'pct', round(100*pairs.in_v1.mean(),1))
print('patents total',pairs.patent_core.nunique(),'patents with any v1 pair',pairs[pairs.in_v1].patent_core.nunique())
v1m=v1[v1.confidence.isin(['high','medium'])]
k2=set(zip(v1m.patent_core,v1m.Target))
pairs['v1_mapped']=[ (p,t) in k2 for p,t in zip(pairs.patent_core,pairs.Target)]
print('patents with >=1 high/medium v1 pair',pairs[pairs.v1_mapped].patent_core.nunique())
print('--- new pairs by endpoint')
nw=pairs[~pairs.in_v1]
print(nw.endpoints.str.split(',').explode().value_counts().to_string())
print('--- new pairs by n_mol bucket')
print(pd.cut(nw.n_mol,[0,1,2,4,9,19,49,10**6]).value_counts().sort_index().to_string())
print('new pairs rows',nw.n_rows.sum(),'patents only in new',len(set(nw.patent_core)-set(pairs[pairs.in_v1].patent_core)))
print('distinct new Target strings',nw.Target.nunique(),'of which also appear in v1 pairs',nw.Target.isin(set(v1.Target)).sum())
