"""Cross-reference each (patent_core, Target) in SAR_merged with BindingDB curation of the same patent,
matching molecules by InChIKey and activity values. Output: per (patent_core, Target, BindingDB target) candidate stats."""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pyarrow.dataset as ds, pandas as pd, numpy as np, re
W=_ws + '/outputs/target-uniprot-map-20261016/work/'
R=_ws + '/outputs/uspto-sar-recovery-20261006/results/'
b=pd.read_parquet(W+'bindingdb_patents_slim.parquet')
def core(p):
    p=p.strip().replace('-','').replace(' ','')
    m=re.match(r'^(US)(RE)?0*(\d+)([A-Z]\d?)?$',p)
    if not m: return p
    return 'US'+(m.group(2) or '')+m.group(3)
# expand to all patents named for the ligand (family members sharing compound numbering)
pat=[]
for i,(p,ln) in enumerate(zip(b.patent.values,b.ligand_name.values)):
    s={core(p)} if p else set()
    for part in ln.split('::'):
        m=re.match(r'\s*(US\s?(?:RE)?\d{5,8})(?:[A-Z]\d)?\s*,',part)
        if m: s.add(core(m.group(1)))
    for x in s: pat.append((i,x))
pp=pd.DataFrame(pat,columns=['bi','patent_core'])
print('bindingdb rows',len(b),'patent links',len(pp),'patents',pp.patent_core.nunique())
long=[]
for at in ['Ki','IC50','Kd','EC50']:
    v=b[at].str.strip()
    ok=v!=''
    t=b.loc[ok,['inchikey','target_name','organism','chain_ids','n_chains','sp_id1','tr_id1','sp_name1','sp_entry1']].copy()
    t['act']=at; raw=v[ok]
    t['b_rel']=raw.str.extract(r'^([<>~=]*)')[0]
    t['b_val']=pd.to_numeric(raw.str.replace(r'^[<>~=\s]+','',regex=True),errors='coerce')
    t['bi']=t.index
    long.append(t)
L=pd.concat(long,ignore_index=True)
L['acc']=np.where(L.chain_ids!='',L.chain_ids,'')
L['ik14']=L.inchikey.str[:14]
L=L.merge(pp,on='bi')
print('long',len(L))
m=ds.dataset(R+'SAR_merged.parquet').to_table(columns=['patent_core','Target','Activity Type','InChIKey','value_nM','Relation']).to_pandas()
m=m[m['Activity Type'].isin(['IC50','Ki','Kd','EC50'])&(m.InChIKey!='')].copy()
m['ik14']=m.InChIKey.str[:14]
pcs=set(L.patent_core)
print('SAR patents',m.patent_core.nunique(),'in bindingdb',len(set(m.patent_core)&pcs))
mm=m[m.patent_core.isin(pcs)]
# molecule-level join on patent + ik14 (connectivity block) ; full-key equality flag kept
J=mm.merge(L[['patent_core','ik14','inchikey','target_name','organism','acc','n_chains','sp_name1','act','b_val','b_rel']],on=['patent_core','ik14'])
J['full_ik']=J.InChIKey==J.inchikey
J['same_act']=J['Activity Type']==J.act
J['vmatch']=J.same_act & (np.abs(np.log10(J.value_nM.clip(lower=1e-6))-np.log10(J.b_val.clip(lower=1e-6)))<=0.02)
key=['patent_core','Target','target_name','organism','acc']
c=J.groupby(key).agg(n_ik_overlap=('ik14','nunique'),n_same_act=('same_act','sum'),n_vmatch=('vmatch','sum'),
    n_chains=('n_chains','first'),sp_name1=('sp_name1','first')).reset_index()
tot=m.groupby(['patent_core','Target']).agg(n_mol=('ik14','nunique'),n_rows=('ik14','size')).reset_index()
c=c.merge(tot,on=['patent_core','Target'])
c.to_parquet(W+'xref_bindingdb_candidates.parquet',index=False)
L.to_parquet(W+'bindingdb_long.parquet',index=False)
print(c.shape)
# best candidate per pair by n_vmatch then overlap
c=c.sort_values(['patent_core','Target','n_vmatch','n_ik_overlap'],ascending=[True,True,False,False])
best=c.groupby(['patent_core','Target']).head(2)
g=pd.read_parquet(W+'ranking_groups.parquet')
rp=g[['patent_core','Target']].drop_duplicates()
top=c.groupby(['patent_core','Target']).first().reset_index()
rt=rp.merge(top,on=['patent_core','Target'],how='left')
print('ranking pairs',len(rp),'with any bindingdb overlap',rt.target_name.notna().sum())
rt['vfrac']=rt.n_vmatch/rt.n_mol
print(rt.vfrac.describe()); print((rt.vfrac>=0.5).sum(),(rt.vfrac>=0.2).sum(),(rt.n_vmatch>=5).sum())
