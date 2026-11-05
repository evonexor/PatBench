"""Build per-(patent_core, Target) cross-reference candidates against curated patent databases
(BindingDB patents, ChEMBL 37 patent sources 37/38). Molecules are matched by InChIKey connectivity block
within the same patent family number; values matched within 5% (|log10 ratio|<=0.02) for the same activity type."""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pyarrow.dataset as ds, pyarrow.parquet as pq, pandas as pd, numpy as np, re, gzip
W=_ws + '/outputs/target-uniprot-map-20261016/work/'
R=_ws + '/outputs/uspto-sar-recovery-20261006/results/'
S=_ws + '/raw/patent-data/'
def core(p):
    p=p.strip().replace('-','').replace(' ','')
    m=re.match(r'^US(RE)?0*(\d+)([A-Z]\d?)?$',p)
    return ('US'+(m.group(1) or '')+m.group(2)) if m else p
# ---- BindingDB
L=pd.read_parquet(W+'bindingdb_long.parquet')
bd=pd.DataFrame({'patent_core':L.patent_core,'ik14':L.ik14,'act':L.act,'ref_val':L.b_val,'ref_rel':L.b_rel,
   'ref_target':L.target_name,'ref_organism':L.organism,'ref_acc':L.acc,'ref_type':np.where(L.n_chains.fillna('1').isin(['1','']),'SINGLE PROTEIN','MULTICHAIN'),
   'ref_desc':'', 'ref_db':'BindingDB_202609'})
# ---- ChEMBL patent sources
c=pq.read_table(S+'standardized/ChEMBL.parquet',columns=['Patent','Target','Activity Type','Relation','Value','Unit','UniProt','molecule_chembl_id','target_chembl_id','target_type','organism','src_id','assay_description']).to_pandas()
c=c[c['Activity Type'].isin(['IC50','Ki','Kd','EC50'])&(c.Unit=='nM')]
ids=set(c.molecule_chembl_id)
ik={}
with gzip.open(S+'chembl/chembl_37/chembl_37_chemreps.txt.gz','rt') as f:
    next(f)
    for line in f:
        a=line.rstrip('\n').split('\t')
        if a[0] in ids: ik[a[0]]=a[3]
c['inchikey']=c.molecule_chembl_id.map(ik).fillna('')
c=c[c.inchikey!='']
ch=pd.DataFrame({'patent_core':c.Patent.map(core),'ik14':c.inchikey.str[:14],'act':c['Activity Type'],'ref_val':pd.to_numeric(c.Value,errors='coerce'),
   'ref_rel':c.Relation,'ref_target':c.Target+' ['+c.target_chembl_id+']','ref_organism':c.organism,'ref_acc':c.UniProt.fillna(''),'ref_type':c.target_type,
   'ref_desc':c.assay_description.fillna(''),'ref_db':'ChEMBL37_src'+c.src_id})
REF=pd.concat([bd,ch],ignore_index=True)
REF.to_parquet(W+'ref_patent_activities.parquet',index=False)
print('ref rows',REF.groupby('ref_db').size().to_dict(),'patents',REF.groupby('ref_db').patent_core.nunique().to_dict())
m=ds.dataset(R+'SAR_merged.parquet').to_table(columns=['patent_core','Target','Activity Type','InChIKey','value_nM']).to_pandas()
m=m[m['Activity Type'].isin(['IC50','Ki','Kd','EC50'])&(m.InChIKey!='')].copy()
m['ik14']=m.InChIKey.str[:14]
tot=m.groupby(['patent_core','Target']).agg(n_mol=('ik14','nunique')).reset_index()
mm=m[m.patent_core.isin(set(REF.patent_core))]
J=mm.merge(REF,on=['patent_core','ik14'])
J['vmatch']=(J['Activity Type']==J.act)&(np.abs(np.log10(J.value_nM.clip(lower=1e-9))-np.log10(J.ref_val.clip(lower=1e-9)))<=0.02)
J['ik_v']=np.where(J.vmatch,J.ik14,None)
key=['patent_core','Target','ref_db','ref_target','ref_organism','ref_acc','ref_type']
C=J.groupby(key).agg(n_overlap=('ik14','nunique'),n_vmatch=('ik_v','nunique'),ref_desc=('ref_desc','first'),
    acts=('act',lambda s:','.join(sorted(set(s))))).reset_index()
C=C.merge(tot,on=['patent_core','Target'])
C['vfrac']=C.n_vmatch/C.n_mol
C.to_parquet(W+'xref_candidates.parquet',index=False)
g=pd.read_parquet(W+'ranking_groups.parquet'); rp=g[['patent_core','Target']].drop_duplicates()
top=C.sort_values('n_vmatch',ascending=False).groupby(['patent_core','Target']).first().reset_index()
rt=rp.merge(top,on=['patent_core','Target'],how='left')
print('ranking pairs',len(rp),'any overlap',rt.ref_db.notna().sum(),'vfrac>=0.5',(rt.vfrac>=0.5).sum(),'vmatch>=5',(rt.n_vmatch>=5).sum())
print(rt.ref_db.value_counts())
print(rt.vfrac.describe())
