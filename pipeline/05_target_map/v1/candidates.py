"""Deterministic layer per ranking (patent_core, Target) pair: stated organism, mutations,
ChEMBL cell-line lexicon hit, dictionary candidate identities, best cross-reference candidate."""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pandas as pd, json, re, sys
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from tmaplib import key, clean_text, strip_noise, stated_organism, mutations
W=_ws + '/outputs/target-uniprot-map-20261016/work/'
D=pd.read_parquet(W+'alias_dictionary.parquet')
k2id=D.groupby('key').identity.agg(lambda s:sorted(set(s))).to_dict()
k2kind=D.groupby('key').kind.agg(lambda s:sorted(set(s))).to_dict()
cells=pd.read_parquet(W+'chembl_cells.parquet')
cellkeys={k for k in cells.cell_name.map(key) if len(k)>=4}
def lookups(t):
    v=[('full',key(t)),('core',key(strip_noise(t)))]
    c=clean_text(strip_noise(t))
    m=re.match(r'^(h|r|m|hu|human|rat|mouse)[\s-]?(?=[A-Z0-9])',c,re.I)
    if m: v.append(('org_stripped',key(c[m.end():])))
    toks=[w for w in re.split(r'[\s/,;]+',c) if len(key(w))>=3]
    for w in toks: v.append(('token',key(w)))
    out=[]
    for kind,k in v:
        if k and k in k2id: out.append({'via':kind,'key':k,'ids':k2id[k][:6],'kinds':k2kind[k]})
    return out
C=pd.read_parquet(W+'xref_candidates.parquet').sort_values(['n_vmatch','n_overlap'],ascending=False)
best={}
for (pc,t),g in C.groupby(['patent_core','Target'],sort=False):
    g=g.head(3)
    best[(pc,t)]=[{'db':r.ref_db,'target':r.ref_target,'acc':r.ref_acc,'org':r.ref_organism,'type':r.ref_type,
                   'n_overlap':int(r.n_overlap),'n_vmatch':int(r.n_vmatch),'n_mol':int(r.n_mol),'desc':r.ref_desc[:200]} for r in g.itertuples()]
rows=[]
with open(W+'pair_context.jsonl') as f:
    for line in f:
        r=json.loads(line); t=r['Target']
        hay=key(t)+' '+key(r['activity_header'])
        cellhit=sorted({ck for ck in cellkeys if len(ck)>=4 and ck in hay})[:5]
        rows.append({**{k:r[k] for k in ['patent_core','Target','acts','n_rows','title','activity_header']},
            'org_stated':stated_organism(t) or stated_organism(r['activity_header']),
            'mutations':mutations(t+' '+r['activity_header']),'cell_lexicon':cellhit,
            'dict':lookups(t),'xref':best.get((r['patent_core'],t),[])})
P=pd.DataFrame(rows)
P.to_pickle(W+'pair_candidates.pkl')
def strong(x): return bool(x) and x[0]['n_vmatch']>=5 and x[0]['n_vmatch']/x[0]['n_mol']>=0.5
P['xref_strong']=P.xref.map(strong)
P['dict_full_unique']=P.dict.map(lambda d:any(h['via'] in('full','core','org_stripped') and len(h['ids'])==1 for h in d))
print(len(P)); print(pd.crosstab(P.xref_strong,P.dict_full_unique))
print('strings',P.Target.nunique(),'strings in >1 patent',(P.groupby('Target').size()>1).sum())
