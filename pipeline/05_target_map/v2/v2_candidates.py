"""Deterministic layer for the new pairs, then build curation items.
Labels that v1 already settled as a clean string-level decision are reused verbatim (same label, same clean
dictionary hit, no cell/mutant cue) instead of being curated again, so the two runs stay consistent."""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pandas as pd, json, re, sys, os
sys.path.insert(0,__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', 'v1'))
from tmaplib import key, clean_text, strip_noise, stated_organism, mutations
V1=_ws + '/outputs/target-uniprot-map-20261016/'
V1W=V1+'work/'
W=_ws + '/outputs/target-uniprot-map-v2-20261019/work/'
D=pd.read_parquet(V1W+'alias_dictionary.parquet')
k2id=D.groupby('key').identity.agg(lambda s:sorted(set(s))).to_dict()
k2kind=D.groupby('key').kind.agg(lambda s:sorted(set(s))).to_dict()
cells=pd.read_parquet(V1W+'chembl_cells.parquet')
cellkeys={k for k in cells.cell_name.map(key) if len(k)>=4}
def lookups(t):
    v=[('full',key(t)),('core',key(strip_noise(t)))]
    c=clean_text(strip_noise(t))
    m=re.match(r'^(h|r|m|hu|human|rat|mouse)[\s-]?(?=[A-Z0-9])',c,re.I)
    if m: v.append(('org_stripped',key(c[m.end():])))
    for w in [w for w in re.split(r'[\s/,;]+',c) if len(key(w))>=3]: v.append(('token',key(w)))
    out=[]
    for kind,k in v:
        if k and k in k2id: out.append({'via':kind,'key':k,'ids':k2id[k][:6],'kinds':k2kind[k]})
    return out
C=pd.read_parquet(W+'xref_candidates_all.parquet').sort_values(['n_vmatch','n_overlap'],ascending=False)
best={}
for (pc,t),g in C.groupby(['patent_core','Target'],sort=False):
    best[(pc,t)]=[{'db':r.ref_db,'target':r.ref_target,'acc':r.ref_acc,'org':r.ref_organism,'type':r.ref_type,
                   'n_overlap':int(r.n_overlap),'n_vmatch':int(r.n_vmatch),'n_mol':int(r.n_mol),'desc':r.ref_desc[:200]} for r in g.head(3).itertuples()]
rows=[]
with open(W+'pair_context_new.jsonl') as f:
    for line in f:
        r=json.loads(line); t=r['Target']
        hay=key(t)+' '+key(r['activity_header'])
        rows.append({**{k:r[k] for k in ['patent_core','Target','endpoints','n_rows','n_mol','title','activity_header','snippets','assay_sections']},
            'org_stated':stated_organism(t) or stated_organism(r['activity_header']),
            'mutations':mutations(t+' '+r['activity_header']),
            'cell_lexicon':sorted({ck for ck in cellkeys if ck in hay})[:5],
            'dict':lookups(t),'xref':best.get((r['patent_core'],t),[])})
P=pd.DataFrame(rows)
P['du']=P.dict.map(lambda d:any(h['via'] in('full','core','org_stripped') and len(h['ids'])==1 for h in d))
P['clean']=P.du & P.cell_lexicon.str.len().eq(0) & P.mutations.str.len().eq(0)
allclean=P.groupby('Target').clean.transform('all')
P['string_item']=allclean
# v1 string-level decisions available for reuse
V=pd.read_parquet(V1+'target_map.parquet')
v1s=V[V.item_id.str.startswith('s')].drop_duplicates('Target')
v1map={r.Target:r for r in v1s.itertuples()}
P['reuse_v1']=P.string_item & P.Target.isin(v1map.keys())
print('new pairs',len(P),'string-level',P.string_item.sum(),'reusable from v1',P.reuse_v1.sum(),
      'pairs needing curation',(~P.reuse_v1).sum())
P.to_pickle(W+'pair_candidates_new.pkl')
print('distinct strings needing curation',P[~P.reuse_v1].Target.nunique())
print('items if deduped:', P[~P.reuse_v1 & P.string_item].Target.nunique() + (~P.string_item & ~P.reuse_v1).sum())
