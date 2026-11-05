"""Build v2 curation items (10-column output: v1's 9 plus an explicit mappable judgment)."""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pandas as pd, json, re, sys, os
sys.path.insert(0,__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', 'v1'))
from tmaplib import key
W=_ws + '/outputs/target-uniprot-map-v2-20261019/work/'
P=pd.read_pickle(W+'pair_candidates_new.pkl')
P=P.reset_index(drop=True); P['pair_id']=['q%05d'%i for i in range(len(P))]
def sh(s,n):
    s=re.sub(r'\s+',' ',s or '').strip(); return s if len(s)<=n else s[:n]+'…'
def dict_str(d):
    ids=[]
    for h in d:
        for i in h['ids']:
            if i not in ids: ids.append(i)
    return ';'.join(ids[:5])
def xref_str(x):
    if not x: return ''
    a=x[0]; s=f"{a['db']}: {sh(a['target'],70)} [{a['acc']}] {a['org']} type={a['type']} overlap={a['n_overlap']} valuematch={a['n_vmatch']}/{a['n_mol']}"
    if a['desc']: s+=f" desc={sh(a['desc'],140)}"
    return s
others=P.groupby('patent_core').Target.agg(list).to_dict()
todo=P[~P.reuse_v1]
items=[]; pair2item={}; done=set()
for r in todo.itertuples():
    generic=not r.dict
    if r.string_item:
        sid='t'+key(r.Target)[:40]
        if sid in done:
            pair2item[r.pair_id]=sid; continue
        sub=todo[todo.Target==r.Target]
        lines=[f"Target: {r.Target}", f"STRING-LEVEL item covering {len(sub)} patents; examples below",
               f"dict-hints: {dict_str(r.dict)}"]
        for q in sub.head(3).itertuples():
            lines.append(f"- {q.patent_core} [{q.endpoints}, {q.n_mol} mols] title: {sh(q.title,110)}"+(f" | header: {sh(q.activity_header,100)}" if q.activity_header else ''))
        for q in sub.itertuples():
            if q.xref: lines.append('xref: '+xref_str(q.xref)); break
        done.add(sid); pair2item[r.pair_id]=sid; items.append((sid,'\n'.join(lines)))
        continue
    lines=[f"Target: {r.Target}",
           f"patent: {r.patent_core} | endpoints: {r.endpoints} | molecules: {r.n_mol} | rows: {r.n_rows}",
           f"title: {sh(r.title,140)}"]
    if r.activity_header: lines.append(f"header: {sh(r.activity_header,160)}")
    if r.snippets: lines.append(f"text: {sh(r.snippets[0],420 if not generic else 520)}")
    if generic or not r.snippets:
        for s in (r.assay_sections or [])[:2]: lines.append(f"assay-section: {sh(s,260)}")
        o=[t for t in others[r.patent_core] if t!=r.Target][:6]
        if o: lines.append('other targets in patent: '+' ; '.join(sh(t,40) for t in o))
    if r.dict: lines.append(f"dict-hints: {dict_str(r.dict)}")
    if r.xref: lines.append(f"xref: {xref_str(r.xref)}")
    pair2item[r.pair_id]=r.pair_id; items.append((r.pair_id,'\n'.join(lines)))
P['item_id']=P.pair_id.map(pair2item).fillna('')
P.to_pickle(W+'pair_candidates_new_ids.pkl')
os.makedirs(W+'curation',exist_ok=True)
B=340; nb=0
for i in range(0,len(items),B):
    nb+=1
    with open(W+f'curation/batch_{nb:02d}.txt','w') as f:
        for iid,txt in items[i:i+B]: f.write(f"### {iid}\n{txt}\n\n")
print('items',len(items),'batches',nb,'chars',sum(len(t) for _,t in items))
