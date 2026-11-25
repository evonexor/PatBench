"""Build compact curation items (one per pair; strings whose every pair is a clean unique dictionary hit are deduped)."""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pandas as pd, json, re, sys
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from tmaplib import key
W=_ws + '/outputs/target-uniprot-map-20261016/work/'
P=pd.read_pickle(W+'pair_candidates.pkl')
ctx={}
with open(W+'pair_context.jsonl') as f:
    for line in f:
        r=json.loads(line); ctx[(r['patent_core'],r['Target'])]=r['snippets']
sec={}
with open(W+'patent_assay_sections.jsonl') as f:
    for line in f:
        r=json.loads(line); sec[r['patent_core']]=r['assay_sections']
others=P.groupby('patent_core').Target.agg(list).to_dict()
def strong(x): return bool(x) and x[0]['n_vmatch']>=5 and x[0]['n_vmatch']/x[0]['n_mol']>=0.5
P['du']=P.dict.map(lambda d:any(h['via'] in('full','core','org_stripped') and len(h['ids'])==1 for h in d))
P['clean']=P.du & P.cell_lexicon.str.len().eq(0) & P.mutations.str.len().eq(0)
allclean=P.groupby('Target').clean.transform('all')
P=P.reset_index(drop=True); P['pair_id']=['p%05d'%i for i in range(len(P))]
def sh(s,n): s=re.sub(r'\s+',' ',s or '').strip(); return s if len(s)<=n else s[:n]+'…'
def dict_str(d):
    ids=[]
    for h in d:
        for i in h['ids']:
            if i not in ids: ids.append(i)
    return ';'.join(ids[:5])
def xref_str(x):
    if not x: return ''
    a=x[0]; s=f"{a['db']}: {sh(a['target'],70)} [{a['acc']}] {a['org']} type={a['type']} valuematch={a['n_vmatch']}/{a['n_mol']}"
    if a['desc']: s+=f" desc={sh(a['desc'],140)}"
    if len(x)>1 and x[1]['n_vmatch']>0: s+=f" | 2nd: {sh(x[1]['target'],50)} [{x[1]['acc']}] vm={x[1]['n_vmatch']}"
    return s
items=[]; pair2item={}
done=set()
for r in P.itertuples():
    snips=ctx.get((r.patent_core,r.Target),[])
    generic=not r.dict
    lines=[f"Target: {r.Target}", f"patent: {r.patent_core} | acts: {r.acts} | rows: {r.n_rows}", f"title: {sh(r.title,140)}"]
    if r.activity_header: lines.append(f"header: {sh(r.activity_header,160)}")
    if snips: lines.append(f"text: {sh(snips[0],420 if not generic else 520)}")
    if generic:
        for s in (sec.get(r.patent_core) or [])[:2]: lines.append(f"assay-section: {sh(s,260)}")
        o=[t for t in others[r.patent_core] if t!=r.Target][:6]
        if o: lines.append('other targets in patent: '+' ; '.join(sh(t,40) for t in o))
    if r.dict: lines.append(f"dict-hints: {dict_str(r.dict)}")
    if r.xref: lines.append(f"xref: {xref_str(r.xref)}")
    if allclean[r.Index]:
        sid='s'+key(r.Target)[:40]
        if sid in done:
            pair2item[r.pair_id]=sid
            continue
        sub=P[P.Target==r.Target]
        lines=[f"Target: {r.Target}", f"STRING-LEVEL item covering {len(sub)} patents; examples below", f"dict-hints: {dict_str(r.dict)}"]
        for q in sub.head(3).itertuples():
            lines.append(f"- {q.patent_core} title: {sh(q.title,110)}" + (f" | header: {sh(q.activity_header,110)}" if q.activity_header else ''))
        xs=[xref_str(q.xref) for q in sub.itertuples() if q.xref][:2]
        for s in xs: lines.append(f"xref: {s}")
        done.add(sid); pair2item[r.pair_id]=sid
        items.append((sid,'\n'.join(lines)))
    else:
        pair2item[r.pair_id]=r.pair_id
        items.append((r.pair_id,'\n'.join(lines)))
P['item_id']=P.pair_id.map(pair2item)
P.to_pickle(W+'pair_candidates_ids.pkl')
B=350; nb=0
import os; os.makedirs(W+'curation',exist_ok=True)
for i in range(0,len(items),B):
    nb+=1
    with open(W+f'curation/batch_{nb:02d}.txt','w') as f:
        for iid,txt in items[i:i+B]: f.write(f"### {iid}\n{txt}\n\n")
print('items',len(items),'batches',nb, 'chars',sum(len(t) for _,t in items))
