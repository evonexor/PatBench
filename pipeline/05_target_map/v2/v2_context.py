"""Patent-text context for the pairs not covered by v1 (same extraction as v1's patent_context.py)."""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pyarrow.dataset as ds, pandas as pd, zipfile, re, json, html
from multiprocessing import Pool
W=_ws + '/outputs/target-uniprot-map-v2-20261019/work/'
R=_ws + '/outputs/uspto-sar-recovery-20261006/results/'
A=_ws + '/outputs/uspto-all-structures-20261004/results/'
F=_ws + '/outputs/uspto-sar-full-20261004/results/'
ROOT=_ws + '/raw/patent-data/'
P=pd.read_parquet(W+'all_pairs.parquet')
v1=pd.read_parquet(_ws + '/outputs/target-uniprot-map-20261016/target_map.parquet',columns=['patent_core','Target'])
key=set(zip(v1.patent_core,v1.Target))
pairs=P[[(p,t) not in key for p,t in zip(P.patent_core,P.Target)]].copy()
print('new pairs',len(pairs),flush=True)
arc=pd.concat([ds.dataset(R+'new_SAR.parquet').to_table(columns=['patent_core','archive']).to_pandas(),
               ds.dataset(A+'SAR_enriched_annotations.parquet').to_table(columns=['patent_core','archive']).to_pandas()]).drop_duplicates('patent_core')
arcmap=dict(zip(arc.patent_core,arc.archive))
hdr=pd.concat([ds.dataset(R+'new_SAR.parquet').to_table(columns=['patent_core','Target','activity_header']).to_pandas(),
   ds.dataset(F+'candidate_SAR.parquet').to_table(columns=['patent_core','Target','activity_header']).to_pandas()])
hdr=hdr[hdr.activity_header!=''].groupby(['patent_core','Target']).activity_header.agg(lambda s:' || '.join(s.value_counts().index[:3])).to_dict()
missing=[p for p in pairs.patent_core.unique() if p not in arcmap]
print('patents without archive path',len(missing),flush=True)
def clean(x):
    x=re.sub(r'<\?.*?\?>|<!--.*?-->','',x,flags=re.S)
    x=re.sub(r'<(entry|td|th)[^>]*>',' | ',x); x=re.sub(r'</(p|heading|row|tr|claim|title|invention-title)>','\n',x)
    x=re.sub(r'<[^>]+>',' ',x); x=html.unescape(x)
    x=re.sub(r'[⁡-⁤​]','',x); x=re.sub(r'[ \t\r\f\v]+',' ',x)
    return re.sub(r'\n\s*','\n',x)
by_pat={p:g for p,g in pairs.groupby('patent_core')}
def work(pc):
    out=[]; a=arcmap.get(pc); text=''; title=''; err=''
    if a is None: err='no_archive_path'
    else:
        try:
            z=zipfile.ZipFile(ROOT+a)
            xs=[n for n in z.namelist() if n.upper().endswith(('.XML','.SGM','.SGML'))]
            xs.sort(key=lambda n:-z.getinfo(n).file_size)
            raw=z.read(xs[0]).decode('utf-8',errors='replace')
            mt=re.search(r'<invention-title[^>]*>(.*?)</invention-title>',raw,re.S) or re.search(r'<B540>(.*?)</B540>',raw,re.S)
            title=clean(mt.group(1)).strip() if mt else ''
            text=clean(re.sub(r'<claims.*?</claims>','',raw,flags=re.S))
        except Exception as ex: err=repr(ex)
    low=text.lower(); heads=[]
    for mo in re.finditer(r'\n([^\n]{3,140})\n',text):
        h=mo.group(1)
        if re.search(r'assay|biolog|inhibit|binding|activity|test',h,re.I) and len(h.split())<=18 and '|' not in h:
            heads.append((mo.start(1),h.strip()))
    assay=[re.sub(r'\s+',' ',text[p:p+700]) for p,_ in heads[:6]]
    for r in by_pat[pc].itertuples():
        t=r.Target.strip(); snips=[]
        if len(t)>=2 and text:
            pat=re.escape(t.lower()) if len(t)>3 else r'(?<![a-z0-9])'+re.escape(t.lower())+r'(?![a-z0-9])'
            hits=[mo.start() for mo in re.finditer(pat,low)]
            def score(p):
                w=text[max(0,p-300):p+300]
                return (w.count('|')<4)*2+bool(re.search(r'assay|inhibit|binding|recombinant|enzyme|receptor|cells?\b',w,re.I))
            hits=sorted(hits[:400],key=lambda p:-score(p))[:3]
            snips=[re.sub(r'\s+',' ',text[max(0,p-350):p+350]) for p in sorted(hits)]
        out.append({'patent_core':pc,'Target':r.Target,'endpoints':r.endpoints,'n_rows':int(r.n_rows),'n_mol':int(r.n_mol),
            'archive':a or '','title':title,'activity_header':hdr.get((pc,r.Target),''),'n_hits':len(snips),
            'snippets':snips,'assay_sections':assay,'err':err})
    return out
if __name__=='__main__':
    pcs=sorted(pairs.patent_core.unique()); allp=[]
    with Pool(22) as p:
        for i,o in enumerate(p.imap_unordered(work,pcs,chunksize=4)):
            allp+=o
            if i%500==0: print(i,len(pcs),flush=True)
    with open(W+'pair_context_new.jsonl','w') as f:
        for r in allp: f.write(json.dumps(r,ensure_ascii=False)+'\n')
    print('pairs',len(allp),'with hits',sum(1 for r in allp if r['n_hits']),'errors',sum(1 for r in allp if r['err']))
