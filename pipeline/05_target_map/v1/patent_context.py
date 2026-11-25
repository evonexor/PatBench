"""Extract per-(patent_core, Target) text context from the original USPTO XML for curation and spot checks.
Writes work/pair_context.jsonl (one record per ranking pair) and work/patent_accessions.parquet (accession-like mentions)."""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pyarrow.dataset as ds, pandas as pd, zipfile, re, json, html, sys
from multiprocessing import Pool
W=_ws + '/outputs/target-uniprot-map-20261016/work/'
R=_ws + '/outputs/uspto-sar-recovery-20261006/results/'
A=_ws + '/outputs/uspto-all-structures-20261004/results/'
ROOT=_ws + '/raw/patent-data/'
g=pd.read_parquet(W+'ranking_groups.parquet')
pairs=g.groupby(['patent_core','Target']).agg(acts=('Activity Type',lambda s:','.join(sorted(set(s)))),n_rows=('n_rows','sum')).reset_index()
arc=pd.concat([ds.dataset(R+'new_SAR.parquet').to_table(columns=['patent_core','archive']).to_pandas(),
               ds.dataset(A+'SAR_enriched_annotations.parquet').to_table(columns=['patent_core','archive']).to_pandas()]).drop_duplicates('patent_core')
arcmap=dict(zip(arc.patent_core,arc.archive))
hdr=pd.concat([ds.dataset(R+'new_SAR.parquet').to_table(columns=['patent_core','Target','activity_header']).to_pandas(),
   ds.dataset(_ws + '/outputs/uspto-sar-full-20261004/results/candidate_SAR.parquet').to_table(columns=['patent_core','Target','activity_header']).to_pandas()])
hdr=hdr[hdr.activity_header!=''].groupby(['patent_core','Target']).activity_header.agg(lambda s:' || '.join(s.value_counts().index[:3])).to_dict()
ACC=re.compile(r'\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:\.\d+)?\b')
ACCCTX=re.compile(r'uniprot|swiss|accession|genbank|ncbi|refseq|sequence ref|acc\.|\bacc\b',re.I)
REFSEQ=re.compile(r'\b([NX][PM]_\d{6,9}(?:\.\d+)?)\b')
def clean(x):
    x=re.sub(r'<\?.*?\?>|<!--.*?-->','',x,flags=re.S)
    x=re.sub(r'<(entry|td|th)[^>]*>',' | ',x)
    x=re.sub(r'</(p|heading|row|tr|claim|title|invention-title)>','\n',x)
    x=re.sub(r'<[^>]+>',' ',x)
    x=html.unescape(x)
    x=re.sub(r'[⁡-⁤​]','',x)
    x=re.sub(r'[ \t\r\f\v]+',' ',x)
    return re.sub(r'\n\s*','\n',x)
def work(pc):
    out=[]; accs=[]
    a=arcmap.get(pc)
    sub=pairs[pairs.patent_core==pc]
    text=''; title=''; err=''
    try:
        z=zipfile.ZipFile(ROOT+a)
        xs=[n for n in z.namelist() if n.upper().endswith(('.XML','.SGM','.SGML'))]
        xs.sort(key=lambda n:-z.getinfo(n).file_size)
        raw=z.read(xs[0]).decode('utf-8',errors='replace')
        m=re.search(r'<invention-title[^>]*>(.*?)</invention-title>',raw,re.S) or re.search(r'<B540>(.*?)</B540>',raw,re.S)
        title=clean(m.group(1)).strip() if m else ''
        # drop claims
        raw=re.sub(r'<claims.*?</claims>','',raw,flags=re.S)
        text=clean(raw)
    except Exception as ex:
        err=repr(ex)
    low=text.lower()
    # assay-section headings
    heads=[]
    for m in re.finditer(r'\n([^\n]{3,140})\n',text):
        h=m.group(1)
        if re.search(r'assay|biolog|inhibit|binding|activity|test',h,re.I) and len(h.split())<=18 and '|' not in h:
            heads.append((m.start(1),h.strip()))
    assay_snips=[]
    for pos,h in heads[:6]:
        assay_snips.append(text[pos:pos+700].replace('\n',' '))
    for m in ACC.finditer(text):
        ctx=text[max(0,m.start()-160):m.end()+80]
        if ACCCTX.search(ctx): accs.append((pc,m.group(1),ctx.replace('\n',' ')))
    for m in REFSEQ.finditer(text):
        accs.append((pc,m.group(1),text[max(0,m.start()-160):m.end()+80].replace('\n',' ')))
    for _,r in sub.iterrows():
        t=r.Target.strip(); snips=[]
        if len(t)>=2 and text:
            pat=re.escape(t.lower()) if len(t)>3 else r'(?<![a-z0-9])'+re.escape(t.lower())+r'(?![a-z0-9])'
            hits=[m.start() for m in re.finditer(pat,low)]
            # prefer prose (no table separators) and assay vocabulary
            def score(p):
                w=text[max(0,p-300):p+300]
                return (w.count('|')<4)*2+bool(re.search(r'assay|inhibit|binding|recombinant|enzyme|receptor|cells?\b',w,re.I))
            hits=sorted(hits[:400],key=lambda p:-score(p))[:3]
            snips=[text[max(0,p-350):p+350].replace('\n',' ') for p in sorted(hits)]
        out.append({'patent_core':pc,'Target':r.Target,'acts':r.acts,'n_rows':int(r.n_rows),'archive':a,'title':title,
            'activity_header':hdr.get((pc,r.Target),''),'n_hits':len(snips),'snippets':snips,'err':err})
    return out,accs,assay_snips,pc,title
if __name__=='__main__':
    pcs=sorted(pairs.patent_core.unique())
    allp=[]; alla=[]; pinfo=[]
    with Pool(22) as p:
        for i,(o,ac,asn,pc,title) in enumerate(p.imap_unordered(work,pcs,chunksize=4)):
            allp+=o; alla+=ac; pinfo.append({'patent_core':pc,'title':title,'assay_sections':asn})
            if i%500==0: print(i,flush=True)
    with open(W+'pair_context.jsonl','w') as f:
        for r in allp: f.write(json.dumps(r,ensure_ascii=False)+'\n')
    with open(W+'patent_assay_sections.jsonl','w') as f:
        for r in pinfo: f.write(json.dumps(r,ensure_ascii=False)+'\n')
    pd.DataFrame(alla,columns=['patent_core','accession','context']).drop_duplicates().to_parquet(W+'patent_accessions.parquet',index=False)
    print('pairs',len(allp),'with hits',sum(1 for r in allp if r['n_hits']),'errors',sum(1 for r in allp if r['err']),'acc mentions',len(alla))
