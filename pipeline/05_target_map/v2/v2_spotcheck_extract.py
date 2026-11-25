"""Evidence windows from the patent XML for a sample of v2 mappings (same extraction as v1's spot check)."""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pandas as pd, zipfile, re, html, pyarrow.dataset as ds
W=_ws + '/outputs/target-uniprot-map-v2-20261019/work/'
ROOT=_ws + '/raw/patent-data/'
arc=pd.concat([ds.dataset(_ws + '/outputs/uspto-sar-recovery-20261006/results/new_SAR.parquet').to_table(columns=['patent_core','archive']).to_pandas(),
   ds.dataset(_ws + '/outputs/uspto-all-structures-20261004/results/SAR_enriched_annotations.parquet').to_table(columns=['patent_core','archive']).to_pandas()]).drop_duplicates('patent_core')
A=dict(zip(arc.patent_core,arc.archive))
def clean(x):
    x=re.sub(r'<\?.*?\?>|<!--.*?-->','',x,flags=re.S)
    x=re.sub(r'<(entry|td|th)[^>]*>',' | ',x); x=re.sub(r'</(p|heading|row|tr|title|invention-title)>','\n',x)
    x=re.sub(r'<[^>]+>',' ',x); x=html.unescape(x)
    x=re.sub(r'[⁡-⁤​]','',x); x=re.sub(r'[ \t]+',' ',x)
    return re.sub(r'\n\s*','\n',x)
S=pd.read_csv(W+'spotcheck_sample.tsv',sep='\t',dtype=str).fillna('')
cache={}
with open(W+'spotcheck_evidence.txt','w') as out:
    for r in S.itertuples():
        pc=r.patent_core
        if pc not in cache:
            try:
                z=zipfile.ZipFile(ROOT+A[pc]); xs=[n for n in z.namelist() if n.upper().endswith(('.XML','.SGM','.SGML'))]
                xs.sort(key=lambda n:-z.getinfo(n).file_size)
                raw=z.read(xs[0]).decode('utf-8',errors='replace')
                m=re.search(r'<invention-title[^>]*>(.*?)</invention-title>',raw,re.S) or re.search(r'<B540>(.*?)</B540>',raw,re.S)
                cache[pc]=(clean(re.sub(r'<claims.*?</claims>','',raw,flags=re.S)), clean(m.group(1)).strip() if m else '')
            except Exception as ex: cache[pc]=('','ERR '+repr(ex))
        text,title=cache[pc]
        low=text.lower(); t=r.Target.strip().lower()
        hits=[m.start() for m in re.finditer(re.escape(t),low)][:200]
        def sc(p):
            w=text[max(0,p-500):p+500]
            return (w.count('|')<5)*2+bool(re.search(r'assay|inhibit|binding|recombinant|enzyme|receptor|express|cells',w,re.I))
        hits=sorted(hits,key=lambda p:-sc(p))[:2]
        out.write(f"=== {r.item_id} | {pc} | TARGET: {r.Target}\nMAPPED: {r.mapped_identity} {r.uniprot_accession} {r.protein_name} | {r.organism} | conf={r.confidence} | class={r.target_class} | mappability={r.mappability} | endpoints={r.activity_types} mols={r.n_mol}\n")
        out.write(f"METHOD: {r.method}\nXREF: {r.xref_db} {r.xref_target} {r.xref_value_match}\nTITLE: {title[:200]}\nHEADER: {r.activity_header[:200]}\nNOTE: {r.curator_note[:160]}\n")
        for p in sorted(hits): out.write('TEXT: '+re.sub(r'\s+',' ',text[max(0,p-600):p+600])+'\n')
        if not hits: out.write('TEXT: <label not found verbatim in patent body>\n')
        out.write('\n')
print('done',len(S))
