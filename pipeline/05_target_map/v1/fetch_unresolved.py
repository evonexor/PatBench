"""Resolve identities the local reference could not cover, via targeted UniProt REST searches.
One query per distinct identity; reviewed entries preferred, then unreviewed; results cached to work/uniprot_fetched.tsv."""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pandas as pd, requests, time, sys, os, re, json
W=_ws + '/outputs/target-uniprot-map-20261016/work/'
CACHE=W+'uniprot_fetched.tsv'
FIELDS='accession,id,protein_name,gene_primary,organism_name,organism_id,length,reviewed,sequence'
def q(query,size=5):
    r=requests.get('https://rest.uniprot.org/uniprotkb/search',params={'query':query,'fields':FIELDS,'format':'tsv','size':size},timeout=60)
    r.raise_for_status()
    lines=r.text.strip().split('\n')
    if len(lines)<2: return []
    hdr=lines[0].split('\t')
    return [dict(zip(hdr,l.split('\t'))) for l in lines[1:]]
def fetch(ident,organism):
    if ident.upper().startswith('ORG:'):
        sp,_,nm=ident[4:].partition('|'); sp=sp.strip(); nm=(nm or sp).strip()
        tries=[f'(organism_name:"{sp}") AND ("{nm}") AND (reviewed:true)',
               f'(organism_name:"{sp}") AND ("{nm}")',
               f'("{sp}") AND ("{nm}") AND (reviewed:true)']
    else:
        org=(organism or 'human').strip() or 'human'
        orgq={'human':'9606','assumed human':'9606','mouse':'10090','rat':'10116'}.get(org.lower())
        oq=f'(organism_id:{orgq})' if orgq else f'(organism_name:"{org}")'
        tries=[f'(gene_exact:{ident}) AND {oq} AND (reviewed:true)',
               f'(gene:{ident}) AND {oq} AND (reviewed:true)',
               f'("{ident}") AND {oq} AND (reviewed:true)']
    for t in tries:
        try: hits=q(t)
        except Exception as ex: print('ERR',t,ex,flush=True); time.sleep(3); continue
        if hits:
            hits.sort(key=lambda h:(h.get('Reviewed')!='reviewed',-int(h.get('Length') or 0)))
            h=hits[0]
            return dict(identity=ident,organism_req=organism,accession=h['Entry'],entry=h.get('Entry Name',''),
                protein_name=h.get('Protein names',''),gene=h.get('Gene Names (primary)',''),organism=h.get('Organism',''),
                tax_id=h.get('Organism (ID)',''),seq_len=h.get('Length',''),reviewed=h.get('Reviewed',''),
                sequence=h.get('Sequence',''),query=t)
        time.sleep(0.2)
    return dict(identity=ident,organism_req=organism,accession='',entry='',protein_name='',gene='',organism='',tax_id='',seq_len='',reviewed='',sequence='',query='NOHIT')
if __name__=='__main__':
    u=pd.read_csv(W+'unresolved_identities.tsv',sep='\t',header=0,names=['identity','organism','n'],dtype=str)
    u=u[u.identity.notna()&(u.identity!='identity')]
    done={}
    if os.path.exists(CACHE):
        d=pd.read_csv(CACHE,sep='\t',dtype=str).fillna('')
        done={(r.identity,r.organism_req) for r in d.itertuples()}
        res=d.to_dict('records')
    else: res=[]
    print('to fetch',len(u))
    for i,r in enumerate(u.itertuples()):
        if (r.identity,r.organism or '') in done: continue
        res.append(fetch(r.identity,r.organism or ''))
        if i%20==0:
            pd.DataFrame(res).to_csv(CACHE,sep='\t',index=False); print(i,r.identity,res[-1]['accession'],flush=True)
        time.sleep(0.15)
    pd.DataFrame(res).to_csv(CACHE,sep='\t',index=False)
    d=pd.DataFrame(res); print('fetched',len(d),'resolved',(d.accession!='').sum())
