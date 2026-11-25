"""Second pass for identities the first UniProt query missed: try the name without parentheses,
the parenthetical alias on its own, and any gene-like token (UL54, leuS, rpoA, CLAG3) as a gene search."""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pandas as pd, re, time, sys
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from fetch_unresolved import q, CACHE, W
d=pd.read_csv(CACHE,sep='\t',dtype=str).fillna('')
todo=d[d.accession=='']
print('retry',len(todo))
def variants(ident):
    sp,_,nm=ident[4:].partition('|') if ident.upper().startswith('ORG:') else ('',None,ident)
    sp=sp.strip(); nm=(nm or ident).strip()
    outer=re.sub(r'\([^)]*\)','',nm).strip()
    inner=[m.group(1).strip() for m in re.finditer(r'\(([^)]*)\)',nm)]
    genes=[t for t in re.findall(r'\b([A-Z]{2}\d{2,3}|[a-z]{3}[A-Z]|[A-Z]{3,5}\d)\b',nm)]
    vs=[]
    if sp:
        oq=f'(organism_name:"{sp}")'
        for n in [outer]+inner:
            if n: vs.append(f'{oq} AND ("{n}") AND (reviewed:true)')
        for g in genes: vs.append(f'{oq} AND (gene:{g}) AND (reviewed:true)')
        for n in [outer]+inner:
            if n: vs.append(f'{oq} AND ("{n}")')
        for g in genes: vs.append(f'{oq} AND (gene:{g})')
    else:
        for n in [outer]+inner:
            if n: vs.append(f'("{n}") AND (reviewed:true)')
    return vs
fix=0
for i,r in todo.iterrows():
    for t in variants(r.identity):
        try: hits=q(t)
        except Exception as ex: print('ERR',ex); time.sleep(2); continue
        if hits:
            hits.sort(key=lambda h:(h.get('Reviewed')!='reviewed',-int(h.get('Length') or 0)))
            h=hits[0]
            d.loc[i,['accession','entry','protein_name','gene','organism','tax_id','seq_len','reviewed','sequence','query']]=[
                h['Entry'],h.get('Entry Name',''),h.get('Protein names',''),h.get('Gene Names (primary)',''),h.get('Organism',''),
                h.get('Organism (ID)',''),h.get('Length',''),h.get('Reviewed',''),h.get('Sequence',''),t]
            print('FIXED',r.identity,'->',h['Entry'],h.get('Organism','')[:40],h.get('Length'),flush=True); fix+=1
            break
        time.sleep(0.2)
    else:
        print('STILL UNRESOLVED',r.identity,flush=True)
d.to_csv(CACHE,sep='\t',index=False)
print('fixed',fix,'remaining',(d.accession=='').sum())
