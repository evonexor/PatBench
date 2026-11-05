"""Assemble target_map.parquet from curation TSVs + deterministic signals + sequence resolution.

Confidence policy (per ranking pair):
  high   = single protein / complex / mutant, curator confidence H, sequence resolved, and corroborated
           (independent same-patent BindingDB/ChEMBL curation agrees on the gene symbol, or the label itself
            is an exact unique alias of that symbol in UniProt/ChEMBL)
  medium = resolved sequence but only one line of evidence (H without corroboration, or M with corroboration)
  low    = weak curator confidence, organism fallback for a non-stated organism, or implied-only target
  none   = phenotypic / unclear / unresolved
"""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pandas as pd, numpy as np, re, sys, glob, json, os
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from tmaplib import key, clean_text
from resolve_lib import Resolver, TAX
W=_ws + '/outputs/target-uniprot-map-20261016/work/'
OUT=_ws + '/outputs/target-uniprot-map-20261016/'
P=pd.read_pickle(W+'pair_candidates_ids.pkl')
from read_curation import read as read_curation
C,bad=read_curation(W+'curation')
C=C.drop_duplicates('item_id',keep='first')
print('curation rows',len(C),'unrepairable lines',len(bad))
for b in bad[:10]: print('BAD',b)
print('curation items',len(C),'expected',P.item_id.nunique())
M=P.merge(C,on='item_id',how='left',suffixes=('',''))
for c in ['class','cell','symbols','implied','organism','mutation','conf','note','batch']:
    M[c]=M[c].fillna('').astype(str)
print('pairs without curation',M['class'].isna().sum())
R=Resolver()
U=R.U; acc2sym=dict(zip(U.accession,U.gene_primary.str.upper()))
cs=pd.read_parquet(W+'chembl_comp_synonyms.parquet')
g=cs[cs.syn_type=='GENE_SYMBOL']; acc2sym.update({a:s.upper() for a,s in zip(g.accession,g.component_synonym) if a not in acc2sym})
D=pd.read_parquet(W+'alias_dictionary.parquet')
alias_ids=D.groupby('key').identity.agg(lambda s:set(s)).to_dict()
def parse_syms(s):
    s=str(s or '').strip()
    if s.lower()=='nan': s=''
    return [x.strip() for x in s.split(';') if x.strip()]
def resolve_one(tok,organism):
    if tok.upper().startswith('ORG:'):
        body=tok[4:]
        sp,_,nm=body.partition('|')
        return ('nonmammal',sp.strip(),nm.strip()), R.nonmammal(sp.strip(),nm.strip() or sp.strip())
    return ('mammal',tok.upper(),''), R.mammal(tok,organism)
unresolved={}
M=M.rename(columns={'class':'cls'})
out=[]
for r in M.itertuples():
    cls=r.cls if isinstance(r.cls,str) else ''
    syms=parse_syms(r.symbols); impl=parse_syms(r.implied)
    use=syms if cls in ('protein','complex','mutant') else []
    kind='curated_target' if use else ('implied_target' if impl else '')
    toks=use or (impl if cls in ('phenotypic',) else [])
    primary=toks[0] if toks else ''
    res=None; how=None
    if primary:
        how,res=resolve_one(primary,r.organism)
        if res is None: unresolved[(primary,r.organism)]=unresolved.get((primary,r.organism),0)+1
    # corroboration
    xr=r.xref[0] if r.xref else None
    xr_sym=acc2sym.get(xr['acc']) if xr and xr.get('acc') and '|' not in (xr.get('acc') or '') else None
    xr_strong=bool(xr) and xr['n_vmatch']>=5 and xr['n_vmatch']/max(xr['n_mol'],1)>=0.5
    agree=bool(xr_sym) and bool(primary) and xr_sym==primary.upper()
    exact_alias=False
    for via in ('full','core','org_stripped'):
        for h in (r.dict or []):
            if h['via']==via and primary and len(h['ids'])==1 and h['ids'][0].upper()==primary.upper(): exact_alias=True
    corrob=(agree and xr_strong) or exact_alias
    conf=(r.conf or '').upper()
    if res is None or kind=='':
        tier='none'
    elif kind=='implied_target':
        tier='low'
    elif conf=='H' and corrob: tier='high'
    elif conf=='H' or (conf=='M' and corrob): tier='medium'
    else: tier='low'
    if tier in ('high','medium') and res and res['organism_used']=='human_fallback' and (r.organism or '').strip().lower() not in ('','human','assumed human','homo sapiens'):
        tier='low'
    meth=[]
    if kind: meth.append(kind+'_from_label_and_patent_text')
    if exact_alias: meth.append('exact_alias_match')
    if agree and xr_strong: meth.append('same_patent_db_agreement')
    elif xr_strong: meth.append('same_patent_db_present_no_symbol_agreement')
    if res: meth.append(res['seq_source'])
    out.append(dict(patent_core=r.patent_core,Target=r.Target,activity_types=r.acts,n_rows=int(r.n_rows),
        item_id=r.item_id,target_class=cls or 'uncurated',measured_in_cells=r.cell if isinstance(r.cell,str) else '',
        target_symbols=';'.join(syms),implied_symbols=';'.join(impl),mapped_identity=primary,mapping_kind=kind,
        mutation=r.mutation if isinstance(r.mutation,str) else '',organism_stated=r.organism if isinstance(r.organism,str) else '',
        organism_detected=r.org_stated,
        uniprot_accession=(res or {}).get('accession',''),uniprot_entry=(res or {}).get('entry',''),
        protein_name=(res or {}).get('protein_name',''),organism=(res or {}).get('organism',''),tax_id=(res or {}).get('tax_id',''),
        sequence=(res or {}).get('sequence',''),seq_len=(res or {}).get('seq_len',0),
        organism_resolution=(res or {}).get('organism_used',''),
        confidence=tier,method='+'.join(meth),curator_confidence=conf,curator_note=r.note if isinstance(r.note,str) else '',
        xref_db=(xr or {}).get('db',''),xref_target=(xr or {}).get('target',''),xref_accession=(xr or {}).get('acc',''),
        xref_value_match=f"{xr['n_vmatch']}/{xr['n_mol']}" if xr else '',xref_symbol=xr_sym or '',
        activity_header=r.activity_header,patent_title=r.title,cell_lexicon=';'.join(r.cell_lexicon)))
T=pd.DataFrame(out)
T.to_pickle(W+'target_map_draft.pkl')
print(T.confidence.value_counts().to_dict())
print(T.target_class.value_counts().to_dict())
pd.Series(unresolved).sort_values(ascending=False).to_csv(W+'unresolved_identities.tsv',sep='\t')
print('unresolved identity count',len(unresolved),'pairs',sum(unresolved.values()))
