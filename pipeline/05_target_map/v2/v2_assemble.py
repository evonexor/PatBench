"""Assemble the v2 (new-pair) mapping rows: curation + reused v1 string decisions, sequence resolution,
confidence tiering. Identical scheme and tier rules to v1; adds the explicit mappable judgment."""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pandas as pd, numpy as np, sys, os
V1=_ws + '/outputs/target-uniprot-map-20261016/'
sys.path.insert(0,__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', 'v1')); sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from tmaplib import key
from resolve_lib import Resolver
from v2_read_curation import read as read_v2
W=_ws + '/outputs/target-uniprot-map-v2-20261019/work/'
P=pd.read_pickle(W+'pair_candidates_new_ids.pkl')
C,bad=read_v2(W+'curation'); C=C.drop_duplicates('item_id',keep='first')
print('curation items',len(C),'unrepairable',len(bad))
for b in bad[:10]: print('BAD',b)
V=pd.read_parquet(V1+'target_map.parquet')
v1s=V[V.item_id.str.startswith('s')].drop_duplicates('Target').set_index('Target')
M=P.merge(C,on='item_id',how='left')
for c in ['class','cell','symbols','implied','organism','mutation','conf','note','mappable']:
    M[c]=M[c].fillna('').astype(str).replace('nan','')
M=M.rename(columns={'class':'cls'})
print('new pairs',len(M),'curated',(M.cls!='').sum(),'reuse_v1',M.reuse_v1.sum(),'uncovered',((M.cls=='')&~M.reuse_v1).sum())
# fill reused rows from the v1 string-level decision
for i,r in M[M.reuse_v1].iterrows():
    v=v1s.loc[r.Target]
    M.at[i,'cls']=v.target_class; M.at[i,'cell']=v.measured_in_cells
    M.at[i,'symbols']=v.target_symbols; M.at[i,'implied']=v.implied_symbols
    M.at[i,'organism']=v.organism_stated; M.at[i,'mutation']=v.mutation
    M.at[i,'conf']=v.curator_confidence; M.at[i,'note']=('reused v1 string decision; '+str(v.curator_note))[:200]
    M.at[i,'mappable']='Y' if v.target_class in ('protein','complex','mutant') or v.implied_symbols else ''
R=Resolver()
U=R.U; acc2sym=dict(zip(U.accession,U.gene_primary.str.upper()))
cs=pd.read_parquet(V1+'work/chembl_comp_synonyms.parquet'); g=cs[cs.syn_type=='GENE_SYMBOL']
acc2sym.update({a:s.upper() for a,s in zip(g.accession,g.component_synonym) if a not in acc2sym})
def parse_syms(s):
    s=str(s or '').strip()
    return [x.strip() for x in s.split(';') if x.strip()] if s.lower()!='nan' else []
def resolve_one(tok,organism):
    if tok.upper().startswith('ORG:'):
        sp,_,nm=tok[4:].partition('|')
        return R.nonmammal(sp.strip(),(nm or sp).strip())
    return R.mammal(tok,organism)
out=[]; unresolved={}
for r in M.itertuples():
    cls=r.cls; syms=parse_syms(r.symbols); impl=parse_syms(r.implied)
    use=syms if cls in ('protein','complex','mutant') else []
    kind='curated_target' if use else ('implied_target' if impl and cls in ('phenotypic','unclear') else '')
    toks=use or (impl if cls in ('phenotypic','unclear') else [])
    primary=toks[0] if toks else ''
    res=resolve_one(primary,r.organism) if primary else None
    if primary and res is None: unresolved[(primary,r.organism)]=unresolved.get((primary,r.organism),0)+1
    xr=r.xref[0] if r.xref else None
    xr_sym=acc2sym.get(xr['acc']) if xr and xr.get('acc') and '|' not in (xr.get('acc') or '') else None
    xr_strong=bool(xr) and xr['n_vmatch']>=5 and xr['n_vmatch']/max(xr['n_mol'],1)>=0.5
    agree=bool(xr_sym) and bool(primary) and xr_sym==primary.upper()
    exact_alias=any(h['via'] in ('full','core','org_stripped') and len(h['ids'])==1 and primary and h['ids'][0].upper()==primary.upper() for h in (r.dict or []))
    corrob=(agree and xr_strong) or exact_alias
    conf=(r.conf or '').upper()
    if res is None or kind=='': tier='none'
    elif kind=='implied_target': tier='low'
    elif conf=='H' and corrob: tier='high'
    elif conf=='H' or (conf=='M' and corrob): tier='medium'
    else: tier='low'
    if tier in ('high','medium') and res and res['organism_used']=='human_fallback' and (r.organism or '').strip().lower() not in ('','human','assumed human','homo sapiens'):
        tier='low'
    meth=[]
    if kind: meth.append(kind+'_from_label_and_patent_text')
    if r.reuse_v1: meth.append('reused_v1_string_decision')
    if exact_alias: meth.append('exact_alias_match')
    if agree and xr_strong: meth.append('same_patent_db_agreement')
    elif xr_strong: meth.append('same_patent_db_present_no_symbol_agreement')
    elif xr: meth.append('same_patent_db_molecule_overlap_only')
    if res: meth.append(res['seq_source'])
    out.append(dict(patent_core=r.patent_core,Target=r.Target,activity_types=r.endpoints,n_rows=int(r.n_rows),n_mol=int(r.n_mol),
        item_id=r.item_id,target_class=cls or 'uncurated',measured_in_cells=r.cell,
        target_symbols=';'.join(syms),implied_symbols=';'.join(impl),mapped_identity=primary,mapping_kind=kind,
        mutation=r.mutation,organism_stated=r.organism,organism_detected=r.org_stated,
        uniprot_accession=(res or {}).get('accession',''),uniprot_entry=(res or {}).get('entry',''),
        protein_name=(res or {}).get('protein_name',''),organism=(res or {}).get('organism',''),tax_id=(res or {}).get('tax_id',''),
        sequence=(res or {}).get('sequence',''),seq_len=(res or {}).get('seq_len',0),
        organism_resolution=(res or {}).get('organism_used',''),confidence=tier,method='+'.join(meth),
        curator_confidence=conf,curator_note=r.note,curator_mappable=r.mappable,
        xref_db=(xr or {}).get('db',''),xref_target=(xr or {}).get('target',''),xref_accession=(xr or {}).get('acc',''),
        xref_symbol=xr_sym or '',xref_value_match=f"{xr['n_vmatch']}/{xr['n_mol']}" if xr else '',
        activity_header=r.activity_header,patent_title=r.title,cell_lexicon=';'.join(r.cell_lexicon)))
T=pd.DataFrame(out); T.to_pickle(W+'v2_rows_draft.pkl')
print(T.confidence.value_counts().to_dict()); print(T.target_class.value_counts().to_dict())
print('mappable', T.curator_mappable.value_counts().to_dict())
pd.Series(unresolved).sort_values(ascending=False).to_csv(W+'unresolved_identities.tsv',sep='\t')
print('unresolved identities',len(unresolved),'pairs',sum(unresolved.values()))
