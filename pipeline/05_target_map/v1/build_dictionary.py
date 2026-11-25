"""Alias dictionary: normalized key -> (identity, accession, organism, source, kind).
identity = uppercase gene symbol for mammalian proteins (orthologs collapse), else accession."""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pandas as pd, re, gzip, sys
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from tmaplib import key, clean_text
W=_ws + '/outputs/target-uniprot-map-20261016/work/'
MAMM={9606,10090,10116,9615,9544,9541,10141,9823,9913,9986,10036}
rows=[]
up=[]
for tax in ['9606','10090','10116']:
    d=pd.read_csv(W+f'uniprot/sprot_{tax}.tsv.gz',sep='\t',dtype=str).fillna('')
    up.append(d)
UP=pd.concat(up,ignore_index=True)
UP.columns=['accession','entry','gene_names','gene_primary','protein_names','tax','organism','length','ec','sequence']
UP.to_parquet(W+'uniprot_sprot_hmr.parquet',index=False)
def parse_names(pn):
    # "Rec (Alt1) (Alt2) [Cleaved into: X (Y)] [Includes: ...]"
    pn=re.sub(r'\[(Cleaved into|Includes):.*?\]$','',pn).strip()
    names=[]; depth=0; buf=''; top=''
    for ch in pn:
        if ch=='(':
            if depth==0: top=top or buf.strip(); buf=''
            else: buf+=ch
            depth+=1
        elif ch==')':
            depth-=1
            if depth==0: names.append(buf.strip()); buf=''
            else: buf+=ch
        else: buf+=ch
    if not top: top=buf.strip()
    alts=[n for n in names if n and not n.startswith('EC ')]
    return top,alts
for r in UP.itertuples():
    tax=int(r.tax); ident=r.gene_primary.upper() or r.accession
    if r.gene_primary: rows.append((r.gene_primary,ident,r.accession,tax,'UniProt','gene_primary'))
    for gname in r.gene_names.split():
        if gname!=r.gene_primary: rows.append((gname,ident,r.accession,tax,'UniProt','gene_synonym'))
    top,alts=parse_names(r.protein_names)
    rows.append((top,ident,r.accession,tax,'UniProt','protein_name'))
    for a in alts: rows.append((a,ident,r.accession,tax,'UniProt','protein_alt'))
# ChEMBL component synonyms (all organisms)
cs=pd.read_parquet(W+'chembl_comp_synonyms.parquet')
sym=cs[cs.syn_type=='GENE_SYMBOL'].groupby('component_id').component_synonym.first().str.upper()
for r in cs.itertuples():
    tax=int(r.tax_id) if pd.notna(r.tax_id) else 0
    ident=(sym.get(r.component_id) if tax in MAMM else None) or r.accession
    kind={'GENE_SYMBOL':'gene_primary','GENE_SYMBOL_OTHER':'gene_synonym','UNIPROT':'protein_alt','EC_NUMBER':None}[r.syn_type]
    if kind: rows.append((r.component_synonym,ident,r.accession,tax,'ChEMBL37',kind))
# ChEMBL single-protein target pref names
ct=pd.read_parquet(W+'chembl_targets.parquet')
sp=ct[(ct.target_type=='SINGLE PROTEIN')&ct.accession.notna()]
csym=dict(zip(cs[cs.syn_type=='GENE_SYMBOL'].accession,cs[cs.syn_type=='GENE_SYMBOL'].component_synonym.str.upper()))
for r in sp.itertuples():
    tax=int(r.tax_id) if pd.notna(r.tax_id) else 0
    ident=(csym.get(r.accession) if tax in MAMM else None) or r.accession
    rows.append((r.pref_name,ident,r.accession,tax,'ChEMBL37','target_pref_name'))
# BindingDB target names (single chain) -> accession
b=pd.read_parquet(W+'bindingdb_patents_slim.parquet',columns=['target_name','organism','chain_ids','n_chains','sp_entry1'])
b=b[(b.chain_ids!='')&(~b.chain_ids.str.contains(r'\|'))]
bt=b.groupby(['target_name','chain_ids','organism']).size().reset_index(name='n')
acc2sym=dict(zip(UP.accession,UP.gene_primary.str.upper())); acc2sym.update({k:v for k,v in csym.items() if k not in acc2sym})
orgtax={'Homo sapiens':9606,'Mus musculus':10090,'Rattus norvegicus':10116}
for r in bt.itertuples():
    tax=orgtax.get(r.organism,0)
    ident=acc2sym.get(r.chain_ids) if tax in MAMM or r.chain_ids in acc2sym else None
    rows.append((r.target_name,ident or r.chain_ids,r.chain_ids,tax,'BindingDB_202609','bindingdb_target_name'))
D=pd.DataFrame(rows,columns=['alias','identity','accession','tax','source','kind'])
D['alias']=D.alias.map(clean_text); D['key']=D.alias.map(key)
D=D[(D.key.str.len()>=2)&(D.identity!='')].drop_duplicates()
D.to_parquet(W+'alias_dictionary.parquet',index=False)
print(D.shape, D.source.value_counts().to_dict())
# identity -> preferred accession per organism
sp_h=UP[UP.tax=='9606'].groupby(UP.gene_primary.str.upper()).accession.apply(list)
print('human identities',len(sp_h), 'multi-acc symbols',(sp_h.str.len()>1).sum())
for k in ['jak2','mglur5','cb1','pi3kdelta','herg','ox1','ox1r','d3','tlr7','mpro','dpp1','dppi','hptpbeta','flap','mcl1','fasn','dgkalpha','5ht2a','mao a','maoa']:
    kk=key(k); x=D[D.key==kk]
    print(k,'->',sorted(set(x.identity))[:8])
