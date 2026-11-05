"""Stream BindingDB_Patents TSV and keep the columns needed for patent cross-referencing."""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import zipfile, io, csv, sys, pandas as pd, re
csv.field_size_limit(sys.maxsize)
Z=_ws + '/raw/patent-data/bindingdb/202609/BindingDB_Patents_202609_tsv.zip'
OUT=_ws + '/outputs/target-uniprot-map-20261016/work/bindingdb_patents_slim.parquet'
zf=zipfile.ZipFile(Z); fh=io.TextIOWrapper(zf.open('BindingDB_Patents.tsv'),encoding='utf-8',errors='replace',newline='')
rd=csv.reader(fh,delimiter='\t',quoting=csv.QUOTE_NONE)
hdr=next(rd); ix={h:i for i,h in enumerate(hdr)}
keep={'rsid':'BindingDB Reactant_set_id','inchikey':'Ligand InChI Key','ligand_name':'BindingDB Ligand Name',
 'target_name':'Target Name','organism':'Target Source Organism According to Curator or DataSource',
 'Ki':'Ki (nM)','IC50':'IC50 (nM)','Kd':'Kd (nM)','EC50':'EC50 (nM)','patent':'Patent Number',
 'n_chains':'Number of Protein Chains in Target (>1 implies a multichain complex)',
 'seq1':'BindingDB Target Chain Sequence 1','sp_name1':'UniProt (SwissProt) Recommended Name of Target Chain 1',
 'sp_entry1':'UniProt (SwissProt) Entry Name of Target Chain 1','sp_id1':'UniProt (SwissProt) Primary ID of Target Chain 1',
 'tr_id1':'UniProt (TrEMBL) Primary ID of Target Chain 1','tr_name1':'UniProt (TrEMBL) Submitted Name of Target Chain 1'}
chain_sp=[(k,ix[f'UniProt (SwissProt) Primary ID of Target Chain {k}'],ix[f'UniProt (TrEMBL) Primary ID of Target Chain {k}']) for k in range(1,51) if f'UniProt (SwissProt) Primary ID of Target Chain {k}' in ix]
rows=[]; n=0
for r in rd:
    n+=1
    if len(r)<len(hdr): r=r+['']*(len(hdr)-len(r))
    d={k:r[ix[v]] for k,v in keep.items()}
    ch=[]
    for k,a,b in chain_sp:
        v=r[a] or r[b]
        if v: ch.append(v)
    d['chain_ids']='|'.join(ch)
    rows.append(d)
    if n%200000==0: print(n,flush=True)
df=pd.DataFrame(rows); print(df.shape)
df.to_parquet(OUT,index=False)
