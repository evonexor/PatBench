_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import sqlite3, pandas as pd
W=_ws + '/outputs/target-uniprot-map-20261016/work/'
con=sqlite3.connect('file:'+W+'chembl/chembl_37/chembl_37_sqlite/chembl_37.db?mode=ro',uri=True)
q={
'chembl_targets':"""select td.tid, td.chembl_id target_chembl_id, td.pref_name, td.target_type, td.organism, td.tax_id,
   cs.component_id, cs.accession, cs.description comp_desc, cs.organism comp_org, cs.tax_id comp_tax
   from target_dictionary td left join target_components tc on td.tid=tc.tid left join component_sequences cs on tc.component_id=cs.component_id""",
'chembl_comp_synonyms':"""select cs.component_id, cs.accession, cs.organism, cs.tax_id, syn.component_synonym, syn.syn_type
   from component_synonyms syn join component_sequences cs on syn.component_id=cs.component_id""",
'chembl_comp_seq':"select component_id, accession, sequence, organism, tax_id, description, db_source from component_sequences",
'chembl_cells':"select cell_id, cell_name, cell_description, cell_source_tissue, cell_source_organism, cellosaurus_id, cell_ontology_id from cell_dictionary",
'chembl_target_synonyms':"select ts.tid, td.chembl_id, ts.target_synonym, ts.syn_type from target_synonyms ts join target_dictionary td on ts.tid=td.tid" ,
}
for k,s in q.items():
    try:
        d=pd.read_sql(s,con); d.to_parquet(W+k+'.parquet',index=False); print(k,d.shape)
    except Exception as ex: print(k,'ERR',ex)
print(pd.read_sql("select name from sqlite_master where type='table'",con).name.tolist())
