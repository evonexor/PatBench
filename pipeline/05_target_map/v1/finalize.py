"""Merge fetched UniProt results into the draft map, write target_map.parquet + coverage stats JSON."""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pandas as pd, numpy as np, os, json, sys
W=_ws + '/outputs/target-uniprot-map-20261016/work/'
OUT=_ws + '/outputs/target-uniprot-map-20261016/'
T=pd.read_pickle(W+'target_map_draft.pkl')
if os.path.exists(W+'uniprot_fetched.tsv'):
    F=pd.read_csv(W+'uniprot_fetched.tsv',sep='\t',dtype=str).fillna('')
    F=F[F.accession!='']
    fm={(r.identity,r.organism_req):r for r in F.itertuples()}
    n=0
    for i,r in T.iterrows():
        if r.uniprot_accession=='' and r.mapped_identity:
            h=fm.get((r.mapped_identity,r.organism_stated)) or fm.get((r.mapped_identity,''))
            if h is not None:
                T.at[i,'uniprot_accession']=h.accession; T.at[i,'uniprot_entry']=h.entry
                T.at[i,'protein_name']=h.protein_name; T.at[i,'organism']=h.organism; T.at[i,'tax_id']=h.tax_id
                T.at[i,'sequence']=h.sequence; T.at[i,'seq_len']=int(h.seq_len or 0)
                T.at[i,'organism_resolution']='uniprot_rest_'+(h.reviewed or '')
                T.at[i,'method']=(r.method+'+UniProt_REST_'+(h.reviewed or '')).strip('+')
                # re-tier: REST-resolved keeps curator-based tier but caps at medium unless corroborated
                corrob='same_patent_db_agreement' in r.method or 'exact_alias_match' in r.method
                T.at[i,'confidence']=('high' if (r.curator_confidence=='H' and corrob) else 'medium' if r.curator_confidence=='H' else 'low') if r.mapping_kind=='curated_target' else 'low'
                n+=1
    print('filled from REST',n)
T['has_sequence']=T.sequence.str.len()>0
T.loc[~T.has_sequence,'confidence']=T.loc[~T.has_sequence,'confidence'].where(T.loc[~T.has_sequence,'confidence']=='none','none')
g=pd.read_parquet(W+'ranking_groups.parquet')
K=g.merge(T[['patent_core','Target','target_class','confidence','uniprot_accession','mapped_identity','mapping_kind','measured_in_cells','seq_len']],on=['patent_core','Target'],how='left')
assert K.confidence.notna().all(), K.confidence.isna().sum()
cols=['patent_core','Target','activity_types','n_rows','target_class','measured_in_cells','mapped_identity','mapping_kind',
 'target_symbols','implied_symbols','mutation','organism_stated','organism_detected','uniprot_accession','uniprot_entry',
 'protein_name','organism','tax_id','seq_len','sequence','organism_resolution','confidence','method','curator_confidence',
 'curator_note','xref_db','xref_target','xref_accession','xref_symbol','xref_value_match','activity_header','patent_title','cell_lexicon','item_id']
T[cols].to_parquet(OUT+'target_map.parquet',index=False)
# group-level table (one row per ranking group)
K.to_parquet(OUT+'ranking_group_targets.parquet',index=False)
hi=K.confidence.isin(['high']); hm=K.confidence.isin(['high','medium'])
prot=K.target_class.isin(['protein','complex','mutant'])
stats=dict(
 ranking_groups=int(len(K)),ranking_rows=int(K.n_rows.sum()),ranking_patents=int(K.patent_core.nunique()),
 distinct_target_strings=int(K.Target.nunique()),distinct_pairs=int(len(T)),
 groups_high=int(hi.sum()),rows_high=int(K.loc[hi,'n_rows'].sum()),patents_high=int(K.loc[hi,'patent_core'].nunique()),
 groups_high_medium=int(hm.sum()),rows_high_medium=int(K.loc[hm,'n_rows'].sum()),patents_high_medium=int(K.loc[hm,'patent_core'].nunique()),
 accessions_high=int(K.loc[hi,'uniprot_accession'].nunique()),accessions_high_medium=int(K.loc[hm,'uniprot_accession'].nunique()),
 by_confidence={k:int(v) for k,v in K.confidence.value_counts().items()},
 rows_by_confidence={k:int(v) for k,v in K.groupby('confidence').n_rows.sum().items()},
 by_class={k:int(v) for k,v in K.target_class.value_counts().items()},
 rows_by_class={k:int(v) for k,v in K.groupby('target_class').n_rows.sum().items()},
 groups_phenotypic=int((K.target_class=='phenotypic').sum()),rows_phenotypic=int(K.loc[K.target_class=='phenotypic','n_rows'].sum()),
 groups_unclear=int((K.target_class=='unclear').sum()),rows_unclear=int(K.loc[K.target_class=='unclear','n_rows'].sum()),
 groups_cell_based=int((K.measured_in_cells=='Y').sum()),rows_cell_based=int(K.loc[K.measured_in_cells=='Y','n_rows'].sum()),
 protein_class_groups=int(prot.sum()),protein_class_rows=int(K.loc[prot,'n_rows'].sum()),
 protein_class_groups_with_sequence=int((prot&K.uniprot_accession.fillna('').ne('')).sum()),
 cell_based_among_protein_mapped=int((prot&(K.measured_in_cells=='Y')).sum()),
)
json.dump(stats,open(OUT+'coverage_stats.json','w'),indent=1)
print(json.dumps(stats,indent=1))
print(T.method.value_counts().head(15))
