"""Merge v1 and v2 rows into one target_map covering every (patent_core, Target) pair in SAR_merged,
add source_run and mappability, write the group-level table and the coverage report (overall, by endpoint,
by group size)."""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pandas as pd, numpy as np, json, os, sys, re
V1=_ws + '/outputs/target-uniprot-map-20261016/'
OUT=_ws + '/outputs/target-uniprot-map-v2-20261019/'
W=OUT+'work/'
T2=pd.read_pickle(W+'v2_rows_draft.pkl')
# fold in UniProt REST results fetched for v2 identities
if os.path.exists(W+'uniprot_fetched.tsv'):
    F=pd.read_csv(W+'uniprot_fetched.tsv',sep='\t',dtype=str).fillna('')
    F=F[F.accession!='']; fm={(r.identity,r.organism_req):r for r in F.itertuples()}
    n=0
    for i,r in T2.iterrows():
        if r.uniprot_accession=='' and r.mapped_identity:
            h=fm.get((r.mapped_identity,r.organism_stated)) or fm.get((r.mapped_identity,''))
            if h is not None:
                T2.at[i,'uniprot_accession']=h.accession; T2.at[i,'uniprot_entry']=h.entry
                T2.at[i,'protein_name']=h.protein_name; T2.at[i,'organism']=h.organism; T2.at[i,'tax_id']=h.tax_id
                T2.at[i,'sequence']=h.sequence; T2.at[i,'seq_len']=int(h.seq_len or 0)
                T2.at[i,'organism_resolution']='uniprot_rest_'+(h.reviewed or '')
                T2.at[i,'method']=(r.method+'+UniProt_REST_'+(h.reviewed or '')).strip('+')
                corrob=('same_patent_db_agreement' in r.method) or ('exact_alias_match' in r.method)
                T2.at[i,'confidence']=('high' if (r.curator_confidence=='H' and corrob) else 'medium' if r.curator_confidence=='H' else 'low') if r.mapping_kind=='curated_target' else 'low'
                n+=1
    print('filled from REST',n)
T2.loc[T2.sequence.str.len()==0,'confidence']='none'
T1=pd.read_parquet(V1+'target_map.parquet')
T1['source_run']='v1_2026-10-16'; T2['source_run']='v2_2026-10-19'
T1['curator_mappable']=''   # v1 has no explicit judgment; derived below
P=pd.read_parquet(W+'all_pairs.parquet')
nm=dict(zip(zip(P.patent_core,P.Target),P.n_mol))
ep=dict(zip(zip(P.patent_core,P.Target),P.endpoints))
T1['n_mol']=[nm[(p,t)] for p,t in zip(T1.patent_core,T1.Target)]
T1['activity_types']=[ep[(p,t)] for p,t in zip(T1.patent_core,T1.Target)]  # all endpoints of the pair, not just the ranking ones
# reused v1 string decisions carry no v2 item id; label them with the v1 item they came from
V1S=pd.read_parquet(V1+'target_map.parquet',columns=['Target','item_id'])
V1S=V1S[V1S.item_id.str.startswith('s')].drop_duplicates('Target').set_index('Target').item_id.to_dict()
m=(T2.item_id=='')&T2.method.str.contains('reused_v1_string_decision')
T2.loc[m,'item_id']=['v1:'+V1S.get(t,'?') for t in T2.loc[m,'Target']]
T=pd.concat([T1,T2],ignore_index=True)
assert T.duplicated(['patent_core','Target']).sum()==0, 'duplicate pairs'
def mapstate(r):
    if r.target_class in ('protein','complex','mutant'):
        return 'target_named' if r.uniprot_accession else 'target_named_sequence_unresolved'
    if r.target_class=='phenotypic':
        if r.implied_symbols: return 'intended_target_named_only'
        return {'N':'no_named_target','?':'indeterminate'}.get(r.curator_mappable,'no_named_target')
    if r.target_class=='unclear':
        return {'N':'no_named_target','Y':'intended_target_named_only'}.get(r.curator_mappable,'label_uninformative')
    return 'uncurated'
T['mappability']=[mapstate(r) for r in T.itertuples()]
POLY=re.compile(r'polyprotein|replicase|genome polyprotein|precursor',re.I)
def caveat(r):
    c=[]
    if r.target_class=='mutant' and r.sequence: c.append('wild_type_sequence_for_mutant')
    if r.target_class=='complex' and r.sequence: c.append('complex_primary_subunit_only')
    if r.sequence and str(r.mapped_identity).upper().startswith('ORG:') and (POLY.search(str(r.protein_name)) or r.seq_len>1500):
        c.append('polyprotein_or_precursor_sequence')
    return ';'.join(c)
T['sequence_caveat']=[caveat(r) for r in T.itertuples()]
cols=['patent_core','Target','activity_types','n_mol','n_rows','source_run','target_class','mappability','measured_in_cells',
 'mapped_identity','mapping_kind','target_symbols','implied_symbols','mutation','organism_stated','organism_detected',
 'uniprot_accession','uniprot_entry','protein_name','organism','tax_id','seq_len','sequence','organism_resolution',
 'confidence','sequence_caveat','method','curator_confidence','curator_mappable','curator_note','xref_db','xref_target','xref_accession',
 'xref_symbol','xref_value_match','activity_header','patent_title','cell_lexicon','item_id']
T=T[cols]
T.to_parquet(OUT+'target_map.parquet',index=False)
G=pd.read_parquet(W+'all_groups.parquet')
K=G.merge(T.drop(columns=['n_mol','n_rows','activity_types']),on=['patent_core','Target'],how='left')
assert K.confidence.notna().all()
K.to_parquet(OUT+'group_targets.parquet',index=False)
hi=K.confidence=='high'; hm=K.confidence.isin(['high','medium'])
def blk(d):
    return dict(groups=int(len(d)),rows=int(d.n_rows.sum()),patents=int(d.patent_core.nunique()),
                accessions=int(d[d.uniprot_accession!=''].uniprot_accession.nunique()))
stats=dict(
 total=blk(K), high=blk(K[hi]), high_medium=blk(K[hm]),
 pairs_total=int(len(T)), pairs_by_run={k:int(v) for k,v in T.source_run.value_counts().items()},
 pairs_high_medium=int(T.confidence.isin(['high','medium']).sum()),
 patents_in_sar=int(K.patent_core.nunique()),
 patents_with_any_high_medium_pair=int(K[hm].patent_core.nunique()),
 by_confidence={k:blk(K[K.confidence==k]) for k in ['high','medium','low','none']},
 by_class={k:blk(K[K.target_class==k]) for k in sorted(K.target_class.unique())},
 by_mappability={k:blk(K[K.mappability==k]) for k in sorted(K.mappability.unique())},
 cell_based=blk(K[K.measured_in_cells=='Y']),
 by_endpoint={e:dict(total=blk(K[K['Activity Type']==e]),high_medium=blk(K[(K['Activity Type']==e)&hm]))
              for e in K['Activity Type'].value_counts().index},
)
bins=[(1,1,'1'),(2,2,'2'),(3,4,'3-4'),(5,9,'5-9'),(10,19,'10-19'),(20,49,'20-49'),(50,10**9,'50+')]
stats['by_group_size']={lab:dict(total=blk(K[(K.n_mol>=lo)&(K.n_mol<=hi_)]),high_medium=blk(K[(K.n_mol>=lo)&(K.n_mol<=hi_)&hm]))
                        for lo,hi_,lab in bins}
json.dump(stats,open(OUT+'coverage_stats.json','w'),indent=1)
print(json.dumps({k:stats[k] for k in ['total','high','high_medium','pairs_total','pairs_by_run','patents_in_sar','patents_with_any_high_medium_pair']},indent=1))
print('--- by endpoint (groups h+m / total, rows h+m / total)')
for e,v in stats['by_endpoint'].items():
    print(f"{e:6s} {v['high_medium']['groups']:6d}/{v['total']['groups']:6d}  {v['high_medium']['rows']:7d}/{v['total']['rows']:7d}  patents {v['high_medium']['patents']}/{v['total']['patents']}")
print('--- by group size')
for s,v in stats['by_group_size'].items():
    print(f"{s:6s} {v['high_medium']['groups']:6d}/{v['total']['groups']:6d}  {v['high_medium']['rows']:7d}/{v['total']['rows']:7d}")
print('--- mappability'); [print(f"{k:34s} {v['groups']:6d} groups {v['rows']:7d} rows {v['patents']:5d} patents") for k,v in stats['by_mappability'].items()]
# agreement between explicit judgment and derivation, on the v2 half only
v2=T[T.source_run.str.startswith('v2')]
x=v2[v2.target_class.isin(['phenotypic','unclear'])]
print('--- sequence caveats'); print(T[T.sequence_caveat!=''].sequence_caveat.value_counts().to_string())
print('--- v2 phenotypic/unclear: explicit mappable vs derived state')
print(pd.crosstab(x.curator_mappable,x.mappability).to_string())
