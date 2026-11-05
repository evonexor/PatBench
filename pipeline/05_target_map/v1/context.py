"""Assemble per-(patent_core, Target) context: title, activity header, assay context text, archive/xml location."""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pyarrow.dataset as ds, pandas as pd, json
W=_ws + '/outputs/target-uniprot-map-20261016/work/'
R=_ws + '/outputs/uspto-sar-recovery-20261006/results/'
F=_ws + '/outputs/uspto-sar-full-20261004/results/'
A=_ws + '/outputs/uspto-all-structures-20261004/results/'
cols=['patent_core','Target','activity_header','target_evidence','archive','xml_file']
c=ds.dataset(F+'candidate_SAR.parquet').to_table(columns=cols+['assay_context','title']).to_pandas()
n=ds.dataset(R+'new_SAR.parquet').to_table(columns=cols+['source_evidence_json','legend_evidence']).to_pandas()
print(n.source_evidence_json.iloc[0][:2000]); print(n.legend_evidence.iloc[0][:500])
e=ds.dataset(A+'SAR_enriched_annotations.parquet').to_table(columns=['patent_core','Target','title','archive','abstract']).to_pandas()
g=pd.read_parquet(W+'ranking_groups.parquet'); rp=g[['patent_core','Target']].drop_duplicates()
def first(df,name):
    return df.drop_duplicates(['patent_core','Target']).set_index(['patent_core','Target'])
ci=first(c,'c'); ni=first(n,'n'); ei=first(e,'e')
rp=rp.set_index(['patent_core','Target'])
print('pairs',len(rp),'in candidate_SAR',rp.index.isin(ci.index).sum(),'in new_SAR',rp.index.isin(ni.index).sum(),'in enriched',rp.index.isin(ei.index).sum())
