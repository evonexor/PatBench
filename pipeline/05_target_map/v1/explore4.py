_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pyarrow.dataset as ds, pandas as pd, pyarrow.parquet as pq
pd.set_option('display.width',250); pd.set_option('display.max_colwidth',300)
A=_ws + '/outputs/uspto-all-structures-20261004/results/'
e=ds.dataset(A+'SAR_enriched_annotations.parquet').to_table(columns=['patent_core','Target','annotation_evidence','activity_locator','archive','title','abstract','technical_field','source_record_id']).to_pandas()
x=e.drop_duplicates(['patent_core','Target']).sample(4,random_state=3)
for _,r in x.iterrows():
    print('=====',r.patent_core,r.Target,'|',r.title[:100],'|',r.archive,'|',r.activity_locator,'|',r.source_record_id)
    print(r.annotation_evidence[:800])
print(e.title.eq('').mean())
F=_ws + '/outputs/uspto-sar-full-20261004/results/'
f=ds.dataset(F+'candidate_SAR.parquet'); print(f.schema); 
print(f.head(2).to_pandas().T.to_string()[:4000])
