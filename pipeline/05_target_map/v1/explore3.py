_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pyarrow.dataset as ds, pandas as pd, json
pd.set_option('display.width',250); pd.set_option('display.max_colwidth',200)
R=_ws + '/outputs/uspto-sar-recovery-20261006/results/'
n=ds.dataset(R+'new_SAR.parquet').to_table(columns=['patent_core','Target','Activity Type','target_evidence','target_scope','UniProt ID','Sequence','activity_header','organism']).to_pandas()
print(n.target_scope.value_counts().head(20))
print((n['UniProt ID']!='').sum(), (n['Sequence']!='').sum())
print(n[n['UniProt ID']!=''][['patent_core','Target','UniProt ID']].drop_duplicates().head(20))
s=n.drop_duplicates(['patent_core','Target']).sample(12,random_state=1)
for _,r in s.iterrows():
    print('====',r.patent_core,'|',r.Target,'|',r.target_scope,'|',r.activity_header[:150])
    print(r.target_evidence[:700])
A=_ws + '/outputs/uspto-all-structures-20261004/results/'
e=ds.dataset(A+'SAR_enriched_annotations.parquet').to_table(columns=['patent_core','Target','UniProt ID','Sequence','annotation_evidence','title']).to_pandas()
print('baseline UniProt ID nonempty', (e['UniProt ID']!='').sum(), (e['Sequence']!='').sum())
print(e.annotation_evidence.iloc[0][:1500])
u=ds.dataset(A+'patent_uniprot_mentions.parquet').to_table().to_pandas()
print(u.link_status.value_counts()); print(u.sample(8,random_state=2)[['Patent','UniProt ID','mention_text']].to_string()[:3000])
