_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pyarrow.dataset as ds, pandas as pd
pd.set_option('display.width',250); pd.set_option('display.max_colwidth',120)
R=_ws + '/outputs/uspto-sar-recovery-20261006/results/'
W=_ws + '/outputs/target-uniprot-map-20261016/work/'
m=ds.dataset(R+'SAR_merged.parquet').to_table(columns=['source','patent_core','Target','Activity Type','InChIKey','organism']).to_pandas()
A=m[m['Activity Type'].isin(['IC50','Ki','Kd','EC50'])]
g=A.groupby(['patent_core','Target','Activity Type']).agg(n_mol=('InChIKey','nunique'),n_rows=('InChIKey','size')).reset_index()
g10=g[g.n_mol>=10].copy()
g10.to_parquet(W+'ranking_groups.parquet')
t=g10.groupby('Target').agg(groups=('n_mol','size'),rows=('n_rows','sum'),patents=('patent_core','nunique')).sort_values('rows',ascending=False)
t.to_csv(W+'ranking_targets.tsv',sep='\t')
print(len(g10), len(t), g10.n_rows.sum(), g10.patent_core.nunique())
print(t.head(150).to_string())
print('empty target groups', (g10.Target=='').sum())
