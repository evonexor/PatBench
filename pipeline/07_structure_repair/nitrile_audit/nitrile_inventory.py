_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import sys,pathlib,json
sys.path.insert(0,_ws + '/outputs/uspto-sar-demo-20261003/runtime/deps')
import pandas as pd
r=pathlib.Path(_ws + '/outputs/patbench-v1.0-20261030-r6');s=pd.read_parquet(r/'sar_all.parquet');a=s[s.structure_repair.fillna('').str.contains('nitrile')];p=pd.read_parquet(_ws + '/outputs/uspto-sar-recovery-20261006/results/patents.parquet',columns=['Patent','archive']);print('patents fields',p.head().to_string(index=False));g=a.groupby(['Patent','patent_core','structure_file'],dropna=False).agg(compound_alias=('compound_alias','first'),n=('source_record_id','size'),smiles_before=('SMILES_before_repair','first'),smiles_after=('SMILES_as_drawn','first'),rule=('structure_repair','first'),groups=('group_id',lambda x:';'.join(map(str,sorted(set(x)))))).reset_index().merge(p,on='Patent',how='left');
test={};
for f in (r/'tasks').glob('rank_*.parquet'):
 d=pd.read_parquet(f);d=d[(d.split_time=='test')&d.structure_repair.fillna('').str.contains('nitrile')]
 for sf in d.structure_file:test.setdefault(sf,set()).add(f.stem)
g['test_tasks']=g.structure_file.map(lambda x:';'.join(sorted(test.get(x,[]))))
out=pathlib.Path(_ws + '/outputs/nitrile-evidence-audit-20261031');out.mkdir(exist_ok=True);g.to_parquet(out/'inventory.parquet',index=False);print('records',len(g),'patents',g.Patent.nunique(),'test files',g.test_tasks.ne('').sum(),'testpatents',g[g.test_tasks.ne('')].Patent.nunique());print(g[g.test_tasks.ne('')].groupby('Patent').agg(files=('structure_file','size'),archive=('archive','first')).to_string());print(g[g.test_tasks.ne('')].head(3).to_string(index=False))
