"""Integration challenges: detect swapped compound IDs, structures and activities."""
import argparse
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import pyarrow as pa
import pyarrow.parquet as pq

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--results',type=Path,required=True);a=p.parse_args()
    source=pq.read_table(a.results/'SAR_with_evidence.parquet').to_pylist()
    first=source[0];second=next(r for r in source if r['Patent']==first['Patent'] and r['compound_alias']!=first['compound_alias'] and r['SMILES']!=first['SMILES'] and r['raw_value']!=first['raw_value'])
    original=[first,second]
    with tempfile.TemporaryDirectory(prefix='alignment-challenge-') as temp:
        base=Path(temp)
        for case,keys in [('unchanged',[]),('swapped_structures',['SMILES','InChIKey']),('swapped_activity',['raw_value','Value','value_nM']),('swapped_compound_ids',['compound_alias'])]:
            rows=copy.deepcopy(original)
            for key in keys:rows[0][key],rows[1][key]=rows[1][key],rows[0][key]
            input_dir=base/case/'input';output=base/case/'output';input_dir.mkdir(parents=True)
            pq.write_table(pa.Table.from_pylist(rows),input_dir/'SAR_with_evidence.parquet')
            run=subprocess.run([sys.executable,str(Path(__file__).with_name('audit_alignment.py')),'--root',str(a.root),'--results',str(input_dir),'--output',str(output)],capture_output=True,text=True)
            report=json.loads((output/'alignment_report.json').read_text())
            assert (run.returncode==0)==(case=='unchanged'),(case,run.stderr)
            assert report['failed_source_checks']==(0 if case=='unchanged' else 2),(case,report)
            print(case+' passed',flush=True)

if __name__=='__main__':main()
