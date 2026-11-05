"""Compare sample exports across worker counts and checkpoint reuse."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--base',type=Path,required=True);a=p.parse_args()
    files=['SAR.csv','SAR.parquet','SAR_enriched.csv','SAR_enriched.parquet','SAR_with_evidence.csv','review.csv','structures.csv']
    files += [p.relative_to(a.base/'results').as_posix() for folder in ['bibliography'] for p in sorted((a.base/'results'/folder).glob('*.json'))]
    def digest(out):return {n:hashlib.sha256((out/n).read_bytes()).hexdigest() for n in files}
    def stats(out):
        r=json.loads((out/'run_report.json').read_text())
        return {k:r[k] for k in ['workers','elapsed_seconds','cache_hits','version']}
    baseline=digest(a.base/'results');runs=[stats(a.base/'results')]
    for label,workers,out in [('1worker-fresh',1,a.base/'benchmark-1/results'),('20worker-cache',20,a.base/'results')]:
        with (a.base/(label+'.log')).open('w') as log:
            subprocess.run([sys.executable,str(Path(__file__).with_name('uspto_sar.py')),'--root',str(a.root),'--manifest',str(a.base/'manifest.json'),'--output',str(out),'--workers',str(workers)],stdout=log,check=True)
        assert digest(out)==baseline,label+' changed scientific, metadata or standalone JSON exports'
        runs.append(stats(out))
    assert runs[-1]['cache_hits']==42
    report=dict(version=runs[0]['version'],runs=runs,sha256=baseline,identical=True,files_checked=files)
    (a.base/'results/performance_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
    print(json.dumps(report,ensure_ascii=False))

if __name__=='__main__':main()
