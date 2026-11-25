"""Prepare a corpus manifest after sample review; does not extract any patents."""
import argparse
import json
from pathlib import Path
from select_samples import core

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--after-sample-review',action='store_true',required=True)
    a=p.parse_args()
    samples=[dict(archive=str(f.relative_to(a.root)),patent_core=core(f.stem.split('-')[0]),selection='Full corpus, classification during extraction') for f in sorted((a.root/'uspto-mol/data').rglob('*.ZIP'))]
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(dict(mode='full',samples=samples),ensure_ascii=False))
    print(f'Prepared {len(samples):,} archive entries. No extraction started.')
