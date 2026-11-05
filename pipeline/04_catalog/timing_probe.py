"""Random-corpus screening probe and <=100-patent timing sample, never a full run."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
import zipfile
from uspto_sar import screening,text,extract_metadata
from xml_source import parse_xml,select_patent_xml

def screen_one(item):
    root,entry=item;start=time.monotonic();result={**entry}
    try:
        with zipfile.ZipFile(Path(root)/entry['archive']) as z:
            names=z.namelist();xml=select_patent_xml(z,entry['archive'])
            result['mol_files']=sum(n.lower().endswith('.mol') for n in names)
            if xml:
                tree=parse_xml(z.read(xml))
                for e in tree.iter():e.tag=e.tag.rsplit('}',1)[-1]
                title_node=tree.find('.//invention-title');abstract_node=tree.find('abstract')
                title=text(title_node) if title_node is not None else ''
                abstract=text(abstract_node) if abstract_node is not None else ''
                cpc=' '.join(text(e) for e in tree.iter() if e.tag in ('classification-cpc','classification-ipcr','main-classification'))
                cpc+=' '+' '.join(e.findtext('section','')+e.findtext('class','')+e.findtext('subclass','') for e in tree.iter('classification-ipcr'))
                screen,why=screening(title,abstract,cpc)
                if tree.find('us-bibliographic-data-grant') is not None:
                    extract_metadata(tree)  # Current extractor also parses metadata before exclusions.
                    result.update(screen=screen,title=title,reason=why,format='xml')
                else:
                    result.update(screen=screen if screen in ('agriculture','materials') else 'legacy_xml_review',title=title,reason=why,format='legacy_xml')
            else:
                sgm=next((n for n in names if n.lower().endswith('.sgm')),None)
                raw=z.read(sgm).decode('utf-8','replace') if sgm else ''
                m=re.search(r'<B540[^>]*>(.*?)</B540>',raw,re.S|re.I)
                title=re.sub('<[^>]+>',' ',m[1]).strip() if m else ''
                screen,why=screening(title,'','')
                result.update(screen=screen if screen in ('agriculture','materials') else 'legacy_sgml_review',title=title,reason=why,format='sgml_or_missing_xml')
    except Exception as e:result.update(screen='error',error=str(e))
    result['seconds']=time.monotonic()-start
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--inventory',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--workers',type=int,default=20);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True);inventory=json.loads(a.inventory.read_text());entries=inventory['random_sample'];start=time.monotonic()
    with ProcessPoolExecutor(max_workers=a.workers) as pool:probes=list(pool.map(screen_one,[(str(a.root),entry) for entry in entries]))
    screen_seconds=time.monotonic()-start
    medical=[r for r in probes if r['screen']=='medical'][:80];other=[r for r in probes if r['screen'] not in ('medical','error')][:20]
    chosen=medical+other;assert len(chosen)<=100
    manifest=a.output/'manifest.json';manifest.write_text(json.dumps(dict(mode='timing_sample',samples=[dict(archive=r['archive'],selection='Random corpus timing probe') for r in chosen])))
    with (a.output/'run.log').open('w') as f:
        subprocess.run([sys.executable,str(Path(__file__).with_name('uspto_sar.py')),'--root',str(a.root),'--manifest',str(manifest),'--output',str(a.output/'results'),'--workers',str(a.workers)],stdout=f,stderr=subprocess.STDOUT,check=True)
    run=json.loads((a.output/'results/run_report.json').read_text());counts=Counter(r['screen'] for r in probes);n=inventory['total_zip'];fraction=counts['medical']/len(probes)
    screen_projection=n*screen_seconds/len(probes)
    medical_projection=n*fraction*run['elapsed_seconds']/len(medical) if medical else None
    report=dict(workers=a.workers,total_zip=n,zip_bytes=inventory['zip_bytes'],random_screen_sample=len(probes),screen_counts=dict(counts),screen_wall_seconds=round(screen_seconds,3),screen_cpu_seconds=sum(r['seconds'] for r in probes),medical_fraction=fraction,timing_sample=len(chosen),timing_medical=len(medical),timing_run=run,screen_only_projection_hours=screen_projection/3600,medical_stage_projection_hours=medical_projection/3600 if medical_projection is not None else None,combined_conservative_projection_hours=(screen_projection+medical_projection)/3600 if medical_projection is not None else None,limitations=['This is a random 600-patent screening sample and at most 100-patent extraction probe, not a full run.','The two-stage projection deliberately includes overlap in screening cost. Small-sample startup, cached reads, output overhead, long-tail patents and corpus disk contention affect extrapolation.','The estimate covers the current parser and its review/skip outcomes, not manual drawing validation or new OCR/legacy-SGML implementations.'],screen_records=probes)
    (a.output/'timing_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
    compact={k:v for k,v in report.items() if k not in ('screen_records','timing_run')};compact['timing_run']={k:run[k] for k in ['statuses','elapsed_seconds','cache_hits','paired_observations','review_observations','structure_observations']}
    print(json.dumps(compact,ensure_ascii=False),flush=True)

if __name__=='__main__':main()
