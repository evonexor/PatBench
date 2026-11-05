"""Deterministic, stratified sample selection. No bioactivity is copied into results."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import xml.etree.ElementTree as ET
import zipfile

MEDICAL = re.compile(r'inhibitor|antagonist|agonist|therapeut|cancer|kinase|pharmaceutical|antiviral|antibacterial|antitumor|treating|treatment', re.I)
EXCLUDE = re.compile(r'plant growth|herbicid|insecticid|agricultur|agrochemical|crop|pesticid|electrochemical|battery|semiconductor|polymer|photovoltaic|hydrogen pump|asphalt|pavement|corrosion|coating|lubricant|construction material', re.I)

def core(p):
    p = re.sub(r'[^A-Z0-9]', '', p.upper())
    p = re.sub(r'[ABCS][0-9]?$', '', p)
    m = re.fullmatch(r'US(\d+)', p)
    return 'US' + str(int(m[1])) if m else p

def metadata(path):
    with zipfile.ZipFile(path) as z:
        from xml_source import select_patent_xml,parse_xml
        name = select_patent_xml(z,path)
        if not name:
            return None
        r = parse_xml(z.read(name))
        return r.findtext('.//invention-title', ''), ' '.join(r.find('abstract').itertext()) if r.find('abstract') is not None else ''

def main():
    a = argparse.ArgumentParser()
    a.add_argument('--root', required=True, type=Path)
    a.add_argument('--output', required=True, type=Path)
    args = a.parse_args()
    root = args.root
    index = {}
    for path in sorted((root / 'uspto-mol/data').rglob('*.ZIP')):
        index.setdefault(core(path.stem.split('-')[0]), path)
    c = sqlite3.connect('file:' + str(root / 'chembl/chembl_37/extracted/chembl_37.db') + '?mode=ro', uri=True)
    docs = c.execute('select doc_id,patent_id,title,year from docs where patent_id is not null and src_id in (37,38)').fetchall()
    candidates = {}
    for doc_id, patent, title, year in docs:
        k = core(patent)
        if k in index and MEDICAL.search(title or '') and not EXCLUDE.search(title or ''):
            candidates.setdefault(k, dict(patent_core=k, archive=str(index[k].relative_to(root)), title=title, year=year, selection='ChEMBL stratified reference', reference_doc_ids=[]))['reference_doc_ids'].append(doc_id)
    anchors = ['US8841312','US8962608','US8889696','US8481536','US8614223','US9102686','US8969341','US9051311','US8957074','US9108978']
    selected = [candidates[k] for k in anchors if k in candidates]
    chosen = {x['patent_core'] for x in selected}
    remaining = sorted((x for k,x in candidates.items() if k not in chosen), key=lambda x:hashlib.sha256(('20261003:'+x['patent_core']).encode()).hexdigest())
    # At most two patents per publication year in the randomized addition.
    counts = {}
    for x in remaining:
        if counts.get(x['year'],0) >= 2:
            continue
        selected.append(x); counts[x['year']] = counts.get(x['year'],0)+1
        if len(selected) == 29:
            break
    p = root / 'uspto-mol/data/2013/I20130101/REISSUE/USRE043899-20130101.ZIP'
    selected.append(dict(patent_core='USRE043899',archive=str(p.relative_to(root)),title=metadata(p)[0],selection='Medicinal no-numeric-assay challenge',reference_doc_ids=[]))
    chosen = {x['patent_core'] for x in selected}
    recent = sorted((p for k,p in index.items() if '/2025/' in str(p) and k not in chosen),key=lambda p:hashlib.sha256(('20261003:'+p.name).encode()).hexdigest())
    screened = 0
    for p in recent[:600]:
        screened += 1
        try:
            t, abstract = metadata(p)
        except (ET.ParseError,ValueError,TypeError,zipfile.BadZipFile):
            continue
        if MEDICAL.search(t) and not EXCLUDE.search(t):
            selected.append(dict(patent_core=core(p.stem.split('-')[0]),archive=str(p.relative_to(root)),title=t,year=2025,selection='Independent 2025 medical title sample',reference_doc_ids=[]))
        if len(selected) >= 40:
            break
    negatives = [
        ('uspto-mol/data/2001/20010102/UTIL06167/US06167721-20010102.ZIP','materials'),
        ('uspto-mol/data/2025/I20250107/UTIL12185/US12185720-20250107.ZIP','agriculture'),
    ]
    for p,reason in negatives:
        selected.append(dict(archive=p,patent_core=core(Path(p).stem.split('-')[0]),selection='Exclusion fixture',expected_exclusion=reason,reference_doc_ids=[]))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(dict(seed=20261003,medical_sample_count=40,recent_screened=screened,archives_indexed=len(index),samples=selected),ensure_ascii=False,indent=2))
    print(json.dumps(dict(samples=len(selected),references=sum(bool(x['reference_doc_ids']) for x in selected),recent_screened=screened,titles=[(x['patent_core'],x.get('title','')) for x in selected]),ensure_ascii=False))

if __name__ == '__main__':
    main()
