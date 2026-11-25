"""Read back original XML/MOL for every paired observation; export MOL unchanged.

This checks source correspondence and data integrity, not drawing-to-MOL accuracy.
Heading/nearby-diagram associations remain candidates, even if all checks pass.
"""
import argparse
from collections import Counter,defaultdict
import csv
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path,PurePosixPath
import re
import xml.etree.ElementTree as ET
import zipfile
import pyarrow.parquet as pq
from rdkit import Chem,RDLogger
from xml_source import parse_xml

RDLogger.DisableLog('rdApp.*')

def plain(e):
    if e is None:return ''
    def walk(n):
        s=n.text or ''
        for c in n:s+=('' if c.tag in ('sub','sup') else ' ')+walk(c)+('' if c.tag in ('sub','sup') else ' ')+(c.tail or '')
        return s
    return ' '.join(walk(e).split())

def number_alias(s):
    s=re.sub(r'^(?:example|compound|cmpd\.?|ex\.?|no\.?|#)\s*#?\s*','',s.strip(),flags=re.I)
    return s.lower().replace('.','-') if re.fullmatch(r'\d+(?:[.-]\d+)?[a-z]?',s,re.I) else None

def cals(table):
    """Resolve the original logical cells, including horizontal/vertical spans."""
    result={};headers=defaultdict(list)
    groups=table.findall('tgroup') or [table]
    width=max(int(g.get('cols','0')) or max(len(r) for r in g.iter('row')) for g in groups)
    for gi,g in enumerate(groups,1):
        names={e.get('colname'):i for i,e in enumerate(g.findall('colspec'))}
        n=int(g.get('cols','0')) or max(len(r) for r in g.iter('row'))
        parents={child:parent for parent in g.iter() for child in parent};spans={};seen_data=False
        for rn,row in enumerate(g.iter('row'),1):
            cells=[None]*n
            for col,(left,entry) in list(spans.items()):
                cells[col]=entry
                if left==1:del spans[col]
                else:spans[col]=(left-1,entry)
            cursor=0
            for entry in row:
                if entry.tag not in ('entry','td','th'):continue
                explicit=entry.get('namest') or entry.get('colname')
                if explicit in names:cursor=names[explicit]
                else:
                    while cursor<n and cells[cursor] is not None:cursor+=1
                last=names.get(entry.get('nameend'),cursor)
                for col in range(cursor,min(last+1,n)):
                    cells[col]=entry
                    if int(entry.get('morerows','0')):spans[col]=(int(entry.get('morerows')),entry)
                cursor=last+1
            result[(gi,rn)]=cells
            data_row=any(e is not None and (e.find('.//img') is not None or re.fullmatch(r'[<>=≤≥]?\s*(?:\d+(?:\.\d+)?|\.\d+)',plain(e))) for e in cells)
            if n==width and (parents[row].tag=='thead' or (not seen_data and not data_row)):
                for col,entry in enumerate(cells):
                    if entry is not None:headers[col].append(plain(entry))
            if parents[row].tag!='thead' and data_row:seen_data=True
    return result,{col:re.sub(r'(?<=[A-Za-z])-\s+(?=[A-Za-z])','',' '.join(dict.fromkeys(h))) for col,h in headers.items()}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--results',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    grouped=defaultdict(list)
    for r in pq.read_table(a.results/'SAR_with_evidence.parquet').to_pylist():grouped[r['archive']].append(r)
    audits=[];mols={};alias_keys=defaultdict(set)
    factors={'pM':Decimal('.001'),'nM':Decimal(1),'uM':Decimal(1000),'mM':Decimal(1000000),'M':Decimal(1000000000)}
    for archive,records in grouped.items():
        with zipfile.ZipFile(a.root/archive) as z:
            root=parse_xml(z.read(records[0]['xml_file']))
            for e in root.iter():e.tag=e.tag.rsplit('}',1)[-1]
            ids={e.get('id'):e for e in root.iter() if e.get('id')};parents={c:p for p in root.iter() for c in p};table_cache={}
            for t in root.iter('table'):
                parent=parents.get(t);loc=(parent.get('id','') if parent is not None else '') or t.get('id','')
                if loc:table_cache[loc]=cals(t)
            def source_row(locator):
                match=re.fullmatch(r'(.*)/tgroup-(\d+)/row-(\d+)(?:/column-(\d+))?',locator)
                if not match:raise ValueError('Invalid source locator: '+locator)
                loc,gi,rn,col=match.groups();cells,headers=table_cache[loc]
                return cells[(int(gi),int(rn))],headers,int(col)-1 if col else None
            for r in records:
                failures=[];activity_cells,activity_headers,column=source_row(r['activity_locator'])
                raw=plain(activity_cells[column]);raw_ok=raw==r['raw_value']
                if not raw_ok:failures.append('raw_activity_cell_mismatch')
                parsed=re.fullmatch(r'\s*([<>=≤≥]+)?\s*(\d+(?:[.,]\d+)?|\.\d+)(?:\s*[a-zA-Zµμ%]+)?\s*',raw)
                numeric_ok=False
                if parsed:
                    relation=(parsed[1] or '=').replace('≤','<=').replace('≥','>=')
                    v=Decimal(parsed[2].replace(',',''));numeric_ok=relation==r['Relation'] and math.isclose(float(v),r['Value'],rel_tol=1e-12,abs_tol=1e-12) and math.isclose(float(v*factors[r['Unit']]),r['value_nM'],rel_tol=1e-12,abs_tol=1e-12)
                if not numeric_ok:failures.append('value_relation_or_conversion_mismatch')
                source_header=activity_headers.get(column,'')
                endpoints={re.sub(r'\s+','',x).lower() for x in re.findall(r'\b(?:IC\s*50|EC\s*50|GI\s*50|AC\s*50|K\s*[id])\b',source_header,re.I)}
                if re.sub(r'\s+','',r['Activity Type']).lower() not in endpoints:failures.append('endpoint_not_in_original_column_header')
                normalize=lambda s:re.sub(r'\s+','',s).lower()
                target_normalize=lambda s:re.sub(r'[\W_]+','',s).lower()
                without_endpoint=re.sub(r'\b(?:IC\s*50|EC\s*50|GI\s*50|AC\s*50|K\s*[id])\b','',source_header,flags=re.I)
                target_ok=target_normalize(r['Target']) in target_normalize(without_endpoint) or any(target_normalize(r['Target']) in target_normalize(plain(e)) for e in activity_cells if e is not None)
                if not target_ok:failures.append('target_literal_not_in_source_header_or_row')
                unit_source=raw if r['unit_evidence']=='activity_cell_unit' else source_header
                unit_source=unit_source.replace('micromolar','uM').replace('nanomolar','nM').replace('millimolar','mM')
                source_units={u.replace('µ','u').replace('μ','u') for u in re.findall(r'(?<![A-Za-z])(?:pM|nM|[uµμ]M|mM|M)(?![A-Za-z])',unit_source)}
                if r['unit_evidence']=='explicit_endpoint_unit_in_context':
                    description=plain(root.find('description'))
                    if normalize(r['assay_context']) not in normalize(description):failures.append('unit_context_not_in_original_description')
                    source_units={u.replace('µ','u').replace('μ','u') for u in re.findall(r'(?:IC\s*50|EC\s*50|Ki|Kd)(?:\s+values?)?[^.;]{0,35}?(?:in|\()\s*([uµμ]M|nM|pM|mM)\b',r['assay_context'],re.I)}
                if source_units!={r['Unit']}:failures.append('reported_unit_not_uniquely_in_source')
                alias_cols=[col for col,h in activity_headers.items() if re.search(r'example|compound|cmpd|\bno\b',h,re.I) and not re.search(r'IC\s*50|EC\s*50|\bK[di]\b|name',h,re.I)]
                row_alias=any(number_alias(plain(activity_cells[col]))==r['compound_alias'] for col in alias_cols if col<len(activity_cells))
                header_example=re.search(r'(?:Example|Ex\.)\s*#?\s*(\d+(?:-\d+)?)',r['activity_header'],re.I)
                context_example=re.search(r'compound of Example\s*#?\s*'+re.escape(r['compound_alias'])+r'\b',r['assay_context'],re.I)
                diagram_alias=r['compound_alias'].startswith('diagram:') and any(r['compound_alias'][8:].upper()==PurePosixPath(im.get('file','')).stem.upper() for e in activity_cells if e is not None for im in e.iter('img'))
                alias_ok=row_alias or bool(header_example and header_example[1]==r['compound_alias']) or bool(context_example) or diagram_alias
                if not alias_ok:failures.append('activity_compound_id_not_found')
                filename=r['structure_file'];key=(archive,filename)
                if key not in mols:
                    data=z.read(filename);mol=Chem.MolFromMolBlock(data.decode('utf-8','replace'),sanitize=True,removeHs=True,strictParsing=True)
                    if mol is None:raise ValueError('Original MOL failed parsing: '+filename)
                    path=PurePosixPath('molfiles')/r['Patent']/PurePosixPath(filename).name
                    destination=a.output/str(path);destination.parent.mkdir(parents=True,exist_ok=True);destination.write_bytes(data)
                    conf=mol.GetConformer() if mol.GetNumConformers() else None
                    coords=[tuple(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())] if conf else []
                    nonzero=bool(coords and any(abs(v)>1e-10 for xyz in coords for v in xyz))
                    nonzero_z=bool(coords and any(abs(xyz[2])>1e-10 for xyz in coords))
                    mols[key]=dict(molfile_path=str(path),mol_sha256=hashlib.sha256(data).hexdigest(),canonical_smiles=Chem.MolToSmiles(mol,canonical=True,isomericSmiles=True),mol_inchikey=Chem.MolToInchiKey(mol),atom_count=mol.GetNumAtoms(),has_coordinates=bool(conf),coordinate_dimension='3D_flag_or_nonzero_z' if conf and (conf.Is3D() or nonzero_z) else '2D' if conf else 'absent',nonzero_coordinates=nonzero)
                m=mols[key];structure_ok=m['canonical_smiles']==r['SMILES'] and m['mol_inchikey']==r['InChIKey']
                if not structure_ok:failures.append('mol_smiles_or_inchikey_mismatch')
                loc=r['structure_locator'];method=r['structure_match']
                if method=='explicit_same_table_row':
                    cells,headers,_=source_row(loc)
                    img_ok=any(PurePosixPath(im.get('file','')).stem.upper()==PurePosixPath(filename).stem.upper() for e in cells if e is not None for im in e.iter('img'))
                    id_ok=diagram_alias or any(number_alias(plain(e))==r['compound_alias'] for e in cells if e is not None)
                    if not(img_ok and id_ok):failures.append('structure_row_id_or_image_mismatch')
                    grade='direct_numbered_structure_row'
                else:
                    heading,diagram=loc.split(' -> ',1) if ' -> ' in loc else (loc,loc)
                    h=ids.get(heading);d=ids.get(diagram)
                    heading_ok=bool(h is not None and re.match(r'^(?:Example|Compound)\s*#?\s*'+re.escape(r['compound_alias'])+r'(?=$|[\s:—–.])',plain(h),re.I))
                    img_ok=bool(d is not None and any(PurePosixPath(im.get('file','')).stem.upper()==PurePosixPath(filename).stem.upper() for im in d.iter('img')))
                    if not(heading_ok and img_ok):failures.append('example_heading_or_image_reference_mismatch')
                    grade='adjacent_example_diagram_candidate'
                alias_keys[(r['Patent'],r['compound_alias'])].add(r['InChIKey'])
                audits.append({**{k:r[k] for k in ['source_record_id','Patent','compound_alias','Target','Activity Type','Relation','Value','Unit','value_nM','SMILES','InChIKey','raw_value','structure_match','structure_locator','activity_locator','archive','xml_file','structure_file']},**m,'alignment_grade':grade,'source_checks_pass':not failures,'drawing_vs_mol_checked':False,'audit_failures':';'.join(failures)})
    inconsistent=[key for key,v in alias_keys.items() if len(v)>1]
    report=dict(paired_records=len(audits),unique_patent_compound_ids=len(alias_keys),original_molfiles=len(mols),passed_source_checks=sum(r['source_checks_pass'] for r in audits),failed_source_checks=sum(not r['source_checks_pass'] for r in audits),conflicting_compound_id_mappings=len(inconsistent),alignment_grades=dict(Counter(r['alignment_grade'] for r in audits)),coordinate_dimensions=dict(Counter(m['coordinate_dimension'] for m in mols.values())),all_molfiles_have_nonzero_coordinates=all(m['has_coordinates'] and m['nonzero_coordinates'] for m in mols.values()),limitations=['Source-row correspondence does not prove that the USPTO MOL exactly matches its original structure drawing.','Adjacent example diagrams remain candidates; human/source-image review is required before asserting complete correctness.','Patent + compound_alias is the compound identity; source_record_id identifies an assay observation.'])
    for name,rows in [('compound_activity_mol_mapping',audits),('SAR_direct_links',[r for r in audits if r['source_checks_pass'] and r['alignment_grade']=='direct_numbered_structure_row']),('SAR_candidate_links',[r for r in audits if not r['source_checks_pass'] or r['alignment_grade']!='direct_numbered_structure_row'])]:
        with (a.output/(name+'.csv')).open('w',encoding='utf-8-sig',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(audits[0]));w.writeheader();w.writerows(rows)
    (a.output/'alignment_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
    print(json.dumps(report,ensure_ascii=False))
    if report['failed_source_checks'] or inconsistent:raise SystemExit(2)

if __name__=='__main__':main()
