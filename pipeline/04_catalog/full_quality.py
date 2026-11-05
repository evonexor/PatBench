"""Per-patent source checks, conservative association policy and original MOL export."""
from collections import defaultdict
from decimal import Decimal
import hashlib
import json
import math
import os
from pathlib import Path,PurePosixPath
import re
import zipfile
from rdkit import Chem
from audit_alignment import cals,number_alias,plain
from uspto_sar import worker,atomic_json
from xml_source import parse_xml

QUALITY_COLUMNS=['compound_id','alignment_grade','source_mapping_checks','drawing_vs_mol_checked','raw_molfile_path','raw_mol_sha256','raw_mol_coordinate_dimension']
FACTORS={'pM':Decimal('.001'),'nM':Decimal(1),'uM':Decimal(1000),'mM':Decimal(1000000),'M':Decimal(1000000000)}

def verify_source(r,tables,ids,description):
    failures=[]
    def row_at(locator):
        m=re.fullmatch(r'(.*)/tgroup-(\d+)/row-(\d+)(?:/column-(\d+))?',locator)
        if not m:raise ValueError('invalid_locator')
        loc,gi,rn,col=m.groups();cells,headers=tables[loc]
        return cells[(int(gi),int(rn))],headers,int(col)-1 if col else None
    normalize=lambda s:re.sub(r'\s+','',s).lower()
    target_norm=lambda s:re.sub(r'[\W_]+','',s).lower()
    try:
        cells,headers,col=row_at(r['activity_locator']);raw=plain(cells[col]);h=headers.get(col,'')
        if raw!=r['raw_value']:failures.append('raw_activity_cell_mismatch')
        numeric=re.fullmatch(r'\s*([<>=≤≥]+)?\s*(\d+(?:[.,]\d+)?|\.\d+)(?:\s*[a-zA-Zµμ%]+)?\s*',raw)
        if numeric and r['Unit'] in FACTORS:
            v=Decimal(numeric[2].replace(',',''));rel=(numeric[1] or '=').replace('≤','<=').replace('≥','>=')
            ok=rel==r['Relation'] and math.isclose(float(v),r['Value'],rel_tol=1e-12,abs_tol=1e-12) and math.isclose(float(v*FACTORS[r['Unit']]),r['value_nM'],rel_tol=1e-12,abs_tol=1e-12)
        else:ok=False
        if not ok:failures.append('value_relation_or_conversion_mismatch')
        pattern=r'\b(?:IC\s*50|EC\s*50|GI\s*50|AC\s*50|K\s*[id])\b'
        endpoints={normalize(x) for x in re.findall(pattern,h,re.I)}
        if normalize(r['Activity Type']) not in endpoints:failures.append('endpoint_not_in_original_column_header')
        without_endpoint=re.sub(pattern,'',h,flags=re.I)
        if not(target_norm(r['Target']) in target_norm(without_endpoint) or any(target_norm(r['Target']) in target_norm(plain(e)) for e in cells if e is not None)):failures.append('target_not_in_source_header_or_row')
        u_source=raw if r['unit_evidence']=='activity_cell_unit' else h
        u_source=u_source.replace('micromolar','uM').replace('nanomolar','nM').replace('millimolar','mM')
        units={u.replace('µ','u').replace('μ','u') for u in re.findall(r'(?<![A-Za-z])(?:pM|nM|[uµμ]M|mM|M)(?![A-Za-z])',u_source)}
        if r['unit_evidence']=='explicit_endpoint_unit_in_context':
            if normalize(r['assay_context']) not in normalize(description):failures.append('unit_context_not_in_source')
            units={u.replace('µ','u').replace('μ','u') for u in re.findall(r'(?:IC\s*50|EC\s*50|Ki|Kd)(?:\s+values?)?[^.;]{0,35}?(?:in|\()\s*([uµμ]M|nM|pM|mM)\b',r['assay_context'],re.I)}
        if units!={r['Unit']}:failures.append('unit_not_uniquely_in_source')
        alias_cols=[i for i,label in headers.items() if re.search(r'example|compound|cmpd|\bno\b',label,re.I) and not re.search(r'IC\s*50|EC\s*50|\bK[di]\b|name',label,re.I)]
        alias_found=any(i<len(cells) and number_alias(plain(cells[i]))==r['compound_alias'] for i in alias_cols)
        column_id=re.search(r'(?:Example|Ex\.)\s*#?\s*(\d+(?:-\d+)?)',r['activity_header'],re.I)
        prose_id=re.search(r'compound of Example\s*#?\s*'+re.escape(r['compound_alias'])+r'\b',r['assay_context'],re.I)
        diagram_id=r['compound_alias'].startswith('diagram:') and any(r['compound_alias'][8:].upper()==PurePosixPath(im.get('file','')).stem.upper() for e in cells if e is not None for im in e.iter('img'))
        if not(alias_found or bool(column_id and column_id[1]==r['compound_alias']) or bool(prose_id) or diagram_id):failures.append('activity_compound_id_not_found')
        file_stem=PurePosixPath(r['structure_file']).stem.upper()
        if r['structure_match']=='explicit_same_table_row':
            structural,_,_=row_at(r['structure_locator'])
            img_ok=any(PurePosixPath(im.get('file','')).stem.upper()==file_stem for e in structural if e is not None for im in e.iter('img'))
            id_ok=diagram_id or any(number_alias(plain(e))==r['compound_alias'] for e in structural if e is not None)
            if not(img_ok and id_ok):failures.append('structure_row_id_or_image_mismatch')
        else:
            loc=r['structure_locator'];heading,diagram=loc.split(' -> ',1) if ' -> ' in loc else (loc,loc)
            hnode=ids.get(heading);dnode=ids.get(diagram)
            heading_ok=hnode is not None and re.match(r'^(?:Example|Compound)\s*#?\s*'+re.escape(r['compound_alias'])+r'(?=$|[\s:—–.])',plain(hnode),re.I)
            img_ok=dnode is not None and any(PurePosixPath(im.get('file','')).stem.upper()==file_stem for im in dnode.iter('img'))
            if not(heading_ok and img_ok):failures.append('example_id_or_image_mismatch')
    except Exception as e:failures.append('source_check_unresolved:'+type(e).__name__)
    return failures

def full_worker(job):
    root,out,spec,code_hash=job;out=Path(out);archive=spec['archive'];archive_hash=hashlib.sha256(archive.encode()).hexdigest()
    cache=out/'checkpoint_cache'/archive_hash[:2];(cache/'checkpoints').mkdir(parents=True,exist_ok=True)
    small,raw_hit=worker((root,str(cache),spec,False))
    if not small['checkpoint']:return dict(checkpoint=None,stats=small['stats'],cached=False)
    original=json.loads(Path(small['checkpoint']).read_text())
    signature=hashlib.sha256((code_hash+'|'+original['signature']).encode()).hexdigest()
    cp=cache/'quality'/ (archive_hash[:20]+'.json');cp.parent.mkdir(exist_ok=True)
    if cp.exists():
        try:
            old=json.loads(cp.read_text())
            if old.get('quality_signature')==signature and all(not r.get('raw_molfile_path') or (out/r['raw_molfile_path']).is_file() for r in old['ready']+old['review']):return dict(checkpoint=str(cp),stats={k:old['stats'].get(k) for k in ['archive','status','Patent','paired_rows','review_rows','candidate_rows']},cached=True)
        except (OSError,ValueError,KeyError):pass
    original['candidates']=[];original['stats']['candidate_rows']=0;original['stats']['quality_check_failures']=0
    records=original['ready']+original['review'];initial_ready={r['source_record_id'] for r in original['ready']};original['ready']=[];original['review']=[]
    if original['stats']['status']=='processed':
        with zipfile.ZipFile(Path(root)/archive) as z:
            tree=parse_xml(z.read(original['stats']['xml_file']))
            for e in tree.iter():e.tag=e.tag.rsplit('}',1)[-1]
            ids={e.get('id'):e for e in tree.iter() if e.get('id')};parents={c:p for p in tree.iter() for c in p};tables={}
            for t in tree.iter('table'):
                parent=parents.get(t);loc=(parent.get('id','') if parent is not None else '') or t.get('id','')
                if loc:
                    try:tables[loc]=cals(t)
                    except (ValueError,IndexError,KeyError):pass
            description=plain(tree.find('description'));mol_info={}
            aliases=defaultdict(set)
            for r in records:
                if r['source_record_id'] in initial_ready:aliases[r['compound_alias']].add(r['InChIKey'])
            for r in records:
                r.update(compound_id=r['Patent']+':'+r['compound_alias'],alignment_grade='direct_numbered_structure_row' if r['structure_match']=='explicit_same_table_row' else 'adjacent_example_diagram_candidate' if r['structure_match'] else 'unresolved',source_mapping_checks='not_checked_in_review',drawing_vs_mol_checked='false',raw_molfile_path='',raw_mol_sha256='',raw_mol_coordinate_dimension='')
                failures=[]
                if r['SMILES'] and r['structure_file']:
                    filename=r['structure_file']
                    if filename not in mol_info:
                        try:
                            data=z.read(filename);mol=Chem.MolFromMolBlock(data.decode('utf-8','replace'),sanitize=True,removeHs=True,strictParsing=True)
                            if mol is None:raise ValueError('mol_parse_failed')
                            key=Chem.MolToInchiKey(mol);smiles=Chem.MolToSmiles(mol,canonical=True,isomericSmiles=True);sha=hashlib.sha256(data).hexdigest()
                            relative=str(PurePosixPath('molfiles')/r['Patent']/archive_hash[:10]/PurePosixPath(filename).name)
                            path=out/relative;path.parent.mkdir(parents=True,exist_ok=True)
                            if not path.exists() or hashlib.sha256(path.read_bytes()).hexdigest()!=sha:
                                tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_bytes(data);tmp.replace(path)
                            if hashlib.sha256(path.read_bytes()).hexdigest()!=sha:raise ValueError('mol_copy_checksum_failed')
                            conf=mol.GetConformer() if mol.GetNumConformers() else None
                            dimension='3D_flag_or_nonzero_z' if conf and (conf.Is3D() or any(abs(conf.GetAtomPosition(i).z)>1e-10 for i in range(mol.GetNumAtoms()))) else '2D' if conf else 'absent'
                            mol_info[filename]=dict(raw_molfile_path=relative,raw_mol_sha256=sha,raw_mol_coordinate_dimension=dimension,key=key,smiles=smiles)
                        except Exception as e:mol_info[filename]=dict(error=type(e).__name__)
                    info=mol_info[filename]
                    if 'error' in info:failures.append('source_mol_verification_failed:'+info['error'])
                    else:
                        r.update({k:info[k] for k in ['raw_molfile_path','raw_mol_sha256','raw_mol_coordinate_dimension']})
                        if info['key']!=r['InChIKey'] or info['smiles']!=r['SMILES']:failures.append('mol_smiles_or_inchikey_mismatch')
                if r['source_record_id'] in initial_ready:
                    failures+=verify_source(r,tables,ids,description)
                    if len(aliases[r['compound_alias']])!=1:failures.append('compound_id_structure_conflict')
                    r['source_mapping_checks']='passed' if not failures else ';'.join(failures)
                    if failures:original['stats']['quality_check_failures']+=1
                    if r['alignment_grade']!='direct_numbered_structure_row':failures.append('example_diagram_association_requires_review')
                    if failures:
                        r['review_reason']=';'.join(dict.fromkeys(filter(None,[r['review_reason'],*failures])));r['extraction_status']='review'
                        original['candidates'].append(r);original['review'].append(r)
                    else:original['ready'].append(r)
                else:
                    if failures:r['review_reason']=';'.join(dict.fromkeys(filter(None,[r['review_reason'],*failures])))
                    original['review'].append(r)
            bib=original.get('documents',{}).get('bibliography')
            if bib:
                relative=str(PurePosixPath('bibliography')/original['stats']['Patent']/(archive_hash[:10]+'.json'))
                original['stats']['bibliography_json_path']=relative;original['bibliography_path']=relative
                for r in records:r['bibliography_json_path']=relative
    original['stats'].update(paired_rows=len(original['ready']),review_rows=len(original['review']),candidate_rows=len(original['candidates']))
    original['quality_signature']=signature;atomic_json(cp,original)
    return dict(checkpoint=str(cp),stats={k:original['stats'].get(k) for k in ['archive','status','Patent','paired_rows','review_rows','candidate_rows']},cached=raw_hit)
