"""Evidence-preserving USPTO MOL/XML to SAR extraction. Python 3.10+, RDKit, PyArrow.

Strict results require explicit compound/structure links, explicit assay units and
an unambiguous assay target. Missing/ambiguous observations remain in review.csv.
Existing ChEMBL activities are used only by a separate validation command.
"""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import csv
from decimal import Decimal, InvalidOperation
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
import sys
import time
import xml.etree.ElementTree as ET
import zipfile
from rdkit import Chem, RDLogger, rdBase
from select_samples import MEDICAL, EXCLUDE, core
from patent_metadata import extract_metadata, METADATA_COLUMNS, date_iso, document
from patent_documents import extract_documents, DOCUMENT_INTEGER_COLUMNS
from xml_source import parse_xml,select_patent_xml

RDLogger.DisableLog('rdApp.*')
VERSION = '0.7.2'
SAR_COLUMNS = ['source','Patent','patent_core','Target','Activity Type','Relation','Value','value_nM','Unit','SMILES','SMILES_original','InChIKey','UniProt','organism','gene','final_label','source_record_id']
EXTRA_COLUMNS = ['compound_alias','raw_value','value_lower','value_upper','activity_header','target_evidence','unit_evidence','assay_context','structure_match','structure_file','structure_name','structure_name_check','structure_locator','activity_locator','archive','xml_file','publication_date','extraction_status','review_reason']
ACTIVITY = re.compile(r'\b(p?IC\s*50|p?EC\s*50|GI\s*50|AC\s*50|K\s*[id])\b', re.I)
ALIAS = re.compile(r'^(?:(?:Example|Ex\.?|Compound|Cmpd\.?|No\.?|#)\s*#?\s*)?(\d+(?:[.-]\d+)?[a-z]?)$', re.I)
EXAMPLE_HEADING = re.compile(r'^(?:Example|Compound)\s*#?\s*(\d+(?:-\d+)?[a-z]?)(?=$|[\s:—–.])(?!\s*(?:and\b|to\b|through\b|&|,))(?:\s*[:.—–]?\s*.*)?$', re.I)
UNIT = re.compile(r'(?<![A-Za-z])(pM|nM|[uµμ]M|mM|M)(?![A-Za-z])')
UNITS = {'pM':Decimal('.001'),'nM':Decimal(1),'uM':Decimal(1000),'mM':Decimal(1000000),'M':Decimal(1000000000)}

def text(e):
    """Keep chemical/endpoint subscripts attached, normalize whitespace elsewhere."""
    if e is None:return ''
    def walk(x):
        s = x.text or ''
        for c in x:
            v = walk(c)
            s += (v if c.tag in ('sub','sup') else ' '+v+' ') + (c.tail or '')
        return s
    return re.sub(r'\s+', ' ', walk(e)).strip()

def alias(s):
    m = ALIAS.fullmatch(re.sub(r'\s+', '',s) if re.fullmatch(r'[\d\s.-]+',s) else s.strip())
    return m[1].lower().replace('.', '-') if m else None

def unit(s):
    s = s.replace('micromolar','uM').replace('nanomolar','nM').replace('millimolar','mM')
    units = {m[1].replace('µ','u').replace('μ','u') for m in UNIT.finditer(s)}
    return next(iter(units)) if len(units)==1 else ''

def measurement(s):
    s = s.strip().replace('−','-').replace('–','-').replace('≤','<=').replace('≥','>=')
    number = r'(?:\d+(?:[.,]\d+)?|\.\d+)(?:[eE][+-]?\d+)?'
    m = re.fullmatch(r'\s*([<>]=?|=|~|≈)?\s*('+number+r')\s*(?:[a-zA-Zµμ%]+)?\s*',s)
    if m:
        v = Decimal(m[2].replace(',',''))
        return m[1] or '=',float(v),None,None
    m = re.fullmatch(r'('+number+r')\s*(?:-|to)\s*('+number+r')\s*(?:[a-zA-Zµμ]+)?',s)
    if m:
        return 'range',None,float(Decimal(m[1].replace(',',''))),float(Decimal(m[2].replace(',','')))
    return 'text',None,None,None

def expanded_rows(table):
    """CALS colspec, namest/nameend, morerows, offset columns, multiple tgroups."""
    groups = list(table.findall('tgroup')) or [table]
    width = max(int(g.get('cols','0')) or max((len(list(r)) for r in g.iter('row')),default=0) for g in groups)
    headers, body = [],[]
    for gi,g in enumerate(groups):
        cols = int(g.get('cols','0')) or width
        if cols != width:
            continue  # separate caption-only tgroup
        names = {c.get('colname'):i for i,c in enumerate(g.findall('colspec'))}
        spans = {}
        parents = {c:p for p in g.iter() for c in p}
        for rn,r in enumerate(g.iter('row')):
            cells = [dict(text='',files=[],loc='',row_locator=f'tgroup-{gi+1}/row-{rn+1}') for _ in range(width)]
            occupied = set()
            for k,(remaining,v) in list(spans.items()):
                cells[k] = v.copy(); occupied.add(k)
                if remaining <= 1:del spans[k]
                else:spans[k]=(remaining-1,v)
            col = 0
            for e in r:
                if e.tag not in ('entry','td','th'):continue
                if e.get('namest') in names:col=names[e.get('namest')]
                elif e.get('colname') in names:col=names[e.get('colname')]
                else:
                    while col in occupied:col+=1
                end = names.get(e.get('nameend'),col)
                v = dict(text=text(e),files=[im.get('file','') for im in e.iter('img')],loc=e.get('id',''),row_locator=f'tgroup-{gi+1}/row-{rn+1}')
                # Title spanning all columns is not a column header.
                full_caption = col==0 and end==width-1 and (not ACTIVITY.search(v['text']) or v['text'].startswith('TABLE'))
                for k in range(col,min(end+1,width)):
                    cells[k] = dict(v, text='' if full_caption else v['text'])
                    occupied.add(k)
                    if int(e.get('morerows','0')):spans[k]=(int(e.get('morerows')),v)
                col=end+1
            p=parents.get(r)
            for c in cells:c['row_locator']=f'tgroup-{gi+1}/row-{rn+1}'
            if p is not None and p.tag=='thead':headers.append(cells)
            else:body.append(cells)
    if not headers:
        # A few documents put ordinary headers in tbody. Accept only an explicit endpoint header.
        split = next((i for i,r in enumerate(body[:12]) if any(ACTIVITY.search(c['text']) for c in r)),None)
        if split is not None:headers,body=body[:split+1],body[split+1:]
    # Many USPTO tables put the second half of a header (units, compound number)
    # at the beginning of tbody. Promote only rows before the first observation.
    while body:
        row=body[0]
        if any(c['files'] or measurement(c['text'])[1] is not None or alias(c['text']) for c in row):break
        if any(len(c['text'])>180 for c in row):break
        if any(c['text'] for c in row):headers.append(row)
        body=body[1:]
    combined = [' '.join(dict.fromkeys(r[i]['text'] for r in headers if r[i]['text'])) for i in range(width)]
    combined = [re.sub(r'(?<=[A-Za-z])-\s+(?=[A-Za-z])','',h) for h in combined]
    return combined,body

def screening(title, abstract, cpc):
    if EXCLUDE.search(title):
        agricultural = re.search(r'agricultur|plant|crop|herbicid|insecticid|agro|pesticid',title,re.I)
        return ('agriculture' if agricultural else 'materials'), 'explicit excluded title term'
    if MEDICAL.search(title) or re.search(r'\bA61[KP]\b',cpc):
        return 'medical','medical title or A61K/A61P classification'
    if re.search(r'\bA01[HN]\b',cpc) and not re.search(r'\bA61[KP]\b',cpc):
        return 'agriculture','agricultural classification without medical evidence'
    return 'uncertain','no sufficiently specific medical evidence'

def name_element_conflicts(name,counts):
    """Conservative presence check, NOT an IUPAC name-to-structure validator."""
    normalized=re.sub(r'[\s-]+','',name.lower())
    required={'F':r'fluoro','Cl':r'chloro','Br':r'bromo','I':r'iodo','S':r'sulfon|sulfin|sulfam|sulph|thio','P':r'phosph','Se':r'selen'}
    return [e for e,pattern in required.items() if re.search(pattern,normalized) and not counts.get(e)]

def load_mol(z,filename,cache):
    if filename in cache:return cache[filename]
    result = dict(structure_file=filename,SMILES='',SMILES_original='',InChIKey='',structure_status='',structure_name='',structure_name_check='not_checked',element_counts={})
    try:
        block=z.read(filename).decode('utf-8','replace')
        mol=Chem.MolFromMolBlock(block,sanitize=True,removeHs=True,strictParsing=True)
        if mol is None:result['structure_status']='invalid_mol'
        elif any(a.GetAtomicNum()==0 or a.HasQuery() for a in mol.GetAtoms()) or any(b.HasQuery() for b in mol.GetBonds()):result['structure_status']='markush_or_query'
        elif len(Chem.GetMolFrags(mol))!=1:result['structure_status']='multiple_fragments_review'
        else:
            result['SMILES_original']=Chem.MolToSmiles(mol,canonical=False,isomericSmiles=True)
            result['SMILES']=Chem.MolToSmiles(mol,canonical=True,isomericSmiles=True)
            result['InChIKey']=Chem.MolToInchiKey(mol)
            result['element_counts']=dict(Counter(a.GetSymbol() for a in mol.GetAtoms()))
            result['structure_status']='concrete' if result['InChIKey'] else 'no_inchikey'
    except Exception as e:result['structure_status']='mol_error:'+type(e).__name__
    cache[filename]=result
    return result

def extract(path, relative, spec):
    start=time.monotonic()
    archive=relative
    stats=dict(archive=archive,selection=spec.get('selection',''),patent_core=spec.get('patent_core',''),title=spec.get('title',''),status='pending')
    ready,review,structure_rows,table_rows=[],[],[],[]
    with zipfile.ZipFile(path) as z:
        names=z.namelist()
        xml=select_patent_xml(z,path)
        mols={Path(n).stem.upper():n for n in names if n.lower().endswith('.mol')}
        if not xml:
            # SGML legacy title is enough for exclusion; do not silently fake an XML parse.
            sgm=next((n for n in names if n.lower().endswith('.sgm')),None)
            raw=z.read(sgm).decode('utf-8','replace') if sgm else ''
            m=re.search(r'<B540[^>]*>(.*?)</B540>',raw,re.S|re.I)
            title=re.sub('<[^>]+>',' ',m[1]).strip() if m else ''
            screen,why=screening(title,'','')
            stats.update(title=title,screen=screen,screen_reason=why,status='excluded' if screen in ('agriculture','materials') else 'legacy_sgml_review',mol_files=len(mols),elapsed_seconds=round(time.monotonic()-start,3))
            return dict(stats=stats,ready=[],review=[],structures=[],tables=[])
        root=parse_xml(z.read(xml))
        # DTDs in this corpus have no namespace; strip any namespace defensively.
        for e in root.iter():e.tag=e.tag.rsplit('}',1)[-1]
        title=text(root.find('.//invention-title'))
        abst=root.find('abstract'); abstract=text(abst) if abst is not None else ''
        cpc=' '.join(text(e) for e in root.iter() if e.tag in ('classification-cpc','classification-ipcr','main-classification'))
        cpc+=' '+' '.join((e.findtext('section','')+e.findtext('class','')+e.findtext('subclass','')) for e in root.iter('classification-ipcr'))
        screen,why=screening(title,abstract,cpc)
        pub=root.find('.//publication-reference/document-id')
        docnum=pub.findtext('doc-number','') if pub is not None else Path(path).stem.split('-')[0][2:]
        kind=pub.findtext('kind','') if pub is not None else ''
        patent_core=core('US'+docnum)
        patent=patent_core+kind
        date=pub.findtext('date','') if pub is not None else ''
        stats.update(title=title,screen=screen,screen_reason=why,Patent=patent,patent_core=patent_core,publication_date=date,mol_files=len(mols),xml_file=xml)
        if root.find('us-bibliographic-data-grant') is None:
            stats.update(status='excluded' if screen in ('agriculture','materials') else 'legacy_xml_review',elapsed_seconds=round(time.monotonic()-start,3))
            return dict(stats=stats,ready=[],review=[],structures=[],tables=[])
        metadata=extract_metadata(root)
        stats.update(metadata)
        if screen!='medical':
            stats.update(status='excluded' if screen in ('agriculture','materials') else 'screen_review',elapsed_seconds=round(time.monotonic()-start,3))
            return dict(stats=stats,ready=[],review=[],structures=[],tables=[])
        doc_metadata,bibliography=extract_documents(root,patent,archive,xml,document)
        metadata.update(doc_metadata);stats.update(doc_metadata)
        d=root.find('description')
        if d is None:raise ValueError('missing description')
        blocks=list(d)
        parent={c:p for p in d.iter() for c in p}
        cache={}
        links=defaultdict(list)
        def files_of(e):
            return [mols[Path(im.get('file','')).stem.upper()] for im in e.iter('img') if Path(im.get('file','')).stem.upper() in mols]
        def add_link(a,files,method,locator,name=''):
            if not a or len(files)!=1:return
            mol=load_mol(z,files[0],cache)
            if name and mol['structure_status']=='concrete':
                conflicts=name_element_conflicts(name,mol['element_counts'])
                mol['structure_name']=name
                # A later unnamed link must not erase a detected source contradiction.
                if conflicts:mol['structure_name_check']='element_mismatch:'+','.join(conflicts)
                elif not mol['structure_name_check'].startswith('element_mismatch'):mol['structure_name_check']='no_element_presence_conflict'
            if mol['structure_status']=='concrete':links[a].append(dict(mol,structure_match=method,structure_locator=locator))
        # Exact same-row compound number/structure associations have highest specificity.
        tables=[]
        for i,b in enumerate(blocks):
            for t in b.iter('table'):
                headers,rows=expanded_rows(t)
                p=parent.get(t)
                loc=(p.get('id','') if p is not None else '') or b.get('id','') or f'table-{len(tables)+1}'
                tables.append(dict(element=t,block=i,headers=headers,rows=rows,locator=loc))
                alias_cols=[j for j,h in enumerate(headers) if re.search(r'example|compound|cmpd|\bno\b',h,re.I) and not ACTIVITY.search(h) and not re.search(r'name|mass',h,re.I)]
                name_cols=[j for j,h in enumerate(headers) if re.search(r'\bname\b',h,re.I)]
                for ri,row in enumerate(rows):
                    compound_name=row[name_cols[0]]['text'] if len(name_cols)==1 else ''
                    for ac in alias_cols:
                        a=alias(row[ac]['text'])
                        refs=[mols[Path(f).stem.upper()] for c in row for f in c['files'] if Path(f).stem.upper() in mols]
                        add_link(a,refs,'explicit_same_table_row',f"{loc}/{row[0]['row_locator']}",compound_name)
                    refs=[mols[Path(f).stem.upper()] for c in row for f in c['files'] if Path(f).stem.upper() in mols]
                    if len(refs)==1:
                        add_link('diagram:'+Path(refs[0]).stem,refs,'explicit_same_table_row',f"{loc}/{row[0]['row_locator']}",compound_name)
        # A product diagram immediately after a single numbered Example heading.
        for i,b in enumerate(blocks):
            s=text(b)
            eligible=b.tag=='heading' or (b.tag=='p' and len(s)<150 and b.find('b') is not None)
            m=EXAMPLE_HEADING.fullmatch(s) if eligible else None
            if not m:continue
            own=files_of(b)
            if own:
                add_link(m[1].lower(),own,'single_example_heading_diagram',b.get('id',''))
                continue
            for nxt in blocks[i+1:i+7]:
                nt=text(nxt)
                if EXAMPLE_HEADING.match(nt) or re.match(r'(?:Step\s+[A-Z0-9]|Intermediate|Scheme|Biolog)',nt,re.I):break
                refs=files_of(nxt)
                if refs:
                    add_link(m[1].lower(),refs,'single_example_heading_first_diagram',b.get('id','')+' -> '+nxt.get('id',''))
                    break
                if nxt.tag=='p' and len(nt)>180:break
        resolved={}
        ambiguous=set()
        for a,entries in links.items():
            keys={e['InChIKey'] for e in entries}
            if len(keys)==1:resolved[a]=sorted(entries,key=lambda e:(e['structure_match']!='explicit_same_table_row',e['structure_file']))[0]
            else:ambiguous.add(a)
        for t in tables:
            h,rows=t['headers'],t['rows']
            actcols=[j for j,s in enumerate(h) if ACTIVITY.search(s)]
            i=t['block']
            preceding=[]
            for b in reversed(blocks[max(0,i-12):i]):
                if b.find('.//table') is not None:break
                preceding.append(text(b))
            ctx=' '.join(reversed(preceding))
            loc=t['locator']
            if not actcols:
                if ACTIVITY.search(text(t['element'])) or re.search(r'biological activity',text(t['element']),re.I):
                    table_rows.append(dict(Patent=patent,archive=archive,table_locator=loc,reason='activity_table_not_supported',headers=h,rows=[[c['text'] for c in r] for r in rows],assay_context=ctx))
                continue
            aliascols=[j for j,s in enumerate(h) if re.search(r'example|compound|cmpd|\bno\b',s,re.I) and j not in actcols and not re.search(r'name',s,re.I)]
            targetcols=[j for j,s in enumerate(h) if re.search(r'kinase|cell line|target|receptor|enzyme',s,re.I) and j not in actcols]
            for ri,row in enumerate(rows):
                if not any(c['text'] for c in row):continue
                for ac in actcols:
                    raw=row[ac]['text']
                    if not raw:continue
                    # A header row in body is not an observation.
                    if ACTIVITY.search(raw):continue
                    a=alias(row[aliascols[0]]['text']) if len(aliascols)==1 else None
                    same_row_files=[mols[Path(f).stem.upper()] for c in row for f in c['files'] if Path(f).stem.upper() in mols]
                    if not a and len(same_row_files)==1:a='diagram:'+Path(same_row_files[0]).stem
                    endpoint=ACTIVITY.search(h[ac])[1].replace(' ','')
                    endpoint=endpoint[:2].upper()+endpoint[2:] if endpoint.lower() not in ('ki','kd') else 'K'+endpoint[1].lower()
                    target=''
                    target_evidence=''
                    if len(targetcols)==1:
                        target=row[targetcols[0]]['text'];target_evidence='explicit_target_row'
                    else:
                        prefix=ACTIVITY.sub('',h[ac])
                        prefix=UNIT.sub('',prefix)
                        prefix=re.sub(r'\b(?:Enzyme|Cell|ATP|inhibition of|micromolar|nanomolar|millimolar)\b',' ',prefix,flags=re.I)
                        prefix=re.sub(r'[()\[\]:%]',' ',prefix).strip()
                        prefix=re.sub(r'\s+',' ',prefix).strip(' /-')
                        if prefix and not re.search(r'compound|example|staurosporine|\bex\.|activity|value',prefix,re.I):
                            target=prefix;target_evidence='explicit_activity_column'
                    # Compound-per-column tables: parse the explicit Example column, never the control column.
                    col_example=re.search(r'(?:Example|Ex\.)\s*#?\s*(\d+(?:-\d+)?)',h[ac],re.I)
                    if col_example:a=col_example[1]
                    reasons=[]
                    if re.search(r'staurosporine|control|reference compound',h[ac],re.I):reasons.append('control_column_not_linked_to_example')
                    if any(re.search(r'\bstaurosporine\b|positive control|negative control|reference compound',c['text'],re.I) for j,c in enumerate(row) if j not in actcols):reasons.append('control_row_not_linked_to_example')
                    if not a:
                        # Single-compound table needs an explicit statement in immediately preceding prose.
                        aliases=set(re.findall(r'compound of Example\s*#?\s*(\d+(?:-\d+)?)',ctx,re.I))
                        if len(aliases)==1 and len(targetcols)==1 and not aliascols and not reasons:a=next(iter(aliases))
                        else:reasons.append('compound_alias_missing_or_ambiguous')
                    mol=resolved.get(a or '')
                    if a in ambiguous:reasons.append('conflicting_structure_for_alias')
                    elif not mol:reasons.append('no_explicit_concrete_structure_link')
                    elif cache[mol['structure_file']]['structure_name_check'].startswith('element_mismatch'):reasons.append('structure_name_element_mismatch')
                    if not target:reasons.append('target_not_explicit')
                    u=unit(h[ac]);ue='activity_column_unit' if u else ''
                    relation,value,lower,upper=measurement(raw)
                    cellunit=unit(raw)
                    if cellunit:
                        if u and u!=cellunit:reasons.append('conflicting_units')
                        u=cellunit;ue='activity_cell_unit'
                    # Only a direct IC50/EC50-unit statement may fill a missing column unit.
                    if not u:
                        statements=re.findall(r'(?:IC\s*50|EC\s*50|Ki|Kd)(?:\s+values?)?[^.;]{0,35}?(?:in|\()\s*([uµμ]M|nM|pM|mM)\b',ctx,re.I)
                        found={v.replace('µ','u').replace('μ','u') for v in statements}
                        if len(found)==1:u=next(iter(found));ue='explicit_endpoint_unit_in_context'
                    if endpoint.lower().startswith('p'):
                        if u:reasons.append('log_endpoint_has_concentration_unit')
                        u='';ue='dimensionless_log_endpoint';reasons.append('log_endpoint_review')
                    elif not u:reasons.append('activity_unit_missing')
                    if relation=='text':reasons.append('qualitative_or_unparsed_value')
                    elif relation=='range':reasons.append('interval_requires_extended_schema')
                    elif relation in ('~','≈'):reasons.append('approximate_value_review')
                    value_nm=float(Decimal(str(value))*UNITS[u]) if value is not None and u in UNITS else None
                    rid='USPTO-'+hashlib.sha256(f'{archive}|{xml}|{loc}|{ri+1}|{ac}|{a}|{raw}'.encode()).hexdigest()[:24]
                    record={k:'' for k in SAR_COLUMNS+EXTRA_COLUMNS}
                    record.update(source='USPTO_MOL',Patent=patent,patent_core=patent_core,Target=target,**{'Activity Type':endpoint},Relation=relation,Value=value,value_nM=value_nm,Unit=u,source_record_id=rid,compound_alias=a or '',raw_value=raw,value_lower=lower,value_upper=upper,activity_header=h[ac],target_evidence=target_evidence,unit_evidence=ue,assay_context=ctx,activity_locator=f"{loc}/{row[0]['row_locator']}/column-{ac+1}",archive=archive,xml_file=xml,publication_date=date,extraction_status='review' if reasons else 'paired',review_reason=';'.join(dict.fromkeys(reasons)))
                    record.update(metadata)
                    if mol:
                        for k in ['SMILES','SMILES_original','InChIKey','structure_file','structure_match','structure_locator']:record[k]=mol[k]
                        for k in ['structure_name','structure_name_check']:record[k]=cache[mol['structure_file']][k]
                    # Species must be adjacent to the target, not just mentioned elsewhere in the context.
                    target_pattern=re.escape(target).replace(r'\ ',r'\s*')
                    organisms={o for adjective,o in [('human','Homo sapiens'),('mouse','Mus musculus'),('rat','Rattus norvegicus')] if target and re.search(r'\b'+adjective+r'\s+(?:c-)?'+target_pattern+r'\b',ctx,re.I)}
                    if len(organisms)==1:record['organism']=next(iter(organisms))
                    (review if reasons else ready).append(record)
        # Preserve all structural observations, including unidentified and generic structures.
        aliases_by_file=defaultdict(list)
        for a,v in resolved.items():
            if not a.startswith('diagram:'):aliases_by_file[v['structure_file']].append(a)
        for filename in sorted(mols.values()):
            row=load_mol(z,filename,cache).copy()
            row.update(Patent=patent,archive=archive,compound_aliases='|'.join(sorted(aliases_by_file[filename])),role='numbered_example' if aliases_by_file[filename] else 'unassigned_diagram')
            structure_rows.append(row)
        stats.update(status='processed',tables=len(tables),explicit_aliases=sum(not a.startswith('diagram:') for a in resolved),ambiguous_aliases=len(ambiguous),paired_rows=len(ready),review_rows=len(review),unparsed_activity_tables=len(table_rows),structure_status_counts=dict(Counter(s['structure_status'] for s in structure_rows)),name_element_mismatches=sum(s['structure_name_check'].startswith('element_mismatch') for s in structure_rows),elapsed_seconds=round(time.monotonic()-start,3))
    return dict(stats=stats,ready=ready,review=review,structures=structure_rows,tables=table_rows,documents=dict(bibliography=bibliography))

def atomic_json(path,obj):
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(obj,ensure_ascii=False,separators=(',',':')))
    tmp.replace(path)

def worker(job):
    root,out,spec,force=job
    path=Path(root)/spec['archive']
    try:stat=path.stat()
    except OSError as exc:return dict(checkpoint=None,stats=dict(archive=spec['archive'],patent_core=spec.get('patent_core',''),status='error',error=repr(exc))),False
    code_hash=hashlib.sha256(b''.join(Path(__file__).with_name(n).read_bytes() for n in ['uspto_sar.py','select_samples.py','patent_metadata.py','patent_documents.py','xml_source.py'])).hexdigest()
    sig=hashlib.sha256((code_hash+'|'+json.dumps(spec,sort_keys=True)+'|'+str(stat.st_size)+'|'+str(stat.st_mtime_ns)).encode()).hexdigest()
    cp=Path(out)/'checkpoints'/ (hashlib.sha256(spec['archive'].encode()).hexdigest()[:20]+'.json')
    if cp.exists() and not force:
        try:
            obj=json.loads(cp.read_text())
            if obj.get('signature')==sig:return dict(checkpoint=str(cp),stats=obj['stats']),True
        except (OSError,json.JSONDecodeError,KeyError):pass
    try:result=extract(path,spec['archive'],spec)
    except Exception as e:result=dict(stats=dict(archive=spec['archive'],patent_core=spec.get('patent_core',''),status='error',error=repr(e)),ready=[],review=[],structures=[],tables=[])
    result['signature']=sig
    # Error checkpoints are reported and retried on a later run.
    if result['stats']['status']!='error':atomic_json(cp,result)
    return dict(checkpoint=str(cp) if result['stats']['status']!='error' else None,stats=result['stats']),False

class RecordWriter:
    """Stream one patent at a time, then atomically replace finished output files."""
    def __init__(self,out,name,columns):
        import pyarrow as pa
        import pyarrow.parquet as pq
        self.pa=pa;self.columns=columns;self.name=name;self.out=out
        self.numeric={'Value','value_nM','value_lower','value_upper','elapsed_seconds'}
        self.integer={'mol_files','tables','explicit_aliases','ambiguous_aliases','paired_rows','review_rows','candidate_rows','quality_check_failures','unparsed_activity_tables'}|DOCUMENT_INTEGER_COLUMNS
        self.schema=pa.schema([pa.field(k,pa.float64() if k in self.numeric else pa.int64() if k in self.integer else pa.string()) for k in columns])
        self.handle=(out/(name+'.csv.tmp')).open('w',encoding='utf-8-sig',newline='')
        self.csv=csv.DictWriter(self.handle,fieldnames=columns);self.csv.writeheader()
        self.parquet=pq.ParquetWriter(out/(name+'.parquet.tmp'),self.schema,compression='zstd')
        self.buffer=[]
    def append(self,rows):
        if not rows:return
        clean=[]
        for row in rows:
            item={k:row.get(k) for k in self.columns}
            item={k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in item.items()}
            item={k:None if k in self.numeric|self.integer and v=='' else v for k,v in item.items()}
            self.csv.writerows([item])
            clean.append({k:v if k in self.numeric|self.integer else str(v) if v is not None else None for k,v in item.items()})
        self.buffer.extend(clean)
        if len(self.buffer)>=2048:self.flush()
    def flush(self):
        if self.buffer:self.parquet.write_table(self.pa.Table.from_pylist(self.buffer,schema=self.schema));self.buffer=[]
    def close(self):
        self.flush()
        self.handle.close();self.parquet.close()
        for ext in ('csv','parquet'):(self.out/(self.name+'.'+ext+'.tmp')).replace(self.out/(self.name+'.'+ext))

def write_table(out,name,rows,columns):
    """Explicit typed Parquet schema; CSV preserves missing as empty, never zero."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    numeric={'Value','value_nM','value_lower','value_upper','elapsed_seconds'}
    integer={'mol_files','tables','explicit_aliases','ambiguous_aliases','paired_rows','review_rows','unparsed_activity_tables'}|DOCUMENT_INTEGER_COLUMNS
    def scalar(v):return json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v
    rows=[{k:scalar(r.get(k)) for k in columns} for r in rows]
    with (out/(name+'.csv')).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(rows)
    fields=[pa.field(k,pa.float64() if k in numeric else pa.int64() if k in integer else pa.string()) for k in columns]
    typed=[{k:(v if k in numeric|integer else str(v) if v is not None else None) for k,v in row.items()} for row in rows]
    pq.write_table(pa.Table.from_pylist(typed,schema=pa.schema(fields)),out/(name+'.parquet'),compression='zstd')

def run(args):
    start=time.monotonic()
    manifest=json.loads(args.manifest.read_text())
    specs=manifest['samples']
    if len(specs)>100 and not args.full:raise SystemExit('More than 100 patents requires --full. Review the sample before a full run.')
    if args.root.resolve() in args.output.resolve().parents:raise SystemExit('Keep output outside the raw dataset tree.')
    args.output.mkdir(parents=True,exist_ok=True);(args.output/'checkpoints').mkdir(exist_ok=True)
    for folder in ['bibliography']:(args.output/folder).mkdir(exist_ok=True)
    if len({s['archive'] for s in specs})!=len(specs):raise SystemExit('Duplicate archive in manifest.')
    results={};hits=0
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        it=iter(specs);pending={}
        def submit_next():
            s=next(it,None)
            if s is not None:pending[pool.submit(worker,(str(args.root),str(args.output),s,args.force))]=s
        for _ in range(min(len(specs),2*args.workers)):submit_next()
        while pending:
            done,_=wait(pending,return_when=FIRST_COMPLETED)
            for fut in done:
                pending.pop(fut);obj,hit=fut.result();results[obj['stats']['archive']]=obj;hits+=hit
                s=obj['stats'];print(json.dumps(dict(done=len(results),total=len(specs),patent=s.get('patent_core'),status=s['status'],paired=s.get('paired_rows',0),review=s.get('review_rows',0),cached=hit),ensure_ascii=False),flush=True)
                submit_next()
    writers={name:RecordWriter(args.output,name,cols) for name,cols in [
        ('SAR',SAR_COLUMNS),('SAR_enriched',SAR_COLUMNS+METADATA_COLUMNS),('SAR_with_evidence',SAR_COLUMNS+EXTRA_COLUMNS+METADATA_COLUMNS),('review',SAR_COLUMNS+EXTRA_COLUMNS+METADATA_COLUMNS),
        ('structures',['Patent','compound_aliases','role','structure_status','SMILES','SMILES_original','InChIKey','structure_file','structure_name','structure_name_check','archive']),
        ('patents',list(dict.fromkeys(['Patent','patent_core','title','publication_date','selection','screen','screen_reason','status','mol_files','tables','explicit_aliases','ambiguous_aliases','paired_rows','review_rows','unparsed_activity_tables','structure_status_counts','elapsed_seconds','archive','xml_file','error']+METADATA_COLUMNS)))]}
    db=sqlite3.connect(args.output/'structure_counts.sqlite')
    db.execute('DROP TABLE IF EXISTS unique_structures');db.execute('CREATE TABLE unique_structures (inchikey TEXT PRIMARY KEY) WITHOUT ROWID')
    reasons=Counter();statuses=Counter();nready=nreview=nstruct=ntables=npat=0
    metadata_coverage=Counter();metadata_reviews=Counter();metadata_statuses=Counter();document_reviews=Counter();sidecar_files=[]
    with (args.output/'unparsed_activity_tables.jsonl.tmp').open('w') as f:
        for s in specs:
            small=results[s['archive']];stats=small['stats'];statuses[stats['status']]+=1
            o=json.loads(Path(small['checkpoint']).read_text()) if small['checkpoint'] else dict(ready=[],review=[],structures=[],tables=[])
            for folder,doc in o.get('documents',{}).items():
                relative=folder+'/'+stats['Patent']+'.json'
                atomic_json(args.output/relative,doc);sidecar_files.append(relative)
            obs=o['ready']+o['review'];assert len({r['source_record_id'] for r in obs})==len(obs)
            for k in ('ready','review'):o[k].sort(key=lambda r:(r['activity_locator'],r['source_record_id']))
            writers['SAR'].append(o['ready']);writers['SAR_enriched'].append(o['ready']);writers['SAR_with_evidence'].append(o['ready']);writers['review'].append(o['review']);writers['structures'].append(o['structures']);writers['patents'].append([stats])
            if stats['status']=='processed':
                metadata_coverage.update(k for k in METADATA_COLUMNS if stats.get(k) not in ('',None,[],{}))
                metadata_reviews.update(k for k in stats.get('metadata_review_reason','').split(';') if k)
                document_reviews.update(k for k in stats.get('document_metadata_review_reason','').split(';') if k)
                metadata_statuses.update([stats.get('is_pct',''),stats.get('technical_field_source','')])
            for t in o['tables']:f.write(json.dumps(t,ensure_ascii=False)+'\n')
            db.executemany('INSERT OR IGNORE INTO unique_structures VALUES (?)',[(r['InChIKey'],) for r in o['structures'] if r['structure_status']=='concrete']);db.commit()
            reasons.update(reason for r in o['review'] for reason in r['review_reason'].split(';'))
            nready+=len(o['ready']);nreview+=len(o['review']);nstruct+=len(o['structures']);ntables+=len(o['tables']);npat+=bool(o['ready'])
    for w in writers.values():w.close()
    (args.output/'unparsed_activity_tables.jsonl.tmp').replace(args.output/'unparsed_activity_tables.jsonl')
    unique=db.execute('SELECT COUNT(*) FROM unique_structures').fetchone()[0];db.close()
    report=dict(version=VERSION,rdkit=rdBase.rdkitVersion,manifest_sha256=hashlib.sha256(args.manifest.read_bytes()).hexdigest(),sample_selection=manifest if len(specs)<=100 else dict(samples_count=len(specs)),attempted=len(specs),statuses=dict(statuses),paired_patents=npat,paired_observations=nready,review_observations=nreview,structure_observations=nstruct,concrete_unique_structures=unique,review_reasons=dict(reasons),unparsed_activity_tables=ntables,cache_hits=hits,elapsed_seconds=round(time.monotonic()-start,3),workers=args.workers,limitations=['No image OCR or IUPAC name-to-structure inference.','Legacy SGML is screened but requires a separate extraction adapter.','Unitless activity tables remain in review even when a likely unit exists in ChEMBL.','UniProt and gene remain empty; no species-specific target registry has been applied.','No salt stripping or tautomer standardization; full stereochemical InChIKey.','Numbers in prose and categorical tables without explicit endpoint columns remain outside strict SAR output.'])
    atomic_json(args.output/'run_report.json',report)
    atomic_json(args.output/'metadata_report.json',dict(version=VERSION,processed_patents=statuses['processed'],field_coverage=dict(metadata_coverage),review_reasons=dict(metadata_reviews),source_statuses=dict(metadata_statuses)))
    atomic_json(args.output/'document_report.json',dict(version=VERSION,processed_patents=statuses['processed'],bibliography_files=sum(s.startswith('bibliography/') for s in sidecar_files),review_reasons=dict(document_reviews),sidecar_files=sidecar_files))
    print(json.dumps({k:v for k,v in report.items() if k!='sample_selection'},ensure_ascii=False,indent=2))
    if statuses['error']:raise SystemExit('One or more archives failed; see patents.csv. Successful checkpoints are retained.')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--workers',type=int,default=4)
    p.add_argument('--force',action='store_true')
    p.add_argument('--full',action='store_true',help='explicit opt-in for >100 patents, after reviewing sample')
    run(p.parse_args())

if __name__=='__main__':main()
