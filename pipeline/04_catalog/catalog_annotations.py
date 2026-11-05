"""Reported aliases/names/ranges and explicitly numbered sequence/UniProt links."""
from collections import defaultdict
import hashlib,json,re
from pathlib import PurePosixPath
import pyarrow as pa
from uspto_sar import expanded_rows,alias,ACTIVITY,text,SAR_COLUMNS
from xml_source import parse_xml

FIELDS=['Sequence','UniProt ID','original_alias','compound_IUPAC_name','original_range','sequence_ids','sequence_link_status','compound_name_status','annotation_evidence']
ACCESSION=re.compile(r'\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?\b')
SEQREF=re.compile(r'\bSEQ\.?\s*ID\s*NO\.?\s*[:#]?\s*(\d+)',re.I)
AA=dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL ASX GLX XAA SEC PYL'.split(),'ARNDCQEGHILKMFPSTWYVBZXUO'))
SEQ_COLUMNS=['sequence_source_id','Patent','archive','source_file','sequence_id','Sequence','sequence_original','sequence_type','organism','declared_length','actual_length','sequence_status','source_format','source_locator','sequence_description','UniProt ID']
SEQ_SCHEMA=pa.schema([pa.field(k,pa.int64() if k in ('declared_length','actual_length') else pa.string()) for k in SEQ_COLUMNS])
ACTIVITY_COLUMNS=list(dict.fromkeys(SAR_COLUMNS+FIELDS+['compound_alias','raw_value','value_lower','value_upper','review_reason','extraction_status','source_structure_id','structure_file','structure_locator','activity_locator','archive','is_formal_SAR']))
ACTIVITY_SCHEMA=pa.schema([pa.field(k,pa.float64() if k in ('Value','value_nM','value_lower','value_upper') else pa.bool_() if k=='is_formal_SAR' else pa.string()) for k in ACTIVITY_COLUMNS])
UNIPROT_SCHEMA=pa.schema([pa.field(k,pa.string()) for k in ['Patent','archive','UniProt ID','source_file','source_locator','mention_text','link_status']])

def explicit_uniprot(s):
    ids=[]
    for label in re.finditer(r'UniProt(?:KB)?|Swiss[ -]?Prot',s,re.I):
        tail=s[label.end():label.end()+100]
        ids.extend(ACCESSION.findall(tail))
    return sorted(set(ids))

def normalize_sequence(raw,kind):
    pieces=[];unknown=False
    for line in raw.splitlines():
        words=re.findall(r'[A-Za-z]+',line)
        if not words:continue
        if 'DNA' in kind.upper() or 'RNA' in kind.upper():
            # ST.25 mixed listings may interleave protein translations; exclude those lines.
            if all(re.fullmatch(r'[A-Z][a-z]{2}',w) and w.upper() in AA for w in words):continue
            if all(re.fullmatch(r'[ACGTURYSWKMBDHVNXacgturyswkmbdhvnx]+',w) for w in words):pieces.extend(w.upper() for w in words)
        else:
            if all(w.upper() in AA for w in words):pieces.extend(AA[w.upper()] for w in words)
            elif all(re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWYBXZUO]+',w.upper()) for w in words):pieces.extend(w.upper() for w in words)
            else:unknown=True
    return ''.join(pieces),unknown

def sequence_row(patent,archive,filename,sid,raw,kind,organism,length,fmt,locator,description=''):
    if fmt=='ST26_XML':
        seq=re.sub(r'\s+','',raw).upper();unknown=not bool(re.fullmatch(r'[A-Z]*',seq))
    else:seq,unknown=normalize_sequence(raw,kind)
    declared=int(length.strip()) if str(length).strip().isdigit() else None
    if fmt=='ST26_XML' and raw.strip()=='000' and declared==0:seq='';unknown=False
    status='nonstandard_sequence_symbols_review' if unknown else 'length_mismatch_review' if declared is not None and len(seq)!=declared else 'parsed' if seq else 'empty_or_skipped_sequence'
    return dict(sequence_source_id='USPTO-SEQUENCE:'+hashlib.sha256(f'{archive}|{filename}|{sid}|{locator}'.encode()).hexdigest(),Patent=patent,archive=archive,source_file=filename,sequence_id=str(sid),Sequence=seq,sequence_original=raw,sequence_type=kind,organism=organism,declared_length=declared,actual_length=len(seq),sequence_status=status,source_format=fmt,source_locator=locator,sequence_description=description,**{'UniProt ID':';'.join(explicit_uniprot(description))})

def parse_sequences(data,filename,patent,archive,root=None):
    source=data.decode('utf-8','replace');rows=[]
    if root is not None or re.search(r'<(?:ST26SequenceListing|SequenceData|sequence-cwu|s400)\b',source,re.I):
        root=root if root is not None else parse_xml(data)
        for e in root.iter():e.tag=e.tag.rsplit('}',1)[-1]
        if list(root.iter('SequenceData')):
            for i,item in enumerate(root.iter('SequenceData')):
                insd=item.find('INSDSeq')
                if insd is None:continue
                kind=text(insd.find('INSDSeq_moltype'));raw=text(insd.find('INSDSeq_sequence'))
                rows.append(sequence_row(patent,archive,filename,item.get('sequenceIDNumber',''),raw,kind,text(insd.find('INSDSeq_organism')),text(insd.find('INSDSeq_length')),'ST26_XML',f'SequenceData-{i+1}',text(insd.find('INSDSeq_feature-table'))))
        else:
            for container in root.iter():
                last=None
                for i,item in enumerate(container):
                    if item.tag=='s200':last=item
                    elif item.tag=='s400' and last is not None:
                        sid=text(last.find('s210'));raw=''.join(item.itertext())
                        rows.append(sequence_row(patent,archive,filename,sid,raw,text(last.find('s212')),text(last.find('s213')),text(last.find('s211')),'USPTO_SEQUENCE_XML',f's400-{i+1}',text(last.find('s220'))))
                        last=None
    elif '<210>' in source and '<400>' in source:
        for i,record in enumerate(re.split(r'(?=<210>\s*\d+)',source)[1:]):
            find=lambda code:(re.search(r'<'+code+r'>\s*([^\r\n<]+)',record).group(1).strip() if re.search(r'<'+code+r'>\s*([^\r\n<]+)',record) else '')
            start=re.search(r'<400>\s*\d+[^\r\n]*[\r\n]+',record)
            if not start:continue
            raw=re.split(r'<\d{3}>',record[start.end():],maxsplit=1)[0]
            rows.append(sequence_row(patent,archive,filename,find('210'),raw,find('212'),find('213'),find('211'),'ST25_TEXT',f'record-{i+1}',find('223')))
    return rows

def extract(z,patent,archive,xml_file,mol_files,old,root=None):
    sequences=[];errors=[]
    if root is None and xml_file and xml_file.lower().endswith('.xml'):
        try:root=parse_xml(z.read(xml_file));[setattr(e,'tag',e.tag.rsplit('}',1)[-1]) for e in root.iter()]
        except Exception as exc:errors.append('annotation_document_xml:'+type(exc).__name__)
    # Only explicit listing formats qualify, never arbitrary letter runs in prose.
    for info in z.infolist():
        if info.is_dir() or not info.filename.lower().endswith(('.xml','.txt','.seq','.sgm','.sgml')):continue
        if info.filename==xml_file:
            if root is None:continue
            if not list(root.iter('s400')) and not list(root.iter('SequenceData')):continue
        else:
            with z.open(info) as f:prefix=f.read(4096)
            if not re.search(br'ST26SequenceListing|sequence-cwu|<s400>|<210>|<400>',prefix,re.I) and 'seqlst' not in info.filename.lower():continue
        try:sequences.extend(parse_sequences(z.read(info),info.filename,patent,archive,root if info.filename==xml_file else None))
        except Exception as exc:errors.append('sequence_file:'+info.filename+':'+type(exc).__name__)
    by_seq=defaultdict(list)
    for r in sequences:by_seq[r['sequence_id']].append(r)
    by_stem=defaultdict(list)
    for filename in mol_files:by_stem[PurePosixPath(filename).stem.upper()].append(filename)
    annotations=defaultdict(list);raw_aliases={};mentions=[]
    if root is not None:
        for i,node in enumerate(e for e in root.iter() if e.tag in ('p','PARA')):
            s=text(node)
            for identifier in explicit_uniprot(s):
                pos=s.find(identifier);mentions.append({'Patent':patent,'archive':archive,'UniProt ID':identifier,'source_file':xml_file,'source_locator':node.get('id') or f'paragraph-{i+1}','mention_text':s[max(0,pos-120):pos+120],'link_status':'patent_mention_not_assigned_to_compound_or_assay'})
        parent={c:p for p in root.iter() for c in p}
        for t in root.iter('table'):
            try:headers,rows=expanded_rows(t)
            except Exception:continue
            p=parent.get(t);loc=(p.get('id','') if p is not None else '') or t.get('id','')
            alias_cols=[i for i,h in enumerate(headers) if re.search(r'example|compound|cmpd|\bno\b',h,re.I) and not re.search(r'name',h,re.I) and not ACTIVITY.search(h)]
            names=[i for i,h in enumerate(headers) if re.search(r'\bname\b|IUPAC',h,re.I)]
            legacy_refs={}
            for gi,group in enumerate(list(t.findall('tgroup')) or [t]):
                for ri,row_node in enumerate(group.iter('row')):
                    refs=[]
                    for element in row_node.iter():
                        if element.tag.upper() in ('CHEMMOL','CHEMCDX','EMI'):
                            for value in element.attrib.values():
                                matched=by_stem[PurePosixPath(value).stem.upper()]
                                if len(matched)==1:refs.append(matched[0])
                    legacy_refs[f'tgroup-{gi+1}/row-{ri+1}']=refs
            for row in rows:
                refs=list(legacy_refs.get(row[0]['row_locator'],[]))
                for cell in row:
                    for file in cell['files']:
                        matched=by_stem[PurePosixPath(file).stem.upper()]
                        if len(matched)==1:refs.append(matched[0])
                refs=set(refs)
                original=row[alias_cols[0]]['text'] if len(alias_cols)==1 else ''
                if loc and original:raw_aliases[loc+'/'+row[0]['row_locator']]=(original,alias(original))
                if len(refs)!=1:continue
                name=row[names[0]]['text'] if len(names)==1 else ''
                sequence_cols=[i for i,h in enumerate(headers) if re.search(r'^(?:compound |peptide )?sequence\b',h,re.I) and not re.search(r'target|receptor|enzyme',h,re.I)]
                # A target's SEQ ID in an assay row is not the compound's sequence.
                context=' '.join(row[i]['text'] for i in names+sequence_cols)
                ids=sorted(set(SEQREF.findall(context)))
                candidate_seqs=[s for sid in ids for s in by_seq[sid] if s['sequence_status']=='parsed']
                seq_values={s['Sequence'] for s in candidate_seqs}
                annotations[next(iter(refs))].append(dict(original_alias=original,compound_IUPAC_name=name,compound_name_status='explicit_IUPAC_column' if names and 'iupac' in headers[names[0]].lower() else 'source_reported_name_IUPAC_unverified' if name else 'name_not_reported',Sequence=next(iter(seq_values)) if len(ids)==1 and len(seq_values)==1 else '',sequence_ids=json.dumps(ids),sequence_link_status='explicit_structure_row_sequence_id' if len(ids)==1 and len(seq_values)==1 else 'sequence_reference_unresolved' if ids else 'not_linked',**{'UniProt ID':''},original_range='',annotation_evidence=json.dumps({'table':loc,'row':row[0]['row_locator'],'name_header':headers[names[0]] if len(names)==1 else ''},ensure_ascii=False)))
    # Keep annotations per activity observation so targets/ranges cannot be mixed across assays.
    activities=[]
    formal={r['source_record_id'] for r in old.get('ready',[])}
    for r in old.get('ready',[])+old.get('review',[]):
        item={k:r.get(k) for k in ACTIVITY_COLUMNS};item.update({k:'' for k in FIELDS})
        item.update(is_formal_SAR=r['source_record_id'] in formal,original_range=r['raw_value'],sequence_link_status='not_linked')
        rowloc=r['activity_locator'].rsplit('/column-',1)[0];original,normalized=raw_aliases.get(rowloc,('',''))
        if normalized==r['compound_alias']:item['original_alias']=original
        else:
            match=re.search(r'(?:Example|Ex\.)\s*#?\s*'+re.escape(r['compound_alias'])+r'\b',r['activity_header'],re.I) if r['compound_alias'] else None
            if match:item['original_alias']=match[0]
        name=r.get('structure_name','')
        if name:item.update(compound_IUPAC_name=name,compound_name_status='source_reported_name_IUPAC_unverified')
        candidates=annotations[r['structure_file']]
        unique={a['compound_IUPAC_name'] for a in candidates if a['compound_IUPAC_name']}
        if len(unique)==1 and not name:item.update(compound_IUPAC_name=next(iter(unique)),compound_name_status='source_reported_name_IUPAC_unverified')
        seqs={a['Sequence'] for a in candidates if a['Sequence']}
        if len(seqs)==1:item.update(Sequence=next(iter(seqs)),sequence_link_status='explicit_structure_row_sequence_id')
        # Accessions in a patent elsewhere are not automatically assigned to every assay.
        context=(r['Target']+' '+r['activity_header']);ids=explicit_uniprot(context)
        item['UniProt ID']=';'.join(ids);item['UniProt']=item['UniProt ID']
        if ids:item['annotation_evidence']=json.dumps({'uniprot_scope':'explicit_target_or_activity_header','text':context},ensure_ascii=False)
        activities.append(item)
    return dict(sequences=sequences,structures=dict(annotations),activities=activities,uniprot_mentions=mentions,errors=errors)
