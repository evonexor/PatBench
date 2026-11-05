"""Bibliography for all domains; metadata failures never suppress chemical files."""
import json
import re
from html import unescape
from pathlib import PurePosixPath
from patent_metadata import extract_metadata, METADATA_COLUMNS, document, date_iso
from patent_documents import extract_documents, DOCUMENT_INTEGER_COLUMNS
from uspto_sar import screening, text
from xml_source import parse_xml, select_patent_xml, _entities, _reserved, _named, _protected

def legacy_value(source, tag):
    found = re.search(r'<'+tag+r'\b[^>]*>(.*?)</'+tag+r'>', source, re.I|re.S)
    if not found:return ''
    return re.sub(r'\s+', ' ', unescape(re.sub(r'<[^>]+>', ' ', found[1]))).strip()

def bibliography(z, archive, fallback):
    meta={k:None if k in DOCUMENT_INTEGER_COLUMNS else '' for k in METADATA_COLUMNS}
    result={'Patent':fallback,'patent_core':fallback,'metadata_format':'unavailable','metadata_status':'unavailable','metadata_unresolved_entities':'','xml_file':'','screen':'uncertain','screen_reason':'No bibliographic classification available',**meta}
    try:xml=select_patent_xml(z, archive)
    except Exception:xml=None
    if xml:
        raw=z.read(xml);result['xml_file']=xml;unknown=[]
        try:root=parse_xml(raw)
        except Exception:
            # Explicit placeholders are confined to bibliography. No assay or MOL text is repaired.
            def replace(match):
                name=match[1].decode('ascii')
                if name in _reserved or name in _entities:return match[0]
                unknown.append(name);return ('[UNRESOLVED_ENTITY:'+name+']').encode('ascii')
            parts=_protected.split(raw)
            rewritten=b''.join(p if i%2 else _named.sub(replace,p) for i,p in enumerate(parts))
            try:root=parse_xml(rewritten)
            except Exception:root=None
        if root is not None:
            for e in root.iter():e.tag=e.tag.rsplit('}',1)[-1]
            if root.find('us-bibliographic-data-grant') is not None:
                result.update(extract_metadata(root))
                pub=document(root.find('us-bibliographic-data-grant/publication-reference/document-id'))
                patent=fallback+pub['kind']
                result.update(Patent=patent,patent_core=fallback,metadata_format='modern_xml',metadata_status='parsed')
                dm,bib=extract_documents(root,patent,archive,xml,document);result.update(dm)
                result['bibliography']=bib
                result['screen'],result['screen_reason']=screening(result['title'],result['abstract'],result['ipc_codes']+' '+result['cpc_codes'])
                if unknown:result.update(metadata_unresolved_entities=json.dumps(sorted(set(unknown))),metadata_status='parsed_with_entity_placeholders')
                result['_root']=root
                return result
            source=raw.decode('utf-8','replace');fmt='legacy_xml'
        else:source=raw.decode('utf-8','replace');fmt='legacy_xml_metadata_unparsed'
        result['metadata_unresolved_entities']=json.dumps(sorted(set(unknown))) if unknown else ''
        result['_root']=root
    else:
        legacy=next((i.filename for i in z.infolist() if i.filename.lower().endswith(('.sgm','.sgml'))),None)
        if not legacy:return result
        source=z.read(legacy).decode('utf-8','replace');fmt='legacy_sgml';result['xml_file']=legacy
    # Literal ST.32 bibliographic blocks, never cited reference dates or organizations.
    source=source.split('</SDOBI>',1)[0] if '</SDOBI>' in source else source
    title=legacy_value(source,'B540');kind=legacy_value(source,'B130')
    grant=date_iso(legacy_value(source,'B140'));filing=date_iso(legacy_value(source,'B220'))
    assignees=[]
    for block in re.findall(r'<B731\b[^>]*>(.*?)</B731>',source,re.S|re.I):
        name=legacy_value(block,'ONM')
        if name:assignees.append(name)
    priorities=[]
    for block in re.findall(r'<B32\d\b[^>]*>(.*?)</B32\d>',source,re.S|re.I):
        date=date_iso(legacy_value(block,'DATE'))
        if date:priorities.append(date)
    result.update(Patent=fallback+kind,metadata_format=fmt,metadata_status='legacy_partial',title=title,
        company_name='；'.join(dict.fromkeys(assignees)),assignee_names='；'.join(dict.fromkeys(assignees)),company_source='ST.32 B731/ONM' if assignees else '',
        grant_date=grant,grant_date_type='本件旧版授权文档日期' if grant else '',document_publication_date=grant,
        application_date=filing,application_number='US'+legacy_value(source,'B210') if legacy_value(source,'B210') else '',
        document_country=legacy_value(source,'B190'),document_kind=kind,
        priority_date=min(priorities) if priorities else '',priority_date_status='ST.32字面优先权日期，未核对完整专利族' if priorities else '旧版书目未提取到完整优先权数据',
        ipc_codes='；'.join(legacy_value(source,t) for t in ('B511','B512') if legacy_value(source,t)),
        abstract=legacy_value(source,'SDOAB'),metadata_review_reason='legacy_partial_metadata;some_fields_not_parsed',
        is_pct='旧版文档尚未完整解析PCT信息')
    # Abstract usually follows SDOBI and must be taken from the full source file.
    full=z.read(result['xml_file']).decode('utf-8','replace')
    result['abstract']=legacy_value(full,'SDOAB')
    result['screen'],result['screen_reason']=screening(title,result['abstract'],result['ipc_codes'])
    return result
