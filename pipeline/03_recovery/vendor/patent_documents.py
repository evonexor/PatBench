"""Patent-level documents stored once; compact fields can be joined to SAR.

Preserve explicitly recorded roles, related applications and citations.
No inferred current ownership or complete family information.
"""
import re
import xml.etree.ElementTree as ET
from pathlib import PurePosixPath

DOCUMENT_COLUMNS = [
    'abstract','abstract_language','application_number','application_date',
    'application_type','document_country','document_kind','applicant_all_names',
    'applicant_address_countries','inventor_names','inventor_source','inventor_count',
    'parent_application_numbers',
    'parent_application_relationships','patent_citation_count','npl_citation_count',
    'bibliography_json_path','document_metadata_review_reason',
]
DOCUMENT_INTEGER_COLUMNS = {
    'inventor_count','patent_citation_count','npl_citation_count',
}

def plain(e):
    return re.sub(r'\s+',' ',''.join(e.itertext())).strip() if e is not None else ''

def unique(values):return list(dict.fromkeys(v for v in values if v))

def person(e):
    a=e.find('addressbook');organization=plain(a.find('orgname')) if a is not None else ''
    personal=' '.join(plain(a.find(k)) for k in ['first-name','middle-name','last-name','suffix'] if a is not None and plain(a.find(k)))
    return dict(name=organization or personal,name_type='organization' if organization else 'person' if personal else 'missing',xml_role=e.tag,attributes=dict(e.attrib),address_country=a.findtext('address/country','') if a is not None else '')

def extract_documents(root,patent,archive,xml_file,document):
    b=root.find('us-bibliographic-data-grant');m={k:None if k in DOCUMENT_INTEGER_COLUMNS else '' for k in DOCUMENT_COLUMNS};notes=[]
    abstract=root.find('abstract');app=b.find('application-reference')
    application=document(app.find('document-id')) if app is not None else document(None)
    pub=document(b.find('publication-reference/document-id'))
    language=(abstract.get('lang') or abstract.get('{http://www.w3.org/XML/1998/namespace}lang') or root.get('lang') or root.get('{http://www.w3.org/XML/1998/namespace}lang') or '').lower() if abstract is not None else ''
    m.update(abstract=plain(abstract),abstract_language=language,application_number=application['number'],application_date=application['date'],application_type=app.get('appl-type','') if app is not None else '',document_country=pub['country'],document_kind=pub['kind'])
    applicants=b.findall('us-parties/us-applicants/us-applicant')+b.findall('parties/applicants/applicant')
    inventors=b.findall('us-parties/inventors/inventor')+b.findall('parties/inventors/inventor')+b.findall('us-parties/us-inventors/us-inventor')
    explicit_inventors=[e for e in applicants if e.get('app-type')=='applicant-inventor']
    inventor_source='XML发明人节点' if inventors else 'XML明确标记的申请人兼发明人' if explicit_inventors else '未记载明确发明人角色'
    inventor_records=[person(e) for e in inventors or explicit_inventors]
    applicant_records=[person(e) for e in applicants]
    assignee_records=[person(e) for e in b.findall('assignees/assignee')]
    inventor_names=unique(p['name'] for p in inventor_records)
    m.update(applicant_all_names='；'.join(unique(p['name'] for p in applicant_records)),applicant_address_countries='；'.join(unique(p['address_country'] for p in applicant_records)),inventor_names='；'.join(inventor_names),inventor_source=inventor_source,inventor_count=len(inventor_names) if inventor_names else None)
    if not inventor_names:notes.append('explicit_inventor_role_missing')
    relations=[];related=b.find('us-related-documents')
    if related is not None:
        for node in related:
            if node.tag in ['related-publication','us-provisional-application']:continue
            for relation in node.findall('relation'):
                for parent in relation.findall('parent-doc'):
                    d=document(parent.find('document-id'))
                    relations.append(dict(relationship_type=node.tag,parent_application=d,parent_grants=[document(e) for e in parent.findall('parent-grant-document/document-id')],child_applications=[document(e) for e in relation.findall('child-doc/document-id')],locator='us-bibliographic-data-grant/us-related-documents/'+node.tag+'/relation',xml=ET.tostring(relation,encoding='unicode')))
    relation_labels={'continuation':'继续申请','continuation-in-part':'部分继续申请','division':'分案','reissue':'再授权关联申请','substitution':'替代申请'}
    m.update(parent_application_numbers='；'.join(unique(e['parent_application']['number'] for e in relations)),parent_application_relationships='；'.join(unique(relation_labels.get(e['relationship_type'],e['relationship_type']) for e in relations)))
    references=b.find('us-references-cited')
    if references is None:references=b.find('references-cited')
    citations=[]
    if references is not None:
        for node in references:
            pc=node.find('patcit');npl=node.find('nplcit')
            citations.append(dict(citation_type='patent' if pc is not None else 'non_patent' if npl is not None else 'other',patent=document(pc.find('document-id')) if pc is not None else None,text=plain(npl) if npl is not None else plain(node),attributes=dict(node.attrib),xml=ET.tostring(node,encoding='unicode')))
    m.update(patent_citation_count=sum(c['citation_type']=='patent' for c in citations) if references is not None else None,npl_citation_count=sum(c['citation_type']=='non_patent' for c in citations) if references is not None else None)
    m.update(bibliography_json_path=str(PurePosixPath('bibliography')/(patent+'.json')),document_metadata_review_reason=';'.join(notes))
    provenance=dict(archive=archive,xml_file=xml_file,application=application,document_reference=pub)
    bibliography=dict(schema_version='1.0',Patent=patent,provenance=provenance,abstract=dict(text=m['abstract'],language=language,locator='abstract'),applicants=applicant_records,inventors=inventor_records,inventor_source=inventor_source,assignees=assignee_records,related_applications=relations,citations=citations,limitations=['Roles and relationships are recorded in this publication, not current ownership or complete family information.','Citations are bibliographic prior-art references, not this application\u2019s filing or priority records.'])
    return m,bibliography
