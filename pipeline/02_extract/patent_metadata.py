"""USPTO grant metadata with explicit provenance; no outside ownership/family inference.

Bibliographic values are confined to this document's bibliographic container.
Indications and document categories are rule-derived patent disclosures, not
approved indications, claim-scope opinions, or a complete family search.
"""
from datetime import datetime
import json
import re
from patent_documents import DOCUMENT_COLUMNS

METADATA_COLUMNS = [
    'title','company_name','company_source','applicant_names','assignee_names',
    'is_pct','pct_number','pct_relationship','pct_filing_date',
    'wo_publication_number','wo_publication_date','pct_evidence',
    'priority_date','priority_date_candidate','priority_date_status','priority_claims',
    'priority_evidence','application_publication_number','application_publication_date',
    'grant_date','grant_date_type','document_publication_date','date_evidence',
    'indication','indication_source','other_disclosed_indications','indication_evidence',
    'ipc_codes','cpc_codes','uspc_codes','patent_category','classification_basis',
    'technical_field','technical_field_source','technical_field_locator',
    'metadata_review_reason',
]
METADATA_COLUMNS += DOCUMENT_COLUMNS
REVIEW_METADATA_COLUMNS = [
    'title','company_name','company_source','is_pct','pct_number','pct_relationship',
    'priority_date','priority_date_candidate','priority_date_status',
    'application_publication_number','application_publication_date','grant_date','grant_date_type',
    'indication','indication_source','other_disclosed_indications','ipc_codes','cpc_codes','patent_category',
    'technical_field_excerpt','technical_field_source','metadata_review_reason',
]

def content(e):
    return re.sub(r'\s+', ' ', ''.join(e.itertext())).strip() if e is not None else ''

def unique(values):
    return list(dict.fromkeys(v for v in values if v))

def date_iso(raw):
    raw=raw.strip()
    try:return datetime.strptime(raw, '%Y%m%d').date().isoformat()
    except ValueError:return ''

def prose_date(raw):
    m=re.fullmatch(r'([A-Za-z]+)\.?\s+(\d{1,2}),?\s+(\d{4})',raw.strip())
    if not m:return ''
    months=['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
    try:return datetime(int(m[3]),months.index(m[1][:3].lower())+1,int(m[2])).date().isoformat()
    except (ValueError,IndexError):return ''

def document(e):
    if e is None:return dict(number='',date='',country='',kind='')
    country=e.findtext('country','').strip();number=e.findtext('doc-number','').strip()
    kind=e.findtext('kind','').strip()
    # PCT/WO identifiers already include an authority. Keep their original number.
    number=number if re.match(r'^(?:PCT/|WO)',number) else country+number if number else ''
    if kind and kind not in ('00',) and number and not number.startswith('PCT/'):number+=kind
    return dict(number=number,date=date_iso(e.findtext('date','')),country=country,kind=kind)

def sections(root):
    current='';out=[];description=root.find('description')
    if description is None:return out
    for i,e in enumerate(description):
        if e.tag=='heading':current=content(e)
        elif e.tag=='p':out.append(dict(section=current,text=content(e),locator='description/'+(e.get('id') or f'p-{i+1}')))
    return out

def code(e):
    parts=[e.findtext(k,'').strip() for k in ['section','class','subclass','main-group','subgroup']]
    return ''.join(parts[:3])+' '+parts[3]+'/'+parts[4] if all(parts) else ''

# Deliberately finite vocabulary. Unmatched therapeutic language stays unclassified.
DISEASES = [
    ('HIV感染',r'\bHIV(?:[- ]?[12])?\b|human immunodeficiency virus'),
    ('乙型肝炎',r'hepatitis B|\bHBV\b'),('丙型肝炎',r'hepatitis C|\bHCV\b'),
    ('细菌感染',r'\bbacterial infections?\b'),('革兰阴性菌感染',r'gram[- ]negative bacterial infections?'),
    ('尿路感染',r'\burinary tract infections?\b'),('膀胱炎',r'\bcystitis\b'),
    ('肾盂肾炎',r'\bpyelonephritis\b'),('尿道炎',r'\burethritis\b'),
    ('膀胱感染',r'\bbladder infections?\b'),('肾脏感染',r'\bkidney infections?\b'),
    ('病毒感染',r'\bviral infections?\b'),('真菌感染',r'\bfungal infections?\b'),
    ('疟疾',r'\bmalaria\b'),('结核病',r'\btuberculosis\b'),
    ('癌症/肿瘤',r'\bcancers?\b|\btumou?rs?\b|\bneoplasms?\b|\bcarcinomas?\b'),
    ('白血病',r'\bleuk[ae]emias?\b'),('淋巴瘤',r'\blymphomas?\b'),
    ('多发性骨髓瘤',r'\bmultiple myeloma\b'),('乳腺癌',r'\bbreast cancer\b'),
    ('肺癌',r'\blung cancer\b'),('前列腺癌',r'\bprostate cancer\b'),
    ('结直肠癌',r'\bcolorectal cancer\b|\bcolon cancer\b'),
    ('黑色素瘤',r'\bmelanoma\b'),('动脉粥样硬化',r'\batherosclerosis\b'),
    ('类风湿关节炎',r'\brheumatoid arthritis\b'),('银屑病',r'\bpsoriasis\b|\bpsoriatic arthritis\b'),
    ('自身免疫性疾病',r'\bautoimmune (?:diseases?|disorders?)\b'),
    ('炎症性疾病',r'\binflammatory (?:diseases?|disorders?|conditions?)\b|\binflammation\b'),
    ('炎症性肠病',r'\binflammatory bowel disease\b|\bCrohn.s disease\b|\bulcerative colitis\b'),
    ('哮喘',r'\basthma\b'),('慢性阻塞性肺病',r'\bCOPD\b|chronic obstructive pulmonary'),
    ('过敏性疾病',r'\ballergic (?:diseases?|disorders?|conditions?)\b|\ballergy\b'),
    ('阿尔茨海默病',r'\bAlzheimer.s'),('帕金森病',r'\bParkinson.s'),
    ('神经退行性疾病',r'\bneurodegenerative (?:diseases?|disorders?)\b'),
    ('多发性硬化',r'\bmultiple sclerosis\b'),('抑郁症',r'\bdepressi(?:on|ve disorder)\b'),
    ('焦虑症',r'\banxiety\b'),('精神分裂症',r'\bschizophrenia\b'),
    ('糖尿病',r'\bdiabetes\b'),('肥胖症',r'\bobesity\b'),
    ('疼痛',r'\bpain\b'),('纤维化',r'\bfibrosis\b|\bfibrotic (?:diseases?|disorders?)\b'),
    ('血栓形成',r'\bthrombosis\b'),('心血管疾病',r'\bcardiovascular (?:diseases?|disorders?)\b'),
    ('心脏病',r'\bheart disease\b'),
    ('非酒精性脂肪性肝病',r'\bNAFLD\b|non[- ]alcoholic fatty liver disease'),
    ('非酒精性脂肪性肝炎',r'\bNASH\b|non[- ]alcohol(?:ic)? steatohepatitis'),
]
THERAPEUTIC=re.compile(r'\btreat(?:ing|ment|ed|s)?\b|\bprevent(?:ion|ing|s)?\b|\btherap(?:y|eutic)\b',re.I)

def indications(root, paragraphs, field):
    # Avoid cited background statements. A patent's use is supported in claims,
    # explicit FIELD/SUMMARY, abstract, or an invention-use paragraph outside BACKGROUND.
    candidates=[]
    therapeutic_claims=set()
    for i,e in enumerate(root.findall('claims/claim')):
        s=content(e)
        refs={int(n) for n in re.findall(r'\bclaim\s+(\d+)\b',s,re.I)}
        if THERAPEUTIC.search(s) or refs & therapeutic_claims:
            num=re.match(r'(\d+)\.',s)
            therapeutic_claims.add(int(num[1]) if num else i+1)
            candidates.append((0,'权利要求', 'claims/'+(e.get('id') or f'claim-{i+1}'), s))
    for p in paragraphs:
        section=p['section']
        eligible=bool(re.search(r'FIELD|SUMMARY',section,re.I)) or (not re.search(r'BACKGROUND|PRIOR ART',section,re.I) and bool(re.search(r'(?:present (?:invention|disclosure)|provided herein|compounds (?:of|according to) (?:the|this|formula)|compounds? (?:described|provided|disclosed) herein)',p['text'],re.I)))
        if eligible:
            for sentence in re.split(r'(?<=[.;])\s+(?=[A-Z])',p['text']):
                if THERAPEUTIC.search(sentence):
                    rank=1 if re.search(r'FIELD',section,re.I) else 3 if re.search(r'SUMMARY',section,re.I) else 4
                    source={1:'技术领域章节治疗用途',3:'发明概要治疗用途',4:'说明书其他治疗用途'}[rank]
                    candidates.append((rank,source,p['locator'],sentence))
    abstract=content(root.find('abstract'))
    for sentence in re.split(r'(?<=[.;])\s+(?=[A-Z])',abstract):
        if THERAPEUTIC.search(sentence):candidates.append((2,'摘要', 'abstract',sentence))
    evidence=[];found=[]
    # Each label's best available source is retained, rather than only the first disease.
    for rank,source,loc,s in sorted(candidates,key=lambda c:c[0]):
        for label,pattern in DISEASES:
            m=re.search(pattern,s,re.I)
            if m and label not in found:
                found.append(label)
                evidence.append(dict(indication=label,matched_text=m[0],source=source,source_rank=rank,locator=loc,text=s))
    best=min((e['source_rank'] for e in evidence),default=None)
    primary=[e for e in evidence if e['source_rank']==best]
    for e in evidence:e['primary']=e['source_rank']==best
    return '；'.join(e['indication'] for e in primary),'；'.join(unique(e['source'] for e in primary)),'；'.join(e['indication'] for e in evidence if not e['primary']),evidence

def extract_metadata(root):
    b=root.find('us-bibliographic-data-grant')
    if b is None:raise ValueError('missing current grant bibliographic container')
    m={k:'' for k in METADATA_COLUMNS};notes=[]
    m['title']=content(b.find('invention-title'))
    assignees=unique(content(e) for e in b.findall('assignees/assignee/addressbook/orgname'))
    applicants=unique(content(e) for e in b.findall('us-parties/us-applicants/us-applicant/addressbook/orgname')+b.findall('parties/applicants/applicant/addressbook/orgname'))
    m.update(company_name='；'.join(assignees or applicants),company_source='授权文献受让人' if assignees else '授权文献申请人（机构）' if applicants else '',applicant_names='；'.join(applicants),assignee_names='；'.join(assignees))
    if not m['company_name']:notes.append('organization_name_missing')
    grant=document(b.find('publication-reference/document-id'))
    published=[document(e) for e in b.findall('us-related-documents/related-publication/document-id')]
    uspub=[p for p in published if p['country']=='US' and p['kind'].startswith('A')]
    uspub.sort(key=lambda p:p['date'] or '9999')
    is_reissue=grant['kind'].startswith('E')
    m.update(document_publication_date=grant['date'],grant_date=grant['date'] if grant['kind'].startswith('B') or is_reissue else '',grant_date_type='再授权文献日期（不是原始授权日）' if is_reissue else '本件授权日期' if grant['kind'].startswith('B') else '')
    if is_reissue:notes.append('reissue_date_not_original_grant_date')
    if uspub:m.update(application_publication_number=uspub[0]['number'],application_publication_date=uspub[0]['date'])
    else:notes.append('us_application_publication_not_recorded')
    m['date_evidence']=dict(grant_reference=grant,related_publications=published)
    paragraphs=sections(root)
    crossrefs=[p for p in paragraphs if re.search(r'CROSS.?REFERENCE|RELATED APPLICATION',p['section'],re.I)]
    current_pct=document(b.find('pct-or-regional-filing-data/document-id'))
    pcte=[]
    if current_pct['number'].startswith('PCT/'):
        pcte.append(dict(**current_pct,relationship='本申请PCT国家阶段',source='us-bibliographic-data-grant/pct-or-regional-filing-data/document-id'))
    related=b.find('us-related-documents')
    if related is not None:
        for relation in related:
            for e in relation.findall('relation/parent-doc/document-id'):
                p=document(e)
                if p['number'].startswith('PCT/'):
                    label={'continuation':'PCT母案的继续申请','division':'PCT母案的分案','continuation-in-part':'PCT母案的部分继续申请'}.get(relation.tag,'PCT关联母案')
                    pcte.append(dict(**p,relationship=label,source='us-bibliographic-data-grant/us-related-documents/'+relation.tag+'/relation/parent-doc/document-id'))
    if not pcte:
        for p in crossrefs:
            numbers=unique(re.findall(r'\bPCT/[A-Z]{2}\s?\d{4}/\d{4,6}\b',p['text']))
            for number in numbers:pcte.append(dict(number=number.replace(' ',''),date='',relationship='说明书记载PCT关联（需核对具体申请关系）',source=p['locator'],text=p['text']))
    m.update(is_pct='是' if pcte else '未识别到PCT（未查完整专利族）',pct_number='；'.join(unique(p['number'] for p in pcte)),pct_relationship='；'.join(unique(p['relationship'] for p in pcte)),pct_filing_date=current_pct['date'] if current_pct['number'].startswith('PCT/') else pcte[0]['date'] if pcte else '',pct_evidence=pcte)
    wo=document(b.find('pct-or-regional-publishing-data/document-id'))
    m.update(wo_publication_number=wo['number'],wo_publication_date=wo['date'])
    if not pcte:notes.append('pct_family_not_verified')
    priorities=[]
    for e in b.findall('priority-claims/priority-claim'):
        p=document(e);priorities.append(dict(**p,source='us-bibliographic-data-grant/priority-claims/priority-claim'))
    for e in b.findall('us-related-documents/us-provisional-application/document-id'):
        p=document(e);priorities.append(dict(**p,source='us-bibliographic-data-grant/us-related-documents/us-provisional-application/document-id'))
    prose=[]
    pattern=r'(\d{2}/\d{3}[,.]?\d{3})[^\n]{0,65}?\bfiled(?:\s+on)?\s+([A-Za-z]+\.?\s+\d{1,2},?\s+\d{4})'
    for p in crossrefs:
        for match in re.finditer(pattern,p['text'],re.I):
            number='US'+re.sub(r'\D','',match[1]);date=prose_date(match[2])
            if date:prose.append(dict(number=number,date=date,source=p['locator'],text=p['text']))
    valid=[p['date'] for p in priorities if p['date']]
    candidate=min(valid) if valid else ''
    conflicting=[dict(bibliographic=p,description=q) for p in priorities for q in prose if p['number']==q['number'] and p['date'] and p['date']!=q['date']]
    m.update(priority_date_candidate=candidate,priority_date=candidate if not conflicting else '',priority_date_status='书目与说明书冲突，待核对' if conflicting else 'XML优先权/临时申请最早记载日期（未核对完整专利族）' if candidate else '本文件未记载可提取优先权日期',priority_claims=priorities,priority_evidence=dict(structured=priorities,cross_reference_dates=prose,conflicts=conflicting))
    if conflicting:notes.append('priority_date_conflict')
    if not candidate:notes.append('priority_date_missing')
    ipc=unique(code(e) for e in b.findall('classifications-ipcr/classification-ipcr'))
    cpc=unique(code(e) for e in b.findall('classifications-cpc/main-cpc/classification-cpc')+b.findall('classifications-cpc/further-cpc/classification-cpc'))
    uspc=unique(content(e) for e in b.findall('classification-national/main-classification')+b.findall('classification-national/further-classification'))
    m.update(ipc_codes='；'.join(ipc),cpc_codes='；'.join(cpc),uspc_codes='；'.join(uspc))
    # Claim-object categories are separate from official assigned classifications.
    categories=[];category_basis=[]
    patterns=[('化合物',r'^\d*\.?\s*(?:A|An|The)?\s*(?:[A-Za-z-]+\s+){0,2}compound\b'),('药物组合物',r'^\d*\.?\s*(?:A|The)?\s*(?:pharmaceutical|therapeutic) composition\b'),('治疗用途/方法',r'^\d*\.?\s*(?:A|The)?\s*method (?:of|for) (?:the )?(?:treating|treatment|preventing|prevention)\b'),('制备方法',r'^\d*\.?\s*(?:A|The)?\s*(?:method|process) (?:of|for) (?:preparing|producing|synthesizing)\b'),('联合用药',r'^\d*\.?\s*(?:A|The)?\s*combination\b'),('靶点调控方法',r'^\d*\.?\s*(?:A|The)?\s*method (?:of|for) (?:inhibiting|modulating)\b'),('抗体/免疫偶联物',r'^\d*\.?\s*(?:An?|The)?\s*(?:immunoconjugate|antibody)\b')]
    for i,e in enumerate(root.findall('claims/claim')):
        s=content(e)
        for label,pattern in patterns:
            if re.search(pattern,s,re.I) and label not in categories:
                categories.append(label);category_basis.append(dict(category=label,locator='claims/'+(e.get('id') or f'claim-{i+1}'),text=s[:600]))
    m.update(patent_category='；'.join(categories) or '未识别（需核对权利要求）',classification_basis=dict(official_source='本件XML已分配IPC/CPC/USPC，不含检索分类范围',derived_categories=category_basis))
    fields=[p for p in paragraphs if re.fullmatch(r'(?:TECHNICAL )?FIELD(?: OF (?:THE )?(?:INVENTION|DISCLOSURE))?',p['section'].strip(' .:'),re.I)]
    if fields:
        m.update(technical_field=' '.join(p['text'] for p in fields),technical_field_source='说明书技术领域章节',technical_field_locator='；'.join(p['locator'] for p in fields))
    else:
        abstract=content(root.find('abstract'))
        m.update(technical_field=abstract or m['title'],technical_field_source='摘要代替（无独立技术领域章节）' if abstract else '标题代替（无技术领域章节和摘要）',technical_field_locator='abstract' if abstract else 'us-bibliographic-data-grant/invention-title')
        notes.append('technical_field_section_missing')
    indication,source,other,evidence=indications(root,paragraphs,fields)
    m.update(indication=indication,indication_source=source or '未识别明确治疗用途',other_disclosed_indications=other,indication_evidence=evidence)
    if not indication:notes.append('indication_not_identified')
    m['metadata_review_reason']=';'.join(notes)
    return m
