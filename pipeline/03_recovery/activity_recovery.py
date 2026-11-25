"""Recover source-linked activities. Never derive experimental IDs from MOL names.

The old accepted SAR is a read-only baseline. Extended ranges/qualitative values
are separate from numeric SAR. All source references stay within one archive.
"""
from collections import Counter, defaultdict
from decimal import Decimal
from html.parser import HTMLParser
from pathlib import Path, PurePosixPath
import hashlib
import json
import re
import sys
import xml.etree.ElementTree as ET
import zipfile

sys.path.insert(0, str(Path(__file__).parent / 'vendor'))
import uspto_sar as cals_adapter
from uspto_sar import expanded_rows, text, measurement, UNIT, SAR_COLUMNS, name_element_conflicts
from xml_source import parse_xml, select_patent_xml, _entities, _reserved
from rdkit import Chem, RDLogger
import pyarrow.parquet as pq

RDLogger.DisableLog('rdApp.*')
CAT = Path('/home/open/patent-work/outputs/uspto-all-structures-20261004/results/checkpoint_cache')
OLD = Path('/home/open/patent-work/outputs/uspto-sar-full-20261004/results/checkpoint_cache')
RAW = Path('/home/open/patent-work/raw/patent-data')
CAT_CODE_SHA256 = '67b554fd4599c2cde756f0689263214bb22fdb97ef46e7631c49970d4f7e0214'
ACTIVITY = re.compile(r'\b(p?IC\s*(?:50|90)|p?EC\s*(?:50|90)|GI\s*50|AC\s*50|CC\s*50|LC\s*50|LD\s*50|ED\s*50|K\s*[id]|MIC|MBC)\b',re.I)
# The vendored row expansion promotes a header only after an explicit endpoint.
# Broaden the recognized endpoints in memory; frozen source files stay unchanged.
cals_adapter.ACTIVITY=ACTIVITY
MASS_UNIT = re.compile(r'(?<![A-Za-z])(?:[uµμn]g|mg|g)\s*/\s*(?:mL|ml|L|l|kg)(?![A-Za-z])')
FACTORS = {'pM':Decimal('.001'), 'nM':Decimal(1), 'uM':Decimal(1000), 'mM':Decimal(1000000), 'M':Decimal(1000000000)}
ALIAS_HEADER = re.compile(r'example|compound|cmpd|comp\.?\b|\bno\b', re.I)
CONTROL = re.compile(r'positive control|negative control|reference compound|staurosporine', re.I)
ASSAY = re.compile(r'\b(?:assay|binding studies|enzyme inhibition|kinase inhibition|biological evaluation)\b', re.I)
HEADING = re.compile(r'^(?:Example|Compound)\s*#?\s*(\d+(?:[.-]\d+)?[a-z]?)(?=$|[\s:—–.])', re.I)
VOID = {'colspec','chemmol','chemcdx','emi','img','br','hr','meta'}
NUMERIC = {'Value','value_nM','value_lower','value_upper','lower_nM','upper_nM'}
INTEGER = {'source_header_offset'}
BOOL = {'drawing_vs_mol_checked'}
EXTRA = ['archive','xml_file','xml_sha256','source_cell_key','activity_locator','activity_header',
         'raw_value','original_range','compound_id','compound_alias','original_alias','alias_evidence',
         'structure_locator','structure_match','source_structure_id','structure_file','raw_mol_sha256',
         'raw_mol_coordinate_dimension','source_header_offset','target_evidence','target_scope',
         'unit_evidence','legend_evidence','value_lower','value_upper','lower_nM','upper_nM',
         'lower_inclusive','upper_inclusive','measurement_kind','extraction_status','review_reason',
         'source_format','source_mapping_checks','source_evidence_json','drawing_vs_mol_checked',
         'Sequence','UniProt ID','compound_IUPAC_name','compound_name_status','phase','prior_status']
RECORD_COLUMNS = list(dict.fromkeys(SAR_COLUMNS + EXTRA))
STRUCT_COLS = ['source_structure_id','structure_file','source_header_offset','source_format',
               'structure_status','SMILES','SMILES_original','InChIKey','source_sha256',
               'raw_mol_coordinate_dimension','Sequence','compound_IUPAC_name','compound_name_status']


def alias(value):
    s = re.sub(r'^(?:Example|Ex\.?|Compound|Cmpd\.?|Comp\.?|No\.?|#)\s*#?\s*','',value.strip(),flags=re.I)
    return s.lower().replace('.','-') if re.fullmatch(r'\d+(?:[.-]\d+)?[a-z]?',s,re.I) else None


def units(value):
    for name, u in [('micromolar','uM'),('nanomolar','nM'),('millimolar','mM'),('picomolar','pM')]:
        value = re.sub(name,u,value,flags=re.I)
    result={m[1].replace('μ','u').replace('µ','u') for m in UNIT.finditer(value)}
    result.update(re.sub(r'\s','',m[0]).replace('μ','u').replace('µ','u').replace('/ml','/mL').replace('/l','/L') for m in MASS_UNIT.finditer(value))
    return result


def source_measurement(value):
    # Preserve mass concentration/dose units without converting them to molarity.
    clean=MASS_UNIT.sub('',value)
    return measurement(clean)


class SGMLTree(HTMLParser):
    """Stack parser with explicit SGML EMPTY tags; mismatched closes are audited."""
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.root=ET.Element('source');self.stack=[self.root];self.mismatches=[]

    def handle_starttag(self, tag, attrs):
        e=ET.SubElement(self.stack[-1],tag,dict((k,v or '') for k,v in attrs))
        if tag not in VOID:self.stack.append(e)

    def handle_startendtag(self,tag,attrs):
        ET.SubElement(self.stack[-1],tag,dict((k,v or '') for k,v in attrs))

    def handle_endtag(self,tag):
        if tag in VOID:return
        found=next((i for i in range(len(self.stack)-1,0,-1) if self.stack[i].tag==tag),None)
        if found is None:self.mismatches.append('orphan_close:'+tag);return
        if found!=len(self.stack)-1:self.mismatches.append('implicit_close:'+tag)
        del self.stack[found:]

    def handle_data(self,s):
        e=self.stack[-1]
        if len(e):e[-1].tail=(e[-1].tail or '')+s
        else:e.text=(e.text or '')+s


def parse_source(raw, sgml=False):
    unknown=[]
    # Preserve unknown named entities literally, never guess symbols/values.
    def entity(m):
        name=m[1]
        if name in _reserved:return m[0]
        if name in _entities:return ''.join('&#%d;'%ord(c) for c in _entities[name])
        unknown.append(name);return '[UNRESOLVED_ENTITY:'+name+']'
    encoding=re.search(br'<\?xml[^>]*encoding=[\"\']([^\"\']+)',raw[:512],re.I)
    s=raw.decode(encoding[1].decode('ascii') if encoding else 'utf-8-sig',errors='strict')
    parts=re.split(r'(<!\[CDATA\[.*?\]\]>|<!--.*?-->)',s,flags=re.S)
    s=''.join(part if i%2 else re.sub(r'&([A-Za-z][A-Za-z0-9]*);',entity,part) for i,part in enumerate(parts))
    s=re.sub(r'(<\?xml[^>]*encoding=)[\"\'][^\"\']+[\"\']',r'\1"utf-8"',s,flags=re.I)
    mismatches=[]
    if sgml:
        # The recovery reads description only. Bibliography SGML has optional
        # end tags and must not cause every otherwise sound assay to be rejected.
        description_match=re.search(r'<SDODE\b[^>]*>.*?</SDODE\s*>',s,re.S|re.I)
        if not description_match:raise ValueError('missing_complete_sgml_description')
        s=description_match[0]
        parser=SGMLTree();parser.feed(s);parser.close();root=parser.root;mismatches=parser.mismatches
        if len(parser.stack)>1:mismatches.append('unclosed_document')
    else:root=parse_xml(s.encode('utf-8'))
    for e in root.iter():
        e.tag=e.tag.rsplit('}',1)[-1].lower()
        e.attrib={k.lower():v for k,v in e.attrib.items()}
        e.tag={'sb':'sub','sp':'sup','h':'heading','para':'p','tr':'row','td':'entry','th':'entry'}.get(e.tag,e.tag)
        if e.tag in ('chemmol','chemcdx','emi'):e.tag='img'
    return root,dict(unresolved_entities=sorted(set(unknown)),sgml_closure_warnings=mismatches)


def blocks_of(description):
    result=[]
    def visit(e):
        if e.tag in ('table','p','heading'):
            result.append(e);return
        for c in e:visit(c)
    visit(description)
    return result


def source_rows(table):
    # HTML tables are normalized into CALS-like row nodes by the adapter.
    return expanded_rows(table)


def table_index(description, blocks):
    parents={c:p for p in description.iter() for c in p}
    result=[]
    for bi,b in enumerate(blocks):
        for t in b.iter('table'):
            p=parents.get(t)
            locator=(p.get('id','') if p is not None else '') or t.get('id','') or b.get('id','') or 'table-%d'%(len(result)+1)
            headers, rows=source_rows(t)
            captions=[]
            groups=list(t.findall('tgroup'))
            width=max((int(g.get('cols','0') or 0) for g in groups),default=0)
            for g in groups:
                if int(g.get('cols','0') or 0)<width:captions.append(text(g))
            result.append(dict(element=t,headers=headers,rows=rows,locator=locator,block=bi,captions=captions))
    if len({t['locator'] for t in result})!=len(result):raise ValueError('duplicate_table_locator')
    return result


def resolve_references(files, structures):
    bystem=defaultdict(list)
    for r in structures:bystem[PurePosixPath(r['structure_file']).stem.upper()].append(r)
    selected={};errors=[]
    for file in set(files):
        stem=PurePosixPath(file.replace('\\','/')).stem.upper()
        choices=bystem.get(stem,[])
        if len(choices)>1:errors.append('ambiguous_duplicate_structure_member:'+stem)
        elif choices:selected[choices[0]['source_structure_id']]=choices[0]
    return list(selected.values()),errors


def refs(e):return [im.get('file','') for im in e.iter('img') if im.get('file')]


def build_links(tables,blocks,structures):
    links=defaultdict(list)
    def put(a, files, method, locator, evidence):
        graph,errs=resolve_references(files,structures)
        if a and len(graph)==1 and not errs:
            links[a].append(dict(graph=graph[0],method=method,locator=locator,evidence=evidence))
    for t in tables:
        acs=[i for i,h in enumerate(t['headers']) if ALIAS_HEADER.search(h) and not ACTIVITY.search(h) and not re.search(r'name|mass',h,re.I)]
        for row in t['rows']:
            identifiers={alias(row[c]['text']) for c in acs}-{None}
            if len(identifiers)==1:
                a=next(iter(identifiers));put(a,[f for c in row for f in c['files']], 'explicit_same_table_row',t['locator']+'/'+row[0]['row_locator'],
                    dict(alias=a,alias_text=[row[c]['text'] for c in acs],row=[c['text'] for c in row]))
    for i,b in enumerate(blocks):
        s=text(b)
        heading=(b.tag=='heading' or (b.tag=='p' and len(s)<180 and b.find('.//b') is not None))
        m=HEADING.match(s) if heading else None
        if m and not re.search(r'\b(?:and|to|through)\b|[,;&]',s[m.end():m.end()+12],re.I):
            # Consider the complete section, not the first adjacent drawing.
            end=next((j for j in range(i+1,len(blocks)) if blocks[j].tag=='heading' or HEADING.match(text(blocks[j]))),len(blocks))
            section=blocks[i:end]
            if len(section)>150:continue
            if any(re.match(r'^(?:Step\b|Intermediate\b|Preparation\b)',text(x),re.I) for x in section[1:]):continue
            # Require a reported product name and a uniquely drawn product at
            # the beginning of the complete example. A single reagent drawing
            # later in a synthesis paragraph is insufficient.
            named_product=s[m.end():].strip(' :.\u2014\u2013-')
            if len(named_product)<12 or not re.search(r'yl|acid|amide|amine|one|ene|ol|nitrile|[()]',named_product,re.I):continue
            first_graph=next((j for j,x in enumerate(section) if refs(x)),None)
            if first_graph is None or first_graph>3:continue
            if any(len(text(x))>180 or re.search(r'\b(?:reaction|reagent|intermediate|starting material|scheme)\b',text(x),re.I) for x in section[1:first_graph+1]):continue
            put(alias(m[1]),[f for x in section for f in refs(x)],'unique_complete_example_section',
                b.get('id','') or 'block-%d'%(i+1),dict(heading=s,reported_product_name=named_product,start_block=i,end_block=end,
                referenced_files=sorted(set(f for x in section for f in refs(x)))))
        # An explicit labeled compound AND its graph in one paragraph.
        if b.tag=='p' and refs(b):
            ids={alias(v) for v in re.findall(r'\b(?:compound|example)\s*#?\s*(\d+(?:[.-]\d+)?[a-z]?)\b',s,re.I)}-{None}
            if len(ids)==1 and not re.search(r'\b(?:intermediate|starting material|reagent|reactant|mixture)\b',s,re.I):
                put(next(iter(ids)),refs(b),'explicit_same_paragraph',b.get('id','') or 'block-%d'%(i+1),dict(paragraph=s))
    resolved={};conflicts=set()
    for a,entries in links.items():
        identities={(r['graph']['SMILES'],r['graph']['InChIKey']) for r in entries}
        if len(identities)==1 and all(r['graph']['InChIKey'] for r in entries):
            resolved[a]=sorted(entries,key=lambda r:r['method']!='explicit_same_table_row')[0]
        else:conflicts.add(a)
    return resolved,conflicts


def contexts(t,blocks):
    before=[];after=[]
    for j in range(t['block']-1,max(-1,t['block']-35),-1):
        b=blocks[j]
        if b.tag=='table' or b.find('.//table') is not None:break
        before.append(dict(locator=b.get('id','') or 'block-%d'%(j+1),text=text(b)))
    for j in range(t['block']+1,min(len(blocks),t['block']+8)):
        b=blocks[j]
        if b.tag=='heading' or b.tag=='table' or b.find('.//table') is not None:break
        after.append(dict(locator=b.get('id','') or 'block-%d'%(j+1),text=text(b)))
    scope=list(reversed(before))+after
    for caption in t['captions']:scope.append(dict(locator=t['locator']+'/caption',text=caption))
    headings=[]
    for j in range(t['block']-1,-1,-1):
        b=blocks[j];s=text(b)
        if (b.tag=='heading' or len(s)<180) and ASSAY.search(s):
            if not re.search(r'general procedures|assays? (?:were|was|are)|synthe',s,re.I):headings=[dict(locator=b.get('id','') or 'block-%d'%(j+1),text=s)]
            break
        if b.tag=='heading' and re.search(r'\b(?:Example\s*\d|synthesis|preparation|Intermediate)\b',s,re.I):break
    return scope,headings


def target_for(header,row,headers,scope,headings):
    # Explicit row target or non-numeric endpoint column identity is strongest.
    cols=[i for i,h in enumerate(headers) if re.search(r'\btarget\b|cell line|receptor|enzyme|kinase',h,re.I) and not ACTIVITY.search(h)]
    if len(cols)==1 and row[cols[0]]['text'] and not alias(row[cols[0]]['text']):
        return row[cols[0]]['text'],'explicit_target_row',dict(column=cols[0],text=row[cols[0]]['text'])
    prefix=MASS_UNIT.sub('',UNIT.sub('',ACTIVITY.sub('',header)))
    prefix=re.sub(r'\b(?:micromolar|nanomolar|millimolar|mean|average|value|values|SD|SEM)\b','',prefix,flags=re.I)
    prefix=re.sub(r'[()\[\]:%]',' ',prefix);prefix=re.sub(r'\s+',' ',prefix).strip(' /-,')
    if valid_target(prefix) and not re.search(r'compound|example|\bex\.|activity|inhibition|concentration|binding affinity|potency',prefix,re.I):
        return prefix,'explicit_activity_column',dict(text=header)
    candidates=[s for s in scope if len(s['text'])<220 and ASSAY.search(s['text']) and not re.search(r'values?|indicates?|determined|following table|protocol|procedure',s['text'],re.I)]
    candidates=headings or candidates
    if candidates:
        distinct={re.sub(r'\W','',s['text']).lower() for s in candidates}
        if len(distinct)==1:
            s=candidates[-1]
            label=re.sub(r'^Example\s+[A-Z0-9]+\s*[:.—-]?\s*','',s['text'],flags=re.I)
            if valid_target(label):return label,'explicit_scoped_assay_heading',s
    return '','',{}


def valid_target(label):
    remaining=re.sub(r'\b(?:biological|biochemical|cellular|activity|activities|binding|enzyme|receptor|kinase|inhibition|assay|assays|study|studies|test|tests|evaluation|human|mouse|rat|proliferation|prol|results|determination|of|the|and|in|vitro|in vivo)\b','',label,flags=re.I)
    return bool(re.search(r'[A-Za-z]',remaining))


def legend_for(raw,endpoint,scope):
    symbol=raw.strip()
    if not re.fullmatch(r'[A-Za-z]|\+{1,4}|-{1,4}',symbol):return None
    matches=[]
    for anchor in scope:
        s=anchor['text'].replace('≦','≤').replace('≧','≥')
        pattern=r'(?<![\w+])(?:[“\"\‘\’\”\”])?'+re.escape(symbol)+r'(?:[”\"\’])?\s*(?:indicates?|means?|denotes?|represents?|[=:])\s*([^;\n]{1,180})'
        for m in re.finditer(pattern,s,re.I):
            definition=re.split(r'\s+[A-Z]\s+(?:indicates|means|denotes)',m[1])[0]
            if not ACTIVITY.search(definition) or re.sub(r'\s','',ACTIVITY.search(definition)[1]).lower()!=endpoint.lower():continue
            u=units(definition)
            if len(u)!=1:continue
            number=r'(?:\d+(?:[.,]\d+)?|\.\d+)(?:[eE][+-]?\d+)?'
            interval=re.search(r'('+number+r')\s*([<≤])\s*=?\s*'+re.escape(endpoint)+r'\s*([<≤])\s*=?\s*('+number+r')',definition,re.I)
            point=re.search(re.escape(endpoint)+r'\s*(?:of|values?\s*(?:of|are)?)?\s*([<>≤≥]=?)\s*('+number+r')',definition,re.I)
            if interval:
                lower=float(Decimal(interval[1].replace(',','')));upper=float(Decimal(interval[4].replace(',','')))
                if lower>upper:continue
                val=dict(Relation='range',Value=None,value_lower=lower,value_upper=upper,lower_inclusive=str(interval[2]=='≤').lower(),upper_inclusive=str(interval[3]=='≤').lower())
            elif point:val=dict(Relation=point[1].replace('≤','<=').replace('≥','>='),Value=float(Decimal(point[2].replace(',',''))),value_lower=None,value_upper=None)
            else:continue
            val.update(Unit=next(iter(u)),legend_evidence=dict(locator=anchor['locator'],text=m[0],symbol=symbol,definition=definition))
            matches.append(val)
    signatures={json.dumps({k:v for k,v in m.items() if k!='legend_evidence'},sort_keys=True) for m in matches}
    return matches[0] if len(signatures)==1 else None


def decode(raw,endpoint,header,scope):
    reasons=[];u='';ue='';legend=legend_for(raw,endpoint,scope)
    cellunits=units(raw);headerunits=units(header)
    if len(cellunits)>1 or len(headerunits)>1 or (cellunits and headerunits and cellunits!=headerunits):reasons.append('conflicting_units')
    if len(cellunits)==1:u=next(iter(cellunits));ue='activity_cell_unit'
    elif len(headerunits)==1:u=next(iter(headerunits));ue='activity_column_unit'
    relation,value,lower,upper=source_measurement(raw)
    # Only a recognized concentration suffix is compatible with a numeric cell.
    # Decimal commas and percentages cannot silently become concentration data.
    if value is not None or lower is not None:
        if re.search(r'\d,\d{1,2}(?!\d)',raw):reasons.append('ambiguous_decimal_comma')
        tail=re.search(r'([A-Za-zµμ%]+)\s*$',MASS_UNIT.sub('',raw))
        if tail and tail[1].replace('µ','u').replace('μ','u') not in FACTORS:reasons.append('unrecognized_numeric_suffix')
    result=dict(Relation=relation,Value=value,value_lower=lower,value_upper=upper,Unit=u,unit_evidence=ue,legend_evidence={},measurement_kind='numeric',lower_inclusive='',upper_inclusive='')
    if legend:
        if u and u!=legend['Unit']:reasons.append('legend_unit_conflict')
        result.update(legend);result['unit_evidence']='explicit_category_legend';result['measurement_kind']='legend_threshold' if legend['Value'] is not None else 'legend_interval'
    if not result['Unit']:
        statements=[]
        for anchor in scope:
            pattern=re.escape(endpoint)+r'(?:\s+values?)?\s*(?:\(([^)]{1,15})\)|(?:are\s+)?(?:expressed|reported|given|measured)\s+in\s+([uµμ]M|nM|pM|mM))'
            for m in re.finditer(pattern,anchor['text'],re.I):
                found=units(m[0])
                if len(found)==1:statements.append((next(iter(found)),anchor,m[0]))
        if len({x[0] for x in statements})==1:
            result.update(Unit=statements[0][0],unit_evidence='explicit_endpoint_unit_statement')
            result['unit_statement']=dict(locator=statements[0][1]['locator'],text=statements[0][2])
    if endpoint.lower().startswith('p'):
        if result['Unit']:reasons.append('log_endpoint_has_concentration_unit')
        result['Unit']='dimensionless';result['measurement_kind']='log_endpoint';result['unit_evidence']='explicit_p_endpoint'
    elif not result['Unit']:reasons.append('activity_unit_missing')
    if result['Relation']=='text':reasons.append('qualitative_value_without_unique_definition')
    if result['Relation'] in ('~','≈'):reasons.append('approximate_value_requires_review')
    for k in ('Value','value_lower','value_upper'):
        if result[k] is not None and result[k]<0:reasons.append('negative_concentration')
    factor=FACTORS.get(result['Unit'])
    result['value_nM']=float(Decimal(str(result['Value']))*factor) if factor and result['Value'] is not None else None
    result['lower_nM']=float(Decimal(str(result['value_lower']))*factor) if factor and result['value_lower'] is not None else None
    result['upper_nM']=float(Decimal(str(result['value_upper']))*factor) if factor and result['value_upper'] is not None else None
    return result,reasons


def verify_graph(z, graph, cache):
    sid=graph['source_structure_id']
    if sid in cache:return cache[sid]
    errors=[]
    choices=[info for info in z.infolist() if info.header_offset==graph['source_header_offset'] and info.filename==graph['structure_file']]
    if len(choices)!=1:errors.append('mol_source_offset_not_unique')
    else:
        raw=z.read(choices[0]);sha=hashlib.sha256(raw).hexdigest()
        if sha!=graph['source_sha256']:errors.append('mol_sha256_mismatch')
        mol=Chem.MolFromMolBlock(raw.decode('utf-8','replace'),sanitize=True,removeHs=True,strictParsing=True)
        if mol is None or any(a.GetAtomicNum()==0 or a.HasQuery() for a in mol.GetAtoms()) or any(b.HasQuery() for b in mol.GetBonds()):errors.append('mol_not_concrete')
        elif Chem.MolToSmiles(mol,canonical=True,isomericSmiles=True)!=graph['SMILES'] or Chem.MolToInchiKey(mol)!=graph['InChIKey']:errors.append('mol_smiles_or_inchikey_mismatch')
    cache[sid]=errors
    return errors


def proof_index(tables,blocks):
    rows={};anchors={b.get('id','') or 'block-%d'%(i+1):text(b) for i,b in enumerate(blocks)}
    for t in tables:
        for row in t['rows']:rows[t['locator']+'/'+row[0]['row_locator']]=(t['headers'],row)
        for caption in t['captions']:anchors[t['locator']+'/caption']=caption
    return rows,anchors


def verify_record_source(r,proof,tables,blocks,structures):
    """Read source cells and evidence again; reject record tampering or mixed scopes.

    This check is separate from the association builder. Raw MOL graph verification
    runs afterwards using the archive member's offset rather than its basename.
    """
    errors=[];row_index,anchors=proof
    loc,number=r['activity_locator'].rsplit('/column-',1);col=int(number)-1
    headers,row=row_index[loc];cell=row[col]['text'];header=headers[col]
    if cell!=r['raw_value'] or header!=r['activity_header']:errors.append('raw_cell_or_header_mismatch')
    ep=ACTIVITY.search(header)
    if not ep or re.sub(r'\s','',ep[1]).lower()!=r['Activity Type'].lower():errors.append('endpoint_source_mismatch')
    a=r['compound_alias'];evidence=json.loads(r['source_evidence_json'])
    if a:
        alias_column=evidence['alias']['column']
        if alias_column is not None:
            if alias(row[alias_column]['text'])!=a or row[alias_column]['text']!=r['original_alias']:errors.append('activity_alias_source_mismatch')
        else:
            colid=re.search(r'\b(?:Example|Compound|Ex\.)\s*#?\s*(\d+(?:[.-]\d+)?[a-z]?)\b',header,re.I)
            # A single compound column to the right of the activity column is
            # accepted only when it is the unique explicit identifier column.
            aliascols=[i for i,h in enumerate(headers) if ALIAS_HEADER.search(h) and not ACTIVITY.search(h) and not re.search(r'name|mass',h,re.I)]
            if not (colid and alias(colid[1])==a) and not (len(aliascols)==1 and alias(row[aliascols[0]]['text'])==a):errors.append('activity_alias_not_in_source')
    normal=lambda s:re.sub(r'[^A-Za-z0-9α-ωΑ-Ω]','',s).lower()
    ta=evidence['target'];te=r['target_evidence']
    if te=='explicit_activity_column':
        if normal(r['Target']) not in normal(ACTIVITY.sub('',header)):errors.append('target_not_in_activity_header')
    elif te=='explicit_target_row':
        if row[ta['column']]['text']!=r['Target']:errors.append('target_row_mismatch')
    elif te=='explicit_scoped_assay_heading':
        if anchors.get(ta['locator'])!=ta['text'] or normal(r['Target']) not in normal(ta['text']):errors.append('target_anchor_mismatch')
        table=next(t for t in tables if r['activity_locator'].startswith(t['locator']+'/'))
        scope,headings=contexts(table,blocks)
        if ta not in (headings or scope):errors.append('target_outside_assay_scope')
    else:errors.append('target_source_missing')
    # Independently inspect the structure row/paragraph or complete section.
    structure_loc=r['structure_locator'];files=[]
    if r['structure_match'] in ('explicit_activity_row_structure','explicit_same_table_row'):
        sh,srow=row_index[structure_loc];files=[f for c in srow for f in c['files']]
        if r['structure_match']=='explicit_same_table_row' and not any(alias(c['text'])==a for c in srow):errors.append('structure_row_alias_mismatch')
    elif r['structure_match']=='explicit_same_paragraph':
        index=next((i for i,b in enumerate(blocks) if (b.get('id','') or 'block-%d'%(i+1))==structure_loc),None)
        if index is None:errors.append('structure_paragraph_not_found')
        else:
            files=refs(blocks[index]);s=text(blocks[index]);ids={alias(v) for v in re.findall(r'\b(?:compound|example)\s*#?\s*(\d+(?:[.-]\d+)?[a-z]?)\b',s,re.I)}-{None}
            if ids!={a}:errors.append('paragraph_compound_not_unique')
    elif r['structure_match']=='unique_complete_example_section':
        index=next((i for i,b in enumerate(blocks) if (b.get('id','') or 'block-%d'%(i+1))==structure_loc),None)
        if index is None:errors.append('example_heading_not_found')
        else:
            heading=HEADING.match(text(blocks[index]))
            if not heading or alias(heading[1])!=a:errors.append('example_heading_identifier_mismatch')
            end=next((j for j in range(index+1,len(blocks)) if blocks[j].tag=='heading' or HEADING.match(text(blocks[j]))),len(blocks))
            files=[f for b in blocks[index:end] for f in refs(b)]
    else:errors.append('unsupported_structure_proof')
    graphs,referrors=resolve_references(files,structures);errors+=referrors
    if len(graphs)!=1 or graphs[0]['source_structure_id']!=r['source_structure_id']:errors.append('structure_source_not_uniquely_linked')
    elif graphs[0]['SMILES']!=r['SMILES'] or graphs[0]['InChIKey']!=r['InChIKey']:errors.append('structure_identity_mismatch')
    # Unit/category evidence must be a literal part of the original anchor.
    legend=evidence['legend'];ue=r['unit_evidence']
    if legend:
        anchor=anchors.get(legend['locator'],'').replace('≦','≤').replace('≧','≥')
        if legend['text'] not in anchor or legend['symbol']!=cell:errors.append('legend_not_in_original_scope')
        checked=legend_for(cell,r['Activity Type'],[dict(locator=legend['locator'],text=anchor)])
        if checked is None or any(checked[k]!=r[k] for k in ('Relation','Value','value_lower','value_upper','Unit')):errors.append('legend_bounds_or_units_mismatch')
    else:
        rel,value,lower,upper=source_measurement(cell)
        if (rel,value,lower,upper)!=(r['Relation'],r['Value'],r['value_lower'],r['value_upper']):errors.append('raw_numeric_value_mismatch')
    if ue=='activity_cell_unit' and units(cell)!={r['Unit']}:errors.append('cell_unit_mismatch')
    if ue=='activity_column_unit' and units(header)!={r['Unit']}:errors.append('header_unit_mismatch')
    if ue=='explicit_endpoint_unit_statement':
        u=evidence['unit_statement'];anchor=anchors.get(u['locator'],'')
        if u['text'] not in anchor or units(u['text'])!={r['Unit']}:errors.append('unit_statement_source_mismatch')
    factor=FACTORS.get(r['Unit'])
    for field,dest in [('Value','value_nM'),('value_lower','lower_nM'),('value_upper','upper_nM')]:
        expected=float(Decimal(str(r[field]))*factor) if factor and r[field] is not None else None
        if r[dest]!=expected:errors.append('nM_conversion_mismatch:'+field)
    return errors


def extract(spec):
    archive=spec['archive'];h=hashlib.sha256(archive.encode()).hexdigest()
    cp=CAT/h[:2]/(h+'.json');catalog=json.loads(cp.read_text());metadata=catalog['stats']
    if metadata['archive']!=archive:raise ValueError('catalog_archive_mismatch')
    stat=(RAW/archive).stat()
    expected=hashlib.sha256(f'{CAT_CODE_SHA256}|{stat.st_size}|{stat.st_mtime_ns}|{archive}'.encode()).hexdigest()
    if catalog['signature']!=expected:raise ValueError('frozen_catalog_or_zip_signature_mismatch')
    structures=pq.read_table(cp.with_suffix('.parquet'),columns=STRUCT_COLS).to_pylist()
    stats=dict(archive=archive,Patent=metadata['Patent'],phase=spec.get('phase','fixed_sample'),prior_status=spec.get('status',''),status='parsed',records=0,numeric_SAR=0,extended_activity=0,review=0,baseline_duplicates=0,tables=0,unparsed_activity_tables=0,source_format='',error='',review_reasons={})
    existing={}
    if spec.get('paired_rows',0):
        qcp=OLD/h[:2]/'quality'/(h[:20]+'.json')
        baseline=json.loads(qcp.read_text())
        for r in baseline['ready']:existing[r['activity_locator']]=r
    records=[]
    with zipfile.ZipFile(RAW/archive) as z:
        source=select_patent_xml(z,archive)
        sgml=False
        if source is None:
            choices=[i.filename for i in z.infolist() if i.filename.lower().endswith(('.sgm','.sgml'))]
            if len(choices)!=1:raise ValueError('no_unique_xml_or_sgml')
            source=choices[0];sgml=True
        raw=z.read(source);xml_sha=hashlib.sha256(raw).hexdigest()
        root,adapter_audit=parse_source(raw,sgml)
        description=root.find('.//description')
        if description is None:description=root.find('.//sdode')
        if description is None:raise ValueError('missing_description')
        blocks=blocks_of(description);tables=table_index(description,blocks)
        stats.update(tables=len(tables),source_format='legacy_sgml' if sgml else 'modern_xml' if root.tag=='us-patent-grant' else 'legacy_xml',adapter_audit=adapter_audit)
        links,conflicts=build_links(tables,blocks,structures);molchecks={};proof=proof_index(tables,blocks)
        for t in tables:
            headers=t['headers'];actcols=[i for i,s in enumerate(headers) if ACTIVITY.search(s)]
            if not actcols:
                if ACTIVITY.search(text(t['element'])):stats['unparsed_activity_tables']+=1
                continue
            aliascols=[i for i,s in enumerate(headers) if ALIAS_HEADER.search(s) and i not in actcols and not re.search(r'name|mass',s,re.I)]
            scope,headings=contexts(t,blocks)
            for row in t['rows']:
                for ac in actcols:
                    rawvalue=row[ac]['text']
                    if not rawvalue or ACTIVITY.search(rawvalue):continue
                    locator=t['locator']+'/'+row[0]['row_locator']+'/column-%d'%(ac+1)
                    reasons=[];matched=None;original_alias='';a=None
                    # Repeated Test Comp / IC50 pairs use the preceding alias column.
                    eligible=[i for i in aliascols if i<ac and not any(i<j<ac for j in actcols)]
                    if len(eligible)==1:original_alias=row[eligible[0]]['text'];a=alias(original_alias)
                    elif len(aliascols)==1:original_alias=row[aliascols[0]]['text'];a=alias(original_alias)
                    colid=re.search(r'\b(?:Example|Compound|Ex\.)\s*#?\s*(\d+(?:[.-]\d+)?[a-z]?)\b',headers[ac],re.I)
                    if colid:
                        if a and a!=alias(colid[1]):reasons.append('row_column_compound_conflict')
                        a=alias(colid[1]);original_alias=colid[0]
                    graph,referrors=resolve_references([f for c in row for f in c['files']],structures)
                    if a in conflicts:reasons.append('conflicting_structure_for_alias')
                    elif a:matched=links.get(a)
                    if len(graph)==1 and not referrors and (not aliascols or len({alias(row[i]['text']) for i in aliascols}-{None})<=1):
                        direct=dict(graph=graph[0],method='explicit_activity_row_structure',locator=t['locator']+'/'+row[0]['row_locator'],evidence=dict(row=[c['text'] for c in row],referenced_files=[f for c in row for f in c['files']]))
                        if matched and matched['graph']['InChIKey']!=graph[0]['InChIKey']:reasons.append('activity_row_structure_conflict')
                        else:matched=direct
                    if referrors:reasons.extend(referrors)
                    if not a and not matched:reasons.append('compound_alias_missing_or_ambiguous')
                    if not matched:reasons.append('no_explicit_concrete_structure_link')
                    elif not matched['graph']['InChIKey'] or matched['graph']['structure_status'] not in ('concrete_single_fragment','concrete_multifragment','concrete'):
                        reasons.append('nonconcrete_structure:'+matched['graph']['structure_status'])
                    target,te,target_anchor=target_for(headers[ac],row,headers,scope,headings)
                    if not target:reasons.append('target_not_explicit')
                    endpoints={re.sub(r'\s','',m[1]).lower() for m in ACTIVITY.finditer(headers[ac])}
                    if len(endpoints)!=1:reasons.append('ambiguous_multiple_endpoints_in_column')
                    if re.search(r'[/;&]|\band\b',target,re.I):reasons.append('multiple_target_identities_require_split')
                    ep=ACTIVITY.search(headers[ac])[1].replace(' ','')
                    ep=('p'+ep[1:].upper() if ep.lower().startswith('p') else ep.upper()) if ep.lower() not in ('ki','kd') else 'K'+ep[-1].lower()
                    decoded,measurement_errors=decode(rawvalue,ep,headers[ac],scope)
                    reasons+=measurement_errors
                    if CONTROL.search(headers[ac]) or any(CONTROL.search(c['text']) for i,c in enumerate(row) if i not in actcols):reasons.append('control_not_experimental_compound')
                    evidence=dict(activity_cell=dict(locator=locator,text=rawvalue,header=headers[ac]),alias=dict(text=original_alias,column=eligible[0] if len(eligible)==1 else None),target=target_anchor,structure=matched['evidence'] if matched else {},unit_statement=decoded.get('unit_statement',{}),legend=decoded['legend_evidence'])
                    if any('[UNRESOLVED_ENTITY:' in s for s in [rawvalue,headers[ac],original_alias,target,json.dumps(evidence,ensure_ascii=False)]):reasons.append('unresolved_entity_in_record_evidence')
                    if sgml and adapter_audit['sgml_closure_warnings']:reasons.append('sgml_closure_requires_review')
                    if matched and not reasons:reasons+=verify_graph(z,matched['graph'],molchecks)
                    if matched and not reasons and matched['evidence'].get('reported_product_name'):
                        g=matched['graph'];mol=Chem.MolFromSmiles(g['SMILES']);counts=Counter(a.GetSymbol() for a in mol.GetAtoms())
                        if name_element_conflicts(matched['evidence']['reported_product_name'],counts):reasons.append('reported_product_name_element_conflict')
                    # A source cell already formally accepted cannot add a second observation.
                    if locator in existing:
                        stats['baseline_duplicates']+=1;continue
                    r={k:None if k in NUMERIC or k in INTEGER else False if k in BOOL else '' for k in RECORD_COLUMNS}
                    cellkey=hashlib.sha256((archive+'|'+source+'|'+locator).encode()).hexdigest()
                    r.update(source='USPTO_MOL_RECOVERY',Patent=metadata['Patent'],patent_core=metadata['patent_core'],archive=archive,xml_file=source,xml_sha256=xml_sha,source_cell_key=cellkey,source_record_id='USPTO-RECOVERY:'+cellkey,Target=target,**{'Activity Type':ep},raw_value=rawvalue,original_range=rawvalue,activity_locator=locator,activity_header=headers[ac],compound_alias=a or '',original_alias=original_alias,compound_id=metadata['Patent']+':'+a if a else '',alias_evidence='explicit_activity_identifier' if a else 'direct_source_graph_without_experimental_id' if matched else '',target_evidence=te,target_scope='reported_assay_or_cell_identity',source_format=stats['source_format'],phase=stats['phase'],prior_status=stats['prior_status'])
                    for k,v in decoded.items():
                        if k in r:r[k]=json.dumps(v,ensure_ascii=False) if k=='legend_evidence' else v
                    if matched:
                        g=matched['graph'];r.update(SMILES=g['SMILES'],SMILES_original=g['SMILES_original'],InChIKey=g['InChIKey'],source_structure_id=g['source_structure_id'],structure_file=g['structure_file'],source_header_offset=g['source_header_offset'],raw_mol_sha256=g['source_sha256'],raw_mol_coordinate_dimension=g['raw_mol_coordinate_dimension'],structure_locator=matched['locator'],structure_match=matched['method'])
                        # Original field annotations may refer to different experiment aliases;
                        # carry a sequence only when the catalog explicitly linked it to this graph.
                        r['Sequence']=g.get('Sequence','') or ''
                        r['compound_IUPAC_name']=g.get('compound_IUPAC_name','') or ''
                        r['compound_name_status']=g.get('compound_name_status','') or ''
                    r['source_evidence_json']=json.dumps(evidence,ensure_ascii=False,separators=(',',':'))
                    if not reasons:reasons+=verify_record_source(r,proof,tables,blocks,structures)
                    r['review_reason']=';'.join(dict.fromkeys(reasons))
                    r['extraction_status']='review' if reasons else 'extended_activity' if decoded['Relation']=='range' or decoded['Unit'] not in FACTORS else 'numeric_SAR'
                    r['source_mapping_checks']='passed_raw_cell_alias_scope_mol_sha_smiles_inchikey_conversion' if not reasons else 'requires_review'
                    records.append(r)
        stats.update(records=len(records),numeric_SAR=sum(r['extraction_status']=='numeric_SAR' for r in records),extended_activity=sum(r['extraction_status']=='extended_activity' for r in records),review=sum(r['extraction_status']=='review' for r in records),review_reasons=dict(Counter(reason for r in records for reason in r['review_reason'].split(';') if reason)))
    return stats,records
