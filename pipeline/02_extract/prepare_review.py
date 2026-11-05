"""Prepare typed review views from the extracted Parquet, without changing results."""
import csv
import json
from pathlib import Path
import pyarrow.parquet as pq
from patent_metadata import REVIEW_METADATA_COLUMNS
from patent_documents import DOCUMENT_COLUMNS

BASE=Path(__file__).resolve().parent.parent
RESULTS=BASE/'results'
reasons={
    'no_explicit_concrete_structure_link':'没有明确的具体结构关联',
    'activity_unit_missing':'活性单位缺失',
    'conflicting_structure_for_alias':'同一编号对应多个结构',
    'qualitative_or_unparsed_value':'定性分级或数值格式需核对',
    'control_column_not_linked_to_example':'阳性对照列',
    'control_row_not_linked_to_example':'对照化合物行',
    'compound_alias_missing_or_ambiguous':'化合物编号缺失或不明确',
    'interval_requires_extended_schema':'区间值需扩展字段',
    'target_not_explicit':'靶点或测量对象不明确',
    'structure_name_element_mismatch':'原始MOL与名称所需元素矛盾',
}
def load(n):return pq.read_table(RESULTS/(n+'.parquet')).to_pylist()

ready=load('SAR_with_evidence');review=load('review');patents=load('patents')
refs=list(csv.DictReader((RESULTS/'reference_validation.csv').open(encoding='utf-8-sig')))
refmap={r['source_record_id']:r['status'] for r in refs}
status_labels={'structure_endpoint_value_match':'结构及活性值一致','same_structure_other_value':'同结构存在其他活性值','no_same_structure_reference':'没有同结构参考记录'}
for r in ready+review:
    r['assay_context_excerpt']=r['assay_context'][:240]+('…' if len(r['assay_context'])>240 else '')
    r['reference_status']=status_labels.get(refmap.get(r['source_record_id'],''),'未作参考比较')
    r['review_reason_cn']='；'.join(reasons.get(k,k) for k in r['review_reason'].split(';') if k)
for r in ready+review+patents:
    field=r.get('technical_field') or ''
    r['technical_field_excerpt']=field[:350]+('…（完整原文见CSV）' if len(field)>350 else '')
    for key in ['indication_evidence','priority_evidence','pct_evidence','classification_basis']:
        evidence=r.get(key) or ''
        r[key+'_excerpt']=evidence[:450]+('…（完整依据见CSV）' if len(evidence)>450 else '')
for r in refs:
    r['extracted_value_nM']=float(r['extracted_value_nM']) if r['extracted_value_nM'] else None
    r['status_cn']=status_labels[r['status']]
for r in patents:
    r['result_note']=('当前规则未提取出可配对的定量活性' if r['status']=='processed' and not r['paired_rows'] else '已产出配对记录' if r['status']=='processed' else '农业或材料，未抽取')
columns=['Patent','compound_alias','Target','Activity Type','Relation','Value','Unit','value_nM','raw_value','SMILES','InChIKey','reference_status','structure_match','activity_header','unit_evidence','structure_file','structure_name','structure_name_check','structure_locator','activity_locator','source_record_id','UniProt','organism','gene','SMILES_original','source','patent_core','final_label','archive','xml_file','publication_date','assay_context_excerpt']
payload=dict(report=json.loads((RESULTS/'run_report.json').read_text()),validation=json.loads((RESULTS/'validation_report.json').read_text()),performance=json.loads((RESULTS/'performance_report.json').read_text()),sheets=[
    dict(name='已配对SAR',columns=columns+REVIEW_METADATA_COLUMNS+DOCUMENT_COLUMNS,rows=ready),
    dict(name='待核对',columns=['Patent','compound_alias','Target','Activity Type','raw_value','Unit','value_nM','review_reason_cn','SMILES','InChIKey','Relation','Value','value_lower','value_upper','activity_header','structure_match','structure_file','structure_name','structure_name_check','activity_locator','structure_locator','source_record_id','review_reason','assay_context_excerpt','archive','xml_file']+REVIEW_METADATA_COLUMNS+DOCUMENT_COLUMNS,rows=review),
    dict(name='专利清单',columns=list(dict.fromkeys(['Patent','title','selection','paired_rows','review_rows','mol_files','explicit_aliases','status','result_note','publication_date','screen','unparsed_activity_tables','structure_status_counts','archive','xml_file']+REVIEW_METADATA_COLUMNS+['indication_evidence_excerpt','priority_evidence_excerpt','pct_evidence_excerpt','classification_basis_excerpt','technical_field_locator','wo_publication_number','wo_publication_date']+DOCUMENT_COLUMNS)),rows=patents),
    dict(name='参考核对',columns=['Patent','compound_alias','extracted_target','extracted_value_nM','Relation','status_cn','reference_targets','reference_values_nM','reference_activity_ids','source_record_id'],rows=refs),
])
payload['metadata']=json.loads((RESULTS/'metadata_report.json').read_text())
payload['documents']=json.loads((RESULTS/'document_report.json').read_text())
(BASE/'review_input.json').write_text(json.dumps(payload,ensure_ascii=False,separators=(',',':')))
# Compact, patent-level review file; the two exclusion fixtures are not included.
sample_columns=['Patent','title','company_name','is_pct','pct_number','pct_relationship','priority_date','priority_date_status','application_publication_number','application_publication_date','grant_date','grant_date_type','indication','indication_source','other_disclosed_indications','ipc_codes','cpc_codes','patent_category','technical_field','technical_field_source','metadata_review_reason','archive','xml_file','wo_publication_number','wo_publication_date']+DOCUMENT_COLUMNS
with (RESULTS/'metadata_sample.csv').open('w',encoding='utf-8-sig',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=sample_columns);writer.writeheader()
    writer.writerows({k:r.get(k) for k in sample_columns} for r in patents if r['status']=='processed')
print(json.dumps(dict(ready=len(ready),review=len(review),patents=len(patents),validation=len(refs)),ensure_ascii=False))
