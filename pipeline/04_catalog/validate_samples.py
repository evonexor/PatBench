"""Independent comparison with existing ChEMBL. Never changes extraction results."""
import argparse
from collections import Counter, defaultdict
import csv
from decimal import Decimal
import json
import math
from pathlib import Path
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger
from uspto_sar import core, UNITS, unit

RDLogger.DisableLog('rdApp.*')

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--results',type=Path,required=True)
    a=p.parse_args()
    observations=pq.read_table(a.results/'SAR_with_evidence.parquet').to_pylist()
    patents={r['patent_core'] for r in observations}
    reference=defaultdict(list);cache={}
    columns=['Patent','Target','Activity Type','Relation','Value','Unit','SMILES','activity_id','UniProt','organism','assay_description']
    relevant=0
    for b in pq.ParquetFile(a.root/'standardized/ChEMBL.parquet').iter_batches(columns=columns,batch_size=16384):
        d=b.to_pydict()
        for i,pat in enumerate(d['Patent']):
            pc=core(pat or '')
            if pc not in patents:continue
            r={k:d[k][i] for k in columns}
            smiles=r['SMILES']
            if smiles not in cache:
                m=Chem.MolFromSmiles(smiles or '')
                cache[smiles]=Chem.MolToInchiKey(m) if m else ''
            k=cache[smiles]
            if not k:continue
            try:v=float(Decimal(str(r['Value']))*UNITS[unit(str(r['Unit']))])
            except (KeyError,ValueError,ArithmeticError):continue
            r['value_nM']=v
            reference[(pc,k,str(r['Activity Type']).upper())].append(r);relevant+=1
    rows=[]
    for r in observations:
        refs=reference.get((r['patent_core'],r['InChIKey'],r['Activity Type'].upper()),[])
        same=[v for v in refs if v['Relation']==r['Relation'] and math.isclose(v['value_nM'],r['value_nM'],rel_tol=1e-7,abs_tol=1e-9)]
        # A structure/value match does not by itself establish target or species equivalence.
        status='structure_endpoint_value_match' if same else 'same_structure_other_value' if refs else 'no_same_structure_reference'
        chosen=same or refs
        rows.append(dict(source_record_id=r['source_record_id'],Patent=r['Patent'],compound_alias=r['compound_alias'],extracted_target=r['Target'],extracted_value_nM=r['value_nM'],Relation=r['Relation'],status=status,reference_targets='|'.join(sorted({v['Target'] or '' for v in chosen})),reference_values_nM='|'.join(sorted({str(v['value_nM']) for v in chosen})),reference_activity_ids='|'.join(str(v['activity_id']) for v in chosen),note='Target and species equivalence requires separate review.'))
    with (a.results/'reference_validation.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['source_record_id','status']);w.writeheader();w.writerows(rows)
    summary=dict(reference_source='standardized/ChEMBL.parquet',reference_rows_selected=relevant,extracted_observations=len(rows),comparison_counts=dict(Counter(r['status'] for r in rows)),matches_by_patent=dict(Counter(r['Patent'] for r in rows if r['status']=='structure_endpoint_value_match')),meaning='Comparison checks full InChIKey, endpoint, relation and nM value. Different values can be different assays; no-reference records are not automatically errors. No reference activities were used to generate the extracted records.')
    (a.results/'validation_report.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
    print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
