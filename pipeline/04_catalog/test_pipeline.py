"""Scientific parsing tests, plus optional checks against actual extracted fixtures."""
import argparse
import json
from pathlib import Path
import unittest
import xml.etree.ElementTree as ET
from uspto_sar import expanded_rows, EXAMPLE_HEADING, measurement, screening, unit, name_element_conflicts

class ParserTests(unittest.TestCase):
    def test_cals_header_split_into_body_and_offset(self):
        t=ET.fromstring('<table><tgroup cols="4"><colspec colname="offset"/><colspec colname="1"/><colspec colname="2"/><colspec colname="3"/><thead><row><entry/><entry>Com-</entry><entry>RET IC<sub>50</sub></entry><entry>Structure</entry></row></thead><tbody><row><entry/><entry>pound No.</entry><entry>(μM)</entry><entry/></row><row><entry/><entry>2</entry><entry>&gt;0.5</entry><entry><img file="X.MOL"/></entry></row></tbody></tgroup></table>')
        h,rows=expanded_rows(t)
        self.assertEqual(h[1],'Compound No.')
        self.assertEqual(h[2],'RET IC50 (μM)')
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0][1]['text'],'2')
    def test_relations_ranges_and_missing(self):
        self.assertEqual(measurement('>1500'),('>',1500.0,None,None))
        self.assertEqual(measurement('0.1–0.5 μM'),('range',None,.1,.5))
        self.assertEqual(measurement('++++')[1],None)
        self.assertEqual(measurement('N/A')[1],None)
        self.assertEqual(measurement('0')[1],0.0)
    def test_units_do_not_guess(self):
        self.assertEqual(unit('IC50 (μM)'),'uM')
        self.assertEqual(unit('ATP 25 μM; enzyme 5 nM'),'')
        self.assertEqual(unit('IC50'),'')
    def test_exclusions(self):
        self.assertEqual(screening('Agricultural oil-based agent for treating asphalt','','')[0],'agriculture')
        self.assertEqual(screening('Electrochemical hydrogen pump','','')[0],'materials')
        self.assertEqual(screening('CDK4 inhibitors for cancer','','')[0],'medical')
    def test_grouped_examples_rejected(self):
        self.assertIsNone(EXAMPLE_HEADING.fullmatch('Example #1-1 and 1-2'))
        self.assertEqual(EXAMPLE_HEADING.fullmatch('Example 2—6-(4-Amino)-derivative')[1],'2')
    def test_source_mol_element_contradiction(self):
        self.assertEqual(name_element_conflicts('2,4-difluorophenylsulfonyl',{'F':2,'C':12}),['S'])
        self.assertEqual(name_element_conflicts('2,4-difluorophenylsulfonyl',{'F':2,'S':1,'C':12}),[])

def integration(results):
    import pyarrow.parquet as pq
    sar=pq.read_table(results/'SAR_with_evidence.parquet').to_pylist()
    review=pq.read_table(results/'review.parquet').to_pylist()
    pats=pq.read_table(results/'patents.parquet').to_pylist()
    assert sum(p['status']=='processed' for p in pats)==40
    assert sum(p['status']=='excluded' for p in pats)==2
    assert len({r['source_record_id'] for r in sar+review})==len(sar)+len(review)
    cdk=[r for r in sar if r['patent_core']=='US8841312' and r['compound_alias']=='1']
    assert {(r['Target'].replace(' ',''),r['value_nM']) for r in cdk}=={('Cdk4',3.2),('Cdk6',4.4)}
    assert not any('staurosporine' in r['activity_header'].lower() for r in sar)
    controls=[r for r in review if 'staurosporine' in r['activity_header'].lower()]
    assert controls and all('control_column_not_linked_to_example' in r['review_reason'] for r in controls)
    jak=[r for r in review if r['patent_core']=='US8962608']
    assert jak and all(r['Unit']=='' and r['value_nM'] is None for r in jak)
    assert all(r['SMILES'] and r['InChIKey'] and r['Target'] and r['Unit'] and r['Value'] is not None and not r['review_reason'] for r in sar)
    assert not any(r['patent_core']=='US8889696' and r['compound_alias']=='39' for r in sar)
    assert any(r['patent_core']=='US8889696' and r['compound_alias']=='39' and 'structure_name_element_mismatch' in r['review_reason'] for r in review)
    print(json.dumps(dict(integration='passed',checks=['40 medical + 2 exclusion fixtures','CDK4/6 Example 1 values and unit conversions','Staurosporine controls excluded','JAK unitless table not silently assigned nM','unique stable observation IDs','all strict SAR prerequisites present']),ensure_ascii=False))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path);a=p.parse_args()
    r=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ParserTests))
    if not r.wasSuccessful():raise SystemExit(1)
    if a.results:integration(a.results)
