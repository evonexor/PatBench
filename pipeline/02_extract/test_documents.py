"""Regression tests for standalone JSON, roles and citation provenance."""
import argparse
import json
from pathlib import Path
import unittest
import xml.etree.ElementTree as ET
from patent_documents import extract_documents
from patent_metadata import document

def fixture(bib='',abstract='Abstract text'):
    return ET.fromstring('<us-patent-grant lang="EN"><us-bibliographic-data-grant><publication-reference><document-id><country>US</country><doc-number>123</doc-number><kind>B2</kind><date>20250107</date></document-id></publication-reference>'+bib+'</us-bibliographic-data-grant><abstract>'+abstract+'</abstract></us-patent-grant>')

def extract(root):return extract_documents(root,'US123B2','fixture.zip','fixture.xml',document)

class DocumentTests(unittest.TestCase):
    def test_applicant_and_inventor_roles_kept_separate(self):
        b='<us-parties><us-applicants><us-applicant app-type="applicant"><addressbook><orgname>Example Co.</orgname></addressbook></us-applicant></us-applicants><inventors><inventor><addressbook><first-name>Min</first-name><last-name>Teng</last-name></addressbook></inventor></inventors></us-parties>'
        m,doc=extract(fixture(b))
        self.assertEqual(m['applicant_all_names'],'Example Co.')
        self.assertEqual(m['inventor_names'],'Min Teng')
        self.assertEqual(doc['applicants'][0]['name_type'],'organization')
    def test_personal_applicant_not_automatically_an_inventor(self):
        b='<parties><applicants><applicant app-type="applicant"><addressbook><first-name>A</first-name><last-name>B</last-name></addressbook></applicant></applicants></parties>'
        m,_=extract(fixture(b))
        self.assertEqual(m['inventor_names'],'')
        self.assertEqual(m['applicant_all_names'],'A B')
    def test_citations_not_mother_applications(self):
        b='<us-references-cited><us-citation><patcit><document-id><country>US</country><doc-number>999</doc-number><kind>B2</kind></document-id></patcit></us-citation><us-citation><nplcit><othercit>Smith et al., Journal.</othercit></nplcit></us-citation></us-references-cited>'
        m,doc=extract(fixture(b))
        self.assertEqual(m['parent_application_numbers'],'')
        self.assertEqual(m['patent_citation_count'],1);self.assertEqual(m['npl_citation_count'],1)
        self.assertEqual(doc['citations'][0]['patent']['number'],'US999B2')
    def test_application_date_and_language_keep_source(self):
        b='<application-reference appl-type="utility"><document-id><country>US</country><doc-number>18616636</doc-number><date>20240326</date></document-id></application-reference>'
        m,doc=extract(fixture(b,abstract='First sentence. Second sentence.'))
        self.assertEqual(m['application_number'],'US18616636')
        self.assertEqual(m['application_date'],'2024-03-26')
        self.assertEqual(m['abstract_language'],'en')
        self.assertEqual(doc['abstract']['text'],'First sentence. Second sentence.')

def integration(out):
    import pyarrow.parquet as pq
    pats=[r for r in pq.read_table(out/'patents.parquet').to_pylist() if r['status']=='processed']
    assert len(pats)==40
    assert all(r['abstract'] and r['inventor_names'] and r['application_date'] for r in pats)
    for r in pats:
        path=out/r['bibliography_json_path']
        assert path.is_file() and path.resolve().is_relative_to(out.resolve())
        doc=json.loads(path.read_text());assert doc['Patent']==r['Patent']
        assert len([c for c in doc['citations'] if c['citation_type']=='patent'])==r['patent_citation_count']
        assert len([c for c in doc['citations'] if c['citation_type']=='non_patent'])==r['npl_citation_count']
    assert sum(bool(r['parent_application_numbers']) for r in pats)==21
    assert sum(r['patent_citation_count'] for r in pats)==2853
    assert sum(r['npl_citation_count'] for r in pats)==2498
    forbidden={'claim_count_declared','claim_count_xml','claim_count_status','independent_claim_candidate_numbers','independent_claim_candidate_count','claim_dependency_status','claims_json_path'}
    for name in ['patents','SAR_enriched','SAR_with_evidence','review']:
        assert not forbidden.intersection(pq.read_schema(out/(name+'.parquet')).names)
        import csv
        with (out/(name+'.csv')).open(encoding='utf-8-sig') as f:
            assert not forbidden.intersection(next(csv.reader(f)))
    assert not (out/'claims').exists()
    example=next(r for r in pats if r['patent_core']=='US12187754')
    assert example['application_number']=='US18616636' and example['application_date']=='2024-03-26'
    assert example['inventor_names']=='Min Teng'
    sar=pq.read_table(out/'SAR_enriched.parquet').to_pylist();patmap={r['Patent']:r for r in pats}
    assert all(r['bibliography_json_path']==patmap[r['Patent']]['bibliography_json_path'] for r in sar)
    assert pq.read_schema(out/'patents.parquet').field('inventor_count').type.bit_width==64
    print('Documents integration passed: 40 bibliography JSON, links, roles, citation counts, SAR joins; claim exports absent.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path);a=p.parse_args()
    r=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(DocumentTests))
    if not r.wasSuccessful():raise SystemExit(1)
    if a.results:integration(a.results)
