"""Regression tests for metadata semantics and real sample joins."""
import argparse
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path
from patent_metadata import extract_metadata, date_iso

def fixture(extra='',description='',claims='',abstract=''):
    return ET.fromstring('<us-patent-grant><us-bibliographic-data-grant><publication-reference><document-id><country>US</country><doc-number>123</doc-number><kind>B2</kind><date>20250107</date></document-id></publication-reference><invention-title>Novel compounds</invention-title>'+extra+'</us-bibliographic-data-grant><abstract>'+abstract+'</abstract><description>'+description+'</description><claims>'+claims+'</claims></us-patent-grant>')

class MetadataTests(unittest.TestCase):
    def test_publication_is_not_grant_or_prior_art(self):
        extra='<us-related-documents><related-publication><document-id><country>US</country><doc-number>20240270773</doc-number><kind>A1</kind><date>20240815</date></document-id></related-publication></us-related-documents><us-references-cited><citation><patcit><document-id><date>19990101</date><doc-number>PCT/US1999/111111</doc-number></document-id></patcit></citation></us-references-cited>'
        m=extract_metadata(fixture(extra))
        self.assertEqual(m['application_publication_date'],'2024-08-15')
        self.assertEqual(m['grant_date'],'2025-01-07')
        self.assertEqual(m['pct_number'],'')
    def test_pct_continuation_not_direct_national_phase(self):
        extra='<us-related-documents><continuation><relation><parent-doc><document-id><country>WO</country><doc-number>PCT/US2022/044710</doc-number><date>20220926</date></document-id></parent-doc></relation></continuation></us-related-documents>'
        m=extract_metadata(fixture(extra))
        self.assertEqual(m['pct_relationship'],'PCT母案的继续申请')
    def test_missing_pct_not_confirmed_no(self):
        self.assertEqual(extract_metadata(fixture())['is_pct'],'未识别到PCT（未查完整专利族）')
    def test_priority_conflict_preserves_both_sources(self):
        extra='<us-related-documents><us-provisional-application><document-id><country>US</country><doc-number>61537978</doc-number><date>20110922</date></document-id></us-provisional-application></us-related-documents>'
        description='<heading>CROSS REFERENCE TO RELATED APPLICATIONS</heading><p id="p-0002">Claims priority from U.S. Provisional Application Ser. No. 61/537,978, filed Sep. 21, 2011.</p>'
        m=extract_metadata(fixture(extra,description))
        self.assertEqual(m['priority_date'],'')
        self.assertEqual(m['priority_date_candidate'],'2011-09-22')
        self.assertEqual(m['priority_evidence']['conflicts'][0]['description']['date'],'2011-09-21')
    def test_indication_does_not_use_background_only(self):
        d='<heading>BACKGROUND</heading><p>Other inhibitors treat lung cancer.</p><heading>FIELD</heading><p>The present compounds treat HIV infection.</p>'
        m=extract_metadata(fixture(description=d))
        self.assertEqual(m['indication'],'HIV感染')
        self.assertEqual(m['technical_field_source'],'说明书技术领域章节')
    def test_classification_excludes_search_scope(self):
        extra='<classifications-ipcr><classification-ipcr><section>A</section><class>61</class><subclass>K</subclass><main-group>31</main-group><subgroup>675</subgroup></classification-ipcr></classifications-ipcr><us-field-of-classification-search><classification-cpc><section>A</section><class>01</class><subclass>N</subclass><main-group>1</main-group><subgroup>2</subgroup></classification-cpc></us-field-of-classification-search>'
        m=extract_metadata(fixture(extra))
        self.assertEqual(m['ipc_codes'],'A61K 31/675')
        self.assertEqual(m['cpc_codes'],'')
    def test_invalid_date_stays_missing(self):
        self.assertEqual(date_iso('20250230'),'')
    def test_claimed_use_separated_from_other_disclosed_uses(self):
        d='<heading>SUMMARY</heading><p>Compounds of the invention treat bacterial infections.</p><heading>DETAILED DESCRIPTION</heading><p>Compounds of the invention treat asthma.</p>'
        c='<claim id="CLM-1">1. A method of treating a bacterial infection.</claim>'
        m=extract_metadata(fixture(description=d,claims=c))
        self.assertEqual(m['indication'],'细菌感染')
        self.assertEqual(m['other_disclosed_indications'],'哮喘')
    def test_reissue_date_labelled(self):
        root=fixture();root.find('.//publication-reference/document-id/kind').text='E1'
        m=extract_metadata(root)
        self.assertEqual(m['grant_date'],'2025-01-07')
        self.assertIn('不是原始授权日',m['grant_date_type'])

def integration(path):
    import pyarrow.parquet as pq
    p=pq.read_table(path/'patents.parquet').to_pylist()
    processed=[r for r in p if r['status']=='processed'];by={r['patent_core']:r for r in processed}
    sar=pq.read_table(path/'SAR_enriched.parquet').to_pylist()
    assert len(processed)==40 and len(sar)==2109
    assert all(r['title'] and r['company_name'] and r['grant_date'] for r in processed)
    assert by['US8841312']['company_name']=='Amgen Inc.'
    assert by['US8841312']['application_publication_date']=='2011-06-16'
    assert by['US8841312']['priority_date']=='2007-12-19'
    assert by['US8962608']['priority_date']==''
    assert 'priority_date_conflict' in by['US8962608']['metadata_review_reason']
    assert by['US12187754']['pct_relationship']=='PCT母案的继续申请'
    assert by['US12187754']['company_name']=='BLACKSMITH MEDICINES, INC.'
    assert 'HIV感染' in by['US12187753']['indication']
    assert all(r['company_name']==by[r['patent_core']]['company_name'] and r['priority_date']==by[r['patent_core']]['priority_date'] for r in sar)
    print('Metadata integration passed: 40 patents, 2109 joined SAR rows, date/PCT/source boundaries checked.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--results',type=Path);a=parser.parse_args()
    r=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(MetadataTests))
    if not r.wasSuccessful():raise SystemExit(1)
    if a.results:integration(a.results)
