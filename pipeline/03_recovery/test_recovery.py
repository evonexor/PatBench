"""Source-association negative controls, thresholds, intervals and legacy layout."""
import unittest
import json
from copy import deepcopy
import xml.etree.ElementTree as ET
from activity_recovery import parse_source, text, blocks_of, table_index, build_links, resolve_references, decode, target_for, proof_index, verify_record_source


def graph(stem,key):
    return dict(source_structure_id=stem,structure_file=stem+'.MOL',InChIKey=key,SMILES=key)


class RecoveryTests(unittest.TestCase):
    def test_category_threshold_not_point(self):
        scope=[dict(locator='p1',text='A indicates an IC50 of >1 μM B indicates an IC50 of <1 μM C indicates an IC50 of <1 nM')]
        r,errors=decode('B','IC50','BZ binding Assay IC50',scope)
        self.assertFalse(errors);self.assertEqual((r['Relation'],r['Value'],r['value_nM']),('<',1.0,1000.0))
        self.assertEqual(r['measurement_kind'],'legend_threshold')

    def test_interval_keeps_limits(self):
        r,errors=decode('++','IC50','IC50 (μM)',[dict(locator='p',text='"++" indicates 1000≦IC50≦10,000 μM')])
        self.assertFalse(errors);self.assertIsNone(r['Value']);self.assertEqual((r['value_lower'],r['value_upper']),(1000.,10000.))

    def test_no_legend_no_unit_guess(self):
        r,errors=decode('B','IC50','IC50',[dict(locator='p',text='Incubated cells with 10 μM ATP.')])
        self.assertIn('activity_unit_missing',errors);self.assertIn('qualitative_value_without_unique_definition',errors)

    def test_conflicting_legend_rejected(self):
        _,errors=decode('A','IC50','IC50',[dict(locator='p',text='A indicates IC50 <1 nM; A indicates IC50 >1 nM')])
        self.assertTrue(errors)

    def test_mixed_units_rejected(self):
        _,errors=decode('1 nM','IC50','IC50 (μM)',[])
        self.assertIn('conflicting_units',errors)

    def test_slash_not_averaged(self):
        r,errors=decode('5.03/2.09','IC50','T gondii IC50 (μM)',[])
        self.assertIsNone(r['Value']);self.assertTrue(errors)

    def test_duplicate_filename_not_selected(self):
        refs,errors=resolve_references(['A.TIF'],[graph('A','X'),graph('A','Y')])
        self.assertTrue(errors);self.assertFalse(refs)

    def test_file_number_not_compound_alias(self):
        root,_=parse_source(b'<description><p><img file="-C00001.TIF"/></p></description>')
        links,conflicts=build_links([],blocks_of(root),[graph('-C00001','X')])
        self.assertFalse(links);self.assertFalse(conflicts)

    def test_full_example_unique_only(self):
        root,_=parse_source(b'<description><heading id="h">Example 1</heading><p><img file="A.TIF"/></p><p><img file="B.TIF"/></p><heading>Biological Assay</heading></description>')
        links,_=build_links([],blocks_of(root),[graph('A','X'),graph('B','Y')]);self.assertNotIn('1',links)
        root,_=parse_source(b'<description><heading id="h">Example 1 N-methylbenzamide</heading><p><img file="A.TIF"/></p><heading>Biological Assay</heading></description>')
        links,_=build_links([],blocks_of(root),[graph('A','X')]);self.assertEqual(links['1']['graph']['InChIKey'],'X')

    def test_target_only_scoped(self):
        target,_,_=target_for('IC50 (nM)',[dict(text='3')],['Compound'],[],[]);self.assertEqual(target,'')
        target,_,_=target_for('TMD8 Prol GI50 [uM]',[],[],[],[]);self.assertEqual(target,'TMD8 Prol')

    def test_legacy_void_and_subscripts(self):
        raw=b'<PATDOC><SDODE><H>IDO Enzyme Assay</H><TABLE-US ID="T"><TABLE><TGROUP COLS="2"><COLSPEC COLNAME="A"><COLSPEC COLNAME="B"><TBODY><ROW><ENTRY>Test Comp.</ENTRY><ENTRY>IC<SB>50</SB> (nM)</ENTRY></ROW><ROW><ENTRY>1</ENTRY><ENTRY>2</ENTRY></ROW></TBODY></TGROUP></TABLE></TABLE-US></SDODE></PATDOC>'
        root,audit=parse_source(raw,True);self.assertFalse(audit['sgml_closure_warnings'])
        d=root.find('.//sdode');ts=table_index(d,blocks_of(d));self.assertEqual(ts[0]['headers'][1],'IC50 (nM)');self.assertEqual(ts[0]['rows'][0][1]['text'],'2')

    def test_unknown_entity_preserved(self):
        root,audit=parse_source(b'<description><p>IC50 &unknownxyz;</p></description>')
        self.assertIn('[UNRESOLVED_ENTITY:unknownxyz]',text(root));self.assertTrue(audit['unresolved_entities'])

    def test_generic_assay_not_target(self):
        target,_,_=target_for('IC50 (nM)',[],[],[],[dict(locator='h',text='Biological Activity Assay')]);self.assertEqual(target,'')

    def test_concentration_suffix_not_percent_or_mass(self):
        for value in ('50%','2 mg','1,23'):
            _,errors=decode(value,'IC50','IDO IC50 (nM)',[]);self.assertTrue(errors)

    def test_mass_concentration_not_molar_conversion(self):
        r,errors=decode('4 ug/mL','MIC','S. aureus MIC (ug/mL)',[])
        self.assertFalse(errors);self.assertEqual(r['Value'],4.);self.assertEqual(r['Unit'],'ug/mL');self.assertIsNone(r['value_nM'])

    def test_legacy_encoding_preserved(self):
        raw='<?xml version="1.0" encoding="iso-8859-1"?><description><p>5 µM</p></description>'.encode('iso-8859-1')
        root,_=parse_source(raw);self.assertEqual(text(root),'5 µM')

    def test_manifest_phase_ids_normalized(self):
        from run_recovery import normalize_specs
        values=normalize_specs([dict(phase=1),dict(phase=2),dict(phase=3)])
        self.assertEqual([s['phase'] for s in values],['review_only','previously_skipped','remaining_processed'])

    def test_source_checker_rejects_mixed_records(self):
        source=b'<description><table id="T"><tgroup cols="2"><thead><row><entry>Compound</entry><entry>IDO IC50 (nM)</entry></row></thead><tbody><row><entry>1<img file="A.TIF"/></entry><entry>2</entry></row></tbody></tgroup></table></description>'
        root,_=parse_source(source);blocks=blocks_of(root);tables=table_index(root,blocks);proof=proof_index(tables,blocks);g=graph('A','X')
        record=dict(activity_locator='T/tgroup-1/row-2/column-2',activity_header='IDO IC50 (nM)',raw_value='2',compound_alias='1',original_alias='1',Target='IDO',target_evidence='explicit_activity_column',structure_locator='T/tgroup-1/row-2',structure_match='explicit_activity_row_structure',source_structure_id='A',SMILES='X',InChIKey='X',**{'Activity Type':'IC50'})
        decoded,_=decode('2','IC50',record['activity_header'],[]);record.update(decoded)
        record['source_evidence_json']=json.dumps(dict(alias=dict(column=0),target=dict(text=record['activity_header']),legend={},unit_statement={}))
        self.assertFalse(verify_record_source(record,proof,tables,blocks,[g]))
        for key,value in [('Target','VEGF'),('compound_alias','2'),('raw_value','3'),('Value',9.),('Unit','uM'),('source_structure_id','B'),('value_nM',3.)]:
            changed=deepcopy(record);changed[key]=value
            self.assertTrue(verify_record_source(changed,proof,tables,blocks,[g]),key)


if __name__=='__main__':unittest.main()
