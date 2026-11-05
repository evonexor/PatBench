import hashlib,json,tempfile,unittest,zipfile
from pathlib import Path
import pyarrow.parquet as pq
from rdkit import Chem
from catalog_annotations import explicit_uniprot,normalize_sequence,parse_sequences,extract
from structure_catalog import worker,selected_sources

class CatalogTests(unittest.TestCase):
    def test_uniprot_requires_label(self):
        self.assertEqual(explicit_uniprot('Compound P12345 and product A2BC19'),[])
        self.assertEqual(explicit_uniprot('UniProt accession P12345; Swiss-Prot A0A023GPI8'),['A0A023GPI8','P12345'])

    def test_dna_excludes_translation_lines(self):
        s,bad=normalize_sequence('atg ccc 6\nMet Pro\n1 2','DNA')
        self.assertEqual(s,'ATGCCC');self.assertFalse(bad)

    def test_st26_amino_acids_are_single_letters(self):
        data=b'<ST26SequenceListing><SequenceData sequenceIDNumber="1"><INSDSeq><INSDSeq_length>3</INSDSeq_length><INSDSeq_moltype>AA</INSDSeq_moltype><INSDSeq_sequence>MET</INSDSeq_sequence></INSDSeq></SequenceData></ST26SequenceListing>'
        row=parse_sequences(data,'seq.xml','US1B2','archive')[0]
        self.assertEqual(row['Sequence'],'MET');self.assertEqual(row['sequence_status'],'parsed')

    def test_sequence_cwu_and_incomplete_listing(self):
        data=b'<sequence-cwu><sequence-list-new-rules><s200><s210>1</s210><s211>4</s211><s212>PRT</s212></s200><s400>1\nVal His Leu Thr\n1</s400><s200><s210>2</s210><s211>100</s211><s212>DNA</s212></s200><s400>2\nacgt 4</s400></sequence-list-new-rules></sequence-cwu>'
        rows=parse_sequences(data,'seq.xml','US1B2','archive')
        self.assertEqual(rows[0]['Sequence'],'VHLT');self.assertEqual(rows[0]['sequence_status'],'parsed')
        self.assertEqual(rows[1]['Sequence'],'ACGT');self.assertEqual(rows[1]['sequence_status'],'length_mismatch_review')

    def test_st25_text(self):
        rows=parse_sequences(b'<210> 7\n<211> 4\n<212> PRT\n<213> human\n<400> 7\nVal His Leu Thr\n1\n','seq.txt','US1B2','archive')
        self.assertEqual(rows[0]['sequence_id'],'7');self.assertEqual(rows[0]['Sequence'],'VHLT')

    def test_embedded_listing_after_large_bibliography(self):
        data=('<us-patent-grant><description>'+(' ' * 12000)+'<sequence-cwu><s200><s210>1</s210><s211>4</s211><s212>DNA</s212></s200><s400>1\nacgt 4</s400></sequence-cwu></description></us-patent-grant>').encode()
        row=parse_sequences(data,'grant.xml','US1B2','archive')[0]
        self.assertEqual(row['Sequence'],'ACGT');self.assertEqual(row['sequence_status'],'parsed')

    def test_target_sequence_is_not_assigned_to_compound(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'test.zip'
            xml=b'<us-patent-grant><description><p id="P1"><table><tgroup cols="4"><thead><row><entry>Compound</entry><entry>Name</entry><entry>Target</entry><entry>Structure</entry></row></thead><tbody><row><entry>1</entry><entry>ethanol</entry><entry>SEQ ID NO: 1</entry><entry><img file="a-C1.TIF"/></entry></row></tbody></tgroup></table></p><p id="P2">UniProt accession P12345 is mentioned here.</p></description></us-patent-grant>'
            seq=b'<sequence-cwu><s200><s210>1</s210><s211>4</s211><s212>PRT</s212></s200><s400>1\nVal His Leu Thr\n1</s400></sequence-cwu>'
            with zipfile.ZipFile(path,'w') as z:z.writestr('grant.xml',xml);z.writestr('SEQLST.XML',seq)
            with zipfile.ZipFile(path) as z:data=extract(z,'US1B2','archive','grant.xml',['a-C1.MOL'],{})
            self.assertEqual(data['structures']['a-C1.MOL'][0]['Sequence'],'')
            self.assertEqual(data['uniprot_mentions'][0]['UniProt ID'],'P12345')
            self.assertEqual(data['structures']['a-C1.MOL'][0]['UniProt ID'],'')

    def test_all_sources_survive_bad_xml_and_bad_mol(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'raw';out=Path(tmp)/'results';root.mkdir();out.mkdir();archive='US1-20250101.ZIP'
            single=Chem.MolToMolBlock(Chem.MolFromSmiles('CCO')).encode()
            salt=Chem.MolToMolBlock(Chem.MolFromSmiles('CC(=O)[O-].[Na+]')).encode()
            generic=Chem.MolToMolBlock(Chem.MolFromSmiles('*CC')).encode()
            entries=[('US1-20250101.XML',b'<PATDOC><bad>&UnknownChar;</bad>'),('same.MOL',single),('same.MOL',salt),('generic.MOL',generic),('invalid.MOL',b'broken MOL'),('extra.CDX',b'CDXDATA'),('patent-C00004.TIF',b'RAWIMAGE'),('patent-D00001.TIF',b'GENERALDRAWING')]
            with zipfile.ZipFile(root/archive,'w') as z:
                for name,data in entries:z.writestr(name,data)
            result=worker((str(root),str(out),str(Path(tmp)/'prior'),{'archive':archive,'patent_core':'US1'},'test-code'))
            rows=pq.read_table(result['part']).to_pylist();stats=json.loads(Path(result['checkpoint']).read_text())['stats']
            self.assertEqual(len(rows),6);self.assertEqual(stats['mol_files'],4);self.assertEqual(stats['source_read_errors'],0)
            self.assertEqual(len({r['source_structure_id'] for r in rows}),6)
            self.assertIn('concrete_multifragment',{r['structure_status'] for r in rows})
            self.assertIn('markush_or_query',{r['structure_status'] for r in rows})
            self.assertIn('invalid_mol',{r['structure_status'] for r in rows})
            for row in rows:self.assertEqual(hashlib.sha256(row['raw_structure_bytes']).hexdigest(),row['source_sha256'])
            self.assertFalse(any(r['has_verified_SAR'] for r in rows))
            cached=worker((str(root),str(out),str(Path(tmp)/'prior'),{'archive':archive,'patent_core':'US1'},'test-code'));self.assertTrue(cached['cached'])

    def test_material_structure_and_numbered_sequence_link(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'raw';out=Path(tmp)/'results';root.mkdir();out.mkdir();archive='US2-20250101.ZIP'
            xml=b'<us-patent-grant><us-bibliographic-data-grant><publication-reference><document-id><country>US</country><doc-number>0000002</doc-number><kind>B2</kind><date>20250101</date></document-id></publication-reference><invention-title>Polymer compositions</invention-title></us-bibliographic-data-grant><description><p id="P1"><table><tgroup cols="3"><thead><row><entry>Compound</entry><entry>IUPAC name</entry><entry>Structure</entry></row></thead><tbody><row><entry>Example 03</entry><entry>ethanol SEQ ID NO: 1</entry><entry><img file="a-C00001.TIF"/></entry></row></tbody></tgroup></table></p></description></us-patent-grant>'
            seq=b'<sequence-cwu><sequence-list-new-rules><s200><s210>1</s210><s211>4</s211><s212>PRT</s212></s200><s400>1\nVal His Leu Thr\n1</s400></sequence-list-new-rules></sequence-cwu>'
            with zipfile.ZipFile(root/archive,'w') as z:
                z.writestr('US2-20250101.XML',xml);z.writestr('a-C00001.MOL',Chem.MolToMolBlock(Chem.MolFromSmiles('CCO')));z.writestr('SEQLST.XML',seq)
            result=worker((str(root),str(out),str(Path(tmp)/'prior'),{'archive':archive,'patent_core':'US2'},'test-code'))
            row=pq.read_table(result['part']).to_pylist()[0];stats=json.loads(Path(result['checkpoint']).read_text())['stats']
            self.assertEqual(stats['screen'],'materials');self.assertEqual(stats['structure_rows'],1)
            self.assertEqual(row['Patent'],'US2B2');self.assertEqual(row['original_alias'],'Example 03')
            self.assertEqual(row['Sequence'],'VHLT');self.assertEqual(row['sequence_link_status'],'explicit_structure_row_sequence_id')

if __name__=='__main__':unittest.main()
