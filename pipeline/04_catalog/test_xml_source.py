import unittest
import xml.etree.ElementTree as ET
from xml_source import parse_xml,select_patent_xml
import io
import zipfile

class SourceXmlTests(unittest.TestCase):
    def test_greek_target_and_unit_preserved(self):
        root=parse_xml(b'<!DOCTYPE patent SYSTEM "not-fetched.dtd"><patent>PI3K&agr; IC50 0.5 &mgr;M &ldquo;test&rdquo;</patent>')
        self.assertEqual(root.text,'PI3Kα IC50 0.5 μM “test”')
    def test_encoding_escaped_markup_and_literal_cdata(self):
        root=parse_xml(b'<?xml version="1.0" encoding="ISO-8859-1"?><patent>&ouml; &lt;b&gt;<![CDATA[&mgr;]]></patent>')
        self.assertEqual(root.text,'ö <b>&mgr;')
    def test_unknown_entity_is_not_guessed(self):
        with self.assertRaises(ET.ParseError):parse_xml(b'<patent>&UNKNOWN_STRUCTURE;</patent>')
    def test_invalid_xml_still_fails(self):
        with self.assertRaises(ET.ParseError):parse_xml(b'<patent><row></patent>')
    def test_sequence_xml_does_not_replace_patent_grant(self):
        data=io.BytesIO()
        with zipfile.ZipFile(data,'w') as z:
            z.writestr('US123-20250101/US123-20250101-S00001.XML','<sequence-cwu/>')
            z.writestr('US123-20250101/US123-20250101.XML','<us-patent-grant/>')
        with zipfile.ZipFile(data) as z:self.assertEqual(select_patent_xml(z,'US123-20250101.ZIP'),'US123-20250101/US123-20250101.XML')
    def test_multiple_unidentified_xmls_are_rejected(self):
        data=io.BytesIO()
        with zipfile.ZipFile(data,'w') as z:
            z.writestr('a.xml','<attachment/>');z.writestr('b.xml','<attachment/>')
        with zipfile.ZipFile(data) as z:
            with self.assertRaises(ValueError):select_patent_xml(z,'US123.ZIP')

if __name__=='__main__':unittest.main(verbosity=2)
