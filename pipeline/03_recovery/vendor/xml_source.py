"""Parse patent XML offline, preserving standard ISO/HTML named characters.

Greek ISO 8879 names and Unicode values follow the W3C entity declarations:
https://www.w3.org/2003/entities/2007/isogrk1.ent
Other standard names use Python's HTML5 entity table. Unknown names still fail.
"""
from html.entities import html5
import re
import xml.etree.ElementTree as ET
from pathlib import PurePosixPath

_greek={'agr':'α','bgr':'β','dgr':'δ','eegr':'η','egr':'ε','ggr':'γ','igr':'ι','kgr':'κ','khgr':'χ','lgr':'λ','mgr':'μ','ngr':'ν','ogr':'ο','ohgr':'ω','pgr':'π','phgr':'φ','psgr':'ψ','rgr':'ρ','sfgr':'ς','sgr':'σ','tgr':'τ','thgr':'θ','ugr':'υ','xgr':'ξ','zgr':'ζ'}
_entities={k.rstrip(';'):v for k,v in html5.items() if k.endswith(';')}
_entities.update(_greek)
_entities.update({k[:-2].upper()+'gr':v.upper() for k,v in _greek.items() if k!='sfgr'})
_reserved={'amp','lt','gt','quot','apos'}
_protected=re.compile(br'(<!\[CDATA\[.*?\]\]>|<!--.*?-->)',re.S)
_named=re.compile(br'&([A-Za-z][A-Za-z0-9]*);')

def select_patent_xml(z,archive):
    names=[n for n in z.namelist() if n.lower().endswith('.xml')]
    if not names:return None
    stem=PurePosixPath(str(archive)).stem.lower()
    exact=[n for n in names if PurePosixPath(n).stem.lower()==stem]
    if len(exact)==1:return exact[0]
    grants=[]
    for n in names:
        with z.open(n) as f:prefix=f.read(2048).lower()
        if b'<us-patent-grant' in prefix or b'<patdoc' in prefix:grants.append(n)
    if len(grants)==1:return grants[0]
    if len(names)==1:return names[0]
    raise ValueError('Cannot identify a unique patent XML; sequence attachments must not substitute for the grant.')

def parse_xml(data):
    try:return ET.fromstring(data)
    except ET.ParseError as exc:
        if 'undefined entity' not in str(exc):raise
    def replace(match):
        name=match[1].decode('ascii')
        if name in _reserved or name not in _entities:return match[0]
        return ''.join('&#'+str(ord(c))+';' for c in _entities[name]).encode('ascii')
    parts=_protected.split(data)
    rewritten=b''.join(part if i%2 else _named.sub(replace,part) for i,part in enumerate(parts))
    return ET.fromstring(rewritten)
