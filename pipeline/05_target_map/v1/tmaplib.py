"""Shared normalization helpers for the Target -> UniProt mapping."""
import re, unicodedata
GREEK={'α':'alpha','β':'beta','γ':'gamma','δ':'delta','ε':'epsilon','ϵ':'epsilon','ζ':'zeta','η':'eta','θ':'theta','ι':'iota',
 'κ':'kappa','λ':'lambda','μ':'mu','µ':'mu','ν':'nu','ξ':'xi','π':'pi','ρ':'rho','σ':'sigma','τ':'tau','φ':'phi','χ':'chi','ψ':'psi','ω':'omega',
 'Α':'alpha','Β':'beta','Γ':'gamma','Δ':'delta','Ε':'epsilon','Ζ':'zeta','Κ':'kappa','Λ':'lambda','Σ':'sigma','Ω':'omega'}
def clean_text(s):
    s=unicodedata.normalize('NFKC',s or '')
    s=re.sub(r'[⁡-⁤​­]','',s)
    s=s.replace('–','-').replace('—','-').replace('‐','-').replace('−','-')
    return re.sub(r'\s+',' ',s).strip()
def greek(s):
    return ''.join(GREEK.get(ch,ch) for ch in s)
def key(s):
    """Aggressive lookup key: NFKC, greek->latin, lowercase, alphanumerics only."""
    return re.sub(r'[^a-z0-9]','',greek(clean_text(s)).lower())
NOISE=r'(assays?|binding|bdg|inhibition|inhibitory|inhibitor|activity|activities|enzymes?|enzymatic|biochemical|biochem|kinase assay|kinetic|htrf|tr-?fret|fret|alphascreen|alphalisa|lanthascreen|z-?lyte|adp-?glo|caliper|spa|elisa|fp|ic50|ic-50|ki|kd|ec50|avg|average|mean|median|stddev|std dev|sd|nm|um|μm|in vitro|protein|functional|affinity|potency|test|data|results?|value|values)'
def strip_noise(s):
    s=greek(clean_text(s))
    s=re.sub(r'\b'+NOISE+r'\b',' ',s,flags=re.I)
    s=re.sub(r'[,;:*()\[\]]',' ',s)
    return re.sub(r'\s+',' ',s).strip()
ORG_WORDS=[(r'\b(human|hu|homo sapiens)\b','Homo sapiens'),(r'\b(rat|rattus)\b','Rattus norvegicus'),(r'\b(mouse|murine|mus musculus)\b','Mus musculus'),
 (r'\b(dog|canine)\b','Canis lupus familiaris'),(r'\b(cyno|cynomolgus|monkey|macaque|rhesus)\b','Macaca'),(r'\b(guinea pig)\b','Cavia porcellus'),
 (r'\b(pig|porcine)\b','Sus scrofa'),(r'\b(bovine)\b','Bos taurus'),(r'\b(rabbit)\b','Oryctolagus cuniculus')]
def stated_organism(s):
    s2=clean_text(s)
    for pat,org in ORG_WORDS:
        if re.search(pat,s2,re.I): return org
    m=re.match(r'^([hrm])(?=[A-Z0-9])',s2)
    if m: return {'h':'Homo sapiens','r':'Rattus norvegicus','m':'Mus musculus'}[m.group(1)]+' (prefix)'
    return ''
MUT=re.compile(r'(?<![A-Za-z0-9])([ACDEFGHIKLMNPQRSTVWY])(\d{1,4})([ACDEFGHIKLMNPQRSTVWY]|del|ins|fs|\*)(?![a-z0-9])')
def mutations(s):
    s2=clean_text(s)
    out=[''.join(m.groups()) for m in MUT.finditer(s2)]
    if re.search(r'(del ?f ?508|ΔF508|deltaF508|F508del)',s2,re.I): out.append('F508del')
    if re.search(r'exon ?20|exon ?19|insertion|deletion|mutant|mutated|\bmut\b',s2,re.I) and not out: out.append('unspecified_mutant')
    return out
