"""Resolve a curated gene symbol (or ORG:species|name) + organism to a UniProt accession and sequence.
Mammalian symbols: UniProt reviewed human/mouse/rat reference, requested organism first then human.
Non-mammalian / ORG: entries: ChEMBL component_sequences name+organism match, else unresolved (needs fetch)."""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pandas as pd, re, sys
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from tmaplib import key, clean_text
W=_ws + '/outputs/target-uniprot-map-20261016/work/'
TAX={'human':'9606','homo sapiens':'9606','mouse':'10090','murine':'10090','mus musculus':'10090',
 'rat':'10116','rattus norvegicus':'10116','assumed human':'9606','':'9606'}
class Resolver:
    def __init__(self):
        U=pd.read_parquet(W+'uniprot_sprot_hmr.parquet')
        self.U=U
        self.by={}
        for r in U.itertuples():
            sym=r.gene_primary.upper()
            if sym: self.by.setdefault((sym,r.tax),[]).append(r)
            for g in r.gene_names.split():
                self.by.setdefault((g.upper(),r.tax),[]).append(r)
        cs=pd.read_parquet(W+'chembl_comp_seq.parquet')
        self.cs=cs[cs.sequence.notna()&(cs.db_source.isin(['SWISS-PROT','TREMBL']))].copy()
        self.cs['namekey']=self.cs.description.fillna('').map(key)
        self.syn=pd.read_parquet(W+'chembl_comp_synonyms.parquet')
        self.syn['namekey']=self.syn.component_synonym.map(key)
        self.seqmap=dict(zip(self.cs.accession,self.cs.sequence))
        self.orgmap=dict(zip(self.cs.accession,self.cs.organism))
        self.descmap=dict(zip(self.cs.accession,self.cs.description))
    def mammal(self,sym,organism):
        tax=TAX.get((organism or '').strip().lower(),'')
        order=[tax] if tax else []
        for t in ['9606','10090','10116']:
            if t not in order: order.append(t)
        for t in order:
            hits=self.by.get((sym.upper(),t))
            if hits:
                exact=[h for h in hits if h.gene_primary.upper()==sym.upper()] or hits
                exact=sorted(exact,key=lambda h:-int(h.length))
                h=exact[0]
                return dict(accession=h.accession,entry=h.entry,protein_name=h.protein_names,organism=h.organism,
                    tax_id=t,sequence=h.sequence,seq_len=int(h.length),seq_source='UniProt_sprot_'+t,
                    matched_symbol=h.gene_primary.upper(),organism_used=('requested' if t==tax else 'human_fallback' if t=='9606' else 'other'))
        return None
    def nonmammal(self,species,name):
        spk=key(species); nk=key(name)
        c=self.cs[self.cs.organism.fillna('').map(key).str.contains(spk,na=False)] if spk else self.cs
        if len(c):
            ex=c[c.namekey==nk]
            if len(ex)==0:
                m=self.syn[(self.syn.namekey==nk)&self.syn.accession.isin(set(c.accession))]
                if len(m): ex=c[c.accession.isin(set(m.accession))]
            if len(ex)==0 and len(nk)>=4:
                ex=c[c.namekey.str.contains(re.escape(nk),na=False)]
            if len(ex):
                r=ex.iloc[0]
                return dict(accession=r.accession,entry='',protein_name=r.description or '',organism=r.organism or species,
                   tax_id=str(r.tax_id),sequence=r.sequence,seq_len=len(r.sequence),seq_source='ChEMBL37_component_sequences',
                   matched_symbol=name,organism_used='chembl_match')
        return None
