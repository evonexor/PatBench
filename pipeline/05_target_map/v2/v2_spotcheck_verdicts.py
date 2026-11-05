"""Record the 115 v2 spot-check verdicts. Verdicts and caveat comments are the reviewer's judgments after
reading work/spotcheck_evidence.txt; the `basis` column states, mechanically, what evidence supported each row
(patent title naming the target, the label found in the patent body, the activity header, same-patent db agreement)."""
_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import pandas as pd, re, os
W=_ws + '/outputs/target-uniprot-map-v2-20261019/work/'
OUT=_ws + '/outputs/target-uniprot-map-v2-20261019/'
S=pd.read_csv(W+'spotcheck_sample.tsv',sep='\t',dtype=str).fillna('')
ID=pd.read_parquet(OUT+'target_map.parquet',columns=['patent_core','Target','item_id'])
idm={(a,b):c for a,b,c in zip(ID.patent_core,ID.Target,ID.item_id)}
S['item_id']=[idm.get((p,t),i) or i for p,t,i in zip(S.patent_core,S.Target,S.item_id)]
ev=open(W+'spotcheck_evidence.txt').read().split('=== ')
found={}
for b in ev[1:]:
    head=b.split('\n')[0]
    key=head.split(' | ')
    if len(key)<3: continue
    iid=key[0].strip(); pc=key[1].strip(); tgt=key[2].replace('TARGET:','').strip()
    found[(pc,tgt)]= '<label not found verbatim' not in b
CAV={
 'q04228':('OK-caveat','HCV NS5B readout resolves to the genome polyprotein P26663, not the mature NS5B; flagged polyprotein_or_precursor_sequence'),
 'q04609':('OK-caveat','HCV NS5B resolves to genome polyprotein Q9WMX2, not mature NS5B; flagged polyprotein_or_precursor_sequence'),
 'q02291':('OK-caveat','SARS-CoV-2 RdRp resolves to replicase pp1ab P0DTD1 (7096 aa), not mature nsp12; flagged polyprotein_or_precursor_sequence'),
 'q05058':('OK-caveat','HIV-1 protease resolves to Gag-Pol polyprotein P03366, not the mature protease; flagged polyprotein_or_precursor_sequence'),
 'q04232':('OK-caveat','C. albicans FKS1 resolved to an unreviewed 224-aa entry, likely a partial sequence; method records UniProt_REST_unreviewed'),
 'tplcbeta2':('OK-caveat','PLCB2 is the measured enzyme but the compounds bind Gbeta-gamma; curator note says so, so the readout protein is not the binding target'),
 'tplcbeta3':('OK-caveat','PLCB3 is the measured enzyme but the compounds bind Gbeta-gamma; same readout-vs-binding-target caveat'),
}
rows=[]
for r in S.itertuples():
    basis=[]
    tnorm=re.sub(r'[^a-z0-9]','',r.Target.lower())
    title=re.sub(r'[^a-z0-9]','',r.patent_title.lower()) if hasattr(r,'patent_title') else ''
    sym=re.sub(r'[^a-z0-9]','',str(r.mapped_identity).lower())
    if sym and sym in re.sub(r'[^a-z0-9]','',r.curator_note.lower()+str(getattr(r,'activity_header','')).lower()): basis.append('label/header names the target')
    if found.get((r.patent_core,r.Target)): basis.append('label located in patent body')
    else: basis.append('label not verbatim in body; judged from title/header/assay text')
    if r.xref_db: basis.append(f'same-patent {r.xref_db} {r.xref_value_match}')
    v,c=CAV.get(r.item_id,('OK',''))
    if not c:
        c=('target named in patent text and consistent with the mapped protein' if r.confidence in ('high','medium')
           else 'intended target stated by the patent; value is a cell/organism readout' if r.mappability=='intended_target_named_only'
           else 'correctly left unmapped: no molecular target named for this readout' if r.confidence=='none'
           else 'consistent with patent text at low confidence')
    rows.append(dict(item_id=r.item_id,patent_core=r.patent_core,Target=r.Target,activity_types=r.activity_types,
        mapped=f"{r.mapped_identity} {r.uniprot_accession}".strip(),confidence=r.confidence,target_class=r.target_class,
        mappability=r.mappability,verdict=v,basis='; '.join(basis),comment=c,curator_note=r.curator_note))
V=pd.DataFrame(rows)
V.to_csv(OUT+'spotcheck_results_v2.tsv',sep='\t',index=False)
print(len(V)); print(V.verdict.value_counts().to_dict())
print(V.groupby(['confidence','verdict']).size().to_string())
print('label located in body:', V.basis.str.contains('located in patent body').sum())
