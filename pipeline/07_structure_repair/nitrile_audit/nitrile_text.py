_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import sys,pathlib,zipfile,re,json,html
sys.path.insert(0,_ws + '/outputs/uspto-sar-demo-20261003/runtime/deps')
import pandas as pd
r=pathlib.Path(_ws + '/outputs/nitrile-evidence-audit-20261031');d=pd.read_parquet(r/'motif_evidence_geometry.parquet');inv=pd.read_parquet(r/'inventory.parquet');arch=inv.set_index('structure_file').archive.to_dict();raw=pathlib.Path(_ws + '/raw/patent-data');texts={};rows=[]
for t in d.itertuples():
 a=arch[t.structure_file]
 if a not in texts:
  with zipfile.ZipFile(raw/a) as z:
   fname=[x for x in z.namelist() if x.lower().endswith('.xml') or x.lower().endswith('.sgm')][0];texts[a]=z.read(fname).decode('utf8',errors='replace')
 text=texts[a];stem=pathlib.Path(t.image_file).name;match=re.search(re.escape(stem),text,re.I)
 if not match:rows.append({'structure_file':t.structure_file,'context':'','near_heading':'','chem_name_flags':''});continue
 before=text[:match.start()];headings=re.findall(r'<heading\b[^>]*>(.*?)</heading>',before,re.S|re.I);h=re.sub('<[^>]*>','',headings[-1]) if headings else ''
 context=html.unescape(re.sub('<[^>]*>',' ',text[max(0,match.start()-1200):match.end()+2500]));context=re.sub(r'\s+',' ',context)
 iso=bool(re.search(r'isocyan(?:id|o)(?!ate)|isonitril',context,re.I));rows.append({'structure_file':t.structure_file,'context':context,'near_heading':h,'chem_name_flags':'isocyanide_text' if iso else ''})
o=pd.DataFrame(rows).drop_duplicates('structure_file');o.to_parquet(r/'text_context.parquet',index=False)
q=d[(d.exact_superatom_labels.isin(['','CN','^NC']))|((d.exact_superatom_labels=='NC')&(d.external_dx<-.15))].merge(o,on='structure_file');q.to_json(r/'exception_context.json',orient='records',force_ascii=False,indent=2)
for t in q.itertuples():print(t.Patent,t.structure_file,t.near_heading,t.context[:1500],flush=True)
print('ISO_FLAGS',o[o.chem_name_flags!=''][['structure_file','near_heading']].to_dict('records'),flush=True)
