_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import sys,pathlib,zipfile,hashlib,json,re,io
sys.path.insert(0,_ws + '/outputs/uspto-sar-demo-20261003/runtime/deps')
import pandas as pd
from rdkit import Chem,RDLogger
from PIL import Image
RDLogger.DisableLog('rdApp.*')
o=pathlib.Path(_ws + '/outputs/nitrile-evidence-audit-20261031');d=pd.read_parquet(o/'inventory.parquet');raw=pathlib.Path(_ws + '/raw/patent-data');pattern=Chem.MolFromSmarts('[#6][N+]#[C-]');rows=[];errors=[]
for arch,g in d.groupby('archive'):
 try:
  with zipfile.ZipFile(raw/arch) as z:
   for t in g.itertuples():
    block=z.read(t.structure_file);m=Chem.MolFromMolBlock(block.decode('latin1'),sanitize=False,removeHs=False)
    if m is None:errors.append([t.structure_file,'mol_parse']);continue
    groups=Chem.GetMolSubstanceGroups(m);matches=m.GetSubstructMatches(pattern);sgs=[]
    for sg in groups:
     p=sg.GetPropsAsDict();sgs.append({'atoms':list(sg.GetAtoms()),'label':str(p.get('LABEL','')),'type':str(p.get('TYPE',''))})
    imagepath=t.structure_file[:-4]+'.TIF';cdxpath=t.structure_file[:-4]+'.CDX'
    imagebytes=z.read(imagepath) if imagepath in z.namelist() else b''
    for j,mt in enumerate(matches):
     pair=set(mt[1:]);cover=[sg for sg in sgs if pair.issubset(set(sg['atoms']))];exact=[sg['label'] for sg in cover if set(sg['atoms'])==pair];conf=m.GetConformer();v=[conf.GetAtomPosition(i) for i in range(m.GetNumAtoms())];pt=conf.GetAtomPosition(mt[1]);row={'Patent':t.Patent,'structure_file':t.structure_file,'compound_alias':str(t.compound_alias),'n_rows':t.n,'test_tasks':t.test_tasks,'motif_index':j,'atom_attachment':mt[0]+1,'atom_n_before':mt[1]+1,'atom_c_before':mt[2]+1,'exact_superatom_labels':';'.join(exact),'covering_groups':json.dumps(cover),'mol_sha256':hashlib.sha256(block).hexdigest(),'image_sha256':hashlib.sha256(imagebytes).hexdigest() if imagebytes else '', 'image_file':imagepath,'x':pt.x,'y':pt.y,'minx':min(a.x for a in v),'maxx':max(a.x for a in v),'miny':min(a.y for a in v),'maxy':max(a.y for a in v)}
     rows.append(row)
    # save all original files for reproducibility; no alterations
    for fname,b in [(t.structure_file,block),(imagepath,imagebytes)]:
     if b:
      p=o/'originals'/fname;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b)
 except Exception as e:errors.append([str(arch),str(e)])
print('errors',errors,flush=True)
pd.DataFrame(rows).to_parquet(o/'motif_evidence.parquet',index=False)
print('motifs',len(rows),flush=True);q=pd.DataFrame(rows);print(q.exact_superatom_labels.value_counts(dropna=False).to_string(),flush=True);print('test',q[q.test_tasks.ne('')].exact_superatom_labels.value_counts(dropna=False).to_string(),flush=True)
(o/'extraction_errors.json').write_text(json.dumps(errors))
