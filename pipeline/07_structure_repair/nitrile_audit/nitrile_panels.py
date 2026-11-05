_ws = (__import__("os").environ.get("PATBENCH_WORKSPACE") or __import__("sys").exit("PATBENCH_WORKSPACE is not set (see docs/REPRODUCE.md)")).rstrip("/")
import sys,pathlib,json
sys.path.insert(0,_ws + '/outputs/uspto-sar-demo-20261003/runtime/deps')
import pandas as pd
from PIL import Image,ImageOps,ImageDraw,ImageFont
r=pathlib.Path(_ws + '/outputs/nitrile-evidence-audit-20261031');d=pd.read_parquet(r/'motif_evidence.parquet');
# test, all non-standard labels, and one non-test N^C example per remaining patent
z=d[(d.test_tasks!='')|(d.exact_superatom_labels!='N^C')].copy();extra=d[(d.test_tasks=='')&(d.exact_superatom_labels=='N^C')&~d.Patent.isin(z.Patent)].groupby('Patent',sort=True).head(1);z=pd.concat([z,extra]).sort_values(['exact_superatom_labels','Patent','structure_file','motif_index']).reset_index(drop=True);z['review_id']=z.index+1
out=r/'panels';out.mkdir(exist_ok=True);font=ImageFont.truetype('DejaVuSans.ttf',13)
for i,t in z.iterrows():
 im=Image.open(r/'originals'/t.image_file).convert('RGB');w,h=im.size
 # original image is an overall drawing, MOL coordinate span gives approximate site position
 px=(t.x-t.minx)/max(t.maxx-t.minx,0.1)*w;py=(t.maxy-t.y)/max(t.maxy-t.miny,0.1)*h
 crop=im.crop((max(0,px-180),max(0,py-125),min(w,px+180),min(h,py+125)))
 crop.thumbnail((360,250));cell=Image.new('RGB',(380,305),'white');cell.paste(crop,((380-crop.width)//2,45));dr=ImageDraw.Draw(cell);dr.text((5,3),f"{t.review_id} {t.Patent} {pathlib.Path(t.structure_file).stem.split('-')[-1]}",font=font,fill='black');dr.text((5,20),f"alias={t.compound_alias[:22]} label={t.exact_superatom_labels} {'TEST' if t.test_tasks else ''}",font=font,fill='black');dr.text((5,287),f"full image {w}x{h}; motif {t.motif_index}",font=font,fill='black')
 if i%20==0:sheet=Image.new('RGB',(1520,1525),'#eeeeee')
 sheet.paste(cell,((i%4)*380,((i%20)//4)*305))
 if i%20==19 or i==len(z)-1:sheet.save(out/f'panel_{i//20+1:03}.jpg',quality=92)
z.to_parquet(r/'visual_review_index.parquet',index=False);print('visual cells',len(z),'panels',(len(z)+19)//20,'test',sum(z.test_tasks.ne('')),'patents',z.Patent.nunique())
