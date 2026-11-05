"""Robust reader for the v2 curation TSVs (10 columns: v1's 9 plus `mappable`).
Repairs lines with a stray/missing tab by anchoring on conf in {H,M,L} and mappable in {Y,N,?}."""
import glob, os, re, pandas as pd
COLS=['item_id','class','cell','symbols','implied','organism','mutation','conf','note','mappable']
CLS={'protein','complex','mutant','phenotypic','unclear'}
HML={'H','M','L'}; MAP={'Y','N','?'}
def read(d):
    rows=[]; bad=[]
    for f in sorted(glob.glob(os.path.join(d,'out_*.tsv'))):
        for ln,line in enumerate(open(f,encoding='utf-8')):
            line=line.rstrip('\n')
            if not line.strip() or line.startswith('item_id\t'): continue
            p=line.split('\t')
            if not re.match(r'^[qt][\w]+$',p[0].strip()): bad.append((f,ln,'id',line[:120])); continue
            ok = len(p)==10 and p[7].strip().upper() in HML and p[9].strip().upper() in MAP
            if not ok:
                cls=p[1].strip().lower() if len(p)>1 else ''
                if cls not in CLS: bad.append((f,ln,'class',line[:160])); continue
                # mappable: last field if it is Y/N/?
                mp=''
                if p and p[-1].strip().upper() in MAP: mp=p[-1].strip().upper(); p=p[:-1]
                idx=[i for i,v in enumerate(p[1:],1) if v.strip().upper() in HML]
                if not idx: bad.append((f,ln,'conf',line[:160])); continue
                ci=idx[-1]
                head=(p[:ci]+['']*(7-ci))[:7] if ci<7 else p[:7]
                note=' '.join(x for x in p[ci+1:] if x)
                p=head+[p[ci],note,mp]
            p=[x.strip() for x in (p+['']*10)[:10]]
            rows.append(dict(zip(COLS,p))|{'batch':os.path.basename(f)})
    return pd.DataFrame(rows), bad
if __name__=='__main__':
    import sys
    C,bad=read(sys.argv[1])
    print('rows',len(C),'unique',C.item_id.nunique(),'unrepairable',len(bad))
    for b in bad[:20]: print('BAD',b)
    print(C['class'].value_counts().to_dict()); print('conf',C.conf.value_counts().to_dict()); print('mappable',C.mappable.value_counts().to_dict())
    print('protein-class with empty symbols',((C['class'].isin(['protein','complex','mutant']))&(C.symbols=='')).sum())
    print('phenotypic marked Y but no implied',((C['class']=='phenotypic')&(C.mappable=='Y')&(C.implied=='')).sum())
