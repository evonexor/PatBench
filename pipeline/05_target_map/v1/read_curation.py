"""Robust reader for the curation TSVs: repairs lines with a missing/extra tab by anchoring on the
confidence field (H/M/L) and the known class vocabulary. Reports anything it cannot repair."""
import glob, os, re, pandas as pd
COLS=['item_id','class','cell','symbols','implied','organism','mutation','conf','note']
CLS={'protein','complex','mutant','phenotypic','unclear'}
def read(d):
    rows=[]; bad=[]
    for f in sorted(glob.glob(os.path.join(d,'out_*.tsv'))):
        for ln,line in enumerate(open(f,encoding='utf-8')):
            line=line.rstrip('\n')
            if not line.strip() or line.startswith('item_id\t'): continue
            p=line.split('\t')
            if not re.match(r'^[ps][\w]+$',p[0].strip()): bad.append((f,ln,'id',line[:120])); continue
            if len(p)==9 and p[7].strip().upper() in ('H','M','L',''): pass
            else:
                idx=[i for i,v in enumerate(p[1:],1) if v.strip().upper() in ('H','M','L')]
                ci=idx[-1] if idx else None
                if ci is None or p[1].strip().lower() not in CLS: bad.append((f,ln,'shape',line[:160])); continue
                head=p[:ci]+['']*(7-ci) if ci<7 else p[:7]
                p=head[:7]+[p[ci]]+[' '.join(x for x in p[ci+1:] if x)]
            p=(p+['']*9)[:9]
            rows.append(dict(zip(COLS,[x.strip() for x in p]))|{'batch':os.path.basename(f)})
    C=pd.DataFrame(rows)
    return C, bad
if __name__=='__main__':
    import sys
    C,bad=read(sys.argv[1])
    print('rows',len(C),'unique',C.item_id.nunique(),'unrepairable',len(bad))
    for b in bad[:20]: print('BAD',b)
    print(C['class'].value_counts().to_dict()); print(C.conf.value_counts().to_dict())
    print('empty symbols among protein-class',((C['class'].isin(['protein','complex','mutant']))&(C.symbols=='')).sum())
