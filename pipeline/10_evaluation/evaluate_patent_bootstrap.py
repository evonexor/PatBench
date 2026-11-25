"""Within-assay ranking evaluation, with patent-cluster intervals and paired comparisons.
python evaluate.py TASK PRED [--reference BASELINE] [--profile PROFILE] [--json-out REPORT]
Higher score means more potent. Groups remain separate; CI resamples whole patents without
averaging different assays. Few-shot repeats are averaged only within the same assay group.
"""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import kendalltau,spearmanr,rankdata
PAIR_MARGIN=.5;NDCG_K=10;TOP_FRACTION=.1;BOOTSTRAP=2000
METRICS=['spearman','kendall','pair_accuracy','ndcg10','top10pct_hit','partial_spearman_heavy_atoms']

def _ndcg(truth,score,k=NDCG_K):
    span=np.ptp(truth)
    if not span:return np.nan
    gain=(truth-truth.min())/span;n=len(truth);k=min(k,n)
    discounts=np.zeros(n);discounts[:k]=1/np.log2(np.arange(2,k+2))
    ideal=np.dot(np.sort(gain)[::-1],discounts)
    order=np.argsort(-score,kind='stable');ss=score[order];gg=gain[order];dcg=0.
    starts=np.r_[0,np.flatnonzero(np.diff(ss)!=0)+1];ends=np.r_[starts[1:],n]
    for a,b in zip(starts,ends):
        if a>=k:break
        dcg+=gg[a:b].mean()*discounts[a:b].sum()
    return float(dcg/ideal) if ideal else np.nan

def _pair_accuracy(truth,score,margin=PAIR_MARGIN):
    dt=truth[:,None]-truth[None,:];ds=score[:,None]-score[None,:];use=np.triu(abs(dt)>=margin,k=1);n=int(use.sum())
    return (float((((np.sign(dt)==np.sign(ds))&use).sum()+.5*((ds==0)&use).sum())/n),n) if n else (np.nan,0)

def _top_hit(truth,score,frac=TOP_FRACTION):
    k=max(1,int(round(len(truth)*frac)))
    def weights(v):
        threshold=np.sort(v)[-k];above=v>threshold;equal=v==threshold
        return above.astype(float)+equal*(k-above.sum())/equal.sum()
    return float(np.dot(weights(truth),weights(score))/k)

def _partial(truth,score,size):
    if size is None:return np.nan
    x=np.c_[np.ones(len(size)),rankdata(size)];a=rankdata(truth);b=rankdata(score)
    a=a-x@np.linalg.lstsq(x,a,rcond=None)[0];b=b-x@np.linalg.lstsq(x,b,rcond=None)[0]
    if np.std(a)<1e-10:return np.nan
    return float(np.corrcoef(a,b)[0,1]) if np.std(b)>1e-10 else 0.

def group_metrics(truth,score,size=None):
    truth=np.asarray(truth,float);score=np.asarray(score,float)
    if len(truth)<3 or np.ptp(truth)==0:return None
    acc,n=_pair_accuracy(truth,score)
    return dict(n_molecules=len(truth),spearman=float(spearmanr(score,truth).statistic) if np.ptp(score)>0 else 0.,kendall=float(kendalltau(score,truth).statistic) if np.ptp(score)>0 else 0.,pair_accuracy=acc,n_pairs=n,ndcg10=_ndcg(truth,score),top10pct_hit=_top_hit(truth,score),partial_spearman_heavy_atoms=_partial(truth,score,size))

def _check_pred(pred,keys):
    if pred.duplicated(keys).any():raise ValueError('Duplicate prediction keys: '+str(keys))
    vals=pd.to_numeric(pred.score,errors='coerce')
    if not np.isfinite(vals).all():raise ValueError('Prediction scores must be finite numbers')

def collapse_stereo(d):
    unit=np.where(d.stereo_set_id.fillna('').ne(''),d.stereo_set_id,d.group_id.astype(str)+':'+d.InChIKey)
    agg={c:'first' for c in d.columns if c not in ['pAct','score','InChIKey']};agg.update(pAct='median',score='median')
    return d.assign(_u=unit).groupby('_u',sort=False).agg(agg).reset_index(drop=True)

def _group_row(g,gid):
    m=group_metrics(g.pAct.to_numpy(float),g.score.to_numpy(float),g.heavy_atoms.to_numpy(float) if 'heavy_atoms' in g else None)
    if m is None:return None
    m.update(group_id=gid,patent_core=g.patent_core.iloc[0],uniprot_accession=g.uniprot_accession.iloc[0],endpoint=g['Activity Type'].iloc[0])
    for c in ['test_tier','target_organism','is_human_target','priority_date']:
        if c in g:m[c]=g[c].iloc[0]
    return m

def per_group(task,pred,collapse=False):
    _check_pred(pred,['group_id','InChIKey'])
    if task.duplicated(['group_id','InChIKey']).any():raise ValueError('Task has duplicate group/molecule keys')
    d=task.merge(pred[['group_id','InChIKey','score']],on=['group_id','InChIKey'],how='left',validate='one_to_one')
    if d.score.isna().any():raise ValueError(f'{d.score.isna().sum()} task rows have no prediction')
    if collapse:d=collapse_stereo(d)
    rows=[_group_row(g,gid) for gid,g in d.groupby('group_id',sort=True)]
    return pd.DataFrame([x for x in rows if x is not None])

def per_group_fewshot(task,pred,protocol,collapse=False):
    keys=['group_id','repeat','InChIKey'];_check_pred(pred,keys)
    if protocol.duplicated(keys).any():raise ValueError('Protocol contains duplicate keys')
    q=protocol[(protocol.role=='query')&protocol.group_id.isin(task.group_id)]
    p=q.merge(pred[keys+['score']],on=keys,how='left',validate='one_to_one')
    if p.score.isna().any():raise ValueError(f'{p.score.isna().sum()} query rows have no prediction')
    d=p.merge(task,on=['group_id','InChIKey'],how='left',validate='many_to_one');rows=[]
    if d.pAct.isna().any():raise ValueError('Protocol contains molecules outside the task slice')
    for (gid,rep),g in d.groupby(['group_id','repeat'],sort=True):
        m=_group_row(collapse_stereo(g) if collapse else g,gid)
        if m is not None:m['repeat']=rep;rows.append(m)
    if not rows:return pd.DataFrame()
    gm=pd.DataFrame(rows);agg={c:'mean' for c in METRICS+['n_molecules','n_pairs']};agg.update({c:'first' for c in gm if c not in agg and c not in ['group_id','repeat']})
    return gm.groupby('group_id',as_index=False).agg(agg)

def _resamples(gm,n=BOOTSTRAP):
    if 'patent_core' not in gm:raise ValueError('Patent id is required for cluster bootstrap')
    clusters=[np.flatnonzero(gm.patent_core.to_numpy()==p) for p in sorted(gm.patent_core.unique())]
    if len(clusters)<2:return
    rng=np.random.default_rng(7)
    for _ in range(n):yield np.concatenate([clusters[i] for i in rng.integers(0,len(clusters),len(clusters))])

def _ci(v,patents=None):
    if patents is None:raise ValueError('Patent labels required; group-only intervals are no longer the default')
    gm=pd.DataFrame({'patent_core':patents,'v':v});gm=gm[np.isfinite(gm.v)].reset_index(drop=True)
    if gm.patent_core.nunique()<2:return np.nan,np.nan
    a=gm.v.to_numpy();values=[np.median(a[ix]) for ix in _resamples(gm)]
    return tuple(map(float,np.percentile(values,[2.5,97.5])))

def summary(gm):
    rec={'groups':len(gm),'patents':int(gm.patent_core.nunique()) if len(gm) else 0,'molecules':int(gm.n_molecules.sum()) if len(gm) else 0,'ci_unit':'patent','bootstrap_replicates':BOOTSTRAP}
    for m in METRICS:
        v=gm[m].to_numpy(float) if len(gm) else np.array([]);good=np.isfinite(v)
        rec[m+'_median']=float(np.median(v[good])) if good.any() else np.nan;rec[m+'_mean']=float(np.mean(v[good])) if good.any() else np.nan
        rec[m+'_ci_lo'],rec[m+'_ci_hi']=_ci(v,gm.patent_core) if len(gm) else (np.nan,np.nan)
    return rec

def paired_comparison(gm,reference,metric='spearman'):
    if set(gm.group_id)!=set(reference.group_id):raise ValueError('Paired comparison requires the same groups')
    p=gm[['group_id','patent_core',metric]].merge(reference[['group_id','patent_core',metric]],on='group_id',suffixes=('_model','_reference'),validate='one_to_one')
    if not p.patent_core_model.equals(p.patent_core_reference):raise ValueError('Patent assignments differ')
    p=p.rename(columns={'patent_core_model':'patent_core'});a=p[metric+'_model'].to_numpy(float);b=p[metric+'_reference'].to_numpy(float);good=np.isfinite(a)&np.isfinite(b);p=p[good].reset_index(drop=True);a,b=a[good],b[good]
    if not len(p):return {'groups':0}
    vals=[float(np.median(a[ix])-np.median(b[ix])) for ix in _resamples(p)]
    lo,hi=np.percentile(vals,[2.5,97.5]) if vals else [np.nan,np.nan]
    return {'metric':metric,'groups':len(p),'patents':p.patent_core.nunique(),'delta_median':float(np.median(a)-np.median(b)),'median_group_delta':float(np.median(a-b)),'mean_group_delta':float(np.mean(a-b)),'share_groups_better':float(np.mean(a>b)),'share_groups_tied':float(np.mean(a==b)),'ci_lo':float(lo),'ci_hi':float(hi),'ci_unit':'patent','delta_definition':'median(model group metric) - median(reference group metric)'}

def cliff_accuracy(task,pred,cliffs,part_groups):
    _check_pred(pred,['group_id','InChIKey']);score=pred.set_index(['group_id','InChIKey']).score;truth=task.set_index(['group_id','InChIKey']).pAct
    c=cliffs[cliffs.group_id.isin(part_groups)].copy()
    if not len(c):return None
    ka=list(zip(c.group_id,c.InChIKey_a));kb=list(zip(c.group_id,c.InChIKey_b))
    if any(k not in score.index or j not in score.index or k not in truth.index or j not in truth.index for k,j in zip(ka,kb)):raise ValueError('Missing prediction/truth for frozen cliff pairs')
    ds=score.loc[ka].to_numpy()-score.loc[kb].to_numpy();dt=truth.loc[ka].to_numpy()-truth.loc[kb].to_numpy();c['correct']=((np.sign(ds)==np.sign(dt)).astype(float)+.5*(ds==0));pc=task.drop_duplicates('group_id').set_index('group_id').patent_core;c['patent_core']=c.group_id.map(pc)
    stat=c.groupby(['group_id','patent_core'],as_index=False).agg(total=('correct','size'),correct=('correct','sum'));nums=stat.correct.to_numpy();dens=stat.total.to_numpy();vals=[nums[ix].sum()/dens[ix].sum() for ix in _resamples(stat)];lo,hi=np.percentile(vals,[2.5,97.5]) if vals else [np.nan,np.nan]
    return {'pairs':len(c),'groups':c.group_id.nunique(),'patents':c.patent_core.nunique(),'accuracy':float(c.correct.mean()),'ci_lo':float(lo),'ci_hi':float(hi),'ci_unit':'patent'}

def stratified_summaries(gm,profile=None):
    d=gm.copy()
    if profile is not None:
        cols=[c for c in ['group_id','target_seen_in_train_same_task','mapping_strict','docking_eligible'] if c in profile];d=d.merge(profile[cols].drop_duplicates('group_id'),on='group_id',how='left',validate='one_to_one')
    out={}
    for col in ['test_tier','endpoint','is_human_target','target_seen_in_train_same_task','mapping_strict','docking_eligible']:
        if col in d:out[col]={str(v):summary(g) for v,g in d.groupby(col,dropna=False)}
    return out

def read(path):
    return pd.read_parquet(path) if str(path).endswith('.parquet') else pd.read_csv(path,sep='\t' if str(path).endswith('.tsv') else ',')

def json_safe(x):
    if isinstance(x,dict):return {str(k):json_safe(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [json_safe(v) for v in x]
    if isinstance(x,(np.integer,)):return int(x)
    if isinstance(x,(float,np.floating)):return float(x) if np.isfinite(x) else None
    if isinstance(x,np.bool_):return bool(x)
    return x

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('task');ap.add_argument('predictions');ap.add_argument('--split',default='split_time');ap.add_argument('--part',default='test');ap.add_argument('--tier');ap.add_argument('--collapse-stereo',action='store_true');ap.add_argument('--cliffs');ap.add_argument('--fewshot');ap.add_argument('--out');ap.add_argument('--reference');ap.add_argument('--profile');ap.add_argument('--subset',choices=['mapping_strict','docking_eligible']);ap.add_argument('--json-out');a=ap.parse_args()
    full=read(a.task);task=full[full[a.split]==a.part]
    scope_rows=len(task);scope_groups=task.group_id.nunique()
    pf=Path(a.profile) if a.profile else (Path(a.task).resolve().parent.parent/'eval_support/test_group_profile.parquet' if a.split=='split_time' and a.part=='test' else None);profile=read(pf) if pf is not None and pf.exists() else None
    if profile is not None and 'task' in profile:profile=profile[profile.task==Path(a.task).stem]
    if a.subset:
        if profile is None:raise ValueError('A profile is required to define a frozen subset')
        ids=set(profile.loc[profile[a.subset].fillna(False).astype(bool),'group_id']);task=task[task.group_id.isin(ids)]
    if a.tier:
        if a.split not in ['split_time','split_all']:raise ValueError('test_tier is defined only for time/all evaluation')
        task=task[task.test_tier==a.tier]
    pred=read(a.predictions);protocol=read(a.fewshot) if a.fewshot else None
    gm=per_group_fewshot(task,pred,protocol,a.collapse_stereo) if a.fewshot else per_group(task,pred,a.collapse_stereo)
    if not len(gm):raise ValueError('No scorable groups in this slice')
    s=summary(gm);report={'summary':s,'coverage':{'task_rows':len(task),'task_groups':task.group_id.nunique(),'scored_groups':len(gm),'score_semantics':'higher means more potent','scope_rows_before_filter':scope_rows,'scope_groups_before_filter':scope_groups,'frozen_subset':a.subset}}
    if a.reference:
        ref=read(a.reference);rg=per_group_fewshot(task,ref,protocol,a.collapse_stereo) if a.fewshot else per_group(task,ref,a.collapse_stereo);report['paired_vs_reference']=paired_comparison(gm,rg)
    if a.cliffs:
        if a.fewshot:raise ValueError('Cliff evaluation requires zero-shot predictions; run separately')
        report['cliffs']=cliff_accuracy(task,pred,read(a.cliffs),set(task.group_id))
    report['strata']=stratified_summaries(gm,profile)
    if a.out:gm.to_csv(a.out,index=False)
    if a.json_out:Path(a.json_out).write_text(json.dumps(json_safe(report),ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(json_safe(report),ensure_ascii=False,indent=2))
if __name__=='__main__':main()
