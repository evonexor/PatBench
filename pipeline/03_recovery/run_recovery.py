"""20-process, ordered-phase SAR recovery with immutable baseline and checkpoints."""
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import argparse
import csv
import fcntl
import hashlib
import json
import os
import shutil
import sys
import time
import traceback
import pyarrow as pa
import pyarrow.parquet as pq
from activity_recovery import extract, RECORD_COLUMNS, NUMERIC, INTEGER, BOOL, SAR_COLUMNS, RAW
from sharded_writer import ShardedWriter

BASELINE=Path('/home/open/patent-work/outputs/uspto-sar-full-20261004/results/SAR.parquet')
SCHEMA=pa.schema([pa.field(k,pa.float64() if k in NUMERIC else pa.int64() if k in INTEGER else pa.bool_() if k in BOOL else pa.string()) for k in RECORD_COLUMNS])


def now():return datetime.now(timezone.utc).isoformat()


def atomic(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2));tmp.replace(path)


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for part in iter(lambda:f.read(4*2**20),b''):h.update(part)
    return h.hexdigest()


def code_digest():
    base=Path(__file__).parent
    files=sorted(base.rglob('*.py'))
    return hashlib.sha256(b''.join(str(p.relative_to(base)).encode()+b'\0'+p.read_bytes() for p in files)).hexdigest()


def normalize_specs(specs):
    labels={1:'review_only',2:'previously_skipped',3:'remaining_processed'}
    normalized=[]
    for original in specs:
        s=dict(original);phase=s.get('phase')
        if phase is None:
            phase=1 if s.get('status')=='processed' and not s.get('paired_rows') and s.get('review_rows',0) else 2 if s.get('status')!='processed' else 3
        s['phase']=labels.get(phase,phase)
        if s['phase'] not in labels.values():raise ValueError('unknown_phase:'+str(s['phase']))
        normalized.append(s)
    return normalized


def checkpoint_worker(job):
    output,spec,code_hash=job;archive=spec['archive'];h=hashlib.sha256(archive.encode()).hexdigest()
    directory=Path(output)/'checkpoints'/h[:2];directory.mkdir(parents=True,exist_ok=True)
    cp=directory/(h+'.json');data=directory/(h+'.parquet')
    st=(RAW/archive).stat()
    sig=hashlib.sha256((code_hash+'|'+json.dumps(spec,sort_keys=True)+'|'+str(st.st_size)+'|'+str(st.st_mtime_ns)).encode()).hexdigest()
    if cp.exists():
        old=json.loads(cp.read_text())
        if old['signature']==sig and old['stats']['status']!='archive_error' and data.is_file() and pq.read_metadata(data).num_rows==old['stats']['records']:
            return dict(checkpoint=str(cp),stats=old['stats'],cached=True)
    started=time.monotonic()
    try:
        stats,records=extract(spec)
        if len({r['source_cell_key'] for r in records})!=len(records):raise ValueError('duplicate_source_activity_cell')
    except (TypeError,KeyError,AssertionError,NameError,AttributeError):
        raise  # Deterministic code failures stop the job rather than hiding as data errors.
    except Exception as e:
        stats=dict(archive=archive,Patent=spec.get('Patent') or '',phase=spec.get('phase','fixed_sample'),prior_status=spec.get('status',''),status='archive_error',records=0,numeric_SAR=0,extended_activity=0,review=0,baseline_duplicates=0,tables=0,unparsed_activity_tables=0,source_format='',error=repr(e),traceback=traceback.format_exc(),review_reasons={});records=[]
    stats['elapsed_seconds']=round(time.monotonic()-started,3)
    table=pa.Table.from_pylist(records,schema=SCHEMA)
    tmp=data.with_suffix('.parquet.tmp');pq.write_table(table,tmp,compression='zstd',row_group_size=32768);tmp.replace(data)
    assert pq.read_metadata(data).num_rows==len(records)
    atomic(cp,dict(signature=sig,code_sha256=code_hash,source_zip_size=st.st_size,source_zip_mtime_ns=st.st_mtime_ns,stats=stats,data_file=str(data)))
    return dict(checkpoint=str(cp),stats=stats,cached=False)


def aggregate(output,specs,status,baseline_hash):
    output=Path(output);results=output/'results';results.mkdir(exist_ok=True)
    destinations={'numeric_SAR':'new_SAR.parquet','extended_activity':'verified_extended_activity.parquet','review':'remaining_review.parquet'}
    writers={}
    for key,name in destinations.items():
        tmp=results/(name+'.tmpdir')
        if tmp.exists():shutil.rmtree(tmp)  # Only this job's interrupted aggregation.
        writers[key]=ShardedWriter(tmp,SCHEMA)
    patents_schema=pa.schema([pa.field(k,pa.int64() if k in {'records','numeric_SAR','extended_activity','review','baseline_duplicates','tables','unparsed_activity_tables'} else pa.string()) for k in ['archive','Patent','phase','prior_status','status','source_format','records','numeric_SAR','extended_activity','review','baseline_duplicates','tables','unparsed_activity_tables','error']])
    pat_tmp=results/'patents.parquet.tmpdir'
    if pat_tmp.exists():shutil.rmtree(pat_tmp)
    patentwriter=ShardedWriter(pat_tmp,patents_schema)
    failures=[];sample={k:[] for k in writers};patentsets={k:set() for k in writers};inchikeys=set();phases=Counter();reasons=Counter();formats=Counter();parsewarnings=Counter();unparsed=0;phase_metrics={};baseline_duplicates=0
    for n,spec in enumerate(specs):
        h=hashlib.sha256(spec['archive'].encode()).hexdigest();cp=output/'checkpoints'/h[:2]/(h+'.json')
        v=json.loads(cp.read_text());stats=v['stats']
        assert v['code_sha256']==status['code_sha256']
        st=(RAW/spec['archive']).stat();assert v['source_zip_size']==st.st_size and v['source_zip_mtime_ns']==st.st_mtime_ns
        assert pq.read_metadata(v['data_file']).num_rows==stats['records']
        assert stats['records']==stats['numeric_SAR']+stats['extended_activity']+stats['review']
        phases[stats['phase']]+=1;formats[stats['source_format']]+=1;reasons.update(stats['review_reasons']);unparsed+=stats['unparsed_activity_tables']
        pm=phase_metrics.setdefault(stats['phase'],Counter());pm['archives']+=1
        for key in ('numeric_SAR','extended_activity','review'):pm[key+'_rows']+=stats[key];pm['patents_with_'+key]+=int(stats[key]>0)
        pm['patents_with_any_new_verified_activity']+=int(stats['numeric_SAR']+stats['extended_activity']>0)
        pm['patents_with_new_activity_and_remaining_review']+=int(stats['numeric_SAR']+stats['extended_activity']>0 and stats['review']>0)
        pm['patents_without_new_verified_activity']+=int(stats['numeric_SAR']+stats['extended_activity']==0)
        baseline_duplicates+=stats['baseline_duplicates']
        audit=stats.get('adapter_audit',{})
        if audit.get('unresolved_entities'):parsewarnings['archives_with_unresolved_entities']+=1
        if audit.get('sgml_closure_warnings'):parsewarnings['archives_with_sgml_closure_warnings']+=1
        if stats['status']=='archive_error':failures.append(stats)
        patentwriter.write_table(pa.Table.from_pylist([{k:stats.get(k) for k in patents_schema.names}],schema=patents_schema))
        for batch in pq.ParquetFile(v['data_file']).iter_batches(batch_size=32768):
            rows=batch.to_pylist()
            for category,writer in writers.items():
                selected=[r for r in rows if r['extraction_status']==category]
                if selected:
                    writer.write_table(pa.Table.from_pylist(selected,schema=SCHEMA))
                    patentsets[category].add(stats['Patent'])
                    if len(sample[category])<100 and len({r['Patent'] for r in sample[category]})<20:sample[category].extend(selected[:min(5,100-len(sample[category]))])
                    if category=='numeric_SAR':inchikeys.update(r['InChIKey'] for r in selected)
        if n%500==0:
            status.update(state='AGGREGATING',aggregated_archives=n+1,updated_at=now());atomic(output/'run_status.json',status)
    dataset_manifest={}
    for category,writer in writers.items():
        writer.close();name=destinations[category];final=results/name
        if final.exists():shutil.rmtree(final)
        writer.path.rename(final);dataset_manifest[name]=dict(rows=writer.total_rows,files=writer.files)
    patentwriter.close();final=results/'patents.parquet'
    if final.exists():shutil.rmtree(final)
    pat_tmp.rename(final);dataset_manifest[final.name]=dict(rows=patentwriter.total_rows,files=patentwriter.files)
    # Core numeric SAR union: stream the unmodified baseline and append new rows.
    mergedtmp=results/'SAR_merged.parquet.tmpdir'
    if mergedtmp.exists():shutil.rmtree(mergedtmp)
    basemeta=pq.read_metadata(BASELINE);base_schema=pq.read_schema(BASELINE)
    mergedwriter=ShardedWriter(mergedtmp,base_schema)
    baseline_rows=0
    for batch in pq.ParquetFile(BASELINE).iter_batches(batch_size=32768):
        table=pa.Table.from_batches([batch]);mergedwriter.write_table(table);baseline_rows+=table.num_rows
    assert baseline_rows==basemeta.num_rows==395143
    if len(specs)==345940:assert baseline_duplicates==baseline_rows,('Original formal cell dedup audit failed',baseline_duplicates,baseline_rows)
    for p in sorted((results/'new_SAR.parquet').glob('part-*.parquet')):
        for batch in pq.ParquetFile(p).iter_batches(batch_size=32768,columns=base_schema.names):
            mergedwriter.write_table(pa.Table.from_batches([batch]).cast(base_schema))
    mergedwriter.close();mergedfinal=results/'SAR_merged.parquet'
    if mergedfinal.exists():shutil.rmtree(mergedfinal)
    mergedtmp.rename(mergedfinal);dataset_manifest[mergedfinal.name]=dict(rows=mergedwriter.total_rows,files=mergedwriter.files)
    assert digest(BASELINE)==baseline_hash
    # Independent read-back equality of every original cell, not just count/SHA.
    left=pq.ParquetFile(BASELINE).read().combine_chunks()
    pieces=[];needed=baseline_rows
    for p in sorted(mergedfinal.glob('part-*.parquet')):
        if needed<=0:break
        table=pq.read_table(p).slice(0,needed);pieces.append(table);needed-=table.num_rows
    assert needed==0 and left.equals(pa.concat_tables(pieces).combine_chunks())
    baseline_patents=set(left['Patent'].to_pylist())
    validation=dict(complete_archive_accounting=patentwriter.total_rows==len(specs),all_checkpoint_source_signatures_verified=True,baseline_rows=baseline_rows,baseline_sha256=baseline_hash,baseline_all_columns_exactly_unchanged=True,baseline_source_cells_deduplicated=baseline_duplicates,original_formal_patents=len(baseline_patents),all_shards_readable_and_counted=True,bounded_parquet_row_groups=True,source_mapping_policy='explicit_identity_and_source_scope_plus_actual_mol_sha_smiles_inchikey',source_mapping_rechecked_for_every_new_accepted_record=True)
    for name,entry in dataset_manifest.items():
        total=0
        for file in entry['files']:
            m=pq.read_metadata(results/name/file['file']);assert m.num_rows==file['rows'] and m.num_row_groups<=32;total+=m.num_rows
        assert total==entry['rows']
    validation['merged_row_count_correct']=mergedwriter.total_rows==baseline_rows+writers['numeric_SAR'].total_rows
    report=dict(output_dir=str(output),code_sha256=status['code_sha256'],manifest_sha256=status['manifest_sha256'],completed_archives=len(specs),total_archives=len(specs),phase_counts=dict(phases),source_formats=dict(formats),new_numeric_SAR_rows=writers['numeric_SAR'].total_rows,new_numeric_SAR_patents=len(patentsets['numeric_SAR']),new_patents_without_prior_formal_SAR=len(patentsets['numeric_SAR']-baseline_patents),merged_formal_SAR_rows=mergedwriter.total_rows,merged_formal_SAR_patents=len(baseline_patents|patentsets['numeric_SAR']),new_unique_InChIKeys=len(inchikeys),verified_extended_activity_rows=writers['extended_activity'].total_rows,verified_extended_activity_patents=len(patentsets['extended_activity']),remaining_review_rows=writers['review'].total_rows,remaining_review_patents=len(patentsets['review']),archive_errors=len(failures),unparsed_activity_tables=unparsed,parser_warning_counts=dict(parsewarnings),review_reason_counts=dict(reasons),elapsed_seconds=round(time.time()-status['start_timestamp'],1),started_at=status['started_at'],completed_at=now())
    report['phase_metrics']={phase:dict(metrics) for phase,metrics in phase_metrics.items()}
    report['baseline_source_cells_deduplicated']=baseline_duplicates
    atomic(results/'dataset_manifest.json',dataset_manifest);atomic(results/'validation_report.json',validation);atomic(results/'run_report.json',report)
    with (results/'failed_archives.jsonl').open('w') as f:
        for v in failures:f.write(json.dumps(v,ensure_ascii=False)+'\n')
    for category,rows in sample.items():
        with (results/(category+'_preview.csv')).open('w',newline='',encoding='utf-8-sig') as f:
            w=csv.DictWriter(f,fieldnames=RECORD_COLUMNS);w.writeheader();w.writerows(rows)
    return report,validation


def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--workers',type=int,default=20);p.add_argument('--preflight',action='store_true');args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=True);lock=(args.output/'run.lock').open('a+')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    manifest=json.loads(args.manifest.read_text());specs=normalize_specs(manifest['samples']);assert len({s['archive'] for s in specs})==len(specs)
    if not args.preflight:
        assert len(specs)==345940
        assert Counter(s['phase'] for s in specs)=={'review_only':25703,'previously_skipped':236178,'remaining_processed':84059}
    code_hash=code_digest();mh=digest(args.manifest)
    if (args.output/'completion.json').exists():
        complete=json.loads((args.output/'completion.json').read_text())
        if complete.get('code_sha256')==code_hash and complete.get('manifest_sha256')==mh and complete.get('state') in ('COMPLETE','COMPLETE_WITH_ERRORS'):print('Already complete; no duplicate work.',flush=True);return
    status=dict(state='RUNNING',workers=args.workers,pid=os.getpid(),started_at=now(),start_timestamp=time.time(),total_archives=len(specs),completed_archives=0,code_sha256=code_hash,manifest_sha256=mh,output_dir=str(args.output),preflight=args.preflight,domain_policy='all_domains',phase_counts={},new_numeric_SAR_rows=0,verified_extended_activity_rows=0,remaining_review_rows=0,archive_errors=0,baseline_duplicates=0,cached_archives=0)
    oldpath=args.output/'run_status.json'
    if oldpath.exists():
        old=json.loads(oldpath.read_text())
        if old.get('code_sha256')==code_hash and old.get('manifest_sha256')==mh:
            status.update(started_at=old['started_at'],start_timestamp=old['start_timestamp'])
    atomic(oldpath,status);baseline_hash=digest(BASELINE)
    phases=list(dict.fromkeys(s.get('phase','fixed_sample') for s in specs))
    phases.sort(key=lambda s:{'review_only':0,'previously_skipped':1,'remaining_processed':2}.get(s,0))
    try:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            for phase in phases:
                jobs=iter([s for s in specs if s.get('phase','fixed_sample')==phase]);futures={};exhausted=False
                status['current_phase']=phase;print('Starting phase:',phase,flush=True)
                lastlog=0;phase_details=Counter();phase_examples=[]
                while not exhausted or futures:
                    while not exhausted and len(futures)<args.workers*2:
                        try:s=next(jobs)
                        except StopIteration:exhausted=True;break
                        future=pool.submit(checkpoint_worker,(str(args.output),s,code_hash));futures[future]=s
                    done,_=wait(futures,timeout=15,return_when=FIRST_COMPLETED)
                    for f in done:
                        item=f.result();futures.pop(f);st=item['stats'];status['completed_archives']+=1
                        status['phase_counts'][phase]=status['phase_counts'].get(phase,0)+1
                        for key,field in [('new_numeric_SAR_rows','numeric_SAR'),('verified_extended_activity_rows','extended_activity'),('remaining_review_rows','review'),('baseline_duplicates','baseline_duplicates')]:status[key]+=st[field]
                        status['archive_errors']+=int(st['status']=='archive_error');status['cached_archives']+=int(item['cached'])
                        phase_details['archives']+=1
                        for field in ('numeric_SAR','extended_activity','review'):
                            phase_details[field+'_rows']+=st[field];phase_details['patents_with_'+field]+=int(st[field]>0)
                        phase_details['patents_with_any_new_verified_activity']+=int(st['numeric_SAR']+st['extended_activity']>0)
                        phase_details['archive_errors']+=int(st['status']=='archive_error')
                        if st['numeric_SAR']+st['extended_activity'] and len(phase_examples)<15:phase_examples.append(item['checkpoint'])
                    if time.time()-lastlog>=30 or exhausted and not futures:
                        status.update(updated_at=now(),elapsed_seconds=round(time.time()-status['start_timestamp'],1),inflight=len(futures));atomic(oldpath,status)
                        print(json.dumps({k:status[k] for k in ['current_phase','completed_archives','total_archives','new_numeric_SAR_rows','verified_extended_activity_rows','remaining_review_rows','archive_errors','elapsed_seconds']},ensure_ascii=False),flush=True);lastlog=time.time()
                # Make phase one reviewable before reading previously unparsed archives.
                atomic(args.output/('phase_'+phase+'_complete.json'),dict(status,phase_completed_at=now()))
                atomic(args.output/('phase_'+phase+'_report.json'),dict(phase=phase,metrics=dict(phase_details),code_sha256=code_hash,manifest_sha256=mh,completed_at=now()))
                with (args.output/('phase_'+phase+'_preview.csv')).open('w',newline='',encoding='utf-8-sig') as handle:
                    preview=csv.DictWriter(handle,fieldnames=RECORD_COLUMNS);preview.writeheader()
                    for cp in phase_examples:
                        data=json.loads(Path(cp).read_text())['data_file']
                        accepted=[r for r in pq.read_table(data).to_pylist() if r['extraction_status']!='review']
                        preview.writerows(accepted[:3])
        status.update(state='AGGREGATING',updated_at=now());atomic(oldpath,status)
        report,validation=aggregate(args.output,specs,status,baseline_hash)
        assert validation['complete_archive_accounting'] and validation['merged_row_count_correct']
        state='COMPLETE_WITH_ERRORS' if report['archive_errors'] or report['unparsed_activity_tables'] or report['parser_warning_counts'] else 'COMPLETE'
        status.update(state=state,completed_archives=len(specs),completed_at=now(),elapsed_seconds=report['elapsed_seconds']);atomic(oldpath,status);atomic(args.output/'completion.json',status)
        print(json.dumps(report,ensure_ascii=False),flush=True)
    except Exception:
        status.update(state='FAILED',error=traceback.format_exc(),updated_at=now());atomic(oldpath,status);raise


if __name__=='__main__':main()
