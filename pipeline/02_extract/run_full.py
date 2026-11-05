"""Bounded 20-process full extraction, immediate streaming exports and durable progress."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
from datetime import datetime,timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import sqlite3
import time
import traceback
from full_quality import full_worker,QUALITY_COLUMNS
from patent_metadata import METADATA_COLUMNS
from uspto_sar import RecordWriter,SAR_COLUMNS,EXTRA_COLUMNS,atomic_json

VERSION='0.8.0-full'

def now():return datetime.now(timezone.utc).isoformat()

def run(args):
    base=args.output.parent;base.mkdir(parents=True,exist_ok=True);args.output.mkdir(parents=True,exist_ok=True)
    lock=(base/'run.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    start=time.monotonic();started=now();manifest=json.loads(args.manifest.read_text());specs=manifest['samples']
    if len(specs)>100 and not args.full:raise ValueError('A full run requires --full.')
    if len({s['archive'] for s in specs})!=len(specs):raise ValueError('Duplicate archive in manifest.')
    if args.root.resolve()==args.output.resolve() or args.root.resolve() in args.output.resolve().parents:raise ValueError('Output must be outside raw data.')
    code_hash=hashlib.sha256(b''.join(p.read_bytes() for p in sorted(Path(__file__).parent.glob('*.py')))).hexdigest()
    manifest_hash=hashlib.sha256(args.manifest.read_bytes()).hexdigest()
    status={'version':VERSION,'state':'RUNNING','started_at':started,'pid':os.getpid(),'workers':args.workers,'total_archives':len(specs),'completed_archives':0,'formal_sar_rows':0,'candidate_rows':0,'review_rows':0,'structure_rows':0,'source_check_failed_rows':0,'statuses':{},'manifest_sha256':manifest_hash,'code_sha256':code_hash,'output_dir':str(args.output),'association_policy':'direct_only_after_source_checks','drawing_vs_mol_checked':False}
    atomic_json(base/'job_config.json',{**status,'root':str(args.root),'manifest':str(args.manifest),'command':'run_full.py --workers 20 --full'})
    writers={};db=None;table_log=None;error_log=None
    try:
        all_cols=SAR_COLUMNS+EXTRA_COLUMNS+METADATA_COLUMNS+QUALITY_COLUMNS
        pat_cols=list(dict.fromkeys(['Patent','patent_core','title','publication_date','selection','screen','screen_reason','status','mol_files','tables','explicit_aliases','ambiguous_aliases','paired_rows','review_rows','candidate_rows','quality_check_failures','unparsed_activity_tables','structure_status_counts','elapsed_seconds','archive','xml_file','error']+METADATA_COLUMNS))
        for name,cols in [('SAR',SAR_COLUMNS),('SAR_enriched',SAR_COLUMNS+METADATA_COLUMNS+QUALITY_COLUMNS),('SAR_with_evidence',all_cols),('candidate_SAR',all_cols),('review',all_cols),('structures',['Patent','compound_aliases','role','structure_status','SMILES','SMILES_original','InChIKey','structure_file','structure_name','structure_name_check','archive']),('patents',pat_cols)]:writers[name]=RecordWriter(args.output,name,cols)
        db=sqlite3.connect(args.output/'counts.sqlite');db.execute('DROP TABLE IF EXISTS structures');db.execute('CREATE TABLE structures (inchikey TEXT PRIMARY KEY) WITHOUT ROWID');db.execute('DROP TABLE IF EXISTS observations');db.execute('CREATE TABLE observations (id TEXT PRIMARY KEY) WITHOUT ROWID')
        table_log=(args.output/'unparsed_activity_tables.jsonl.tmp').open('w');error_log=(args.output/'failed_archives.jsonl.tmp').open('w')
        statuses=Counter();reasons=Counter();coverage=Counter();paired_patents=0;bib_files=0;cached=0;last_status=0;last_print=0
        def publish(inflight=0):
            elapsed=time.monotonic()-start;rate=status['completed_archives']/elapsed if elapsed else 0
            status.update(updated_at=now(),elapsed_seconds=round(elapsed,1),archives_per_second=round(rate,3),estimated_remaining_seconds=round((len(specs)-status['completed_archives'])/rate) if rate else None,statuses=dict(statuses),inflight=inflight,cache_hits=cached,free_disk_GiB=round(shutil.disk_usage(args.output).free/2**30,1))
            atomic_json(base/'run_status.json',status)
        publish()
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            iterator=iter(specs);pending={}
            def submit():
                spec=next(iterator,None)
                if spec is not None:pending[pool.submit(full_worker,(str(args.root),str(args.output),spec,code_hash))]=spec
            for _ in range(min(len(specs),args.workers*2)):submit()
            while pending:
                completed,_=wait(pending,timeout=15,return_when=FIRST_COMPLETED)
                for future in completed:
                    spec=pending.pop(future)
                    try:small=future.result()
                    except Exception as exc:small={'checkpoint':None,'stats':{'archive':spec['archive'],'status':'error','error':repr(exc)},'cached':False}
                    if small['checkpoint']:obj=json.loads(Path(small['checkpoint']).read_text());stats=obj['stats']
                    else:stats=small['stats'];obj={'ready':[],'review':[],'candidates':[],'structures':[],'tables':[]}
                    statuses[stats['status']]+=1;cached+=bool(small['cached'])
                    if stats['status']=='error':error_log.write(json.dumps(stats,ensure_ascii=False)+'\n');error_log.flush()
                    for r in obj['ready']:
                        assert r['alignment_grade']=='direct_numbered_structure_row' and r['source_mapping_checks']=='passed' and not r['review_reason'] and r['raw_molfile_path'],r['source_record_id']
                    ids=[(r['source_record_id'],) for r in obj['ready']+obj['review']]
                    db.executemany('INSERT INTO observations VALUES (?)',ids)
                    for name,key in [('SAR','ready'),('SAR_enriched','ready'),('SAR_with_evidence','ready'),('candidate_SAR','candidates'),('review','review'),('structures','structures')]:writers[name].append(obj[key])
                    writers['patents'].append([stats])
                    for t in obj['tables']:table_log.write(json.dumps(t,ensure_ascii=False)+'\n')
                    db.executemany('INSERT OR IGNORE INTO structures VALUES (?)',[(r['InChIKey'],) for r in obj['structures'] if r['structure_status']=='concrete'])
                    if status['completed_archives']%100==0:db.commit()
                    bib=obj.get('documents',{}).get('bibliography')
                    if bib:
                        p=args.output/obj['bibliography_path'];p.parent.mkdir(parents=True,exist_ok=True);atomic_json(p,bib);bib_files+=1
                    coverage.update(k for k in METADATA_COLUMNS if stats.get(k) not in ('',None,[],{}))
                    reasons.update(reason for r in obj['review'] for reason in r['review_reason'].split(';') if reason)
                    status['completed_archives']+=1;status['formal_sar_rows']+=len(obj['ready']);status['candidate_rows']+=len(obj['candidates']);status['review_rows']+=len(obj['review']);status['structure_rows']+=len(obj['structures']);status['source_check_failed_rows']+=stats.get('quality_check_failures',0);paired_patents+=bool(obj['ready'])
                    submit()
                elapsed=time.monotonic()-start
                if elapsed-last_status>=5 or not completed:
                    publish(len(pending));last_status=elapsed
                if elapsed-last_print>=30:
                    print(json.dumps({k:status[k] for k in ['completed_archives','total_archives','formal_sar_rows','candidate_rows','review_rows','elapsed_seconds']},ensure_ascii=False),flush=True);last_print=elapsed
                if shutil.disk_usage(args.output).free<20*2**30:raise RuntimeError('Less than 20 GiB disk space available; checkpoints retained.')
        status['state']='FINALIZING';publish()
        db.commit();unique=db.execute('SELECT COUNT(*) FROM structures').fetchone()[0];measurements=db.execute('SELECT COUNT(*) FROM observations').fetchone()[0]
        for writer in writers.values():writer.close()
        table_log.close();error_log.close();(args.output/'unparsed_activity_tables.jsonl.tmp').replace(args.output/'unparsed_activity_tables.jsonl');(args.output/'failed_archives.jsonl.tmp').replace(args.output/'failed_archives.jsonl')
        metadata={'version':VERSION,'field_coverage':dict(coverage),'review_reasons':dict(reasons),'bibliography_files':bib_files}
        atomic_json(args.output/'metadata_report.json',metadata)
        report={**status,'paired_patents':paired_patents,'unique_concrete_structures':unique,'unique_measurements':measurements,'bibliography_files':bib_files,'review_reasons':dict(reasons),'row_order':'completion order; identities are stable','limitations':['Original MOL graphs have not been checked against every original structure drawing.','Only direct numbered structure associations passing independent source checks enter formal SAR.','Adjacent example diagrams and any failed source checks remain in candidate_SAR and review.','Legacy XML/SGML and unsupported layouts remain explicitly listed for further parsing.','UniProt/gene mapping, image OCR and inferred assay units are not performed.']}
        atomic_json(args.output/'run_report.json',report)
        status['state']='VERIFYING_EXPORTS';publish()
        import pyarrow.parquet as pq
        counts={name:pq.read_metadata(args.output/(name+'.parquet')).num_rows for name in writers}
        assert counts['SAR']==counts['SAR_enriched']==counts['SAR_with_evidence']==status['formal_sar_rows']
        assert counts['candidate_SAR']==status['candidate_rows'] and counts['review']==status['review_rows'] and counts['patents']==len(specs)
        assert measurements==counts['SAR']+counts['review']
        atomic_json(args.output/'validation_report.json',{'parquet_row_counts':counts,'unique_measurement_ids':measurements,'source_checks_in_workers':'passed or isolated in review','complete_archive_accounting':True,'drawing_vs_mol_checked':False})
        status['state']='COMPLETE_WITH_ERRORS' if statuses['error'] else 'COMPLETE';status.update(finished_at=now(),unique_concrete_structures=unique,bibliography_files=bib_files);publish()
        report.update(status);atomic_json(args.output/'run_report.json',report)
        atomic_json(base/'completion.json',{**status,'validation_report':str(args.output/'validation_report.json'),'run_report':str(args.output/'run_report.json')})
        print(json.dumps(status,ensure_ascii=False),flush=True)
    except BaseException as exc:
        status.update(state='FAILED',updated_at=now(),error=repr(exc),traceback=traceback.format_exc(),elapsed_seconds=round(time.monotonic()-start,1));atomic_json(base/'run_status.json',status);raise
    finally:
        if db:db.close()
        if table_log and not table_log.closed:table_log.close()
        if error_log and not error_log.closed:error_log.close()

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--workers',type=int,default=20);p.add_argument('--full',action='store_true');run(p.parse_args())
