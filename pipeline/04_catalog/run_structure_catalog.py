"""20-process extraction of all original chemical sources, all domains, no activity gate."""
import argparse,fcntl,hashlib,json,os,shutil,sqlite3,time,traceback
from collections import Counter
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
from datetime import datetime,timezone
from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq
from structure_catalog import worker,SCHEMA,STRUCTURE_COLUMNS
from catalog_annotations import SEQ_SCHEMA,ACTIVITY_SCHEMA,UNIPROT_SCHEMA,FIELDS
from patent_metadata import METADATA_COLUMNS
from patent_documents import DOCUMENT_INTEGER_COLUMNS
from uspto_sar import atomic_json

VERSION='0.9.0-all-structures'
def now():return datetime.now(timezone.utc).isoformat()
def scalar(v):return json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v
def extras(schema):return pa.schema(list(schema)+[pa.field(k,pa.int64() if k in DOCUMENT_INTEGER_COLUMNS else pa.string()) for k in METADATA_COLUMNS if k not in schema.names])
def run(a):
    base=a.output.parent;base.mkdir(parents=True,exist_ok=True);a.output.mkdir(parents=True,exist_ok=True)
    lock=(base/'run.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    specs=json.loads(a.manifest.read_text())['samples']
    if len(specs)>100 and not a.full:raise ValueError('--full is required for a corpus run')
    assert len({s['archive'] for s in specs})==len(specs)
    assert a.root.resolve() not in a.output.resolve().parents
    started=time.monotonic();code_hash=hashlib.sha256(b''.join(p.read_bytes() for p in sorted(Path(__file__).parent.glob('*.py')))).hexdigest()
    status=dict(version=VERSION,state='RUNNING',pid=os.getpid(),workers=a.workers,started_at=now(),total_archives=len(specs),completed_archives=0,structure_rows=0,mol_files=0,non_mol_sources=0,sequences_rows=0,activity_annotation_rows=0,formal_SAR_annotation_rows=0,patents_with_sources=0,previous_verified_SAR_expected=0,previous_verified_SAR_linked=0,source_read_errors=0,archive_errors=0,annotation_files_with_errors=0,code_sha256=code_hash,manifest_sha256=hashlib.sha256(a.manifest.read_bytes()).hexdigest(),output_dir=str(a.output),domain_policy='all_domains_including_agriculture_materials',structure_gate='every_mol_regardless_metadata_or_activity',drawing_vs_mol_checked=False)
    atomic_json(base/'job_config.json',{**status,'raw_root':str(a.root),'prior_results':str(a.prior),'manifest':str(a.manifest)})
    status['uniprot_mentions_rows']=0
    pat_cols=list(dict.fromkeys(['archive','Patent','patent_core','status','prior_status','metadata_format','metadata_status','metadata_unresolved_entities','xml_file','screen','screen_reason','structure_rows','mol_files','non_mol_sources','image_files','source_read_errors','sequences_rows','activities_rows','uniprot_mentions_rows','previous_verified_SAR_expected','previous_verified_SAR_linked','elapsed_seconds','structure_status_counts','source_format_counts','annotation_errors','error']+METADATA_COLUMNS))
    ints={'structure_rows','mol_files','non_mol_sources','image_files','source_read_errors','sequences_rows','activities_rows','uniprot_mentions_rows','previous_verified_SAR_expected','previous_verified_SAR_linked'}|DOCUMENT_INTEGER_COLUMNS
    pat_schema=pa.schema([pa.field(k,pa.int64() if k in ints else pa.float64() if k=='elapsed_seconds' else pa.string()) for k in pat_cols])
    schemas={'structures':extras(SCHEMA),'patents':pat_schema,'patent_sequences':SEQ_SCHEMA,'activity_annotations':ACTIVITY_SCHEMA,'SAR_enriched_annotations':extras(ACTIVITY_SCHEMA),'patent_uniprot_mentions':UNIPROT_SCHEMA}
    writers={};db=None;fail=None
    try:
        for name,schema in schemas.items():writers[name]=pq.ParquetWriter(a.output/(name+'.parquet.tmp'),schema,compression='zstd',use_dictionary=[k for k in schema.names if k not in ('raw_structure_bytes','source_structure_id','Sequence','sequence_original','source_record_id','source_sha256')])
        db=sqlite3.connect(a.output/'counts.sqlite');db.execute('PRAGMA journal_mode=WAL');db.execute('PRAGMA synchronous=NORMAL');db.execute('PRAGMA cache_size=-65536');db.execute('DROP TABLE IF EXISTS structures');db.execute('CREATE TABLE structures (inchikey TEXT PRIMARY KEY) WITHOUT ROWID')
        fail=(a.output/'failed_archives.jsonl.tmp').open('w');states=Counter();formats=Counter();structure_states=Counter();prior_states=Counter();metadata_states=Counter();coverage=Counter();cache_hits=0;last_publish=0;last_log=0
        def publish(inflight=0):
            elapsed=time.monotonic()-started
            status.update(updated_at=now(),elapsed_seconds=round(elapsed,1),inflight=inflight,statuses=dict(states),source_formats=dict(formats),structure_status_counts=dict(structure_states),prior_statuses=dict(prior_states),metadata_statuses=dict(metadata_states),cache_hits=cache_hits,free_disk_GiB=round(shutil.disk_usage(a.output).free/2**30,1))
            atomic_json(base/'run_status.json',status)
        def enriched(batch,schema,stats):
            fields={k:batch.column(k) for k in batch.schema.names}
            for k in schema.names:
                if k not in fields:
                    v=scalar(stats.get(k));v=str(v) if v is not None and pa.types.is_string(schema.field(k).type) else v
                    fields[k]=pa.array([v]*batch.num_rows,type=schema.field(k).type)
            return pa.Table.from_arrays([fields[k] for k in schema.names],schema=schema)
        publish()
        with ProcessPoolExecutor(max_workers=a.workers) as pool:
            pending={};iterator=iter(specs)
            def submit():
                spec=next(iterator,None)
                if spec is not None:pending[pool.submit(worker,(str(a.root),str(a.output),str(a.prior),spec,code_hash))]=spec
            for _ in range(min(len(specs),a.workers*2)):submit()
            while pending:
                completed,_=wait(pending,timeout=15,return_when=FIRST_COMPLETED)
                for future in completed:
                    spec=pending.pop(future)
                    result=None
                    try:
                        result=future.result();obj=json.loads(Path(result['checkpoint']).read_text());stats=obj['stats'];part=Path(result['part']);cache_hits+=result['cached']
                        assert pq.read_metadata(part).num_rows==stats['mol_files']+stats['non_mol_sources']==stats['structure_rows']
                        assert stats['previous_verified_SAR_expected']==stats['previous_verified_SAR_linked'],('Previous SAR correspondence changed',spec['archive'])
                        for batch in pq.ParquetFile(part).iter_batches(batch_size=2048):
                            writers['structures'].write_table(enriched(batch,schemas['structures'],stats))
                            keys=[(s,) for s in batch.column('InChIKey').to_pylist() if s];db.executemany('INSERT OR IGNORE INTO structures VALUES (?)',keys)
                            for k in FIELDS:coverage['structures:'+k]+=sum(bool(v) and v!='[]' for v in batch.column(k).to_pylist())
                        key=part.stem
                        for name,destination in [('sequences','patent_sequences'),('activities','activity_annotations'),('uniprot_mentions','patent_uniprot_mentions')]:
                            if not stats.get(name+'_rows'):continue
                            source=part.with_name(key+'.'+name+'.parquet')
                            assert pq.read_metadata(source).num_rows==stats[name+'_rows']
                            for batch in pq.ParquetFile(source).iter_batches(batch_size=2048):
                                writers[destination].write_batch(batch)
                                if name=='activities':
                                    selected=pa.Table.from_batches([batch]).filter(batch.column('is_formal_SAR'))
                                    if selected.num_rows:
                                        assert all(selected.column('source_structure_id').to_pylist())
                                        writers['SAR_enriched_annotations'].write_table(enriched(selected,schemas['SAR_enriched_annotations'],stats));status['formal_SAR_annotation_rows']+=selected.num_rows
                                    for k in FIELDS:coverage['activities:'+k]+=sum(bool(v) for v in batch.column(k).to_pylist())
                        status['patents_with_sources']+=bool(stats['structure_rows']);status['annotation_files_with_errors']+=bool(stats['annotation_errors'])
                        for dest,src in [('structure_rows','structure_rows'),('mol_files','mol_files'),('non_mol_sources','non_mol_sources'),('sequences_rows','sequences_rows'),('activity_annotation_rows','activities_rows'),('uniprot_mentions_rows','uniprot_mentions_rows'),('source_read_errors','source_read_errors'),('previous_verified_SAR_expected','previous_verified_SAR_expected'),('previous_verified_SAR_linked','previous_verified_SAR_linked')]:status[dest]+=stats[src]
                        structure_states.update(stats['structure_status_counts']);formats.update(stats['source_format_counts']);prior_states[stats['prior_status']]+=1;metadata_states[stats.get('metadata_status','unavailable')]+=1
                    except Exception as exc:
                        if result is not None:raise
                        stats=dict(archive=spec['archive'],patent_core=spec['patent_core'],status='archive_error',error=repr(exc));status['archive_errors']+=1;fail.write(json.dumps(stats,ensure_ascii=False)+'\n');fail.flush()
                    clean={k:scalar(stats.get(k)) for k in pat_cols}
                    clean={k:str(v) if v is not None and pa.types.is_string(pat_schema.field(k).type) else v for k,v in clean.items()}
                    writers['patents'].write_table(pa.Table.from_pylist([clean],schema=pat_schema));states[stats['status']]+=1;status['completed_archives']+=1
                    if status['completed_archives']%100==0:db.commit()
                    submit()
                elapsed=time.monotonic()-started
                if elapsed-last_publish>=5 or not completed:publish(len(pending));last_publish=elapsed
                if elapsed-last_log>=30:print(json.dumps({k:status[k] for k in ('completed_archives','total_archives','structure_rows','mol_files','sequences_rows','archive_errors')},ensure_ascii=False),flush=True);last_log=elapsed
                if shutil.disk_usage(a.output).free<30*2**30:raise RuntimeError('Less than 30 GiB available; checkpoints retained')
        status['state']='FINALIZING';publish();db.commit()
        for name,w in writers.items():w.close();(a.output/(name+'.parquet.tmp')).replace(a.output/(name+'.parquet'))
        fail.close();(a.output/'failed_archives.jsonl.tmp').replace(a.output/'failed_archives.jsonl')
        counts={name:pq.read_metadata(a.output/(name+'.parquet')).num_rows for name in writers}
        assert counts['patents']==len(specs) and counts['structures']==status['mol_files']+status['non_mol_sources']==status['structure_rows']
        assert counts['patent_sequences']==status['sequences_rows'] and counts['activity_annotations']==status['activity_annotation_rows'] and counts['SAR_enriched_annotations']==status['formal_SAR_annotation_rows']==status['previous_verified_SAR_linked']
        assert counts['patent_uniprot_mentions']==status['uniprot_mentions_rows']
        unique=db.execute('SELECT COUNT(*) FROM structures').fetchone()[0]
        atomic_json(a.output/'validation_report.json',{'parquet_row_counts':counts,'complete_archive_accounting':True,'complete_readable_archive_source_accounting':status['archive_errors']==0,'all_domains_included':True,'every_mol_in_readable_archives_recorded':True,'original_bytes_retained':True,'previous_verified_SAR_correspondence_preserved':status['previous_verified_SAR_expected']==status['previous_verified_SAR_linked'],'drawing_vs_mol_checked':False})
        status.update(state='COMPLETE_WITH_ERRORS' if status['archive_errors'] or status['source_read_errors'] or status['annotation_files_with_errors'] else 'COMPLETE',finished_at=now(),unique_valid_inchikeys=unique);publish()
        atomic_json(a.output/'run_report.json',{**status,'field_coverage':dict(coverage),'limitations':['Chemical images/CDX without MOL are retained as sources requiring conversion/OCR.','Unidentified images in drawing pages are not classified as chemical structures automatically.','Source-reported chemical names are not independently validated IUPAC nomenclature.','Sequence or UniProt data are not assigned to every compound in the same patent.','Original MOL graphs still require comparison with original drawings for image-recognition QA.','Additional SAR extraction in newly included domains is not inferred; previously verified SAR is preserved.']})
        atomic_json(base/'completion.json',{**status,'validation_report':str(a.output/'validation_report.json'),'run_report':str(a.output/'run_report.json')});print(json.dumps(status,ensure_ascii=False),flush=True)
    except BaseException as exc:
        status.update(state='FAILED',updated_at=now(),error=repr(exc),traceback=traceback.format_exc());atomic_json(base/'run_status.json',status);raise
    finally:
        if db:db.close()
        if fail and not fail.closed:fail.close()

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ('root','manifest','output','prior'):p.add_argument('--'+flag,type=Path,required=True)
    p.add_argument('--workers',type=int,default=20);p.add_argument('--full',action='store_true');run(p.parse_args())
