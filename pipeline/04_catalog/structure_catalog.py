"""Every MOL and unmatched CDX/chemical-image source, independent of domain or assays."""
from collections import Counter,defaultdict
import hashlib,json,re,time,zipfile
from pathlib import Path,PurePosixPath
from rdkit import Chem,RDLogger
import pyarrow as pa
import pyarrow.parquet as pq
from catalog_metadata import bibliography
from catalog_annotations import extract as extract_annotations,FIELDS,SEQ_SCHEMA,ACTIVITY_SCHEMA,UNIPROT_SCHEMA
from patent_metadata import METADATA_COLUMNS
from uspto_sar import atomic_json

RDLogger.DisableLog('rdApp.*')
STRUCTURE_COLUMNS=['source_structure_id','Patent','archive','source_archive_sha256_key','structure_file','source_header_offset','source_format','structure_status','SMILES','SMILES_original','InChIKey','query_SMILES','source_sha256','source_bytes_count','raw_structure_bytes','raw_mol_coordinate_dimension','atom_count','fragment_count','compound_id','compound_aliases','alias_association_status','verified_SAR_record_ids','has_verified_SAR','activity_status','drawing_vs_mol_checked','structure_parse_error','related_source_files']+FIELDS
INTEGERS={'source_header_offset','source_bytes_count','atom_count','fragment_count'}
BOOLEANS={'has_verified_SAR','drawing_vs_mol_checked'}
SCHEMA=pa.schema([pa.field(k,pa.binary() if k=='raw_structure_bytes' else pa.int64() if k in INTEGERS else pa.bool_() if k in BOOLEANS else pa.string()) for k in STRUCTURE_COLUMNS])

def selected_sources(z):
    infos=[i for i in z.infolist() if not i.is_dir()]
    stems={str(PurePosixPath(i.filename).with_suffix('')).casefold() for i in infos if i.filename.lower().endswith('.mol')}
    chosen=[(i,'MOL') for i in infos if i.filename.lower().endswith('.mol')]
    fallback_stems=set()
    for i in infos:
        ext=PurePosixPath(i.filename).suffix.lower();stem=str(PurePosixPath(i.filename).with_suffix('')).casefold()
        if ext in ('.cdx','.cdxml','.sdf','.cml','.smi','.smiles','.inchi') and stem not in stems:
            chosen.append((i,ext[1:].upper()));fallback_stems.add(stem)
    for i in infos:
        ext=PurePosixPath(i.filename).suffix.lower();stem=str(PurePosixPath(i.filename).with_suffix('')).casefold()
        if ext in ('.tif','.tiff','.png','.jpg','.jpeg') and stem not in stems|fallback_stems and re.search(r'-C\d+$',PurePosixPath(i.filename).stem,re.I):chosen.append((i,'CHEMICAL_IMAGE_CANDIDATE'))
    return sorted(chosen,key=lambda pair:pair[0].header_offset)

def mol_graph(data):
    r=dict(SMILES='',SMILES_original='',InChIKey='',query_SMILES='',raw_mol_coordinate_dimension='',atom_count=None,fragment_count=None,structure_status='invalid_mol',structure_parse_error='')
    try:
        block=data.decode('utf-8','replace');raw=Chem.MolFromMolBlock(block,sanitize=False,removeHs=False,strictParsing=True)
        if raw is None:raise ValueError('mol_parse_failed')
        r['atom_count']=raw.GetNumAtoms();r['fragment_count']=len(Chem.GetMolFrags(raw))
        conf=raw.GetConformer() if raw.GetNumConformers() else None
        if conf:
            points=[conf.GetAtomPosition(i) for i in range(raw.GetNumAtoms())]
            r['raw_mol_coordinate_dimension']='3D_flag_or_nonzero_z' if conf.Is3D() or any(abs(p.z)>1e-10 for p in points) else '2D' if any(abs(p.x)+abs(p.y)>1e-10 for p in points) else 'zero_coordinates'
        else:r['raw_mol_coordinate_dimension']='absent'
        if any(a.GetAtomicNum()==0 or a.HasQuery() for a in raw.GetAtoms()) or any(b.HasQuery() for b in raw.GetBonds()):
            r['structure_status']='markush_or_query'
            try:r['query_SMILES']=Chem.MolToSmiles(raw,canonical=True,isomericSmiles=True)
            except Exception:pass
            return r
        mol=Chem.MolFromMolBlock(block,sanitize=True,removeHs=True,strictParsing=True)
        if mol is None:raise ValueError('mol_sanitization_failed')
        r['SMILES_original']=Chem.MolToSmiles(mol,canonical=False,isomericSmiles=True);r['SMILES']=Chem.MolToSmiles(mol,canonical=True,isomericSmiles=True)
        r['InChIKey']=Chem.MolToInchiKey(mol)
        r['structure_status']='concrete_multifragment' if r['fragment_count']>1 else 'concrete' if r['InChIKey'] else 'no_inchikey'
    except Exception as exc:r['structure_parse_error']=type(exc).__name__+':'+str(exc)
    return r

def worker(job):
    root,out,prior,spec,code_hash=job;root=Path(root);out=Path(out);prior=Path(prior)
    archive=spec['archive'];key=hashlib.sha256(archive.encode()).hexdigest();path=root/archive;started=time.monotonic()
    folder=out/'checkpoint_cache'/key[:2];folder.mkdir(parents=True,exist_ok=True)
    checkpoint=folder/(key+'.json');part=folder/(key+'.parquet');st=path.stat()
    signature=hashlib.sha256(f'{code_hash}|{st.st_size}|{st.st_mtime_ns}|{archive}'.encode()).hexdigest()
    if checkpoint.is_file() and part.is_file():
        obj=json.loads(checkpoint.read_text())
        if obj.get('signature')==signature and pq.read_metadata(part).num_rows==obj['stats']['structure_rows'] and all(not obj['stats'].get(name+'_rows') or part.with_name(key+'.'+name+'.parquet').is_file() for name in ('sequences','activities','uniprot_mentions')):return {'checkpoint':str(checkpoint),'part':str(part),'cached':True}
    fallback=spec['patent_core'];stats=dict(archive=archive,Patent=fallback,patent_core=fallback,status='pending',structure_rows=0,mol_files=0,non_mol_sources=0,image_files=0,source_read_errors=0,prior_status='not_seen',elapsed_seconds=0)
    old={};alias_map={};sar=defaultdict(list)
    cp=prior/'checkpoint_cache'/key[:2]/'quality'/(key[:20]+'.json')
    if cp.is_file():
        old=json.loads(cp.read_text());stats['prior_status']=old['stats']['status'];alias_map={r['structure_file']:r.get('compound_aliases','') for r in old.get('structures',[])}
        for r in old.get('ready',[]):sar[r['structure_file']].append(r)
    tmp=part.with_suffix('.parquet.tmp');writer=pq.ParquetWriter(tmp,SCHEMA,compression='zstd',use_dictionary=[k for k in STRUCTURE_COLUMNS if k not in ('raw_structure_bytes','source_structure_id','structure_file','source_sha256')]);buffer=[];counts=Counter();formats=Counter();ids=set();seen_sar=set()
    try:
        with zipfile.ZipFile(path) as z:
            infos=z.infolist();stats['mol_files']=sum(i.filename.lower().endswith('.mol') and not i.is_dir() for i in infos);stats['image_files']=sum(i.filename.lower().endswith(('.tif','.tiff','.png','.jpg','.jpeg')) for i in infos)
            try:meta=bibliography(z,archive,fallback)
            except Exception as exc:meta={'Patent':fallback,'patent_core':fallback,'metadata_status':'metadata_error','metadata_format':'unavailable','metadata_review_reason':repr(exc),'screen':'uncertain'}
            relative=str(PurePosixPath('bibliography')/key[:2]/(key+'.json'))
            bpath=out/relative;bpath.parent.mkdir(parents=True,exist_ok=True)
            bib=meta.get('bibliography') or {'Patent':meta['Patent'],'archive':archive,'source_file':meta.get('xml_file',''),'metadata':{k:v for k,v in meta.items() if not k.startswith('_')},'limitations':['Legacy metadata is partial; blank fields are not inferred.']}
            atomic_json(bpath,bib);meta['bibliography_json_path']=relative
            stats.update({k:v for k,v in meta.items() if k!='bibliography' and not k.startswith('_')});patent=stats['Patent']
            annotation=extract_annotations(z,patent,archive,stats.get('xml_file',''),[i.filename for i in infos if i.filename.lower().endswith('.mol')],old,meta.get('_root'))
            stats.update(annotation_errors=annotation['errors'],sequences_rows=len(annotation['sequences']),activities_rows=len(annotation['activities']),uniprot_mentions_rows=len(annotation['uniprot_mentions']))
            activity_by_id={r['source_record_id']:r for r in annotation['activities']}
            sources=selected_sources(z);stats['non_mol_sources']=len(sources)-stats['mol_files']
            related=defaultdict(list)
            for i in infos:related[str(PurePosixPath(i.filename).with_suffix('')).casefold()].append(i.filename)
            for info,fmt in sources:
                rid='USPTO-STRUCTURE:'+key+':'+str(info.header_offset);assert rid not in ids;ids.add(rid)
                r={k:None if k in INTEGERS else False if k in BOOLEANS else b'' if k=='raw_structure_bytes' else '' for k in STRUCTURE_COLUMNS}
                r.update(source_structure_id=rid,Patent=patent,archive=archive,source_archive_sha256_key=key,structure_file=info.filename,source_header_offset=info.header_offset,source_format=fmt,compound_aliases=alias_map.get(info.filename,''),alias_association_status='previous_alias_requires_review' if alias_map.get(info.filename) else 'experimental_number_unresolved',activity_status='no_verified_activity_link',related_source_files=json.dumps(related[str(PurePosixPath(info.filename).with_suffix('')).casefold()]))
                try:
                    data=z.read(info);r.update(raw_structure_bytes=data,source_bytes_count=len(data),source_sha256=hashlib.sha256(data).hexdigest())
                    if fmt=='MOL':r.update(mol_graph(data))
                    else:r.update(structure_status='source_preserved_requires_conversion' if fmt!='CHEMICAL_IMAGE_CANDIDATE' else 'chemical_image_candidate_requires_ocr')
                    authoritative=info.header_offset==z.getinfo(info.filename).header_offset
                    matching=[a for a in sar[info.filename] if authoritative and a['SMILES']==r['SMILES'] and a['InChIKey']==r['InChIKey'] and r['InChIKey'] and a['raw_mol_sha256']==r['source_sha256']]
                    if matching:
                        aliases={a['compound_alias'] for a in matching};r.update(has_verified_SAR=True,activity_status='linked_to_previous_verified_SAR',verified_SAR_record_ids=json.dumps(sorted(a['source_record_id'] for a in matching)),alias_association_status='verified_direct_structure_activity_link',compound_id=patent+':'+next(iter(aliases)) if len(aliases)==1 else '')
                        seen_sar.update(a['source_record_id'] for a in matching)
                        matched_annotations=[activity_by_id[a['source_record_id']] for a in matching]
                        for field in ('UniProt ID','original_range'):
                            r[field]=json.dumps([{'source_record_id':a['source_record_id'],'Target':a['Target'],'value':a[field]} for a in matched_annotations if a.get(field)],ensure_ascii=False)
                        for a in matched_annotations:a['source_structure_id']=rid
                    records=annotation['structures'].get(info.filename,[]) if authoritative else []
                    for field in ('original_alias','compound_IUPAC_name','Sequence','compound_name_status','sequence_ids','sequence_link_status'):
                        values={a[field] for a in records if a.get(field)}
                        if len(values)==1:r[field]=next(iter(values))
                    if records:r['annotation_evidence']=json.dumps([a['annotation_evidence'] for a in records],ensure_ascii=False)
                    if not r['compound_IUPAC_name']:
                        names={a.get('structure_name','') for a in sar[info.filename] if a.get('structure_name')}
                        if len(names)==1:r.update(compound_IUPAC_name=next(iter(names)),compound_name_status='source_reported_name_IUPAC_unverified')
                except Exception as exc:r.update(structure_status='source_read_error',structure_parse_error=repr(exc));stats['source_read_errors']+=1
                counts[r['structure_status']]+=1;formats[fmt]+=1;buffer.append(r)
                if len(buffer)>=512:writer.write_table(pa.Table.from_pylist(buffer,schema=SCHEMA));buffer=[]
            if buffer:writer.write_table(pa.Table.from_pylist(buffer,schema=SCHEMA))
            assert len(sources)==stats['mol_files']+stats['non_mol_sources']
        writer.close();tmp.replace(part)
        for name,rows,schema in [('sequences',annotation['sequences'],SEQ_SCHEMA),('activities',annotation['activities'],ACTIVITY_SCHEMA),('uniprot_mentions',annotation['uniprot_mentions'],UNIPROT_SCHEMA)]:
            if rows:
                target=part.with_name(key+'.'+name+'.parquet');temporary=target.with_suffix('.parquet.tmp')
                clean=[{k:(None if schema.field(k).type==pa.float64() and row.get(k)=='' else row.get(k)) for k in schema.names} for row in rows]
                pq.write_table(pa.Table.from_pylist(clean,schema=schema),temporary,compression='zstd');temporary.replace(target)
        stats.update(status='source_read_errors' if stats['source_read_errors'] else 'extracted' if sources else 'no_identified_chemical_source',structure_rows=len(sources),structure_status_counts=dict(counts),source_format_counts=dict(formats),previous_verified_SAR_expected=len(old.get('ready',[])),previous_verified_SAR_linked=len(seen_sar),elapsed_seconds=round(time.monotonic()-started,3))
        atomic_json(checkpoint,{'signature':signature,'stats':stats})
        return {'checkpoint':str(checkpoint),'part':str(part),'cached':False}
    except BaseException:
        writer.close();raise
