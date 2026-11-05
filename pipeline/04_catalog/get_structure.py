"""Export the exact source MOL/CDX/image from its stable structure ID."""
import argparse,hashlib
from pathlib import Path
import pyarrow.parquet as pq
import pyarrow.compute as pc

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,required=True);p.add_argument('--id',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    prefix,key,offset=a.id.split(':');assert prefix=='USPTO-STRUCTURE' and len(key)==64
    path=a.results/'checkpoint_cache'/key[:2]/(key+'.parquet')
    table=pq.read_table(path,columns=['source_structure_id','raw_structure_bytes','source_sha256','structure_file','source_format'])
    found=table.filter(pc.equal(table['source_structure_id'],a.id)).to_pylist();assert len(found)==1
    row=found[0];data=row['raw_structure_bytes'];assert data and hashlib.sha256(data).hexdigest()==row['source_sha256']
    a.output.parent.mkdir(parents=True,exist_ok=True);assert not a.output.exists(),'Output already exists';a.output.write_bytes(data)
    print(str(a.output));print(row['source_format'],row['source_sha256'])
