"""Bounded row groups and files; preserve Arrow schemas and binary source bytes."""
from pathlib import Path
import struct
import pyarrow as pa
import pyarrow.parquet as pq

STATS_COLUMNS = {
    'Patent', 'patent_core', 'InChIKey', 'source_structure_id', 'source_record_id',
    'source_format', 'structure_status', 'has_verified_SAR', 'is_formal_SAR',
    'Unit', 'Activity Type', 'Relation', 'Value', 'value_nM',
    'document_publication_date', 'grant_date', 'application_date',
    'metadata_status', 'status', 'prior_status', 'sequence_status',
}


class ShardedWriter:
    def __init__(self, path, schema, row_group_rows=32768,
                 buffer_bytes=64 * 2**20, max_row_groups=32,
                 shard_rows=262144, shard_bytes=512 * 2**20):
        self.path = Path(path)
        self.path.mkdir(parents=True, exist_ok=False)
        self.schema = schema
        self.row_group_rows = row_group_rows
        self.buffer_bytes = buffer_bytes
        self.max_row_groups = max_row_groups
        self.shard_rows = shard_rows
        self.shard_bytes = shard_bytes
        self.buffer = []
        self.buffer_rows = self.buffer_size = 0
        self.writer = None
        self.file_rows = self.file_groups = 0
        self.total_rows = 0
        self.files = []
        self.closed = False

    def write_batch(self, batch):
        self.write_table(pa.Table.from_batches([batch]))

    def write_table(self, table):
        assert not self.closed and table.schema.equals(self.schema)
        offset = 0
        while offset < table.num_rows:
            rows = min(table.num_rows-offset, self.row_group_rows-self.buffer_rows)
            piece = table.slice(offset, rows)
            self.buffer.append(piece)
            self.buffer_rows += rows
            self.buffer_size += piece.nbytes
            offset += rows
            if self.buffer_rows >= self.row_group_rows or self.buffer_size >= self.buffer_bytes:
                self._flush()

    def _open(self):
        self.current = self.path / ('part-%05d.parquet' % len(self.files))
        self.writer = pq.ParquetWriter(
            self.current, self.schema, compression='zstd',
            write_statistics=[k for k in self.schema.names if k in STATS_COLUMNS],
            use_dictionary=[k for k in self.schema.names if k not in (
                'raw_structure_bytes', 'Sequence', 'sequence_original',
                'source_structure_id', 'source_record_id', 'source_sha256',
                'SMILES', 'SMILES_original', 'query_SMILES',
            )],
        )

    def _flush(self):
        if not self.buffer_rows:
            return
        if self.writer is None:
            self._open()
        table = pa.concat_tables(self.buffer)
        self.writer.write_table(table, row_group_size=self.row_group_rows)
        self.file_rows += table.num_rows
        self.total_rows += table.num_rows
        self.file_groups += 1
        self.buffer.clear()
        self.buffer_rows = self.buffer_size = 0
        if (self.file_groups >= self.max_row_groups or
                self.file_rows >= self.shard_rows or
                self.current.stat().st_size >= self.shard_bytes):
            self._close_file()

    def _close_file(self):
        if self.writer is None:
            return
        self.writer.close()
        self.writer = None
        metadata = pq.read_metadata(self.current)
        assert metadata.num_rows == self.file_rows
        assert metadata.num_row_groups <= self.max_row_groups
        with self.current.open('rb') as handle:
            handle.seek(-8, 2)
            trailer = handle.read(8)
        assert trailer[4:] == b'PAR1'
        footer_bytes = struct.unpack('<I', trailer[:4])[0]
        assert footer_bytes < 64 * 2**20, 'Unexpectedly large Parquet footer'
        self.files.append(dict(
            file=self.current.name, rows=self.file_rows,
            row_groups=metadata.num_row_groups, bytes=self.current.stat().st_size,
            footer_bytes=footer_bytes,
        ))
        self.file_rows = self.file_groups = 0

    def close(self):
        if self.closed:
            return
        self._flush()
        if self.writer is None and not self.files:
            self._open()
        self._close_file()
        self.closed = True
        assert sum(f['rows'] for f in self.files) == self.total_rows


def dataset_row_count(path):
    files = sorted(Path(path).glob('part-*.parquet'))
    assert files, 'Dataset has no Parquet files'
    return sum(pq.read_metadata(p).num_rows for p in files)
