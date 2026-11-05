#!/usr/bin/env python3
"""Copy the NAS patent dataset onto the open server's dedicated 8 TB SSD."""
import concurrent.futures as cf
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

SOURCE = Path('/home/open/nas/syn-data/patent-data')
DEST = Path('/mnt/ssd8tb/open-home/patent-work/raw/patent-data')
RUN = Path(sys.argv[1])
RUN.mkdir(parents=True, exist_ok=True)
STATE = {'source': str(SOURCE), 'destination': str(DEST), 'started_at': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'phase': 'preflight'}
BASE = ['rsync', '-a', '--no-owner', '--no-group', '--omit-dir-times']


def save(**fields):
    STATE.update(fields)
    STATE['updated_at'] = time.strftime('%Y-%m-%dT%H:%M:%S%z')
    tmp = RUN / 'status.json.tmp'
    tmp.write_text(json.dumps(STATE, ensure_ascii=False, indent=2) + '\n')
    tmp.replace(RUN / 'status.json')


def copy_one(task):
    name, src, dst, excludes = task
    dst.mkdir(parents=True, exist_ok=True)
    cmd = BASE + ['--partial-dir=.rsync-partial', '--info=progress2', '--stats']
    for pattern in excludes:
        cmd += ['--exclude', pattern]
    cmd += [str(src) + '/', str(dst) + '/']
    log = RUN / (name + '.rsync.log')
    for attempt in range(1, 4):
        with log.open('ab') as out:
            out.write(('\nAttempt %d\n' % attempt).encode())
            out.flush()
            proc = subprocess.run(cmd, stdout=out, stderr=subprocess.STDOUT)
        if proc.returncode == 0:
            return {'task': name, 'exit_code': 0, 'attempts': attempt}
        if attempt < 3:
            time.sleep(5)
    return {'task': name, 'exit_code': proc.returncode, 'attempts': 3, 'log': str(log)}


def safe_file(relative):
    p = Path(relative)
    if p.is_absolute() or '..' in p.parts:
        raise ValueError('Unsafe manifest path: ' + str(p))
    result = DEST / p
    if result.is_symlink() or not result.is_file():
        raise ValueError('Missing or non-regular manifest file: ' + str(p))
    return result


def check_hash(item):
    relative, algorithm, expected, expected_size = item
    try:
        path = safe_file(relative)
        size = path.stat().st_size
        if expected_size is not None and size != expected_size:
            return {'path': relative, 'error': 'size mismatch', 'actual_size': size, 'expected_size': expected_size}
        if algorithm == 'sha256':
            digest = hashlib.sha256()
        elif algorithm == 'git-sha1':
            digest = hashlib.sha1()
            digest.update(('blob %d\0' % size).encode())
        else:
            raise ValueError('Unsupported hash algorithm: ' + algorithm)
        with path.open('rb') as f:
            for block in iter(lambda: f.read(4 * 1024 * 1024), b''):
                digest.update(block)
        actual = digest.hexdigest()
        if actual != expected:
            return {'path': relative, 'error': 'hash mismatch', 'expected': expected, 'actual': actual}
        return None
    except Exception as exc:
        return {'path': relative, 'error': str(exc)}


def main():
    source_fs = subprocess.check_output(['findmnt', '-n', '-o', 'FSTYPE', '-T', str(SOURCE)], text=True).strip()
    disk = subprocess.check_output(['findmnt', '-n', '-o', 'SOURCE,FSTYPE', '-T', str(DEST.parent)], text=True).strip()
    if source_fs != 'cifs' or disk.split() != ['/dev/nvme0n1', 'ext4']:
        raise RuntimeError('Unexpected source/destination mount: ' + source_fs + ' / ' + disk)
    space = os.statvfs('/mnt/ssd8tb')
    if space.f_bavail * space.f_frsize < 1024 ** 4:
        raise RuntimeError('Less than 1 TiB free on destination SSD')
    DEST.mkdir(parents=True, exist_ok=True)
    tasks = [('root-files', SOURCE, DEST, ['/*/'])]
    for entry in sorted(os.scandir(SOURCE), key=lambda e: e.name):
        if not entry.is_dir(follow_symlinks=False):
            continue
        src, dst = Path(entry.path), DEST / entry.name
        if entry.name != 'uspto-mol' or not (src / 'data').is_dir():
            tasks.append((entry.name, src, dst, []))
            continue
        tasks.append(('uspto-headers', src, dst, ['/data/***']))
        tasks.append(('uspto-data-files', src / 'data', dst / 'data', ['/*/']))
        for year in sorted(os.scandir(src / 'data'), key=lambda e: e.name, reverse=True):
            if year.is_dir(follow_symlinks=False):
                year_src, year_dst = Path(year.path), dst / 'data' / year.name
                tasks.append(('uspto-' + year.name + '-files', year_src, year_dst, ['/*/']))
                for week in sorted(os.scandir(year_src), key=lambda e: e.name):
                    if week.is_dir(follow_symlinks=False):
                        tasks.append(('uspto-' + year.name + '-' + week.name, Path(week.path), year_dst / week.name, []))
    (RUN / 'tasks.json').write_text(json.dumps([{'name': n, 'source': str(s), 'destination': str(d), 'excludes': x} for n, s, d, x in tasks], indent=2) + '\n')
    save(phase='copying', tasks_total=len(tasks), tasks_completed=0, completed=[], pending=[t[0] for t in tasks])
    finished = []
    with cf.ThreadPoolExecutor(max_workers=32) as pool:
        futures = {pool.submit(copy_one, task): task[0] for task in tasks}
        while futures:
            done, _ = cf.wait(futures, timeout=10, return_when=cf.FIRST_COMPLETED)
            for future in done:
                finished.append(future.result())
                del futures[future]
            space = os.statvfs('/mnt/ssd8tb')
            save(tasks_completed=len(finished), completed=finished, pending=list(futures.values()), ssd_used_bytes=(space.f_blocks - space.f_bfree) * space.f_frsize)
    failures = [r for r in finished if r['exit_code'] != 0]
    if failures:
        raise RuntimeError('Copy tasks failed: ' + json.dumps(failures))

    save(phase='verifying-file-list')
    with (RUN / 'file-list-differences.log').open('wb') as out, (RUN / 'file-list-errors.log').open('wb') as err:
        result = subprocess.run(BASE + ['--dry-run', '--delete', '--itemize-changes', str(SOURCE) + '/', str(DEST) + '/'], stdout=out, stderr=err)
    differences = (RUN / 'file-list-differences.log').read_text(errors='replace')
    if result.returncode != 0 or differences.strip():
        save(file_list_exit_code=result.returncode, file_list_differences=differences[:12000])
        raise RuntimeError('Source and copied file list differ; see file-list logs')
    save(file_list_verified=True, phase='verifying-checksums')
    checks = []
    for line in (DEST / 'SHA256SUMS').read_text().splitlines():
        if line.strip():
            digest, relative = line.split(maxsplit=1)
            checks.append((relative.lstrip('*'), 'sha256', digest, None))
    base_checks = len(checks)
    manifest = json.loads((DEST / '_metadata/uspto_mol_verified/manifest.json').read_text())
    for item in manifest['files']:
        checks.append(('uspto-mol/' + item['path'], item['algorithm'], item['digest'], item['size']))
    save(checksum_total=len(checks), checksum_checked=0, checksum_failed=0, base_checksum_files=base_checks, uspto_checksum_files=len(manifest['files']))
    failures = []
    failed_count, last_save = 0, time.monotonic()
    with cf.ThreadPoolExecutor(max_workers=4) as pool:
        for count, failure in enumerate(pool.map(check_hash, checks), 1):
            if failure:
                failed_count += 1
                if len(failures) < 100:
                    failures.append(failure)
            if time.monotonic() - last_save >= 10 or count == len(checks):
                save(checksum_checked=count, checksum_failed=failed_count, checksum_errors=failures)
                last_save = time.monotonic()
    if failed_count:
        raise RuntimeError('%d checksum failures' % failed_count)
    regular_files = total_bytes = symlinks = 0
    for root, dirs, files in os.walk(DEST, followlinks=False):
        for name in dirs:
            if (Path(root) / name).is_symlink():
                symlinks += 1
        for name in files:
            p = Path(root) / name
            if p.is_symlink():
                symlinks += 1
            else:
                regular_files += 1
                total_bytes += p.stat().st_size
    save(phase='complete', completed_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'), regular_files=regular_files, symlinks=symlinks, total_file_bytes=total_bytes)
    (RUN / 'verification-report.json').write_text(json.dumps(STATE, ensure_ascii=False, indent=2) + '\n')


if __name__ == '__main__':
    try:
        main()
    except BaseException as exc:
        save(phase='failed', error=str(exc))
        raise
