"""Launch only the validated recovery in its own output directory."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone

B=Path('/home/open/patent-work/outputs/uspto-sar-recovery-20261006')
UNIT='uspto-sar-recovery-20261006.service'
PY='/home/open/miniconda3/envs/molscribe/bin/python'
DEPS='/home/open/patent-work/outputs/uspto-sar-demo-20261003/runtime/deps'
SCRIPT=str(B/'code/run_recovery.py')
sha=hashlib.sha256(b''.join(str(p.relative_to(B/'code')).encode()+b'\0'+p.read_bytes() for p in sorted((B/'code').rglob('*.py')))).hexdigest()
assert sha=='81fd68912e6f5bb92bb3629ba3fbc8f1ba68e40ebd9431f91598d1b8a3b89fd2'
completion=json.loads((B/'preflight-signed/completion.json').read_text())
validation=json.loads((B/'preflight-signed/results/validation_report.json').read_text())
report=json.loads((B/'preflight-signed/results/run_report.json').read_text())
assert completion['code_sha256']==sha and completion['completed_archives']==completion['total_archives']==50
assert report['archive_errors']==0 and report['new_numeric_SAR_rows']==95 and report['verified_extended_activity_rows']==26
for field in ('complete_archive_accounting','all_checkpoint_source_signatures_verified','baseline_all_columns_exactly_unchanged','all_shards_readable_and_counted','bounded_parquet_row_groups','merged_row_count_correct'):
    assert validation[field] is True
manifest=json.loads((B/'manifest.json').read_text());assert len(manifest['samples'])==345940
assert manifest['counts']=={'review_only':25703,'previously_skipped':236178,'remaining_processed':84059}
if (B/'completion.json').exists():raise RuntimeError('Full recovery completion already exists; inspect before any restart.')
unit=subprocess.run(['systemctl','--user','show',UNIT,'--property=LoadState,ActiveState,SubState,MainPID'],capture_output=True,text=True)
props=dict(line.split('=',1) for line in unit.stdout.splitlines() if '=' in line)
if props.get('ActiveState') in ('active','activating') or int(props.get('MainPID','0')):raise RuntimeError('Same service already running')
for directory in Path('/proc').iterdir():
    if not directory.name.isdigit() or int(directory.name)==os.getpid():continue
    try:argv=(directory/'cmdline').read_bytes().split(b'\0')
    except (OSError,PermissionError):continue
    if SCRIPT.encode() in argv:raise RuntimeError('Same extraction process already exists: '+directory.name)
tests=subprocess.run([PY,'-m','unittest','-v','test_recovery'],cwd=B/'code',capture_output=True,text=True)
(B/'audit/unit_tests.txt').write_text(tests.stdout+tests.stderr)
assert tests.returncode==0 and 'Ran 18 tests' in tests.stderr
command=['systemd-run','--user','--unit='+UNIT,'--description=USPTO prioritized SAR recovery',
         '--working-directory='+str(B),'--property=Type=exec','--property=Restart=on-failure',
         '--property=RestartSec=30','--property=StartLimitIntervalSec=3600',
         '--property=StartLimitBurst=3','--property=Nice=5','--property=TasksMax=512',
         '--property=StandardOutput=append:'+str(B/'run.log'),
         '--property=StandardError=append:'+str(B/'run.log'),
         '--setenv=PYTHONPATH='+DEPS,'--setenv=OMP_NUM_THREADS=1',
         '--setenv=OPENBLAS_NUM_THREADS=1','--setenv=MKL_NUM_THREADS=1',
         '--setenv=NUMEXPR_NUM_THREADS=1',PY,SCRIPT,'--manifest',str(B/'manifest.json'),
         '--output',str(B),'--workers','20']
config=dict(unit=UNIT,systemd_run_argv=command,code_sha256=sha,
            manifest_sha256=hashlib.sha256((B/'manifest.json').read_bytes()).hexdigest(),
            preflight_manifest_sha256=completion['manifest_sha256'],baseline_sha256=validation['baseline_sha256'],
            preflight_output=str(B/'preflight-signed'),workers=20,
            launch_at=datetime.now(timezone.utc).isoformat())
(B/'audit/launch_config.json').write_text(json.dumps(config,ensure_ascii=False,indent=2))
result=subprocess.run(command,capture_output=True,text=True,check=True)
config['launch_result']=result.stdout+result.stderr
(B/'audit/launch_config.json').write_text(json.dumps(config,ensure_ascii=False,indent=2))
print(config['launch_result']);print(json.dumps(config,ensure_ascii=False))
