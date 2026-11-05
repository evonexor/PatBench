#!/usr/bin/env bash
set -euo pipefail
task_script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
task_base="$(dirname -- "$task_script_dir")"
task_root="${1:-/home/open/patent-work/raw/patent-data}"
task_python="${PATENT_PYTHON:-/home/open/miniconda3/envs/molscribe/bin/python}"
if [[ -d "$task_base/runtime/deps" ]]; then
  export PYTHONPATH="$task_base/runtime/deps${PYTHONPATH:+:$PYTHONPATH}"
fi
"$task_python" "$task_script_dir/select_samples.py" --root "$task_root" --output "$task_base/manifest.json" > "$task_base/selection.log"
"$task_python" "$task_script_dir/uspto_sar.py" --root "$task_root" --manifest "$task_base/manifest.json" --output "$task_base/results" --workers "${PATENT_WORKERS:-20}" | tee "$task_base/run.log"
"$task_python" "$task_script_dir/test_pipeline.py" --results "$task_base/results" 2>&1 | tee "$task_base/tests.log"
"$task_python" "$task_script_dir/test_metadata.py" --results "$task_base/results" 2>&1 | tee "$task_base/metadata_tests.log"
"$task_python" "$task_script_dir/test_documents.py" --results "$task_base/results" 2>&1 | tee "$task_base/document_tests.log"
"$task_python" "$task_script_dir/test_xml_source.py" 2>&1 | tee "$task_base/xml_tests.log"
"$task_python" "$task_script_dir/validate_samples.py" --root "$task_root" --results "$task_base/results" | tee "$task_base/validation.log"
