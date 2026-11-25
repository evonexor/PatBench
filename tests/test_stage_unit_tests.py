"""Run the unit tests that ship with the extraction, recovery and catalog code (they need RDKit)."""
import os
import subprocess
import sys

import pytest

from patbench import PIPELINE

pytest.importorskip("rdkit")

STAGES = {
    "02_extract": ["test_pipeline", "test_documents", "test_metadata", "test_xml_source", "test_alignment"],
    "03_recovery": ["test_recovery"],
    "04_catalog": ["test_catalog"],
}
CASES = [(d, m) for d, ms in STAGES.items() for m in ms]


@pytest.mark.rdkit
@pytest.mark.parametrize("directory,module", CASES)
def test_stage_unit_tests(directory, module):
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join([str(PIPELINE / "02_extract"), env.get("PYTHONPATH", "")])
    r = subprocess.run([sys.executable, "-m", "unittest", module], cwd=PIPELINE / directory, env=env,
                       capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-3000:]
