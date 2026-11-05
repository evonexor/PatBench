"""PatBench data-processing pipeline (see README.md)."""
from pathlib import Path

__version__ = "1.0.0.dev0"

# the code runs from a source checkout (`pip install -e .`): pipeline/ and provenance/ sit beside src/
REPO = Path(__file__).resolve().parents[2]
PIPELINE = REPO / "pipeline"
PROVENANCE = REPO / "provenance"
CONFIGS = REPO / "configs"
