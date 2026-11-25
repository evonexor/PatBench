"""PatBench data-processing and within-assay evaluation toolkit."""
from pathlib import Path

__version__ = "1.0.0"

# Editable source checkout layout: pipeline/, provenance/ and configs/ sit beside src/
REPO = Path(__file__).resolve().parents[2]
PIPELINE = REPO / "pipeline"
PROVENANCE = REPO / "provenance"
CONFIGS = REPO / "configs"

from patbench.evaluate import evaluate_predictions, unpack_figshare_release  # noqa: E402
from patbench.verify import verify_release  # noqa: E402

__all__ = [
    "__version__",
    "REPO",
    "PIPELINE",
    "PROVENANCE",
    "CONFIGS",
    "evaluate_predictions",
    "unpack_figshare_release",
    "verify_release",
]
