"""Multi-omics integration methods."""

from immunomics.integration.multivi import run_multivi
from immunomics.integration.wnn import run_wnn
from immunomics.integration.mofa import run_mofa
from immunomics.integration.benchmark import compare_methods

__all__ = ["run_multivi", "run_wnn", "run_mofa", "compare_methods"]
