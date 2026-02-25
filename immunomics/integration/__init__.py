#!/usr/bin/env python3
"""
ImmunOmics v0.1.0

Integration subpackage — multi-omics integration methods.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

from .benchmark import compare_methods
from .mofa import run_mofa
from .multivi import run_multivi
from .wnn import run_wnn

__all__ = ["run_multivi", "run_wnn", "run_mofa", "compare_methods"]

# ImmunOmics v0.1.0
# Any usage is subject to this software's license.
