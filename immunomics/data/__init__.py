#!/usr/bin/env python3
"""
ImmunOmics v0.1.0

Data subpackage — loaders and preprocessors.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

from .datasets import load_multiome_pbmc
from .preprocess_atac import preprocess_atac
from .preprocess_rna import preprocess_rna

__all__ = ["load_multiome_pbmc", "preprocess_rna", "preprocess_atac"]

# ImmunOmics v0.1.0
# Any usage is subject to this software's license.
