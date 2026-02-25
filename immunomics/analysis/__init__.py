#!/usr/bin/env python3
"""
ImmunOmics v0.1.0

Analysis subpackage — downstream multi-omics analyses.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

from .differential import differential_genes, differential_peaks
from .peak_gene_links import link_peaks_to_genes
from .tf_activity import infer_tf_activity

__all__ = [
    "link_peaks_to_genes",
    "infer_tf_activity",
    "differential_peaks",
    "differential_genes",
]

# ImmunOmics v0.1.0
# Any usage is subject to this software's license.
