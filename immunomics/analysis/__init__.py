"""Downstream analysis modules for multi-omics immune data."""

from immunomics.analysis.peak_gene_links import link_peaks_to_genes
from immunomics.analysis.tf_activity import infer_tf_activity
from immunomics.analysis.differential import differential_peaks, differential_genes

__all__ = [
    "link_peaks_to_genes",
    "infer_tf_activity",
    "differential_peaks",
    "differential_genes",
]
