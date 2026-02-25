#!/usr/bin/env python3
"""
ImmunOmics v0.1.0

scRNA-seq preprocessing pipeline.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import logging

import anndata as ad
import scanpy as sc

logger = logging.getLogger(__name__)


def preprocess_rna(
    adata: ad.AnnData,
    min_genes: int = 200,
    min_cells: int = 3,
    max_pct_mito: float = 20.0,
    n_top_genes: int = 3000,
    target_sum: float = 1e4,
    n_pcs: int = 50,
    copy: bool = True,
) -> ad.AnnData:
    """
    Preprocess scRNA-seq data for multi-omics integration.

    Steps:
    1. QC filtering
    2. Normalization + log-transform
    3. HVG selection
    4. Scaling
    5. PCA

    Parameters
    ----------
    adata
        AnnData with raw RNA counts
    min_genes, min_cells
        QC filtering thresholds
    max_pct_mito
        Maximum mitochondrial percentage
    n_top_genes
        Number of highly variable genes
    target_sum
        Library size normalization target
    n_pcs
        Number of principal components
    copy
        Whether to operate on a copy

    Returns
    -------
    Preprocessed AnnData
    """
    if copy:
        adata = adata.copy()

    logger.info(f"RNA preprocessing: {adata.n_obs} cells, {adata.n_vars} genes")

    # Store raw counts
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    # QC
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata = adata[adata.obs["pct_counts_mt"] < max_pct_mito].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells)
    logger.info(f"QC: {n_before} -> {adata.n_obs} cells")

    # Normalize
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.raw = adata

    # HVG — use seurat_v3 if scikit-misc is available, else fall back to seurat
    try:
        sc.pp.highly_variable_genes(
            adata, n_top_genes=n_top_genes, flavor="seurat_v3", layer="counts"
        )
    except ImportError:
        logger.warning(
            "scikit-misc not installed; falling back to flavor='seurat'. "
            "Install scikit-misc for seurat_v3 variance-stabilization."
        )
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="seurat")
    logger.info(f"Selected {adata.var['highly_variable'].sum()} HVGs")

    # Scale and PCA (on HVGs only for PCA, but keep all genes)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack")

    logger.info(f"RNA preprocessing complete: {adata.n_obs} cells, {n_pcs} PCs")
    return adata


# ImmunOmics v0.1.0
# Any usage is subject to this software's license.
