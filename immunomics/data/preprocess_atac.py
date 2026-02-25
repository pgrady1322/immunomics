#!/usr/bin/env python3
"""
ImmunOmics v0.1.0

scATAC-seq preprocessing (TF-IDF, LSI, peak filtering).

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import logging

import anndata as ad
import numpy as np
import scanpy as sc

logger = logging.getLogger(__name__)


def preprocess_atac(
    adata: ad.AnnData,
    min_peaks: int = 500,
    min_cells: int = 5,
    n_top_peaks: int = 50000,
    n_components: int = 50,
    copy: bool = True,
) -> ad.AnnData:
    """
    Preprocess scATAC-seq peak-by-cell matrix.

    Steps:
    1. QC filtering (cells with too few/many peaks)
    2. Peak selection (by variability)
    3. TF-IDF normalization
    4. LSI (Latent Semantic Indexing) dimensionality reduction

    Parameters
    ----------
    adata
        AnnData with raw ATAC counts (peaks x cells)
    min_peaks
        Minimum number of peaks per cell
    min_cells
        Minimum number of cells per peak
    n_top_peaks
        Number of top variable peaks to retain
    n_components
        Number of LSI components
    copy
        Whether to operate on a copy

    Returns
    -------
    Preprocessed AnnData with LSI coordinates in obsm['X_lsi']
    """
    if copy:
        adata = adata.copy()

    logger.info(f"ATAC preprocessing: {adata.n_obs} cells, {adata.n_vars} peaks")

    # Store raw counts
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    # QC
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=min_peaks)  # min_genes works for peaks too
    sc.pp.filter_genes(adata, min_cells=min_cells)
    logger.info(f"QC: {n_before} -> {adata.n_obs} cells, {adata.n_vars} peaks")

    # Select top variable peaks
    if adata.n_vars > n_top_peaks:
        # Use mean accessibility as proxy for variability
        import scipy.sparse as sp

        X = adata.X
        if sp.issparse(X):
            peak_means = np.array(X.mean(axis=0)).flatten()
            peak_vars = np.array(X.power(2).mean(axis=0)).flatten() - peak_means**2
        else:
            peak_means = X.mean(axis=0)
            peak_vars = X.var(axis=0)

        # Select peaks with highest variance-to-mean ratio (overdispersed)
        ratios = peak_vars / (peak_means + 1e-8)
        top_peaks = np.argsort(ratios)[-n_top_peaks:]
        adata = adata[:, top_peaks].copy()
        logger.info(f"Selected top {n_top_peaks} variable peaks")

    # TF-IDF normalization
    adata = _tfidf_normalize(adata)

    # LSI (SVD on TF-IDF matrix)
    adata = _run_lsi(adata, n_components=n_components)

    logger.info(f"ATAC preprocessing complete: {adata.n_obs} cells, {n_components} LSI components")
    return adata


def _tfidf_normalize(adata: ad.AnnData) -> ad.AnnData:
    """
    TF-IDF normalization for scATAC-seq data.

    Term Frequency: normalized by total fragments per cell
    Inverse Document Frequency: log(1 + n_cells / n_cells_with_peak)
    """
    import scipy.sparse as sp

    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)

    # Term frequency (per-cell normalization)
    tf = X.multiply(1.0 / X.sum(axis=1))

    # Inverse document frequency
    n_cells = X.shape[0]
    n_cells_per_peak = np.array((X > 0).sum(axis=0)).flatten()
    idf = np.log1p(n_cells / (n_cells_per_peak + 1))

    # TF-IDF
    tfidf = tf.multiply(idf).tocsr()

    adata.X = tfidf
    adata.layers["tfidf"] = tfidf.copy()

    logger.info("Applied TF-IDF normalization")
    return adata


def _run_lsi(adata: ad.AnnData, n_components: int = 50) -> ad.AnnData:
    """
    Latent Semantic Indexing (truncated SVD on TF-IDF matrix).

    The first component typically correlates with sequencing depth
    and is excluded from downstream analysis.
    """
    from sklearn.decomposition import TruncatedSVD

    X = adata.X
    X_dense = X

    svd = TruncatedSVD(n_components=n_components + 1, random_state=42)
    lsi = svd.fit_transform(X_dense)

    # Drop first component (correlated with depth)
    adata.obsm["X_lsi"] = lsi[:, 1:]
    adata.uns["lsi_variance_ratio"] = svd.explained_variance_ratio_[1:]

    logger.info(
        f"LSI: {n_components} components, "
        f"explained variance: {svd.explained_variance_ratio_[1:].sum():.3f}"
    )

    return adata


def parse_peak_coordinates(adata: ad.AnnData) -> ad.AnnData:
    """
    Parse peak names (chr:start-end) into structured annotation.

    Adds var columns: chrom, start, end
    """
    import re

    chroms, starts, ends = [], [], []
    for peak_name in adata.var_names:
        match = re.match(r"(chr\w+)[:\-](\d+)[:\-](\d+)", peak_name)
        if match:
            chroms.append(match.group(1))
            starts.append(int(match.group(2)))
            ends.append(int(match.group(3)))
        else:
            chroms.append("unknown")
            starts.append(0)
            ends.append(0)

    adata.var["chrom"] = chroms
    adata.var["start"] = starts
    adata.var["end"] = ends

    return adata


# ImmunOmics v0.1.0
# Any usage is subject to this software's license.
