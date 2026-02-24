#!/usr/bin/env python3
"""
ImmunOmics v0.1.0

Differential peak and gene expression analysis.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import logging

import anndata as ad
import pandas as pd
import scanpy as sc

logger = logging.getLogger(__name__)


def differential_genes(
    adata: ad.AnnData,
    groupby: str = "cell_type",
    method: str = "wilcoxon",
    n_genes: int = 100,
) -> pd.DataFrame:
    """
    Find differentially expressed genes per cell type.

    Parameters
    ----------
    adata
        RNA AnnData
    groupby
        Column to group by (cell type or condition)
    method
        Statistical test: 'wilcoxon', 't-test', 'logreg'
    n_genes
        Number of top DEGs per group

    Returns
    -------
    DataFrame with DE results
    """
    logger.info(f"Finding DE genes grouped by '{groupby}' using {method}")

    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method, n_genes=n_genes)

    result = sc.get.rank_genes_groups_df(adata, group=None)
    logger.info(f"Found {len(result)} DE results across {adata.obs[groupby].nunique()} groups")

    return result


def differential_peaks(
    adata_atac: ad.AnnData,
    groupby: str = "cell_type",
    method: str = "wilcoxon",
    n_peaks: int = 500,
) -> pd.DataFrame:
    """
    Find differentially accessible peaks per cell type.

    Parameters
    ----------
    adata_atac
        ATAC AnnData (binarized or TF-IDF normalized)
    groupby
        Column to group by
    method
        Statistical test
    n_peaks
        Number of top DA peaks per group

    Returns
    -------
    DataFrame with DA results
    """
    logger.info(f"Finding DA peaks grouped by '{groupby}'")

    sc.tl.rank_genes_groups(adata_atac, groupby=groupby, method=method, n_genes=n_peaks)

    result = sc.get.rank_genes_groups_df(adata_atac, group=None)
    logger.info(f"Found {len(result)} DA results across {adata_atac.obs[groupby].nunique()} groups")

    return result


# ImmunOmics v0.1.0
# Any usage is subject to this software's license.
