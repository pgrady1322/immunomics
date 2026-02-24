#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ImmunOmics v0.1.0

MOFA+ multi-omics factor analysis via muon.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import logging
from typing import Optional

import numpy as np
import anndata as ad

logger = logging.getLogger(__name__)


def run_mofa(
    adata_rna: ad.AnnData,
    adata_atac: ad.AnnData,
    n_factors: int = 15,
    max_iterations: int = 1000,
    convergence_mode: str = "medium",
    compute_umap: bool = True,
    n_neighbors: int = 20,
) -> ad.AnnData:
    """
    Integrate RNA + ATAC using MOFA+.

    Parameters
    ----------
    adata_rna
        Preprocessed RNA AnnData
    adata_atac
        Preprocessed ATAC AnnData
    n_factors
        Number of latent factors to learn
    max_iterations
        Maximum MOFA training iterations
    convergence_mode
        'fast', 'medium', or 'slow'
    compute_umap
        Whether to compute UMAP on factor space
    n_neighbors
        Number of neighbors for UMAP

    Returns
    -------
    Integrated AnnData with MOFA factors in obsm['X_mofa']
    """
    import muon as mu
    import scanpy as sc

    logger.info(f"Running MOFA+ integration with {n_factors} factors...")

    # Ensure shared barcodes
    shared = adata_rna.obs_names.intersection(adata_atac.obs_names)
    logger.info(f"Shared barcodes: {len(shared)}")

    adata_rna = adata_rna[shared].copy()
    adata_atac = adata_atac[shared].copy()

    # Create MuData
    mdata = mu.MuData({"rna": adata_rna, "atac": adata_atac})

    # Run MOFA
    mu.tl.mofa(
        mdata,
        n_factors=n_factors,
        convergence_mode=convergence_mode,
        maxiter=max_iterations,
        outfile=None,  # Don't save intermediate HDF5
        use_obs="intersection",
    )

    logger.info(f"MOFA+ training complete: {n_factors} factors")

    # Build integrated AnnData
    adata_integrated = adata_rna.copy()

    if "X_mofa" in mdata.obsm:
        adata_integrated.obsm["X_mofa"] = mdata.obsm["X_mofa"]
    elif "X_mofa" in mdata.mod["rna"].obsm:
        adata_integrated.obsm["X_mofa"] = mdata.mod["rna"].obsm["X_mofa"]

    adata_integrated.obsm["X_lsi"] = adata_atac.obsm.get("X_lsi", None)
    adata_integrated.uns["integration_method"] = "MOFA+"
    adata_integrated.uns["mofa_variance_explained"] = mdata.uns.get(
        "mofa_variance_explained", {}
    )

    # Compute UMAP on MOFA factors
    if compute_umap and "X_mofa" in adata_integrated.obsm:
        sc.pp.neighbors(adata_integrated, use_rep="X_mofa", n_neighbors=n_neighbors)
        sc.tl.umap(adata_integrated)
        sc.tl.leiden(adata_integrated, resolution=1.0)

    logger.info("MOFA+ integration complete")
    return adata_integrated

# ImmunOmics v0.1.0
# Any usage is subject to this software's license.
