#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ImmunOmics v0.1.0

MultiVI joint RNA+ATAC integration via scvi-tools.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import logging
from typing import Optional, Dict, Any

import numpy as np
import anndata as ad

logger = logging.getLogger(__name__)


def run_multivi(
    adata_rna: ad.AnnData,
    adata_atac: ad.AnnData,
    n_latent: int = 20,
    n_layers: int = 2,
    n_hidden: int = 256,
    max_epochs: int = 100,
    batch_size: int = 256,
    batch_key: Optional[str] = None,
    n_neighbors: int = 20,
    compute_umap: bool = True,
) -> ad.AnnData:
    """
    Integrate RNA + ATAC using MultiVI.

    Parameters
    ----------
    adata_rna
        Preprocessed scRNA-seq AnnData (with raw counts in layers['counts'])
    adata_atac
        Preprocessed scATAC-seq AnnData (with raw counts in layers['counts'])
    n_latent
        Latent space dimensionality
    n_layers
        Number of hidden layers in encoder/decoder
    n_hidden
        Hidden layer width
    max_epochs
        Maximum training epochs
    batch_size
        Training batch size
    batch_key
        Batch covariate for batch correction
    n_neighbors
        Neighbors for UMAP computation
    compute_umap
        Whether to compute UMAP on latent space

    Returns
    -------
    Integrated AnnData with shared latent space in obsm['X_multivi']
    """
    import scvi
    import muon as mu
    import scanpy as sc

    logger.info(
        f"Running MultiVI integration: "
        f"RNA ({adata_rna.n_obs} cells, {adata_rna.n_vars} genes) + "
        f"ATAC ({adata_atac.n_obs} cells, {adata_atac.n_vars} peaks)"
    )

    # Ensure barcodes match
    shared_barcodes = adata_rna.obs_names.intersection(adata_atac.obs_names)
    if len(shared_barcodes) == 0:
        raise ValueError("No shared barcodes between RNA and ATAC data.")

    logger.info(f"Shared barcodes: {len(shared_barcodes)}")

    adata_rna = adata_rna[shared_barcodes].copy()
    adata_atac = adata_atac[shared_barcodes].copy()

    # Concatenate for MultiVI (RNA features + ATAC peaks)
    # MultiVI expects a single AnnData with both modalities
    adata_multi = ad.concat(
        [adata_rna, adata_atac],
        axis=1,
        merge="unique",
    )

    # Use raw counts
    if "counts" in adata_rna.layers and "counts" in adata_atac.layers:
        import scipy.sparse as sp

        rna_counts = adata_rna.layers["counts"]
        atac_counts = adata_atac.layers["counts"]
        adata_multi.layers["counts"] = sp.hstack([rna_counts, atac_counts]).tocsr()

    # Mark modality for each feature
    adata_multi.var["modality"] = (
        ["Gene Expression"] * adata_rna.n_vars + ["Peaks"] * adata_atac.n_vars
    )

    # Setup MultiVI
    scvi.model.MULTIVI.setup_anndata(
        adata_multi,
        batch_key=batch_key,
        layer="counts" if "counts" in adata_multi.layers else None,
    )

    # Train
    model = scvi.model.MULTIVI(
        adata_multi,
        n_latent=n_latent,
        n_layers_encoder=n_layers,
        n_hidden=n_hidden,
    )

    model.train(
        max_epochs=max_epochs,
        batch_size=batch_size,
        early_stopping=True,
    )

    # Get latent representation
    latent = model.get_latent_representation()
    adata_multi.obsm["X_multivi"] = latent
    logger.info(f"MultiVI latent shape: {latent.shape}")

    # Also store back in original modalities
    adata_rna.obsm["X_multivi"] = latent
    adata_atac.obsm["X_multivi"] = latent

    # Compute neighborhood graph and UMAP
    if compute_umap:
        sc.pp.neighbors(adata_multi, use_rep="X_multivi", n_neighbors=n_neighbors)
        sc.tl.umap(adata_multi)
        sc.tl.leiden(adata_multi, resolution=1.0)
        logger.info(f"UMAP + Leiden complete: {adata_multi.obs['leiden'].nunique()} clusters")

    # Store model reference
    adata_multi.uns["multivi_model"] = {"n_latent": n_latent, "method": "MultiVI"}

    return adata_multi

# ImmunOmics v0.1.0
# Any usage is subject to this software's license.
