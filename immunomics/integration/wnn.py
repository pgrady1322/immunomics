#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ImmunOmics v0.1.0

Weighted Nearest Neighbor integration via rpy2 + Seurat v5.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import logging
from typing import Optional

import numpy as np
import anndata as ad

logger = logging.getLogger(__name__)


def run_wnn(
    adata_rna: ad.AnnData,
    adata_atac: ad.AnnData,
    n_neighbors: int = 20,
    compute_umap: bool = True,
) -> ad.AnnData:
    """
    Integrate RNA + ATAC using Seurat v5 Weighted Nearest Neighbors.

    Requires R with Seurat v5 and Signac installed.
    Uses rpy2 for R interoperability.

    Parameters
    ----------
    adata_rna
        Preprocessed RNA AnnData (needs PCA in obsm['X_pca'])
    adata_atac
        Preprocessed ATAC AnnData (needs LSI in obsm['X_lsi'])
    n_neighbors
        Number of neighbors for KNN graph
    compute_umap
        Whether to compute WNN UMAP

    Returns
    -------
    Integrated AnnData with WNN graph and UMAP
    """
    logger.info("Running Seurat v5 WNN integration via rpy2...")

    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri, numpy2ri
        from rpy2.robjects.packages import importr

        pandas2ri.activate()
        numpy2ri.activate()
    except ImportError:
        raise ImportError(
            "WNN integration requires rpy2 and R packages: Seurat (v5+), Signac. "
            "Install with: conda install -c conda-forge rpy2 r-seurat r-signac"
        )

    # Import R packages
    seurat = importr("Seurat")
    base = importr("base")

    # Ensure shared barcodes
    shared = adata_rna.obs_names.intersection(adata_atac.obs_names)
    logger.info(f"Shared barcodes: {len(shared)}")

    adata_rna = adata_rna[shared].copy()
    adata_atac = adata_atac[shared].copy()

    # Convert to R matrices
    rna_pca = ro.r.matrix(
        ro.FloatVector(adata_rna.obsm["X_pca"].flatten()),
        nrow=len(shared),
        ncol=adata_rna.obsm["X_pca"].shape[1],
    )

    atac_lsi = ro.r.matrix(
        ro.FloatVector(adata_atac.obsm["X_lsi"].flatten()),
        nrow=len(shared),
        ncol=adata_atac.obsm["X_lsi"].shape[1],
    )

    # Run WNN in R
    ro.r("""
    run_wnn <- function(rna_pca, atac_lsi, n_neighbors) {
        library(Seurat)

        # Create dummy Seurat object
        dummy_counts <- matrix(0, nrow=10, ncol=nrow(rna_pca))
        colnames(dummy_counts) <- paste0("cell_", 1:nrow(rna_pca))
        rownames(dummy_counts) <- paste0("gene_", 1:10)

        obj <- CreateSeuratObject(counts=dummy_counts)

        # Add reductions
        obj[["pca"]] <- CreateDimReducObject(
            embeddings=rna_pca,
            key="PC_",
            assay="RNA"
        )
        obj[["lsi"]] <- CreateDimReducObject(
            embeddings=atac_lsi,
            key="LSI_",
            assay="RNA"
        )

        # Find multimodal neighbors
        obj <- FindMultiModalNeighbors(
            obj,
            reduction.list=list("pca", "lsi"),
            dims.list=list(1:ncol(rna_pca), 1:ncol(atac_lsi)),
            k.nn=n_neighbors
        )

        # Run WNN UMAP
        obj <- RunUMAP(
            obj,
            nn.name="weighted.nn",
            reduction.name="wnn.umap",
            reduction.key="wnnUMAP_"
        )

        # Cluster on WNN graph
        obj <- FindClusters(obj, graph.name="wsnn", resolution=1.0, algorithm=3)

        # Extract results
        list(
            umap=Embeddings(obj, "wnn.umap"),
            clusters=obj$seurat_clusters,
            modality_weights=obj@neighbors$weighted.nn@nn.dist
        )
    }
    """)

    wnn_results = ro.r["run_wnn"](rna_pca, atac_lsi, n_neighbors)

    # Build result AnnData
    adata_integrated = adata_rna.copy()
    adata_integrated.obsm["X_umap_wnn"] = np.array(wnn_results[0])
    adata_integrated.obs["wnn_clusters"] = np.array(wnn_results[1]).astype(str)
    adata_integrated.obsm["X_lsi"] = adata_atac.obsm["X_lsi"]
    adata_integrated.uns["integration_method"] = "WNN"

    logger.info(
        f"WNN integration complete: {len(np.unique(adata_integrated.obs['wnn_clusters']))} clusters"
    )

    return adata_integrated

# ImmunOmics v0.1.0
# Any usage is subject to this software's license.
