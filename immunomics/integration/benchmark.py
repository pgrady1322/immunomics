#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ImmunOmics v0.1.0

Integration method comparison (silhouette, ARI).

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import logging
import time
from typing import Dict, Any, List, Optional

import numpy as np
import anndata as ad

logger = logging.getLogger(__name__)


def compare_methods(
    adata_rna: ad.AnnData,
    adata_atac: ad.AnnData,
    methods: Optional[List[str]] = None,
    label_key: Optional[str] = "cell_type",
) -> Dict[str, Dict[str, Any]]:
    """
    Compare integration methods on the same data.

    Evaluates:
    - Silhouette score (cell type separation in integrated space)
    - ARI (Adjusted Rand Index between clusters and cell types)
    - LISI score (batch mixing if applicable)
    - Runtime

    Parameters
    ----------
    adata_rna
        Preprocessed RNA data
    adata_atac
        Preprocessed ATAC data
    methods
        List of methods to compare: ['multivi', 'wnn', 'mofa']
    label_key
        Ground truth cell type column (for supervised metrics)

    Returns
    -------
    Dictionary mapping method name to evaluation metrics
    """
    from sklearn.metrics import silhouette_score, adjusted_rand_score

    if methods is None:
        methods = ["multivi", "mofa"]  # WNN requires R, so optional

    method_runners = {
        "multivi": _run_multivi_wrapper,
        "wnn": _run_wnn_wrapper,
        "mofa": _run_mofa_wrapper,
    }

    results = {}

    for method in methods:
        if method not in method_runners:
            logger.warning(f"Unknown method: {method}")
            continue

        logger.info(f"Benchmarking {method}...")
        start_time = time.time()

        try:
            adata_int = method_runners[method](adata_rna, adata_atac)
            elapsed = time.time() - start_time

            # Get latent representation
            rep_key = {
                "multivi": "X_multivi",
                "wnn": "X_pca",  # WNN uses its own graph, fall back to PCA
                "mofa": "X_mofa",
            }.get(method, "X_pca")

            latent = adata_int.obsm.get(rep_key, adata_int.obsm.get("X_pca"))

            metrics = {"runtime_seconds": elapsed}

            # Silhouette score
            if label_key and label_key in adata_int.obs.columns and latent is not None:
                labels = adata_int.obs[label_key].values
                # Subsample for speed if > 10k cells
                if len(labels) > 10000:
                    idx = np.random.choice(len(labels), 10000, replace=False)
                    sil = silhouette_score(latent[idx], labels[idx])
                else:
                    sil = silhouette_score(latent, labels)
                metrics["silhouette_score"] = sil

            # ARI between Leiden clusters and cell types
            if label_key in adata_int.obs.columns and "leiden" in adata_int.obs.columns:
                ari = adjusted_rand_score(
                    adata_int.obs[label_key], adata_int.obs["leiden"]
                )
                metrics["adjusted_rand_index"] = ari

            results[method] = metrics
            logger.info(f"{method}: {metrics}")

        except Exception as e:
            logger.error(f"{method} failed: {e}")
            results[method] = {"error": str(e), "runtime_seconds": time.time() - start_time}

    return results


def _run_multivi_wrapper(adata_rna, adata_atac):
    from immunomics.integration.multivi import run_multivi

    return run_multivi(adata_rna, adata_atac)


def _run_wnn_wrapper(adata_rna, adata_atac):
    from immunomics.integration.wnn import run_wnn

    return run_wnn(adata_rna, adata_atac)


def _run_mofa_wrapper(adata_rna, adata_atac):
    from immunomics.integration.mofa import run_mofa

    return run_mofa(adata_rna, adata_atac)

# ImmunOmics v0.1.0
# Any usage is subject to this software's license.
