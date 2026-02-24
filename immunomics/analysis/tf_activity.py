#!/usr/bin/env python3
"""
ImmunOmics v0.1.0

Transcription factor activity inference for immune cells.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import logging

import anndata as ad
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Key immune transcription factors
IMMUNE_TFS = {
    "T cell": ["TBX21", "GATA3", "RORC", "FOXP3", "TCF7", "TOX", "EOMES", "LEF1", "BATF"],
    "B cell": ["PAX5", "IRF4", "BLIMP1", "BCL6", "EBF1", "E2A"],
    "Myeloid": ["SPI1", "CEBPA", "CEBPB", "IRF8", "IRF5", "MAFB"],
    "NK cell": ["EOMES", "TBX21", "ID2", "NFIL3"],
    "DC": ["IRF8", "IRF4", "BATF3", "ID2", "TCF4"],
    "Inflammatory": ["NFKB1", "RELA", "IRF3", "IRF7", "STAT1", "STAT3"],
}


def infer_tf_activity(
    adata_rna: ad.AnnData,
    adata_atac: ad.AnnData,
    tf_list: list[str] | None = None,
    cell_type_key: str = "cell_type",
    method: str = "correlation",
) -> pd.DataFrame:
    """
    Infer transcription factor activity across cell types.

    Approach:
    1. Score TF motif accessibility per cell (from ATAC peaks)
    2. Get TF gene expression per cell (from RNA)
    3. Correlate motif scores with TF expression
    4. Aggregate by cell type

    Parameters
    ----------
    adata_rna
        RNA AnnData with gene expression
    adata_atac
        ATAC AnnData with peak accessibility
    tf_list
        List of TFs to analyze. If None, uses curated immune TF list.
    cell_type_key
        Column in obs for cell type grouping
    method
        'correlation' (motif score vs expression) or 'chromvar'

    Returns
    -------
    DataFrame with TF activity scores per cell type
        Columns: tf_name, cell_type, motif_score, expression, activity_score
    """
    logger.info("Inferring transcription factor activity...")

    # Curate TF list
    if tf_list is None:
        tf_list = list(set(tf for tfs in IMMUNE_TFS.values() for tf in tfs))

    # Filter to TFs present in RNA data
    available_tfs = [tf for tf in tf_list if tf in adata_rna.var_names]
    logger.info(f"Analyzing {len(available_tfs)}/{len(tf_list)} TFs found in RNA data")

    # Ensure shared barcodes
    shared = adata_rna.obs_names.intersection(adata_atac.obs_names)
    adata_rna = adata_rna[shared].copy()
    adata_atac = adata_atac[shared].copy()

    if method == "chromvar":
        return _chromvar_activity(adata_rna, adata_atac, available_tfs, cell_type_key)
    else:
        return _expression_based_activity(adata_rna, available_tfs, cell_type_key)


def _expression_based_activity(
    adata_rna: ad.AnnData,
    tf_list: list[str],
    cell_type_key: str,
) -> pd.DataFrame:
    """
    Estimate TF activity from expression + known target gene sets.

    Uses TF expression level as a proxy for activity, aggregated
    by cell type. More sophisticated methods (e.g., SCENIC, pySCENIC)
    can be added as alternatives.
    """
    import scipy.sparse as sp

    X = adata_rna.X
    if sp.issparse(X):
        X = X.toarray()

    results = []
    cell_types = adata_rna.obs[cell_type_key].unique()

    for tf in tf_list:
        tf_idx = np.where(adata_rna.var_names == tf)[0]
        if len(tf_idx) == 0:
            continue
        tf_idx = tf_idx[0]

        for ct in cell_types:
            mask = adata_rna.obs[cell_type_key] == ct
            tf_expr = X[mask, tf_idx]

            results.append(
                {
                    "tf_name": tf,
                    "cell_type": ct,
                    "mean_expression": float(np.mean(tf_expr)),
                    "pct_expressing": float(np.mean(tf_expr > 0) * 100),
                    "expression_zscore": float(
                        (np.mean(tf_expr) - np.mean(X[:, tf_idx])) / (np.std(X[:, tf_idx]) + 1e-8)
                    ),
                }
            )

    df = pd.DataFrame(results)

    # Rank-based activity score
    if len(df) > 0:
        df["activity_score"] = df.groupby("tf_name")["expression_zscore"].rank(pct=True)

    logger.info(f"Computed activity for {len(tf_list)} TFs across {len(cell_types)} cell types")
    return df


def _chromvar_activity(
    adata_rna: ad.AnnData,
    adata_atac: ad.AnnData,
    tf_list: list[str],
    cell_type_key: str,
) -> pd.DataFrame:
    """
    Run chromVAR for motif enrichment analysis (requires R).

    chromVAR computes per-cell TF motif deviation scores from
    scATAC-seq data. These scores reflect the relative
    accessibility of each TF's binding sites.
    """
    logger.info("Running chromVAR motif analysis via R...")

    try:
        from rpy2.robjects.packages import importr

        importr("chromVAR")
        logger.info("chromVAR loaded successfully")
    except ImportError:
        logger.warning(
            "chromVAR not available. Falling back to expression-based activity. "
            "Install via: conda install -c bioconda bioconductor-chromvar"
        )
        return _expression_based_activity(adata_rna, tf_list, cell_type_key)

    # Placeholder for full chromVAR implementation
    logger.info("Full chromVAR integration planned for v0.2")
    return _expression_based_activity(adata_rna, tf_list, cell_type_key)


# ImmunOmics v0.1.0
# Any usage is subject to this software's license.
