#!/usr/bin/env python3
"""
ImmunOmics v0.1.0

Dataset loader for 10x Multiome PBMC.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import logging
from pathlib import Path

import anndata as ad

logger = logging.getLogger(__name__)

CACHE_DIR = Path.home() / ".cache" / "immunomics" / "data"

DATASETS = {
    "multiome_pbmc_10k": {
        "description": "10x Genomics Multiome PBMC 10k (matched RNA + ATAC)",
        "url_rna": "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5",
        "url_atac_fragments": "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz",
        "n_cells": "~10,000",
        "reference": "10x Genomics (2022)",
    },
}


def _ensure_cache_dir() -> Path:
    """Create cache directory if it doesn't exist."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return CACHE_DIR


def load_multiome_pbmc(
    cache_dir: str | None = None,
    force_download: bool = False,
) -> tuple[ad.AnnData, ad.AnnData]:
    """
    Load 10x Genomics Multiome PBMC dataset.

    Returns matched scRNA-seq and scATAC-seq AnnData objects
    from the same cells (shared barcodes).

    Parameters
    ----------
    cache_dir
        Directory to cache downloaded data
    force_download
        If True, re-download even if cached

    Returns
    -------
    Tuple of (adata_rna, adata_atac)
    """
    import muon as mu

    cache = Path(cache_dir) if cache_dir else _ensure_cache_dir()
    mdata_path = cache / "multiome_pbmc.h5mu"

    if mdata_path.exists() and not force_download:
        logger.info(f"Loading cached Multiome PBMC from {mdata_path}")
        mdata = mu.read(str(mdata_path))
        return mdata.mod["rna"], mdata.mod["atac"]

    logger.info("Downloading 10x Multiome PBMC dataset...")
    logger.info("This includes both gene expression and ATAC-seq data.")

    # Download RNA (gene expression + peaks in same h5)
    import urllib.request

    h5_path = cache / "pbmc_multiome_filtered.h5"
    dataset = DATASETS["multiome_pbmc_10k"]

    if not h5_path.exists() or force_download:
        urllib.request.urlretrieve(dataset["url_rna"], h5_path)
        logger.info(f"Downloaded filtered matrix to {h5_path}")

    # Read as MuData (muon handles the multiome format)
    import scanpy as sc

    adata = sc.read_10x_h5(h5_path, gex_only=False)
    adata.var_names_make_unique()

    # Split into RNA and ATAC modalities
    # In 10x Multiome h5, features are labeled as "Gene Expression" or "Peaks"
    if "feature_types" in adata.var.columns:
        rna_mask = adata.var["feature_types"] == "Gene Expression"
        atac_mask = adata.var["feature_types"] == "Peaks"

        adata_rna = adata[:, rna_mask].copy()
        adata_atac = adata[:, atac_mask].copy()
    else:
        # Fallback: try genome annotation
        logger.warning("Could not find feature_types column. Attempting heuristic split.")
        # Peaks are formatted as chr:start-end
        is_peak = adata.var_names.str.match(r"^chr\d+[:-]")
        adata_rna = adata[:, ~is_peak].copy()
        adata_atac = adata[:, is_peak].copy()

    # Store raw counts
    adata_rna.layers["counts"] = adata_rna.X.copy()
    adata_atac.layers["counts"] = adata_atac.X.copy()

    logger.info(
        f"RNA: {adata_rna.n_obs} cells, {adata_rna.n_vars} genes | "
        f"ATAC: {adata_atac.n_obs} cells, {adata_atac.n_vars} peaks"
    )

    # Save as MuData
    mdata = mu.MuData({"rna": adata_rna, "atac": adata_atac})
    mdata.write(str(mdata_path))
    logger.info(f"Saved to {mdata_path}")

    return adata_rna, adata_atac


# ImmunOmics v0.1.0
# Any usage is subject to this software's license.
