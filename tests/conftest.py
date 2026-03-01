#!/usr/bin/env python3
"""
ImmunOmics v0.1.0

Shared fixtures for the test suite.

Author: Patrick Grady
License: MIT License - See LICENSE
"""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix


@pytest.fixture
def mock_rna():
    """Create mock RNA AnnData with realistic structure."""
    np.random.seed(42)
    n_cells, n_genes = 200, 500
    X = csr_matrix(np.random.rand(n_cells, n_genes).astype(np.float32))

    adata = ad.AnnData(X=X)
    adata.var_names = pd.Index(
        [f"MT-Gene_{i}" if i < 5 else f"Gene_{i}" for i in range(n_genes)]
    )
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.obs["cell_type"] = np.random.choice(["T cell", "B cell", "Monocyte"], n_cells)
    adata.layers["counts"] = X.copy()

    return adata


@pytest.fixture
def mock_atac():
    """Create mock ATAC AnnData with peak-format var names."""
    np.random.seed(42)
    n_cells, n_peaks = 100, 200
    X = csr_matrix((np.random.rand(n_cells, n_peaks) > 0.9).astype(np.float32))

    adata = ad.AnnData(X=X)
    adata.var_names = [f"chr1:{i * 1000}-{i * 1000 + 500}" for i in range(n_peaks)]
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.obs["cell_type"] = np.random.choice(["T cell", "B cell", "Monocyte"], n_cells)
    adata.layers["counts"] = X.copy()

    return adata


@pytest.fixture
def mock_rna_with_tfs():
    """Create RNA AnnData containing known immune TF genes."""
    np.random.seed(42)
    n_cells = 80
    tfs = ["TBX21", "GATA3", "PAX5", "SPI1", "FOXP3", "IRF4"]
    other = [f"Gene_{i}" for i in range(20)]
    genes = tfs + other
    X = csr_matrix(np.random.rand(n_cells, len(genes)).astype(np.float32))

    adata = ad.AnnData(X=X)
    adata.var_names = genes
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.obs["cell_type"] = np.random.choice(["T cell", "B cell", "Monocyte"], n_cells)
    adata.layers["counts"] = X.copy()

    return adata


@pytest.fixture
def mock_tf_activity_df():
    """Create a mock TF activity DataFrame for plot testing."""
    rows = []
    for tf in ["TBX21", "GATA3", "PAX5", "SPI1"]:
        for ct in ["T cell", "B cell", "Monocyte"]:
            rows.append(
                {
                    "tf_name": tf,
                    "cell_type": ct,
                    "expression_zscore": np.random.randn(),
                    "activity_score": np.random.rand(),
                }
            )
    return pd.DataFrame(rows)


@pytest.fixture
def mock_peak_gene_links_df():
    """Create a mock peak-gene links DataFrame for plot testing."""
    return pd.DataFrame(
        {
            "peak": [f"chr1:{i * 1000}-{i * 1000 + 500}" for i in range(10)],
            "gene": [f"Gene_{i}" for i in range(10)],
            "correlation": np.random.uniform(0.1, 0.9, 10),
            "pvalue": np.random.uniform(0.001, 0.05, 10),
            "distance": np.random.randint(1000, 500000, 10),
            "chrom": ["chr1"] * 10,
        }
    )
