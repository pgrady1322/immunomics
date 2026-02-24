"""
Tests for ImmunOmics.
"""

import pytest
import numpy as np
import anndata as ad
from scipy.sparse import csr_matrix


@pytest.fixture
def mock_rna():
    """Create mock RNA AnnData."""
    np.random.seed(42)
    n_cells, n_genes = 100, 50
    X = csr_matrix(np.random.rand(n_cells, n_genes).astype(np.float32))

    adata = ad.AnnData(X=X)
    adata.var_names = [f"Gene_{i}" for i in range(n_genes)]
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.obs["cell_type"] = np.random.choice(["T cell", "B cell", "Monocyte"], n_cells)
    adata.layers["counts"] = X.copy()

    return adata


@pytest.fixture
def mock_atac():
    """Create mock ATAC AnnData."""
    np.random.seed(42)
    n_cells, n_peaks = 100, 200
    X = csr_matrix((np.random.rand(n_cells, n_peaks) > 0.9).astype(np.float32))

    adata = ad.AnnData(X=X)
    adata.var_names = [f"chr1:{i*1000}-{i*1000+500}" for i in range(n_peaks)]
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.obs["cell_type"] = np.random.choice(["T cell", "B cell", "Monocyte"], n_cells)
    adata.layers["counts"] = X.copy()

    return adata


class TestPreprocessRNA:
    def test_basic_preprocessing(self, mock_rna):
        from immunomics.data.preprocess_rna import preprocess_rna

        result = preprocess_rna(mock_rna, n_top_genes=20, n_pcs=10)
        assert "X_pca" in result.obsm
        assert result.obsm["X_pca"].shape[1] == 10


class TestPreprocessATAC:
    def test_tfidf_and_lsi(self, mock_atac):
        from immunomics.data.preprocess_atac import preprocess_atac

        result = preprocess_atac(mock_atac, min_peaks=1, min_cells=1, n_components=10)
        assert "X_lsi" in result.obsm
        assert result.obsm["X_lsi"].shape[1] == 10

    def test_parse_peak_coordinates(self, mock_atac):
        from immunomics.data.preprocess_atac import parse_peak_coordinates

        result = parse_peak_coordinates(mock_atac)
        assert "chrom" in result.var.columns
        assert "start" in result.var.columns
        assert "end" in result.var.columns


class TestDifferential:
    def test_differential_genes(self, mock_rna):
        import scanpy as sc
        from immunomics.analysis.differential import differential_genes

        sc.pp.normalize_total(mock_rna)
        sc.pp.log1p(mock_rna)

        result = differential_genes(mock_rna, groupby="cell_type", n_genes=5)
        assert len(result) > 0
        assert "names" in result.columns or "group" in result.columns
