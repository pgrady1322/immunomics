#!/usr/bin/env python3
"""
ImmunOmics v0.1.0

Extended preprocessing tests — RNA and ATAC edge cases, internal helpers.

Author: Patrick Grady
License: MIT License - See LICENSE
"""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix


class TestPreprocessRNAExtended:
    """Extended tests for preprocess_rna edge cases."""

    def test_inplace_modifies_input(self, mock_rna):
        from immunomics.data.preprocess_rna import preprocess_rna

        result = preprocess_rna(mock_rna, n_top_genes=20, n_pcs=10, copy=False)
        # When copy=False the same object is mutated
        assert "X_pca" in mock_rna.obsm

    def test_custom_target_sum(self, mock_rna):
        from immunomics.data.preprocess_rna import preprocess_rna

        result = preprocess_rna(mock_rna, n_top_genes=20, n_pcs=10, target_sum=1e3)
        assert "X_pca" in result.obsm

    def test_custom_max_pct_mito(self, mock_rna):
        from immunomics.data.preprocess_rna import preprocess_rna

        # Very strict filter
        result = preprocess_rna(mock_rna, n_top_genes=20, n_pcs=10, max_pct_mito=0.01)
        # Should still return valid AnnData even if filtered heavily
        assert isinstance(result, ad.AnnData)

    def test_no_mt_genes(self):
        """Data with no MT genes should still preprocess cleanly."""
        from immunomics.data.preprocess_rna import preprocess_rna

        np.random.seed(42)
        n_cells, n_genes = 100, 200
        X = csr_matrix(np.random.rand(n_cells, n_genes).astype(np.float32))
        adata = ad.AnnData(X=X)
        adata.var_names = [f"Gene_{i}" for i in range(n_genes)]
        adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
        adata.obs["cell_type"] = "T cell"
        adata.layers["counts"] = X.copy()

        result = preprocess_rna(adata, n_top_genes=20, n_pcs=10, min_genes=1, min_cells=1)
        assert "X_pca" in result.obsm

    def test_dense_matrix_input(self):
        """Dense matrix input should work."""
        from immunomics.data.preprocess_rna import preprocess_rna

        np.random.seed(42)
        n_cells, n_genes = 100, 200
        X = np.random.rand(n_cells, n_genes).astype(np.float32)
        adata = ad.AnnData(X=X)
        adata.var_names = [f"Gene_{i}" for i in range(n_genes)]
        adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
        adata.obs["cell_type"] = "T cell"
        adata.layers["counts"] = X.copy()

        result = preprocess_rna(adata, n_top_genes=20, n_pcs=10, min_genes=1, min_cells=1)
        assert "X_pca" in result.obsm

    def test_all_qc_metrics_present(self, mock_rna):
        from immunomics.data.preprocess_rna import preprocess_rna

        result = preprocess_rna(mock_rna, n_top_genes=20, n_pcs=10)
        assert "mt" in result.var.columns
        assert "pct_counts_mt" in result.obs.columns

    def test_raw_stored(self, mock_rna):
        from immunomics.data.preprocess_rna import preprocess_rna

        result = preprocess_rna(mock_rna, n_top_genes=20, n_pcs=10)
        assert result.raw is not None


class TestPreprocessATACExtended:
    """Extended tests for preprocess_atac edge cases and internal helpers."""

    def test_inplace_modifies_input(self, mock_atac):
        from immunomics.data.preprocess_atac import preprocess_atac

        preprocess_atac(mock_atac, min_peaks=1, min_cells=1, n_components=10, copy=False)
        assert "X_lsi" in mock_atac.obsm

    def test_dense_matrix_input(self):
        """Dense matrix should be handled by TF-IDF."""
        from immunomics.data.preprocess_atac import preprocess_atac

        np.random.seed(42)
        n_cells, n_peaks = 50, 100
        X = (np.random.rand(n_cells, n_peaks) > 0.9).astype(np.float32)
        adata = ad.AnnData(X=X)
        adata.var_names = [f"chr1:{i * 1000}-{i * 1000 + 500}" for i in range(n_peaks)]
        adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
        adata.layers["counts"] = X.copy()

        result = preprocess_atac(adata, min_peaks=1, min_cells=1, n_components=5)
        assert "X_lsi" in result.obsm
        assert "tfidf" in result.layers

    def test_tfidf_normalize_directly(self, mock_atac):
        from immunomics.data.preprocess_atac import _tfidf_normalize

        result = _tfidf_normalize(mock_atac.copy())
        assert "tfidf" in result.layers
        # TF-IDF should not have negative values
        import scipy.sparse as sp

        X = result.layers["tfidf"]
        if sp.issparse(X):
            assert X.data.min() >= 0
        else:
            assert X.min() >= 0

    def test_run_lsi_directly(self, mock_atac):
        from immunomics.data.preprocess_atac import _run_lsi, _tfidf_normalize

        adata = _tfidf_normalize(mock_atac.copy())
        result = _run_lsi(adata, n_components=5)
        assert "X_lsi" in result.obsm
        assert result.obsm["X_lsi"].shape[1] == 5
        assert "lsi_variance_ratio" in result.uns
        assert len(result.uns["lsi_variance_ratio"]) == 5

    def test_parse_peak_coordinates_malformed(self):
        """Malformed peak names get 'unknown' chromosome."""
        from immunomics.data.preprocess_atac import parse_peak_coordinates

        adata = ad.AnnData(
            X=np.zeros((3, 3)),
            var=pd.DataFrame(index=["bad_peak", "chr1:1000-2000", "also_bad"]),
        )
        result = parse_peak_coordinates(adata)
        assert result.var.loc["bad_peak", "chrom"] == "unknown"
        assert result.var.loc["bad_peak", "start"] == 0
        assert result.var.loc["chr1:1000-2000", "chrom"] == "chr1"
        assert result.var.loc["chr1:1000-2000", "start"] == 1000

    def test_n_components_affects_output_dimension(self, mock_atac):
        from immunomics.data.preprocess_atac import preprocess_atac

        for n_comp in [5, 10]:
            result = preprocess_atac(mock_atac, min_peaks=1, min_cells=1, n_components=n_comp)
            assert result.obsm["X_lsi"].shape[1] == n_comp

    def test_counts_layer_auto_stored(self):
        """If no counts layer exists, preprocess_atac should store one."""
        from immunomics.data.preprocess_atac import preprocess_atac

        np.random.seed(42)
        n_cells, n_peaks = 50, 100
        X = csr_matrix((np.random.rand(n_cells, n_peaks) > 0.9).astype(np.float32))
        adata = ad.AnnData(X=X)
        adata.var_names = [f"chr1:{i * 1000}-{i * 1000 + 500}" for i in range(n_peaks)]
        adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
        # Deliberately no counts layer
        assert "counts" not in adata.layers

        result = preprocess_atac(adata, min_peaks=1, min_cells=1, n_components=5)
        assert "counts" in result.layers
