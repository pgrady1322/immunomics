#!/usr/bin/env python3
"""
ImmunOmics v0.1.0

Extended integration benchmark tests — mock-based method comparison.

Author: Patrick Grady
License: MIT License - See LICENSE
"""

import anndata as ad
import numpy as np
import pytest
from scipy.sparse import csr_matrix
from unittest.mock import patch, MagicMock


@pytest.fixture
def preprocessed_pair():
    """Create a preprocessed RNA + ATAC pair with shared barcodes."""
    np.random.seed(42)
    n_cells = 80

    # RNA
    rna_X = csr_matrix(np.random.rand(n_cells, 50).astype(np.float32))
    adata_rna = ad.AnnData(X=rna_X)
    adata_rna.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata_rna.var_names = [f"Gene_{i}" for i in range(50)]
    adata_rna.obs["cell_type"] = np.random.choice(["T cell", "B cell", "Monocyte"], n_cells)
    adata_rna.obsm["X_pca"] = np.random.randn(n_cells, 10).astype(np.float32)

    # ATAC
    atac_X = csr_matrix((np.random.rand(n_cells, 40) > 0.8).astype(np.float32))
    adata_atac = ad.AnnData(X=atac_X)
    adata_atac.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata_atac.var_names = [f"chr1:{i * 1000}-{i * 1000 + 500}" for i in range(40)]
    adata_atac.obs["cell_type"] = adata_rna.obs["cell_type"].values
    adata_atac.obsm["X_lsi"] = np.random.randn(n_cells, 10).astype(np.float32)

    return adata_rna, adata_atac


class TestBenchmarkCompareMethods:
    """Tests for compare_methods benchmark function."""

    def test_import(self):
        from immunomics.integration.benchmark import compare_methods

        assert callable(compare_methods)

    def test_unknown_method_skipped(self, preprocessed_pair):
        from immunomics.integration.benchmark import compare_methods

        adata_rna, adata_atac = preprocessed_pair
        results = compare_methods(adata_rna, adata_atac, methods=["nonexistent"])
        assert "nonexistent" not in results

    def test_empty_methods_list(self, preprocessed_pair):
        from immunomics.integration.benchmark import compare_methods

        adata_rna, adata_atac = preprocessed_pair
        results = compare_methods(adata_rna, adata_atac, methods=[])
        assert results == {}

    def test_method_failure_captured(self, preprocessed_pair):
        """When a method fails, the error should be captured, not raised."""
        from immunomics.integration.benchmark import compare_methods

        adata_rna, adata_atac = preprocessed_pair

        # Mock multivi to raise
        with patch(
            "immunomics.integration.benchmark._run_multivi_wrapper",
            side_effect=RuntimeError("scvi not installed"),
        ):
            results = compare_methods(adata_rna, adata_atac, methods=["multivi"])
            assert "multivi" in results
            assert "error" in results["multivi"]
            assert "scvi not installed" in results["multivi"]["error"]

    def test_successful_method_has_metrics(self, preprocessed_pair):
        """A successful integration run should produce metrics."""
        from immunomics.integration.benchmark import compare_methods

        adata_rna, adata_atac = preprocessed_pair

        # Mock multivi to return a fake integrated result
        mock_integrated = adata_rna.copy()
        mock_integrated.obsm["X_multivi"] = np.random.randn(
            mock_integrated.n_obs, 10
        ).astype(np.float32)

        with patch(
            "immunomics.integration.benchmark._run_multivi_wrapper",
            return_value=mock_integrated,
        ):
            results = compare_methods(adata_rna, adata_atac, methods=["multivi"])
            assert "multivi" in results
            assert "runtime_seconds" in results["multivi"]
            assert "silhouette_score" in results["multivi"]
            assert "error" not in results["multivi"]

    def test_default_methods(self, preprocessed_pair):
        """Default methods should be multivi + mofa."""
        from immunomics.integration.benchmark import compare_methods

        adata_rna, adata_atac = preprocessed_pair

        # Mock both to fail quickly (we just want to verify which methods run)
        with patch(
            "immunomics.integration.benchmark._run_multivi_wrapper",
            side_effect=ImportError("no scvi"),
        ), patch(
            "immunomics.integration.benchmark._run_mofa_wrapper",
            side_effect=ImportError("no muon"),
        ):
            results = compare_methods(adata_rna, adata_atac, methods=None)
            # Default is ["multivi", "mofa"]
            assert "multivi" in results
            assert "mofa" in results

    def test_no_label_key(self, preprocessed_pair):
        """If label_key is None, silhouette score should not be computed."""
        from immunomics.integration.benchmark import compare_methods

        adata_rna, adata_atac = preprocessed_pair

        mock_integrated = adata_rna.copy()
        mock_integrated.obsm["X_multivi"] = np.random.randn(
            mock_integrated.n_obs, 10
        ).astype(np.float32)

        with patch(
            "immunomics.integration.benchmark._run_multivi_wrapper",
            return_value=mock_integrated,
        ):
            results = compare_methods(
                adata_rna, adata_atac, methods=["multivi"], label_key=None
            )
            assert "multivi" in results
            assert "silhouette_score" not in results["multivi"]


class TestIntegrationModuleImportsExtended:
    """Extended import and callable tests."""

    def test_run_multivi_signature(self):
        from immunomics.integration.multivi import run_multivi
        import inspect

        sig = inspect.signature(run_multivi)
        assert "adata_rna" in sig.parameters
        assert "adata_atac" in sig.parameters

    def test_run_mofa_signature(self):
        from immunomics.integration.mofa import run_mofa
        import inspect

        sig = inspect.signature(run_mofa)
        assert "adata_rna" in sig.parameters
        assert "adata_atac" in sig.parameters

    def test_run_wnn_signature(self):
        from immunomics.integration.wnn import run_wnn
        import inspect

        sig = inspect.signature(run_wnn)
        assert "adata_rna" in sig.parameters
        assert "adata_atac" in sig.parameters
