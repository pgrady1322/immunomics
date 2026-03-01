#!/usr/bin/env python3
"""
ImmunOmics v0.1.0

Extended visualization tests — joint UMAP, save paths, edge cases.

Author: Patrick Grady
License: MIT License - See LICENSE
"""

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")

import matplotlib.pyplot as plt


@pytest.fixture
def adata_with_umap():
    """Create AnnData with X_umap coordinates for UMAP plotting."""
    np.random.seed(42)
    n_cells = 80
    X = np.random.rand(n_cells, 20).astype(np.float32)
    adata = ad.AnnData(X=X)
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.obs["cell_type"] = np.random.choice(["T cell", "B cell", "Monocyte"], n_cells)
    adata.obs["leiden"] = np.random.choice(["0", "1", "2"], n_cells)
    adata.obsm["X_umap"] = np.random.randn(n_cells, 2).astype(np.float32)
    return adata


class TestPlotJointUMAP:
    """Tests for plot_joint_umap function."""

    def test_basic_umap(self, adata_with_umap):
        from immunomics.visualization.plots import plot_joint_umap

        fig = plot_joint_umap(adata_with_umap)
        assert fig is not None
        plt.close(fig)

    def test_custom_color_keys(self, adata_with_umap):
        from immunomics.visualization.plots import plot_joint_umap

        fig = plot_joint_umap(adata_with_umap, color_keys=["cell_type"])
        assert fig is not None
        plt.close(fig)

    def test_save_path(self, adata_with_umap, tmp_path):
        import os

        from immunomics.visualization.plots import plot_joint_umap

        path = str(tmp_path / "umap.png")
        fig = plot_joint_umap(adata_with_umap, save_path=path)
        assert os.path.exists(path)
        plt.close(fig)

    def test_figsize(self, adata_with_umap):
        from immunomics.visualization.plots import plot_joint_umap

        fig = plot_joint_umap(adata_with_umap, figsize=(10, 4))
        assert fig is not None
        plt.close(fig)


class TestPlotTFHeatmapExtended:
    """Extended TF heatmap tests."""

    def test_few_tfs(self, mock_tf_activity_df):
        from immunomics.visualization.plots import plot_tf_heatmap

        fig = plot_tf_heatmap(mock_tf_activity_df, top_n=2)
        assert fig is not None
        plt.close(fig)

    def test_expression_zscore_fallback(self):
        """When activity_score is missing, should use expression_zscore."""
        from immunomics.visualization.plots import plot_tf_heatmap

        rows = []
        for tf in ["TBX21", "GATA3"]:
            for ct in ["T cell", "B cell"]:
                rows.append(
                    {
                        "tf_name": tf,
                        "cell_type": ct,
                        "expression_zscore": np.random.randn(),
                    }
                )
        df = pd.DataFrame(rows)
        fig = plot_tf_heatmap(df, top_n=5)
        assert fig is not None
        plt.close(fig)

    def test_save_path(self, mock_tf_activity_df, tmp_path):
        import os

        from immunomics.visualization.plots import plot_tf_heatmap

        path = str(tmp_path / "heatmap.png")
        fig = plot_tf_heatmap(mock_tf_activity_df, save_path=path)
        assert os.path.exists(path)
        plt.close(fig)


class TestPlotPeakGeneLinksExtended:
    """Extended peak-gene link plot tests."""

    def test_save_path(self, mock_peak_gene_links_df, tmp_path):
        import os

        from immunomics.visualization.plots import plot_peak_gene_links

        path = str(tmp_path / "links.png")
        fig = plot_peak_gene_links(mock_peak_gene_links_df, save_path=path)
        assert os.path.exists(path)
        plt.close(fig)

    def test_top_n_filter(self, mock_peak_gene_links_df):
        from immunomics.visualization.plots import plot_peak_gene_links

        fig = plot_peak_gene_links(mock_peak_gene_links_df, top_n=3)
        assert fig is not None
        plt.close(fig)


class TestPlotIntegrationComparisonExtended:
    """Extended integration comparison plot tests."""

    def test_single_method(self):
        from immunomics.visualization.plots import plot_integration_comparison

        results = {"multivi": {"silhouette_score": 0.4, "adjusted_rand_index": 0.5}}
        fig = plot_integration_comparison(results)
        assert fig is not None
        plt.close(fig)

    def test_all_errors(self):
        """All methods failed — should still return a figure."""
        from immunomics.visualization.plots import plot_integration_comparison

        results = {
            "multivi": {"error": "scvi not installed"},
            "mofa": {"error": "muon not installed"},
        }
        fig = plot_integration_comparison(results)
        assert fig is not None
        plt.close(fig)

    def test_custom_metrics(self):
        from immunomics.visualization.plots import plot_integration_comparison

        results = {"multivi": {"silhouette_score": 0.4, "runtime_seconds": 120.5}}
        fig = plot_integration_comparison(
            results, metrics=["silhouette_score", "runtime_seconds"]
        )
        assert fig is not None
        plt.close(fig)

    def test_save_path(self, tmp_path):
        import os

        from immunomics.visualization.plots import plot_integration_comparison

        results = {"multivi": {"silhouette_score": 0.4, "adjusted_rand_index": 0.5}}
        path = str(tmp_path / "comparison.png")
        fig = plot_integration_comparison(results, save_path=path)
        assert os.path.exists(path)
        plt.close(fig)
