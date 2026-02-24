#!/usr/bin/env python3
"""
ImmunOmics v0.1.0

Comprehensive test suite for all modules.

Author: Patrick Grady
License: MIT License - See LICENSE
"""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner
from scipy.sparse import csr_matrix

# ─── Fixtures ────────────────────────────────────────────────────────────────


@pytest.fixture
def mock_rna():
    """Create mock RNA AnnData with realistic structure."""
    np.random.seed(42)
    n_cells, n_genes = 200, 500
    X = csr_matrix(np.random.rand(n_cells, n_genes).astype(np.float32))

    adata = ad.AnnData(X=X)
    adata.var_names = [f"Gene_{i}" for i in range(n_genes)]
    # Mark a few as mitochondrial so QC filter works
    adata.var_names = pd.Index([f"MT-Gene_{i}" if i < 5 else f"Gene_{i}" for i in range(n_genes)])
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


# ─── Package Metadata ───────────────────────────────────────────────────────


class TestPackageMetadata:
    def test_version(self):
        from immunomics import __version__

        assert __version__ == "0.1.0"

    def test_subpackage_imports(self):
        """All subpackages should be importable."""
        import immunomics

        assert hasattr(immunomics, "data")
        assert hasattr(immunomics, "integration")
        assert hasattr(immunomics, "analysis")
        assert hasattr(immunomics, "visualization")


# ─── Data Preprocessing ─────────────────────────────────────────────────────


class TestPreprocessRNA:
    def test_basic_preprocessing(self, mock_rna):
        from immunomics.data.preprocess_rna import preprocess_rna

        result = preprocess_rna(mock_rna, n_top_genes=20, n_pcs=10)
        assert "X_pca" in result.obsm
        assert result.obsm["X_pca"].shape[1] == 10

    def test_preserves_counts_layer(self, mock_rna):
        from immunomics.data.preprocess_rna import preprocess_rna

        result = preprocess_rna(mock_rna, n_top_genes=20, n_pcs=10)
        assert "counts" in result.layers

    def test_copy_default(self, mock_rna):
        from immunomics.data.preprocess_rna import preprocess_rna

        original_shape = mock_rna.shape
        _ = preprocess_rna(mock_rna, n_top_genes=20, n_pcs=10)
        assert mock_rna.shape == original_shape

    def test_hvg_selection(self, mock_rna):
        from immunomics.data.preprocess_rna import preprocess_rna

        result = preprocess_rna(mock_rna, n_top_genes=20, n_pcs=10)
        assert "highly_variable" in result.var.columns


class TestPreprocessATAC:
    def test_tfidf_and_lsi(self, mock_atac):
        from immunomics.data.preprocess_atac import preprocess_atac

        result = preprocess_atac(mock_atac, min_peaks=1, min_cells=1, n_components=10)
        assert "X_lsi" in result.obsm
        assert result.obsm["X_lsi"].shape[1] == 10

    def test_tfidf_layer_stored(self, mock_atac):
        from immunomics.data.preprocess_atac import preprocess_atac

        result = preprocess_atac(mock_atac, min_peaks=1, min_cells=1, n_components=10)
        assert "tfidf" in result.layers

    def test_lsi_variance_stored(self, mock_atac):
        from immunomics.data.preprocess_atac import preprocess_atac

        result = preprocess_atac(mock_atac, min_peaks=1, min_cells=1, n_components=10)
        assert "lsi_variance_ratio" in result.uns

    def test_parse_peak_coordinates(self, mock_atac):
        from immunomics.data.preprocess_atac import parse_peak_coordinates

        result = parse_peak_coordinates(mock_atac)
        assert "chrom" in result.var.columns
        assert "start" in result.var.columns
        assert "end" in result.var.columns
        assert all(result.var["chrom"] == "chr1")

    def test_copy_default(self, mock_atac):
        from immunomics.data.preprocess_atac import preprocess_atac

        original_shape = mock_atac.shape
        _ = preprocess_atac(mock_atac, min_peaks=1, min_cells=1, n_components=10)
        assert mock_atac.shape == original_shape


class TestDatasets:
    def test_datasets_registry(self):
        from immunomics.data.datasets import DATASETS

        assert "multiome_pbmc_10k" in DATASETS
        assert "description" in DATASETS["multiome_pbmc_10k"]
        assert "url_rna" in DATASETS["multiome_pbmc_10k"]

    def test_cache_dir_creation(self):
        from immunomics.data.datasets import CACHE_DIR

        assert CACHE_DIR.name == "data"
        assert "immunomics" in str(CACHE_DIR)


# ─── Differential Analysis ──────────────────────────────────────────────────


class TestDifferential:
    def test_differential_genes(self, mock_rna):
        import scanpy as sc

        from immunomics.analysis.differential import differential_genes

        sc.pp.normalize_total(mock_rna)
        sc.pp.log1p(mock_rna)

        result = differential_genes(mock_rna, groupby="cell_type", n_genes=5)
        assert len(result) > 0
        assert "names" in result.columns or "group" in result.columns

    def test_differential_peaks(self, mock_atac):
        import scanpy as sc

        from immunomics.analysis.differential import differential_peaks

        sc.pp.normalize_total(mock_atac)
        sc.pp.log1p(mock_atac)

        result = differential_peaks(mock_atac, groupby="cell_type", n_peaks=5)
        assert len(result) > 0

    def test_differential_custom_method(self, mock_rna):
        import scanpy as sc

        from immunomics.analysis.differential import differential_genes

        sc.pp.normalize_total(mock_rna)
        sc.pp.log1p(mock_rna)

        result = differential_genes(mock_rna, groupby="cell_type", method="t-test", n_genes=5)
        assert len(result) > 0


# ─── TF Activity ────────────────────────────────────────────────────────────


class TestTFActivity:
    def test_immune_tfs_defined(self):
        from immunomics.analysis.tf_activity import IMMUNE_TFS

        assert "T cell" in IMMUNE_TFS
        assert "B cell" in IMMUNE_TFS
        assert "Myeloid" in IMMUNE_TFS
        assert len(IMMUNE_TFS) >= 5

    def test_expression_based_activity(self, mock_rna_with_tfs, mock_atac):
        from immunomics.analysis.tf_activity import infer_tf_activity

        shared = mock_rna_with_tfs.obs_names[:80]
        mock_atac_sub = mock_atac[:80].copy()
        mock_atac_sub.obs_names = shared

        result = infer_tf_activity(
            mock_rna_with_tfs,
            mock_atac_sub,
            tf_list=["TBX21", "GATA3", "PAX5"],
        )
        assert len(result) > 0
        assert "tf_name" in result.columns
        assert "cell_type" in result.columns
        assert "expression_zscore" in result.columns

    def test_activity_score_computed(self, mock_rna_with_tfs, mock_atac):
        from immunomics.analysis.tf_activity import infer_tf_activity

        shared = mock_rna_with_tfs.obs_names[:80]
        mock_atac_sub = mock_atac[:80].copy()
        mock_atac_sub.obs_names = shared

        result = infer_tf_activity(
            mock_rna_with_tfs,
            mock_atac_sub,
            tf_list=["TBX21", "GATA3"],
        )
        assert "activity_score" in result.columns
        assert result["activity_score"].between(0, 1).all()


# ─── Peak-Gene Links ────────────────────────────────────────────────────────


class TestPeakGeneLinks:
    def test_parse_peaks(self):
        from immunomics.analysis.peak_gene_links import _parse_peaks

        peaks = ["chr1:1000-2000", "chr2:5000-6000", "chrX:100-200"]
        df = _parse_peaks(peaks)
        assert len(df) == 3
        assert "chrom" in df.columns
        assert "start" in df.columns
        assert "center" in df.columns
        assert df.iloc[0]["chrom"] == "chr1"

    def test_find_nearby_pairs(self):
        from immunomics.analysis.peak_gene_links import _find_nearby_pairs

        peaks = pd.DataFrame(
            {
                "peak": ["chr1:1000-2000", "chr1:50000-51000"],
                "chrom": ["chr1", "chr1"],
                "start": [1000, 50000],
                "end": [2000, 51000],
                "center": [1500, 50500],
            }
        )
        genes = pd.DataFrame(
            {
                "gene_name": ["GeneA"],
                "chrom": ["chr1"],
                "tss": [2000],
            }
        )
        pairs = _find_nearby_pairs(peaks, genes, max_distance=10000)
        assert len(pairs) == 1
        assert pairs.iloc[0]["gene"] == "GeneA"

    def test_link_peaks_empty_annotation(self, mock_rna, mock_atac):
        from immunomics.analysis.peak_gene_links import link_peaks_to_genes

        empty_annot = pd.DataFrame(columns=["gene_name", "chrom", "tss"])
        result = link_peaks_to_genes(mock_rna, mock_atac, gene_annotation=empty_annot)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0


# ─── Visualization ───────────────────────────────────────────────────────────


class TestVisualization:
    def test_plot_tf_heatmap(self, mock_tf_activity_df):
        import matplotlib

        matplotlib.use("Agg")
        from immunomics.visualization.plots import plot_tf_heatmap

        fig = plot_tf_heatmap(mock_tf_activity_df, top_n=4)
        assert fig is not None
        import matplotlib.pyplot as plt

        plt.close(fig)

    def test_plot_peak_gene_links(self, mock_peak_gene_links_df):
        import matplotlib

        matplotlib.use("Agg")
        from immunomics.visualization.plots import plot_peak_gene_links

        fig = plot_peak_gene_links(mock_peak_gene_links_df, top_n=5)
        assert fig is not None
        import matplotlib.pyplot as plt

        plt.close(fig)

    def test_plot_integration_comparison(self):
        import matplotlib

        matplotlib.use("Agg")
        from immunomics.visualization.plots import plot_integration_comparison

        results = {
            "multivi": {"silhouette_score": 0.35, "adjusted_rand_index": 0.42},
            "mofa": {"silhouette_score": 0.30, "adjusted_rand_index": 0.38},
        }
        fig = plot_integration_comparison(results)
        assert fig is not None
        import matplotlib.pyplot as plt

        plt.close(fig)

    def test_plot_integration_comparison_with_error(self):
        import matplotlib

        matplotlib.use("Agg")
        from immunomics.visualization.plots import plot_integration_comparison

        results = {
            "multivi": {"silhouette_score": 0.35, "adjusted_rand_index": 0.42},
            "wnn": {"error": "R not available"},
        }
        fig = plot_integration_comparison(results)
        assert fig is not None
        import matplotlib.pyplot as plt

        plt.close(fig)

    def test_plot_tf_heatmap_save(self, mock_tf_activity_df, tmp_path):
        import matplotlib

        matplotlib.use("Agg")
        from immunomics.visualization.plots import plot_tf_heatmap

        path = str(tmp_path / "tf_heatmap.png")
        fig = plot_tf_heatmap(mock_tf_activity_df, save_path=path)
        import os

        import matplotlib.pyplot as plt

        assert os.path.exists(path)
        plt.close(fig)


# ─── Config ──────────────────────────────────────────────────────────────────


class TestConfig:
    def test_load_config(self):
        from immunomics.utils.config import load_config

        config = load_config("configs/integration.yaml")
        assert "data" in config
        assert "rna_preprocessing" in config
        assert "atac_preprocessing" in config
        assert "integration" in config

    def test_load_config_missing_file(self):
        from immunomics.utils.config import load_config

        with pytest.raises(FileNotFoundError):
            load_config("nonexistent.yaml")

    def test_config_integration_methods(self):
        from immunomics.utils.config import load_config

        config = load_config("configs/integration.yaml")
        assert "multivi" in config["integration"]
        assert "mofa" in config["integration"]
        assert "wnn" in config["integration"]


# ─── CLI ─────────────────────────────────────────────────────────────────────


class TestCLI:
    def test_main_help(self):
        from immunomics.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "ImmunOmics" in result.output

    def test_version_flag(self):
        from immunomics.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["--version"])
        assert result.exit_code == 0
        assert "0.1.0" in result.output

    def test_download_help(self):
        from immunomics.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["download", "--help"])
        assert result.exit_code == 0
        assert "Download" in result.output

    def test_integrate_help(self):
        from immunomics.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["integrate", "--help"])
        assert result.exit_code == 0
        assert "integration" in result.output.lower()

    def test_benchmark_help(self):
        from immunomics.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["benchmark", "--help"])
        assert result.exit_code == 0
        assert "Compare" in result.output

    def test_subcommands_registered(self):
        from immunomics.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["--help"])
        assert "download" in result.output
        assert "integrate" in result.output
        assert "benchmark" in result.output


# ─── Integration Module Imports ──────────────────────────────────────────────


class TestIntegrationImports:
    def test_benchmark_module_importable(self):
        from immunomics.integration.benchmark import compare_methods

        assert callable(compare_methods)

    def test_multivi_module_importable(self):
        from immunomics.integration.multivi import run_multivi

        assert callable(run_multivi)

    def test_mofa_module_importable(self):
        from immunomics.integration.mofa import run_mofa

        assert callable(run_mofa)

    def test_wnn_module_importable(self):
        from immunomics.integration.wnn import run_wnn

        assert callable(run_wnn)

    def test_subpackage_exports(self):
        from immunomics.integration import compare_methods, run_mofa, run_multivi, run_wnn

        assert all(callable(f) for f in [run_multivi, run_mofa, run_wnn, compare_methods])


# ImmunOmics v0.1.0
# Any usage is subject to this software's license.
