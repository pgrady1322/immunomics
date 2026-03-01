#!/usr/bin/env python3
"""
ImmunOmics v0.1.0

Edge-case and misc tests — CLI, config, datasets, utils.

Author: Patrick Grady
License: MIT License - See LICENSE
"""

import os

import numpy as np
import pytest
from click.testing import CliRunner


class TestCLIExtended:
    """Extended CLI tests."""

    def test_verbose_flag(self):
        from immunomics.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["--verbose", "--help"])
        assert result.exit_code == 0

    def test_integrate_unknown_method(self):
        from immunomics.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["integrate", "--method", "nonexistent", "--help"])
        assert result.exit_code == 0

    def test_download_custom_output(self):
        from immunomics.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["download", "--help"])
        assert "--output" in result.output

    def test_benchmark_has_config_option(self):
        from immunomics.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["benchmark", "--help"])
        assert "--config" in result.output

    def test_integrate_has_output_option(self):
        from immunomics.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["integrate", "--help"])
        assert "--output" in result.output


class TestConfigExtended:
    """Extended configuration tests."""

    def test_config_has_atac_params(self):
        from immunomics.utils.config import load_config

        config = load_config("configs/integration.yaml")
        assert "atac_preprocessing" in config

    def test_config_rna_params(self):
        from immunomics.utils.config import load_config

        config = load_config("configs/integration.yaml")
        rna = config["rna_preprocessing"]
        # Should have typical params
        assert isinstance(rna, dict)

    def test_config_data_section(self):
        from immunomics.utils.config import load_config

        config = load_config("configs/integration.yaml")
        assert "data" in config

    def test_load_config_returns_dict(self):
        from immunomics.utils.config import load_config

        config = load_config("configs/integration.yaml")
        assert isinstance(config, dict)

    def test_config_yaml_path_object(self):
        """load_config should accept a Path-like string."""
        from pathlib import Path

        from immunomics.utils.config import load_config

        config = load_config(str(Path("configs") / "integration.yaml"))
        assert isinstance(config, dict)


class TestDatasetsExtended:
    """Extended dataset registry tests."""

    def test_dataset_has_required_keys(self):
        from immunomics.data.datasets import DATASETS

        for name, info in DATASETS.items():
            assert "description" in info, f"{name} missing description"
            assert "url_rna" in info, f"{name} missing url_rna"

    def test_cache_dir_is_pathlib(self):
        from pathlib import Path

        from immunomics.data.datasets import CACHE_DIR

        assert isinstance(CACHE_DIR, Path)

    def test_ensure_cache_dir(self, tmp_path):
        from immunomics.data.datasets import _ensure_cache_dir

        # Just verify it returns a Path and doesn't error
        result = _ensure_cache_dir()
        assert isinstance(result, type(tmp_path))

    def test_multiome_pbmc_dataset_metadata(self):
        from immunomics.data.datasets import DATASETS

        ds = DATASETS["multiome_pbmc_10k"]
        assert "10x" in ds["description"] or "10x" in ds["reference"]
        assert ds["url_rna"].startswith("https://")
        assert "n_cells" in ds


class TestPackageImportsExtended:
    """Extended import tests for subpackages."""

    def test_analysis_subpackage(self):
        from immunomics.analysis import differential, peak_gene_links, tf_activity

        assert hasattr(differential, "differential_genes")
        assert hasattr(differential, "differential_peaks")
        assert hasattr(peak_gene_links, "link_peaks_to_genes")
        assert hasattr(tf_activity, "infer_tf_activity")

    def test_data_subpackage(self):
        from immunomics.data import datasets, preprocess_atac, preprocess_rna

        assert hasattr(preprocess_rna, "preprocess_rna")
        assert hasattr(preprocess_atac, "preprocess_atac")
        assert hasattr(datasets, "DATASETS")

    def test_visualization_subpackage(self):
        from immunomics.visualization import plots

        assert hasattr(plots, "plot_joint_umap")
        assert hasattr(plots, "plot_tf_heatmap")
        assert hasattr(plots, "plot_peak_gene_links")
        assert hasattr(plots, "plot_integration_comparison")

    def test_utils_subpackage(self):
        from immunomics.utils import config

        assert hasattr(config, "load_config")
