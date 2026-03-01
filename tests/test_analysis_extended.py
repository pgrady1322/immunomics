#!/usr/bin/env python3
"""
ImmunOmics v0.1.0

Extended analysis tests — differential edge cases, TF activity, peak-gene links.

Author: Patrick Grady
License: MIT License - See LICENSE
"""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix


class TestDifferentialExtended:
    """Extended differential expression / accessibility tests."""

    def test_logreg_method(self, mock_rna):
        import scanpy as sc

        from immunomics.analysis.differential import differential_genes

        sc.pp.normalize_total(mock_rna)
        sc.pp.log1p(mock_rna)

        result = differential_genes(mock_rna, groupby="cell_type", method="logreg", n_genes=5)
        assert len(result) > 0

    def test_two_group_differential(self):
        """Differential expr with exactly 2 groups."""
        import scanpy as sc

        from immunomics.analysis.differential import differential_genes

        np.random.seed(42)
        n_cells, n_genes = 100, 50
        X = csr_matrix(np.random.rand(n_cells, n_genes).astype(np.float32))
        adata = ad.AnnData(X=X)
        adata.var_names = [f"Gene_{i}" for i in range(n_genes)]
        adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
        adata.obs["group"] = ["A"] * 50 + ["B"] * 50
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)

        result = differential_genes(adata, groupby="group", n_genes=10)
        assert len(result) > 0
        assert "group" in result.columns

    def test_differential_peaks_custom_method(self, mock_atac):
        import scanpy as sc

        from immunomics.analysis.differential import differential_peaks

        sc.pp.normalize_total(mock_atac)
        sc.pp.log1p(mock_atac)

        result = differential_peaks(mock_atac, groupby="cell_type", method="t-test", n_peaks=5)
        assert len(result) > 0

    def test_result_has_expected_columns(self, mock_rna):
        import scanpy as sc

        from immunomics.analysis.differential import differential_genes

        sc.pp.normalize_total(mock_rna)
        sc.pp.log1p(mock_rna)

        result = differential_genes(mock_rna, groupby="cell_type", n_genes=5)
        assert "names" in result.columns
        assert "pvals" in result.columns or "pvals_adj" in result.columns


class TestTFActivityExtended:
    """Extended TF activity inference tests."""

    def test_immune_tfs_all_categories(self):
        from immunomics.analysis.tf_activity import IMMUNE_TFS

        expected = {"T cell", "B cell", "Myeloid", "NK cell", "DC", "Inflammatory"}
        assert set(IMMUNE_TFS.keys()) == expected

    def test_immune_tfs_all_lists_nonempty(self):
        from immunomics.analysis.tf_activity import IMMUNE_TFS

        for category, tfs in IMMUNE_TFS.items():
            assert len(tfs) > 0, f"Category {category} has no TFs"

    def test_tf_not_in_rna(self, mock_rna_with_tfs, mock_atac):
        """TFs not present in RNA data should be silently skipped."""
        from immunomics.analysis.tf_activity import infer_tf_activity

        shared = mock_rna_with_tfs.obs_names[:80]
        mock_atac_sub = mock_atac[:80].copy()
        mock_atac_sub.obs_names = shared

        result = infer_tf_activity(
            mock_rna_with_tfs,
            mock_atac_sub,
            tf_list=["NONEXISTENT_TF", "TBX21"],
        )
        # Should still get results for TBX21
        assert len(result) > 0
        assert "TBX21" in result["tf_name"].values
        assert "NONEXISTENT_TF" not in result["tf_name"].values

    def test_default_tf_list(self, mock_rna_with_tfs, mock_atac):
        """Passing tf_list=None should use the curated IMMUNE_TFS."""
        from immunomics.analysis.tf_activity import infer_tf_activity

        shared = mock_rna_with_tfs.obs_names[:80]
        mock_atac_sub = mock_atac[:80].copy()
        mock_atac_sub.obs_names = shared

        result = infer_tf_activity(
            mock_rna_with_tfs,
            mock_atac_sub,
            tf_list=None,
        )
        assert len(result) > 0

    def test_mean_expression_present(self, mock_rna_with_tfs, mock_atac):
        from immunomics.analysis.tf_activity import infer_tf_activity

        shared = mock_rna_with_tfs.obs_names[:80]
        mock_atac_sub = mock_atac[:80].copy()
        mock_atac_sub.obs_names = shared

        result = infer_tf_activity(
            mock_rna_with_tfs,
            mock_atac_sub,
            tf_list=["TBX21", "GATA3"],
        )
        assert "mean_expression" in result.columns
        assert "pct_expressing" in result.columns

    def test_chromvar_fallback(self, mock_rna_with_tfs, mock_atac):
        """chromVAR method should fall back to expression-based when R is not available."""
        from immunomics.analysis.tf_activity import infer_tf_activity

        shared = mock_rna_with_tfs.obs_names[:80]
        mock_atac_sub = mock_atac[:80].copy()
        mock_atac_sub.obs_names = shared

        result = infer_tf_activity(
            mock_rna_with_tfs,
            mock_atac_sub,
            tf_list=["TBX21"],
            method="chromvar",
        )
        # Should still return results via fallback
        assert len(result) > 0


class TestPeakGeneLinksExtended:
    """Extended peak-gene linkage tests."""

    def test_parse_peaks_empty(self):
        from immunomics.analysis.peak_gene_links import _parse_peaks

        df = _parse_peaks([])
        assert len(df) == 0

    def test_parse_peaks_malformed(self):
        from immunomics.analysis.peak_gene_links import _parse_peaks

        df = _parse_peaks(["bad_name", "no_format", "123"])
        assert len(df) == 0

    def test_parse_peaks_mixed(self):
        from immunomics.analysis.peak_gene_links import _parse_peaks

        df = _parse_peaks(["chr1:100-200", "bad", "chrX:999-1500"])
        assert len(df) == 2
        assert df.iloc[0]["chrom"] == "chr1"
        assert df.iloc[1]["chrom"] == "chrX"

    def test_find_nearby_pairs_no_match(self):
        from immunomics.analysis.peak_gene_links import _find_nearby_pairs

        peaks = pd.DataFrame(
            {
                "peak": ["chr1:1000-2000"],
                "chrom": ["chr1"],
                "start": [1000],
                "end": [2000],
                "center": [1500],
            }
        )
        genes = pd.DataFrame(
            {
                "gene_name": ["GeneA"],
                "chrom": ["chr2"],  # Different chromosome
                "tss": [1500],
            }
        )
        pairs = _find_nearby_pairs(peaks, genes, max_distance=1000000)
        assert len(pairs) == 0

    def test_find_nearby_pairs_multiple_genes(self):
        from immunomics.analysis.peak_gene_links import _find_nearby_pairs

        peaks = pd.DataFrame(
            {
                "peak": ["chr1:1000-2000", "chr1:3000-4000"],
                "chrom": ["chr1", "chr1"],
                "start": [1000, 3000],
                "end": [2000, 4000],
                "center": [1500, 3500],
            }
        )
        genes = pd.DataFrame(
            {
                "gene_name": ["GeneA", "GeneB"],
                "chrom": ["chr1", "chr1"],
                "tss": [1600, 3600],
            }
        )
        pairs = _find_nearby_pairs(peaks, genes, max_distance=5000)
        # GeneA within 5kb of both peaks, GeneB within 5kb of both peaks
        assert len(pairs) >= 2

    def test_link_peaks_shared_barcodes(self, mock_rna, mock_atac):
        """link_peaks_to_genes should intersect on shared barcodes."""
        from immunomics.analysis.peak_gene_links import link_peaks_to_genes

        # Ensure some shared barcodes
        shared = min(mock_rna.n_obs, mock_atac.n_obs)
        mock_atac_sub = mock_atac[:shared].copy()
        mock_atac_sub.obs_names = mock_rna.obs_names[:shared]

        empty_annot = pd.DataFrame(columns=["gene_name", "chrom", "tss"])
        result = link_peaks_to_genes(
            mock_rna[:shared],
            mock_atac_sub,
            gene_annotation=empty_annot,
        )
        assert isinstance(result, pd.DataFrame)

    def test_get_gene_tss_fallback(self):
        """_get_gene_tss should return empty DataFrame on network failure."""
        from immunomics.analysis.peak_gene_links import _get_gene_tss

        # This will fail (pybiomart likely not installed, or network unavailable)
        # but should not raise — it returns an empty DataFrame
        result = _get_gene_tss(["NONEXISTENT_GENE_XYZ123"])
        assert isinstance(result, pd.DataFrame)
        assert "gene_name" in result.columns
