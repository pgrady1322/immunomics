"""
Peak-to-gene linkage analysis.

Links scATAC-seq peaks to putative target genes by correlating
peak accessibility with gene expression across cells. This
identifies cell-type-specific regulatory elements.
"""

import logging
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats

logger = logging.getLogger(__name__)


def link_peaks_to_genes(
    adata_rna: ad.AnnData,
    adata_atac: ad.AnnData,
    max_distance: int = 500_000,
    min_correlation: float = 0.1,
    pvalue_threshold: float = 0.05,
    gene_annotation: Optional[pd.DataFrame] = None,
    n_jobs: int = -1,
) -> pd.DataFrame:
    """
    Link ATAC peaks to target genes by correlation analysis.

    For each peak within `max_distance` of a gene's TSS, computes
    the Pearson correlation between peak accessibility and gene
    expression across cells. Significant positive correlations
    suggest regulatory relationships.

    Parameters
    ----------
    adata_rna
        RNA AnnData (log-normalized expression)
    adata_atac
        ATAC AnnData (raw or TF-IDF counts)
    max_distance
        Maximum distance (bp) from peak to gene TSS
    min_correlation
        Minimum absolute correlation to report
    pvalue_threshold
        P-value cutoff after BH correction
    gene_annotation
        DataFrame with columns: gene_name, chrom, tss_position.
        If None, fetched from pybiomart.
    n_jobs
        Number of parallel jobs (-1 = all cores)

    Returns
    -------
    DataFrame with columns:
        peak, gene, correlation, pvalue, pvalue_adj, distance, chrom
    """
    logger.info("Linking peaks to genes by correlation analysis...")

    # Ensure shared barcodes
    shared = adata_rna.obs_names.intersection(adata_atac.obs_names)
    adata_rna = adata_rna[shared].copy()
    adata_atac = adata_atac[shared].copy()

    logger.info(f"Using {len(shared)} shared cells")

    # Get gene TSS positions
    if gene_annotation is None:
        gene_annotation = _get_gene_tss(adata_rna.var_names.tolist())

    # Parse peak coordinates
    peak_coords = _parse_peaks(adata_atac.var_names.tolist())

    # Find peak-gene pairs within distance
    pairs = _find_nearby_pairs(peak_coords, gene_annotation, max_distance)
    logger.info(f"Found {len(pairs)} peak-gene pairs within {max_distance/1000:.0f}kb")

    if len(pairs) == 0:
        return pd.DataFrame(columns=["peak", "gene", "correlation", "pvalue", "distance", "chrom"])

    # Compute correlations
    import scipy.sparse as sp

    rna_X = adata_rna.X
    atac_X = adata_atac.X
    if sp.issparse(rna_X):
        rna_X = rna_X.toarray()
    if sp.issparse(atac_X):
        atac_X = atac_X.toarray()

    results = []
    gene_idx_map = {g: i for i, g in enumerate(adata_rna.var_names)}
    peak_idx_map = {p: i for i, p in enumerate(adata_atac.var_names)}

    for _, row in pairs.iterrows():
        gene = row["gene"]
        peak = row["peak"]

        if gene not in gene_idx_map or peak not in peak_idx_map:
            continue

        gene_expr = rna_X[:, gene_idx_map[gene]]
        peak_acc = atac_X[:, peak_idx_map[peak]]

        # Skip if zero variance
        if np.std(gene_expr) == 0 or np.std(peak_acc) == 0:
            continue

        r, p = stats.pearsonr(gene_expr.flatten(), peak_acc.flatten())

        if abs(r) >= min_correlation:
            results.append({
                "peak": peak,
                "gene": gene,
                "correlation": r,
                "pvalue": p,
                "distance": row["distance"],
                "chrom": row["chrom"],
            })

    df = pd.DataFrame(results)

    if len(df) > 0:
        # Multiple testing correction (Benjamini-Hochberg)
        from statsmodels.stats.multitest import multipletests

        _, pvals_adj, _, _ = multipletests(df["pvalue"].values, method="fdr_bh")
        df["pvalue_adj"] = pvals_adj
        df = df[df["pvalue_adj"] < pvalue_threshold].copy()
        df = df.sort_values("correlation", ascending=False)

    logger.info(f"Significant peak-gene links: {len(df)}")
    return df


def _parse_peaks(peak_names: list) -> pd.DataFrame:
    """Parse peak names (chr:start-end) to coordinates."""
    import re

    rows = []
    for name in peak_names:
        match = re.match(r"(chr\w+)[:\-](\d+)[:\-](\d+)", name)
        if match:
            rows.append({
                "peak": name,
                "chrom": match.group(1),
                "start": int(match.group(2)),
                "end": int(match.group(3)),
                "center": (int(match.group(2)) + int(match.group(3))) // 2,
            })
    return pd.DataFrame(rows)


def _get_gene_tss(gene_names: list) -> pd.DataFrame:
    """Get TSS positions for genes from pybiomart."""
    try:
        from pybiomart import Server

        server = Server(host="http://www.ensembl.org")
        dataset = server.marts["ENSEMBL_MART_ENSEMBL"].datasets["hsapiens_gene_ensembl"]

        result = dataset.query(
            attributes=[
                "external_gene_name",
                "chromosome_name",
                "transcription_start_site",
            ],
            filters={"external_gene_name": gene_names},
        )

        result.columns = ["gene_name", "chrom", "tss"]
        result["chrom"] = "chr" + result["chrom"].astype(str)
        result = result.drop_duplicates(subset=["gene_name"])

        return result

    except Exception as e:
        logger.warning(f"Could not fetch gene annotations: {e}")
        logger.info("Returning empty annotation. Provide gene_annotation manually.")
        return pd.DataFrame(columns=["gene_name", "chrom", "tss"])


def _find_nearby_pairs(
    peaks: pd.DataFrame,
    genes: pd.DataFrame,
    max_distance: int,
) -> pd.DataFrame:
    """Find peak-gene pairs within max_distance."""
    pairs = []

    for _, gene_row in genes.iterrows():
        chrom = gene_row["chrom"]
        tss = gene_row["tss"]
        gene_name = gene_row["gene_name"]

        # Filter peaks on same chromosome
        chrom_peaks = peaks[peaks["chrom"] == chrom]

        # Filter by distance
        distances = np.abs(chrom_peaks["center"].values - tss)
        nearby = chrom_peaks[distances <= max_distance].copy()
        nearby["gene"] = gene_name
        nearby["distance"] = np.abs(nearby["center"] - tss)

        pairs.append(nearby[["peak", "gene", "distance", "chrom"]])

    if pairs:
        return pd.concat(pairs, ignore_index=True)
    return pd.DataFrame(columns=["peak", "gene", "distance", "chrom"])
