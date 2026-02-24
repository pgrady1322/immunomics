"""
Publication-quality visualizations for multi-omics immune analysis.

Includes joint UMAP plots, TF activity heatmaps, peak-gene link
visualization, and integration method comparison.
"""

import logging
from typing import Optional, List, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

logger = logging.getLogger(__name__)


def plot_joint_umap(
    adata,
    color_keys: Optional[List[str]] = None,
    figsize: tuple = (18, 5),
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    Plot joint UMAP colored by cell type, modality weights, etc.

    Parameters
    ----------
    adata
        Integrated AnnData with UMAP coordinates
    color_keys
        Columns in obs to color by. Default: ['cell_type', 'leiden']
    figsize
        Figure size
    save_path
        If provided, save figure

    Returns
    -------
    matplotlib Figure
    """
    import scanpy as sc

    if color_keys is None:
        color_keys = [c for c in ["cell_type", "leiden"] if c in adata.obs.columns]

    n_panels = len(color_keys)
    fig, axes = plt.subplots(1, n_panels, figsize=figsize)
    if n_panels == 1:
        axes = [axes]

    for ax, key in zip(axes, color_keys):
        sc.pl.umap(adata, color=key, ax=ax, show=False, title=key)

    plt.suptitle("Multi-Omics Integrated UMAP", fontsize=14, y=1.02)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig


def plot_tf_heatmap(
    tf_activity: pd.DataFrame,
    top_n: int = 20,
    figsize: tuple = (12, 8),
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    Plot TF activity heatmap across cell types.

    Parameters
    ----------
    tf_activity
        DataFrame from infer_tf_activity()
    top_n
        Number of top TFs to show (by variance across cell types)
    figsize
        Figure size
    save_path
        If provided, save figure

    Returns
    -------
    matplotlib Figure
    """
    # Pivot to matrix
    score_col = "activity_score" if "activity_score" in tf_activity.columns else "expression_zscore"

    matrix = tf_activity.pivot_table(
        index="tf_name",
        columns="cell_type",
        values=score_col,
        aggfunc="mean",
    )

    # Select top TFs by variance
    if len(matrix) > top_n:
        variances = matrix.var(axis=1)
        top_tfs = variances.nlargest(top_n).index
        matrix = matrix.loc[top_tfs]

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        matrix,
        cmap="RdBu_r",
        center=0,
        annot=True,
        fmt=".2f",
        ax=ax,
        cbar_kws={"label": "Activity Score"},
    )
    ax.set_title("Transcription Factor Activity Across Immune Cell Types", fontsize=14)
    ax.set_xlabel("Cell Type", fontsize=12)
    ax.set_ylabel("Transcription Factor", fontsize=12)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig


def plot_peak_gene_links(
    links: pd.DataFrame,
    top_n: int = 30,
    figsize: tuple = (10, 8),
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    Plot top peak-to-gene linkages.

    Parameters
    ----------
    links
        DataFrame from link_peaks_to_genes()
    top_n
        Number of top links to display
    figsize
        Figure size
    save_path
        If provided, save figure

    Returns
    -------
    matplotlib Figure
    """
    top_links = links.nlargest(top_n, "correlation")

    fig, ax = plt.subplots(figsize=figsize)

    colors = top_links["correlation"].values
    scatter = ax.scatter(
        range(len(top_links)),
        top_links["correlation"],
        c=colors,
        cmap="viridis",
        s=100,
        edgecolors="black",
        linewidth=0.5,
    )

    ax.set_xticks(range(len(top_links)))
    ax.set_xticklabels(
        [f"{row['gene']}\n{row['peak'][:20]}" for _, row in top_links.iterrows()],
        rotation=90,
        fontsize=8,
    )
    ax.set_ylabel("Correlation (peak accessibility ↔ gene expression)", fontsize=12)
    ax.set_title("Top Peak-to-Gene Regulatory Links", fontsize=14)
    plt.colorbar(scatter, label="Correlation")
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig


def plot_integration_comparison(
    results: Dict[str, Dict],
    metrics: Optional[List[str]] = None,
    figsize: tuple = (10, 6),
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    Compare integration methods side by side.

    Parameters
    ----------
    results
        Output from compare_methods()
    metrics
        Metrics to compare
    figsize
        Figure size
    save_path
        If provided, save figure

    Returns
    -------
    matplotlib Figure
    """
    if metrics is None:
        metrics = ["silhouette_score", "adjusted_rand_index"]

    rows = []
    for method, method_results in results.items():
        if "error" in method_results:
            continue
        for metric in metrics:
            if metric in method_results:
                rows.append({
                    "Method": method.upper(),
                    "Metric": metric.replace("_", " ").title(),
                    "Score": method_results[metric],
                })

    df = pd.DataFrame(rows)

    fig, ax = plt.subplots(figsize=figsize)
    sns.barplot(data=df, x="Method", y="Score", hue="Metric", ax=ax, palette="Set2")
    ax.set_title("Integration Method Comparison", fontsize=14)
    ax.set_ylabel("Score", fontsize=12)
    ax.legend(title="Metric", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig
