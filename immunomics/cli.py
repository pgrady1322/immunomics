#!/usr/bin/env python3
"""
ImmunOmics v0.1.0

Click CLI — preprocess, integrate, analyze, and visualize commands.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import logging

import click

from . import __version__

logger = logging.getLogger(__name__)


@click.group()
@click.version_option(version=__version__)
@click.option("--verbose", "-v", is_flag=True, help="Enable verbose logging")
def main(verbose):
    """ImmunOmics — Multi-omics integration of immune cell states."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )


@main.command()
@click.option("--dataset", "-d", default="multiome_pbmc_10k", help="Dataset to download")
@click.option("--output", "-o", default=None, help="Output directory")
def download(dataset, output):
    """Download multi-omics dataset."""
    from .data.datasets import load_multiome_pbmc

    click.echo(f"Downloading {dataset}...")
    adata_rna, adata_atac = load_multiome_pbmc(cache_dir=output)
    click.echo(
        f"RNA: {adata_rna.n_obs} cells, {adata_rna.n_vars} genes | "
        f"ATAC: {adata_atac.n_obs} cells, {adata_atac.n_vars} peaks"
    )


@main.command()
@click.option("--config", "-c", default="configs/integration.yaml", help="Config YAML")
@click.option("--method", "-m", default="multivi", help="Integration method: multivi, wnn, mofa")
@click.option("--output", "-o", default="results", help="Output directory")
def integrate(config, method, output):
    """Run multi-omics integration."""
    from pathlib import Path

    from .data.datasets import load_multiome_pbmc
    from .data.preprocess_atac import preprocess_atac
    from .data.preprocess_rna import preprocess_rna
    from .utils.config import load_config

    cfg = load_config(config)
    click.echo(f"Running {method} integration...")

    # Load data
    adata_rna, adata_atac = load_multiome_pbmc()
    adata_rna = preprocess_rna(adata_rna, **cfg.get("rna_preprocessing", {}))
    adata_atac = preprocess_atac(adata_atac, **cfg.get("atac_preprocessing", {}))

    # Run integration
    runners = {
        "multivi": lambda: __import__(
            "immunomics.integration.multivi", fromlist=["run_multivi"]
        ).run_multivi(adata_rna, adata_atac, **cfg.get("integration", {}).get("multivi", {})),
        "mofa": lambda: __import__("immunomics.integration.mofa", fromlist=["run_mofa"]).run_mofa(
            adata_rna, adata_atac, **cfg.get("integration", {}).get("mofa", {})
        ),
        "wnn": lambda: __import__("immunomics.integration.wnn", fromlist=["run_wnn"]).run_wnn(
            adata_rna, adata_atac
        ),
    }
    if method not in runners:
        click.echo(f"Unknown method: {method}. Choose from: {list(runners.keys())}")
        return

    adata_int = runners[method]()

    # Save results
    out_dir = Path(output)
    out_dir.mkdir(parents=True, exist_ok=True)
    adata_int.write(out_dir / f"{method}_integrated.h5ad")
    click.echo(f"Saved integrated result to {out_dir / f'{method}_integrated.h5ad'}")


@main.command()
@click.option("--config", "-c", default="configs/integration.yaml", help="Config YAML")
@click.option("--output", "-o", default="results", help="Output directory")
def benchmark(config, output):
    """Compare all integration methods."""
    from pathlib import Path

    from .data.datasets import load_multiome_pbmc
    from .data.preprocess_atac import preprocess_atac
    from .data.preprocess_rna import preprocess_rna
    from .integration.benchmark import compare_methods
    from .utils.config import load_config

    cfg = load_config(config)
    click.echo("Running integration benchmark...")

    adata_rna, adata_atac = load_multiome_pbmc()
    adata_rna = preprocess_rna(adata_rna, **cfg.get("rna_preprocessing", {}))
    adata_atac = preprocess_atac(adata_atac, **cfg.get("atac_preprocessing", {}))

    methods = [k for k, v in cfg.get("integration", {}).items() if v.get("enabled", True)]
    results = compare_methods(adata_rna, adata_atac, methods=methods)

    # Display results
    for method, metrics in results.items():
        click.echo(f"\n{method.upper()}:")
        for k, v in metrics.items():
            click.echo(f"  {k}: {v}")

    # Save
    out_dir = Path(output)
    out_dir.mkdir(parents=True, exist_ok=True)
    import json

    with open(out_dir / "benchmark_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    click.echo(f"\nSaved benchmark results to {out_dir / 'benchmark_results.json'}")


if __name__ == "__main__":
    main()

# ImmunOmics v0.1.0
# Any usage is subject to this software's license.
