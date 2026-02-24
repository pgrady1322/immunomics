"""
ImmunOmics CLI — command-line interface for multi-omics integration.
"""

import logging
import click

from immunomics import __version__

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
    from immunomics.data.datasets import load_multiome_pbmc

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
    click.echo(f"Running {method} integration...")
    click.echo(f"Results will be saved to {output}/")


@main.command()
@click.option("--config", "-c", default="configs/integration.yaml", help="Config YAML")
@click.option("--output", "-o", default="results", help="Output directory")
def benchmark(config, output):
    """Compare all integration methods."""
    click.echo("Running integration benchmark...")


if __name__ == "__main__":
    main()
