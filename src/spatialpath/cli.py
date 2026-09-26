from pathlib import Path

import typer
import yaml

from .cluster import cluster_spots, find_marker_genes, reduce_dimensions
from .io import load_visium, save_h5ad
from .normalize import normalize, select_hvg, subset_hvg
from .plot import plot_clusters, plot_qc
from .qc import compute_qc_metrics, filter_spots, qc_summary
from .report import build_report, save_report

app = typer.Typer(add_completion=False, help="SpatialPath pipeline")


def _load_config(path):
    with open(path, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


@app.command("run")
def run(
    data_dir: Path = typer.Argument(..., exists=True, file_okay=False),
    config: Path = typer.Option(Path("configs/default.yaml"), exists=True),
    output_dir: Path = typer.Option(Path("results")),
):
    cfg = _load_config(config)
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_dir / "figures"

    typer.echo(f"loading data from {data_dir}")
    adata = load_visium(data_dir)

    typer.echo("running qc")
    adata = compute_qc_metrics(adata, mito_prefix=cfg["qc"]["mito_prefix"])
    adata = filter_spots(
        adata,
        min_genes=cfg["qc"]["min_genes"],
        min_counts=cfg["qc"]["min_counts"],
        max_mito=cfg["qc"]["max_mito"],
    )
    stats = qc_summary(adata)
    typer.echo(f"spots={stats['n_spots']} genes={stats['n_genes']}")

    typer.echo("normalizing")
    adata = normalize(adata, target_sum=cfg["normalize"]["target_sum"])
    adata = select_hvg(adata, n_top_genes=cfg["normalize"]["n_top_genes"])
    adata_hvg = subset_hvg(adata)

    typer.echo("clustering")
    adata_hvg = reduce_dimensions(
        adata_hvg,
        n_comps=cfg["cluster"]["n_comps"],
        n_neighbors=cfg["cluster"]["n_neighbors"],
    )
    adata_hvg = cluster_spots(
        adata_hvg,
        resolution=cfg["cluster"]["resolution"],
        method=cfg["cluster"]["method"],
    )
    adata_hvg = find_marker_genes(adata_hvg, groupby="cluster")

    adata.obs["cluster"] = adata_hvg.obs["cluster"].values
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]

    typer.echo("saving outputs")
    save_h5ad(adata, output_dir / "processed.h5ad")
    plot_qc(adata, fig_dir / "qc.png")
    plot_clusters(adata, output_path=fig_dir / "clusters.png")

    report = build_report(stats, adata)
    save_report(report, output_dir / "report.json")

    typer.echo(f"done: {output_dir}")


if __name__ == "__main__":
    app()
