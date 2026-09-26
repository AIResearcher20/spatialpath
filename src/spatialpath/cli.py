import argparse
from pathlib import Path

import yaml

from .cluster import cluster_spots, find_marker_genes, reduce_dimensions
from .io import load_visium, save_h5ad
from .normalize import normalize, select_hvg, subset_hvg
from .plot import plot_clusters, plot_qc
from .qc import compute_qc_metrics, filter_spots, qc_summary
from .report import build_report, save_report


def _load_config(path):
    with open(path, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def run(data_dir, config, output_dir):
    cfg = _load_config(config)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_dir / "figures"

    print(f"loading data from {data_dir}")
    adata = load_visium(data_dir)

    print("running qc")
    adata = compute_qc_metrics(adata, mito_prefix=cfg["qc"]["mito_prefix"])
    adata = filter_spots(
        adata,
        min_genes=cfg["qc"]["min_genes"],
        min_counts=cfg["qc"]["min_counts"],
        max_mito=cfg["qc"]["max_mito"],
    )
    stats = qc_summary(adata)
    print(f"spots={stats['n_spots']} genes={stats['n_genes']}")

    print("normalizing")
    adata = normalize(adata, target_sum=cfg["normalize"]["target_sum"])
    adata = select_hvg(adata, n_top_genes=cfg["normalize"]["n_top_genes"])
    adata_hvg = subset_hvg(adata)

    print("clustering")
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

    print("saving outputs")
    save_h5ad(adata, output_dir / "processed.h5ad")
    plot_qc(adata, fig_dir / "qc.png")
    plot_clusters(adata, output_path=fig_dir / "clusters.png")

    report = build_report(stats, adata)
    save_report(report, output_dir / "report.json")

    print(f"done: {output_dir}")


def main():
    parser = argparse.ArgumentParser(description="SpatialPath pipeline")
    parser.add_argument("data_dir", type=str, help="Visium data directory")
    parser.add_argument(
        "--config",
        type=str,
        default="configs/default.yaml",
        help="YAML configuration file",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results",
        help="Output directory",
    )
    args = parser.parse_args()

    run(args.data_dir, args.config, args.output_dir)


if __name__ == "__main__":
    main()
