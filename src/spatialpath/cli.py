import argparse
import json
from pathlib import Path

import yaml

from .cluster import (
    cluster_spots,
    find_marker_genes,
    marker_table,
    reduce_dimensions,
    top_marker_genes,
)
from .io import load_visium, save_h5ad
from .normalize import normalize, select_hvg, subset_hvg
from .plot import (
    plot_clusters,
    plot_gene,
    plot_marker_dotplot,
    plot_qc,
)
from .qc import compute_qc_metrics, filter_spots, qc_summary
from .report import build_report, save_report


def _load_config(path):
    with open(path, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def run(data_dir, config, output_dir):
    data_path = Path(data_dir)
    if not data_path.is_dir():
        raise NotADirectoryError(f"Data directory not found: {data_path}")

    cfg = _load_config(config)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    print(f"loading data from {data_path}")
    adata = load_visium(data_path)

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
    adata.uns["rank_genes_groups"] = adata_hvg.uns["rank_genes_groups"]

    print("marker genes")
    markers = top_marker_genes(adata, n=5)
    marker_df = marker_table(adata, n=10)
    marker_df.to_csv(output_dir / "marker_genes.csv", index=False)

    top_genes = []
    for cluster, gene_list in markers.items():
        top_genes.append(gene_list[0][0])

    print("saving outputs")
    save_h5ad(adata, output_dir / "processed.h5ad")

    plot_qc(adata, fig_dir / "qc.png")
    plot_clusters(adata, output_path=fig_dir / "clusters.png")
    plot_marker_dotplot(
        adata,
        top_genes,
        output_path=fig_dir / "markers.png",
    )
    if top_genes:
        plot_gene(adata, top_genes[0], output_path=fig_dir / "top_gene.png")

    qc_report = {
        "n_spots": stats["n_spots"],
        "n_genes": stats["n_genes"],
        "median_genes_per_spot": stats["median_genes_per_spot"],
        "median_counts_per_spot": stats["median_counts_per_spot"],
        "median_mito_pct": stats["median_mito_pct"],
    }
    with open(output_dir / "qc_metrics.json", "w") as f:
        json.dump(qc_report, f, indent=2)

    cluster_sizes = adata.obs["cluster"].value_counts().sort_index()
    clustering_report = {
        "n_clusters": int(cluster_sizes.shape[0]),
        "sizes": {str(k): int(v) for k, v in cluster_sizes.items()},
    }
    with open(output_dir / "clustering_metrics.json", "w") as f:
        json.dump(clustering_report, f, indent=2)

    marker_json = {
        str(cluster): [{"gene": g, "score": s} for g, s in genes]
        for cluster, genes in markers.items()
    }
    with open(output_dir / "marker_genes.json", "w") as f:
        json.dump(marker_json, f, indent=2)

    summary = {
        "n_spots": stats["n_spots"],
        "n_genes": stats["n_genes"],
        "n_clusters": int(cluster_sizes.shape[0]),
        "top_markers": {c: genes[0][0] for c, genes in markers.items()},
    }
    with open(output_dir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    report = build_report(stats, adata)
    save_report(report, output_dir / "report.json")

    print(f"done: {output_dir}")


def main():
    parser = argparse.ArgumentParser(description="SpatialPath pipeline")
    parser.add_argument("data_dir", type=str)
    parser.add_argument("--config", type=str, default="configs/default.yaml")
    parser.add_argument("--output-dir", type=str, default="results")
    args = parser.parse_args()

    run(args.data_dir, args.config, args.output_dir)


if __name__ == "__main__":
    main()
