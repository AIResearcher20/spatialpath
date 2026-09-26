from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc


def _save(path):
    if path is None:
        return
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_qc(adata, output_path=None):
    sc.pl.spatial(
        adata,
        color=["total_counts", "n_genes_by_counts", "pct_counts_mt"],
        show=False,
    )
    _save(output_path)


def plot_clusters(adata, cluster_key="cluster", output_path=None):
    sc.pl.spatial(adata, color=cluster_key, show=False)
    _save(output_path)


def plot_gene(adata, gene, output_path=None):
    if gene not in adata.var_names:
        raise KeyError(f"gene not found: {gene}")
    sc.pl.spatial(adata, color=gene, show=False)
    _save(output_path)


def plot_marker_dotplot(adata, markers, cluster_key="cluster", output_path=None):
    sc.pl.dotplot(
        adata,
        markers,
        groupby=cluster_key,
        show=False,
    )
    _save(output_path)


def plot_marker_heatmap(adata, markers, cluster_key="cluster", output_path=None):
    sc.pl.heatmap(
        adata,
        markers,
        groupby=cluster_key,
        show=False,
    )
    _save(output_path)
