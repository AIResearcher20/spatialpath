import pandas as pd
import scanpy as sc


def reduce_dimensions(adata, n_comps=50, n_neighbors=15):
    n_comps = min(n_comps, adata.n_obs - 1, adata.n_vars - 1)
    sc.tl.pca(adata, n_comps=n_comps)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_pca")
    return adata


def cluster_spots(adata, resolution=1.0, method="leiden", key_added="cluster"):
    if method == "leiden":
        sc.tl.leiden(adata, resolution=resolution, key_added=key_added)
    elif method == "louvain":
        sc.tl.louvain(adata, resolution=resolution, key_added=key_added)
    else:
        raise ValueError(f"unknown method: {method}")
    return adata


def find_marker_genes(adata, groupby="cluster", n_genes=25, method="wilcoxon"):
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        n_genes=n_genes,
        method=method,
    )
    return adata


def top_marker_genes(adata, n=5, groupby="cluster"):
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    markers = {}
    for g in groups:
        names = result["names"][g][:n].tolist()
        scores = result["scores"][g][:n].tolist()
        markers[str(g)] = list(zip(names, scores))
    return markers


def marker_table(adata, n=10, groupby="cluster"):
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    rows = []
    for g in groups:
        names = result["names"][g][:n]
        scores = result["scores"][g][:n]
        for name, score in zip(names, scores):
            rows.append({"cluster": str(g), "gene": name, "score": float(score)})
    return pd.DataFrame(rows)
