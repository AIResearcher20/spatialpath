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
