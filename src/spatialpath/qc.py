import scanpy as sc


def compute_qc_metrics(adata, mito_prefix="MT-"):
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(mito_prefix)
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    return adata


def filter_spots(adata, min_genes=200, min_counts=500, max_mito=20.0):
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs["pct_counts_mt"] < max_mito].copy()
    return adata


def qc_summary(adata):
    return {
        "n_spots": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "median_genes_per_spot": float(adata.obs["n_genes_by_counts"].median()),
        "median_counts_per_spot": float(adata.obs["total_counts"].median()),
        "median_mito_pct": float(adata.obs["pct_counts_mt"].median()),
    }
