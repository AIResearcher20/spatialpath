import scanpy as sc


def normalize(adata, target_sum=1e4):
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.raw = adata
    return adata


def select_hvg(adata, n_top_genes=2000, flavor="seurat"):
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor=flavor)
    return adata


def subset_hvg(adata):
    if "highly_variable" not in adata.var:
        raise KeyError("highly_variable column not found")
    return adata[:, adata.var["highly_variable"]].copy()


def scale(adata, max_value=10):
    sc.pp.scale(adata, max_value=max_value)
    return adata
