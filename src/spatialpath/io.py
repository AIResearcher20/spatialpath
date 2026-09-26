from pathlib import Path

import anndata as ad
import squidpy as sq


def load_visium(data_dir, library_id="visium"):
    data_dir = Path(data_dir)
    adata = sq.read.visium(data_dir, library_id=library_id)
    adata.var_names_make_unique()
    return adata


def load_h5ad(path):
    return ad.read_h5ad(path)


def save_h5ad(adata, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path)
