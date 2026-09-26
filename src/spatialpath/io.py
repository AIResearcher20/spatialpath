from pathlib import Path

import anndata as ad
import scanpy as sc


def load_visium(data_dir, library_id="visium"):
    data_dir = Path(data_dir)
    count_file = data_dir / "filtered_feature_bc_matrix.h5"

    if not count_file.exists():
        raise FileNotFoundError(f"missing file: {count_file}")

    adata = sc.read_visium(
        data_dir,
        count_file=count_file.name,
        library_id=library_id,
    )
    adata.var_names_make_unique()
    return adata


def load_h5ad(path):
    return ad.read_h5ad(path)


def save_h5ad(adata, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path)
