from spatialpath.io import load_h5ad, save_h5ad


def test_roundtrip(adata, tmp_path):
    path = tmp_path / "x.h5ad"
    save_h5ad(adata, path)
    loaded = load_h5ad(path)
    assert loaded.shape == adata.shape
    assert list(loaded.var_names) == list(adata.var_names)


def test_save_creates_parent_dirs(adata, tmp_path):
    path = tmp_path / "nested" / "deep" / "x.h5ad"
    save_h5ad(adata, path)
    assert path.exists()
