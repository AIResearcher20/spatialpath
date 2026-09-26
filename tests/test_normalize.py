import numpy as np

from spatialpath.normalize import normalize, scale, select_hvg, subset_hvg


def test_normalize_shape(adata):
    out = normalize(adata)
    assert out.shape == adata.shape
    assert out.raw is not None


def test_select_hvg_marks(adata):
    out = normalize(adata)
    out = select_hvg(out, n_top_genes=30)
    assert "highly_variable" in out.var
    assert out.var["highly_variable"].sum() > 0


def test_subset_hvg(adata):
    out = normalize(adata)
    out = select_hvg(out, n_top_genes=20)
    sub = subset_hvg(out)
    assert sub.n_vars <= out.n_vars
    assert sub.n_vars > 0


def test_scale_bounds(adata):
    out = normalize(adata)
    out = select_hvg(out, n_top_genes=50)
    sub = subset_hvg(out)
    sub = scale(sub, max_value=10)
    values = np.asarray(sub.X)
    assert values.max() <= 10
