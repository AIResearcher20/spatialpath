from spatialpath.cluster import cluster_spots, reduce_dimensions
from spatialpath.normalize import normalize, select_hvg, subset_hvg


def _prep(adata):
    out = normalize(adata)
    out = select_hvg(out, n_top_genes=30)
    out = subset_hvg(out)
    return out


def test_reduce_dimensions(adata):
    out = _prep(adata)
    out = reduce_dimensions(out, n_comps=10, n_neighbors=5)
    assert "X_pca" in out.obsm


def test_cluster_adds_key(adata):
    out = _prep(adata)
    out = reduce_dimensions(out, n_comps=10, n_neighbors=5)
    out = cluster_spots(out, resolution=0.5)
    assert "cluster" in out.obs


def test_invalid_method_raises(adata):
    out = _prep(adata)
    out = reduce_dimensions(out, n_comps=10, n_neighbors=5)
    try:
        cluster_spots(out, method="kmeans")
    except ValueError:
        return
    raise AssertionError("expected ValueError")
