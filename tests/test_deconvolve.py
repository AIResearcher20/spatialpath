import numpy as np

from spatialpath.deconvolve import build_signatures, deconvolve
from spatialpath.normalize import normalize


def test_signature_shape(reference):
    sig = build_signatures(reference, cell_type_key="cell_type")
    assert sig.shape[1] == reference.obs["cell_type"].nunique()
    assert sig.shape[0] == reference.n_vars


def test_proportions_sum_to_one(adata, reference):
    a = normalize(adata)
    r = normalize(reference)
    props = deconvolve(a, r, cell_type_key="cell_type")
    sums = props.values.sum(axis=1)
    assert np.allclose(sums, 1.0, atol=1e-6)


def test_proportions_non_negative(adata, reference):
    a = normalize(adata)
    r = normalize(reference)
    props = deconvolve(a, r, cell_type_key="cell_type")
    assert (props.values >= 0).all()
