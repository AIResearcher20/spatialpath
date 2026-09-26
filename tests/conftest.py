import numpy as np
import pandas as pd
import pytest
import anndata as ad


def _make_adata(n_spots=80, n_genes=120, n_types=3, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.poisson(4, size=(n_spots, n_genes)).astype(np.float32)
    obs = pd.DataFrame(index=[f"spot_{i}" for i in range(n_spots)])
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs["cell_type"] = rng.integers(0, n_types, size=n_spots).astype(str)
    return adata


@pytest.fixture
def adata():
    return _make_adata()


@pytest.fixture
def reference():
    return _make_adata(n_spots=60, n_genes=120, n_types=3, seed=1)
