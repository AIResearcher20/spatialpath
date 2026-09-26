import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge


def _dense(matrix):
    return matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)


def build_signatures(reference, cell_type_key="cell_type"):
    X = _dense(reference.X)
    types = reference.obs[cell_type_key].astype("category")
    genes = pd.Index(reference.var_names)

    signatures = {}
    for name in types.cat.categories:
        mask = (types == name).values
        if mask.sum() == 0:
            continue
        signatures[name] = X[mask].mean(axis=0)

    return pd.DataFrame(signatures, index=genes)


def deconvolve(adata, reference, cell_type_key="cell_type", alpha=1.0):
    signatures = build_signatures(reference, cell_type_key=cell_type_key)
    common = adata.var_names.intersection(signatures.index)
    if len(common) == 0:
        raise ValueError("no shared genes between spots and reference")

    X_spot = _dense(adata[:, common].X)
    S = signatures.loc[common].values

    model = Ridge(alpha=alpha, positive=True, fit_intercept=False)
    model.fit(S, X_spot.T)
    proportions = model.coef_.T
    proportions = np.clip(proportions, 0, None)

    row_sums = proportions.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    proportions = proportions / row_sums

    return pd.DataFrame(
        proportions,
        index=adata.obs_names,
        columns=signatures.columns,
    )
