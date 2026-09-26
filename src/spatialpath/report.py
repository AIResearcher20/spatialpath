import json
from pathlib import Path


def build_report(qc_stats, adata, cluster_key="cluster"):
    cluster_sizes = (
        adata.obs[cluster_key].value_counts().sort_index()
        if cluster_key in adata.obs
        else None
    )

    report = {
        "data": {
            "n_spots": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
        },
        "qc": qc_stats,
    }

    if cluster_sizes is not None:
        report["clustering"] = {
            "n_clusters": int(cluster_sizes.shape[0]),
            "sizes": {str(k): int(v) for k, v in cluster_sizes.items()},
        }

    return report


def save_report(report, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
