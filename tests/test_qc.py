from spatialpath.qc import compute_qc_metrics, filter_spots, qc_summary


def test_metrics_added(adata):
    out = compute_qc_metrics(adata)
    for col in ("total_counts", "n_genes_by_counts", "pct_counts_mt"):
        assert col in out.obs


def test_filter_reduces_or_keeps(adata):
    out = compute_qc_metrics(adata)
    before = out.n_obs
    out = filter_spots(out, min_genes=1, min_counts=1, max_mito=100.0)
    assert out.n_obs <= before


def test_summary_keys(adata):
    out = compute_qc_metrics(adata)
    summary = qc_summary(out)
    assert set(summary) == {
        "n_spots",
        "n_genes",
        "median_genes_per_spot",
        "median_counts_per_spot",
        "median_mito_pct",
    }
