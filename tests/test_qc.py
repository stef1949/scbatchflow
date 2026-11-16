import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from scbatchflow.qc import call_cells, cell_qc, mad_upper_bound


def make_adata(n_cells=120, n_empty=300, n_genes=60, seed=0):
    rng = np.random.default_rng(seed)
    cells = rng.poisson(8, size=(n_cells, n_genes))
    empty = rng.poisson(0.05, size=(n_empty, n_genes))
    genes = [f"G{i}" for i in range(n_genes - 5)] + [f"MT-{i}" for i in range(5)]
    x = np.vstack([cells, empty]).astype("float32")
    obs = pd.DataFrame({"sample": "s1", "batch": "b1"},
                       index=[f"s1_bc{i}" for i in range(n_cells + n_empty)])
    return ad.AnnData(sp.csr_matrix(x), obs=obs, var=pd.DataFrame(index=genes))


def test_call_cells_separates_cells_from_empty_droplets():
    totals = np.concatenate([np.full(100, 1000.0), np.full(900, 5.0)])
    threshold = call_cells(totals, expected_cells=100)
    assert 5 < threshold < 1000
    assert (totals >= threshold).sum() == 100


def test_mad_upper_bound_ignores_outliers():
    values = np.array([1.0] * 50 + [1.2] * 50 + [80.0])
    assert mad_upper_bound(values, 5) < 80


def test_cell_qc_removes_empties_and_dying_cells():
    adata = make_adata()
    # make ten real cells mostly mitochondrial
    x = adata.X.toarray()
    x[:10, -5:] += 400
    adata.X = sp.csr_matrix(x)
    filtered, summary = cell_qc(adata, expected_cells=120, min_genes=20, doublets=False)
    assert summary["n_barcodes"] == 420
    assert summary["n_cells_called"] == 120
    assert summary["n_removed_high_mito"] == 10
    assert filtered.n_obs == summary["n_cells_pass_qc"] == 110
    assert {"n_umis", "n_genes", "pct_mito", "doublet_score"} <= set(filtered.obs.columns)
