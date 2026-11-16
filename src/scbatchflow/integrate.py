"""Batch-integration methods behind one interface.

Every method takes the merged, normalised AnnData and returns a cells x dims
embedding. To add a method (for example a trained HarmonizeNN encoder), write
a function with the same signature and register it in METHODS; the pipeline,
the metrics and the database pick it up from the name alone.
"""

from __future__ import annotations

from collections.abc import Callable

import anndata as ad
import numpy as np
import scanpy as sc


def _pca(adata: ad.AnnData, batch_key: str, seed: int) -> np.ndarray:
    """No correction: the baseline every method has to beat."""
    return np.asarray(adata.obsm["X_pca"])


def _harmony(adata: ad.AnnData, batch_key: str, seed: int) -> np.ndarray:
    import harmonypy

    res = harmonypy.run_harmony(
        np.asarray(adata.obsm["X_pca"], dtype=np.float64), adata.obs, [batch_key],
        random_state=seed, verbose=False,
    )
    z = np.asarray(res.Z_corr)
    return z if z.shape[0] == adata.n_obs else z.T


def _combat(adata: ad.AnnData, batch_key: str, seed: int) -> np.ndarray:
    hvg = adata[:, adata.var["highly_variable"]].copy()
    if hasattr(hvg.X, "toarray"):
        hvg.X = hvg.X.toarray()
    # drop genes that are constant inside any batch: ComBat cannot scale them
    keep = np.ones(hvg.n_vars, dtype=bool)
    for b in hvg.obs[batch_key].unique():
        keep &= np.asarray(hvg[hvg.obs[batch_key] == b].X.var(axis=0)).ravel() > 0
    hvg = hvg[:, keep].copy()
    sc.pp.combat(hvg, key=batch_key)
    hvg.X = np.array(hvg.X, dtype=np.float32, copy=True)  # combat can hand back a read-only view
    sc.pp.scale(hvg, max_value=10)
    n = int(min(adata.obsm["X_pca"].shape[1], hvg.n_vars - 1, hvg.n_obs - 1))
    sc.pp.pca(hvg, n_comps=n, random_state=seed)
    return np.asarray(hvg.obsm["X_pca"])


METHODS: dict[str, Callable[[ad.AnnData, str, int], np.ndarray]] = {
    "pca": _pca,
    "harmony": _harmony,
    "combat": _combat,
}


def run_method(adata: ad.AnnData, method: str, batch_key: str = "sample",
               seed: int = 0) -> np.ndarray:
    if method not in METHODS:
        raise ValueError(f"unknown method '{method}'; available: {sorted(METHODS)}")
    if adata.obs[batch_key].nunique() < 2 and method != "pca":
        raise ValueError(f"'{method}' needs at least two levels of '{batch_key}'")
    emb = np.asarray(METHODS[method](adata, batch_key, seed), dtype=np.float32)
    if emb.shape[0] != adata.n_obs or not np.isfinite(emb).all():
        raise RuntimeError(f"{method} returned an invalid embedding with shape {emb.shape}")
    return emb
