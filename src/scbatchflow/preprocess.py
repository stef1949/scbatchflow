"""Merge samples, normalise, pick variable genes, PCA and marker-based labels."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def merge_and_normalise(
    adatas: list[ad.AnnData], n_hvg: int = 2000, n_pcs: int = 30, seed: int = 0
) -> ad.AnnData:
    adata = ad.concat(adatas, join="inner", merge="same")
    adata.obs_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=3)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    n_hvg = min(n_hvg, adata.n_vars)
    batch_key = "sample" if adata.obs["sample"].nunique() > 1 else None
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat", batch_key=batch_key)
    if adata.var["highly_variable"].sum() < 10:  # tiny references: keep everything
        adata.var["highly_variable"] = True

    n_pcs = int(min(n_pcs, adata.var["highly_variable"].sum() - 1, adata.n_obs - 1))
    hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(hvg, max_value=10)
    sc.pp.pca(hvg, n_comps=n_pcs, random_state=seed)
    adata.obsm["X_pca"] = hvg.obsm["X_pca"]
    adata.uns["n_pcs"] = n_pcs
    return adata


def load_markers(path: str | Path) -> dict[str, list[str]]:
    markers = json.loads(Path(path).read_text())
    if not isinstance(markers, dict) or not markers:
        raise ValueError(f"{path} must be a JSON object of cell type -> list of gene symbols")
    return {k: list(v) for k, v in markers.items()}


def label_by_markers(
    adata: ad.AnnData, markers: dict[str, list[str]], min_margin: float = 0.1, seed: int = 0
) -> pd.Series:
    """Label each cell by its best marker-set score, computed within its own sample.

    Scoring inside each sample keeps the labels independent of any batch
    correction, so they can serve as the biological reference when the
    integration methods are compared. Cells whose best score does not beat the
    runner-up by `min_margin` are labelled 'Unassigned'.
    """
    usable = {ct: [g for g in genes if g in adata.var_names] for ct, genes in markers.items()}
    usable = {ct: genes for ct, genes in usable.items() if len(genes) >= 2}
    labels = pd.Series("Unassigned", index=adata.obs_names, dtype=object)
    if len(usable) < 2:
        return labels
    for sample in adata.obs["sample"].unique():
        sub = adata[adata.obs["sample"] == sample].copy()
        scores = {}
        for ct, genes in usable.items():
            ctrl = max(5, min(50, sub.n_vars // 4))
            try:
                sc.tl.score_genes(sub, genes, ctrl_size=ctrl, score_name="_s", random_state=seed)
                scores[ct] = sub.obs["_s"].to_numpy()
            except Exception:
                scores[ct] = np.asarray(sub[:, genes].X.mean(axis=1)).ravel()
        mat = pd.DataFrame(scores, index=sub.obs_names)
        ordered = np.sort(mat.to_numpy(), axis=1)
        margin = ordered[:, -1] - ordered[:, -2]
        best = mat.idxmax(axis=1)
        labels.loc[sub.obs_names] = best.where((margin >= min_margin) & (ordered[:, -1] > 0),
                                               "Unassigned")
    return labels
