"""Final clustering, UMAP, cluster labels and marker genes on the chosen embedding."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def finalise(adata: ad.AnnData, embedding: np.ndarray, method: str, resolution: float = 0.8,
             n_markers: int = 25, seed: int = 0) -> tuple[ad.AnnData, pd.DataFrame]:
    adata = adata.copy()
    adata.obsm["X_integrated"] = embedding
    adata.uns["integration_method"] = method
    sc.pp.neighbors(adata, use_rep="X_integrated", random_state=seed,
                    n_neighbors=min(15, adata.n_obs - 1))
    sc.tl.leiden(adata, resolution=resolution, random_state=seed, flavor="igraph",
                 n_iterations=2, directed=False, key_added="cluster")
    sc.tl.umap(adata, random_state=seed)

    # name each cluster after the majority marker label of its cells
    majority = (
        adata.obs.groupby("cluster", observed=True)["marker_label"]
        .agg(lambda s: s.value_counts().idxmax())
    )
    adata.obs["cell_type"] = adata.obs["cluster"].map(majority).astype(str)

    markers = pd.DataFrame(columns=["cluster", "gene", "rank", "logfc", "pval_adj", "score"])
    if adata.obs["cluster"].nunique() > 1:
        sc.tl.rank_genes_groups(adata, "cluster", method="wilcoxon", n_genes=n_markers)
        frames = []
        for cl in adata.obs["cluster"].cat.categories:
            df = sc.get.rank_genes_groups_df(adata, group=cl).head(n_markers)
            frames.append(pd.DataFrame({
                "cluster": cl, "gene": df["names"], "rank": np.arange(1, len(df) + 1),
                "logfc": df["logfoldchanges"], "pval_adj": df["pvals_adj"], "score": df["scores"],
            }))
        markers = pd.concat(frames, ignore_index=True)
        del adata.uns["rank_genes_groups"]  # recarrays do not round-trip cleanly to h5ad
    return adata, markers
