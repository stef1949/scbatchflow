"""Cell calling and per-cell quality control for one sample."""

from __future__ import annotations

import anndata as ad
import numpy as np
import scanpy as sc


def call_cells(total_umis: np.ndarray, expected_cells: int = 3000) -> float:
    """Order-of-magnitude cell calling (the Cell Ranger 2 rule).

    Take the 99th percentile of UMI totals among the top `expected_cells`
    barcodes and call every barcode with at least a tenth of it a cell.
    """
    totals = np.sort(np.asarray(total_umis, dtype=float))[::-1]
    top = totals[: max(1, min(expected_cells, totals.size))]
    return float(np.percentile(top, 99) / 10.0)


def mad_upper_bound(values: np.ndarray, n_mads: float) -> float:
    med = float(np.median(values))
    mad = float(np.median(np.abs(values - med)))
    return med + n_mads * 1.4826 * mad


def cell_qc(
    adata: ad.AnnData,
    expected_cells: int = 3000,
    min_genes: int = 200,
    max_pct_mito: float = 20.0,
    mito_mads: float | None = 5.0,
    doublets: bool = True,
    seed: int = 0,
) -> tuple[ad.AnnData, dict]:
    """Return the filtered AnnData and a summary of what was removed and why."""
    adata = adata.copy()
    adata.var["mito"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mito"], percent_top=None, log1p=False, inplace=True
    )
    adata.obs = adata.obs.rename(
        columns={"total_counts": "n_umis", "n_genes_by_counts": "n_genes",
                 "pct_counts_mito": "pct_mito"}
    )
    summary: dict = {"n_barcodes": int(adata.n_obs)}

    threshold = call_cells(adata.obs["n_umis"].to_numpy(), expected_cells)
    is_cell = adata.obs["n_umis"].to_numpy() >= threshold
    summary["umi_threshold"] = threshold
    summary["n_cells_called"] = int(is_cell.sum())
    adata = adata[is_cell].copy()

    mito_cut = max_pct_mito
    if mito_mads is not None and adata.n_obs > 10:
        mito_cut = min(max_pct_mito, mad_upper_bound(adata.obs["pct_mito"].to_numpy(), mito_mads))
        # never let a tight distribution push the cut below a sane floor
        mito_cut = max(mito_cut, 5.0)
    low_genes = adata.obs["n_genes"].to_numpy() < min_genes
    high_mito = adata.obs["pct_mito"].to_numpy() > mito_cut
    summary["pct_mito_cut"] = float(mito_cut)
    summary["n_removed_low_genes"] = int(low_genes.sum())
    summary["n_removed_high_mito"] = int((high_mito & ~low_genes).sum())
    adata = adata[~(low_genes | high_mito)].copy()

    adata.obs["doublet_score"] = np.nan
    adata.obs["predicted_doublet"] = False
    summary["n_removed_doublets"] = 0
    if doublets and adata.n_obs >= 100:
        try:
            sc.pp.scrublet(adata, random_state=seed, verbose=False)
            adata.obs["predicted_doublet"] = adata.obs["predicted_doublet"].fillna(False)
            summary["n_removed_doublets"] = int(adata.obs["predicted_doublet"].sum())
            adata = adata[~adata.obs["predicted_doublet"].to_numpy()].copy()
        except Exception as err:  # scrublet can fail on very small or degenerate inputs
            summary["doublet_warning"] = f"scrublet skipped: {err}"

    summary["n_cells_pass_qc"] = int(adata.n_obs)
    if adata.n_obs:
        summary["median_umis"] = float(np.median(adata.obs["n_umis"]))
        summary["median_genes"] = float(np.median(adata.obs["n_genes"]))
        summary["median_pct_mito"] = float(np.median(adata.obs["pct_mito"]))
    return adata, summary
