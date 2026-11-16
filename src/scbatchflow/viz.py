"""Figures for the run report.

Conventions used throughout:
- batch/sample and cell type are identities, so they take fixed categorical colours
  assigned in sorted order (the same entity keeps its colour in every figure);
- magnitudes (scores, expression) use one blue ramp, light to dark;
- 'Unassigned' is always neutral grey.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6",
               "#e34948"]
INK, SECONDARY, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURFACE, NEUTRAL = "#e1e0d9", "#c3c2b7", "#fcfcfb", "#b9b8b2"
BAD = "#d03b3b"
SEQUENTIAL = LinearSegmentedColormap.from_list(
    "scbf_blue", ["#eef4fd", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": AXIS, "axes.labelcolor": SECONDARY, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "font.size": 9, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False,
})


def colour_map(values) -> dict[str, str]:
    """Fixed colour per category, by sorted name. More than 8 categories reuse hues, so
    figures that can exceed 8 (clusters) also label the groups on the plot itself."""
    cats = sorted({str(v) for v in values} - {"Unassigned"}, key=_natural)
    out = {c: CATEGORICAL[i % len(CATEGORICAL)] for i, c in enumerate(cats)}
    out["Unassigned"] = NEUTRAL
    return out


def _natural(s: str):
    return (0, int(s)) if s.isdigit() else (1, s)


def _save(fig, path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def umap_coords(embedding: np.ndarray, seed: int = 0) -> np.ndarray:
    tmp = ad.AnnData(X=np.zeros((embedding.shape[0], 1), dtype=np.float32))
    tmp.obsm["X_emb"] = np.asarray(embedding)
    sc.pp.neighbors(tmp, use_rep="X_emb", n_neighbors=min(15, tmp.n_obs - 1), random_state=seed)
    sc.tl.umap(tmp, random_state=seed)
    return np.asarray(tmp.obsm["X_umap"])


def _scatter(ax, xy, labels, title, legend=True, centroids=False, seed=0):
    labels = np.asarray([str(v) for v in labels])
    colours = colour_map(labels)
    order = np.random.default_rng(seed).permutation(len(labels))  # no group drawn on top
    size = float(np.clip(6000 / max(len(labels), 1), 1.5, 14))
    ax.scatter(xy[order, 0], xy[order, 1], c=[colours[v] for v in labels[order]], s=size,
               linewidths=0, alpha=0.85, rasterized=True)
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)
    cats = sorted(set(labels) - {"Unassigned"}, key=_natural)
    cats += ["Unassigned"] if "Unassigned" in labels else []
    if centroids:
        for c in cats:
            pts = xy[labels == c]
            ax.text(np.median(pts[:, 0]), np.median(pts[:, 1]), c, ha="center", va="center",
                    fontsize=8, fontweight="bold", color=INK,
                    bbox={"boxstyle": "round,pad=0.15", "fc": SURFACE, "ec": "none", "alpha": 0.8})
    if legend and len(cats) <= 12:
        handles = [plt.Line2D([], [], marker="o", ls="", ms=6, color=colours[c], label=c)
                   for c in cats]
        ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8,
                  handletextpad=0.2, labelcolor=SECONDARY)


def plot_final_umaps(adata: ad.AnnData, outdir, batch_key: str = "sample") -> list[Path]:
    """One file per colouring, plus a three-panel overview for the report."""
    out = Path(outdir)
    xy = np.asarray(adata.obsm["X_umap"])
    method = adata.uns.get("integration_method", "")
    panels = [(batch_key, f"By {batch_key}", False), ("cell_type", "By cell type", False),
              ("cluster", "By cluster", True)]
    paths = []
    for key, title, cen in panels:
        fig, ax = plt.subplots(figsize=(5.6, 4.4))
        _scatter(ax, xy, adata.obs[key], f"{title} ({method})", centroids=cen, legend=not cen)
        paths.append(_save(fig, out / f"umap_{'sample' if key == batch_key else key}.png"))
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    for ax, (key, title, cen) in zip(axes, panels, strict=True):
        _scatter(ax, xy, adata.obs[key], title, centroids=cen, legend=not cen)
    fig.subplots_adjust(wspace=0.45)
    paths.append(_save(fig, out / "umap_overview.png"))
    return paths


def plot_method_grid(adata: ad.AnnData, embeddings: dict[str, np.ndarray], path,
                     batch_key: str = "sample", max_cells: int = 20000, seed: int = 0) -> Path:
    """Before and after: one row per method, coloured by batch and by marker label."""
    idx = np.arange(adata.n_obs)
    if adata.n_obs > max_cells:
        idx = np.sort(np.random.default_rng(seed).choice(adata.n_obs, max_cells, replace=False))
    methods = sorted(embeddings, key=lambda m: (m != "pca", m))  # uncorrected baseline first
    fig, axes = plt.subplots(len(methods), 2, figsize=(10.5, 3.9 * len(methods)), squeeze=False)
    for row, method in enumerate(methods):
        xy = umap_coords(embeddings[method][idx], seed=seed)
        name = "pca (no correction)" if method == "pca" else method
        _scatter(axes[row, 0], xy, adata.obs[batch_key].to_numpy()[idx],
                 f"{name}: by {batch_key}", seed=seed)
        _scatter(axes[row, 1], xy, adata.obs["marker_label"].to_numpy()[idx],
                 f"{name}: by marker label", seed=seed)
    fig.subplots_adjust(wspace=0.5, hspace=0.25)
    return _save(fig, path)


METRIC_LABELS = {
    "ilisi": "iLISI", "kbet_acceptance": "kBET", "batch_asw": "Batch ASW", "clisi": "cLISI",
    "label_asw": "Label ASW", "ari": "ARI", "nmi": "NMI", "batch_score": "Batch removal",
    "bio_score": "Biology kept", "overall_score": "Overall",
}


def plot_metric_heatmap(table: pd.DataFrame, path) -> Path:
    cols = [c for c in METRIC_LABELS if c in table.columns]
    data = table[cols].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(1.0 + 0.95 * len(cols), 0.9 + 0.55 * len(table)))
    ax.imshow(np.nan_to_num(data, nan=0.0), cmap=SEQUENTIAL, vmin=0, vmax=1, aspect="auto")
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            v = data[i, j]
            ax.text(j, i, "n/a" if np.isnan(v) else f"{v:.2f}", ha="center", va="center",
                    fontsize=9, color="#ffffff" if (not np.isnan(v) and v > 0.55) else INK)
    ax.set_xticks(range(len(cols)), [METRIC_LABELS[c] for c in cols], rotation=30, ha="right")
    ax.set_yticks(range(len(table)), [str(m) for m in table.index])
    ax.tick_params(length=0, colors=SECONDARY)
    for s in ax.spines.values():
        s.set_visible(False)
    # white rules separate the raw metrics from the three summary scores
    ax.set_xticks(np.arange(-0.5, len(cols)), minor=True)
    ax.set_yticks(np.arange(-0.5, len(table)), minor=True)
    ax.grid(which="minor", color=SURFACE, linewidth=2)
    ax.tick_params(which="minor", length=0)
    ax.set_title("Integration scores (0 to 1, higher is better)")
    return _save(fig, path)


def plot_tradeoff(table: pd.DataFrame, path) -> Path:
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    ax.grid(color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    best = table["overall_score"].idxmax()
    for method, row in table.iterrows():
        top = method == best
        ax.scatter(row["batch_score"], row["bio_score"], s=90, color=CATEGORICAL[0] if top
                   else NEUTRAL, edgecolors=SURFACE, linewidths=2, zorder=3)
        ax.annotate(f"{method}" + ("  (chosen)" if top else ""),
                    (row["batch_score"], row["bio_score"]), textcoords="offset points",
                    xytext=(8, 7), color=INK if top else SECONDARY,
                    fontweight="bold" if top else "normal")
    ax.set_xlabel("Batch removal")
    ax.set_ylabel("Biology kept")
    ax.set_xlim(0, 1.08)
    ax.set_ylim(0, 1.08)
    ax.set_title("Trade-off between mixing batches and keeping biology")
    return _save(fig, path)


def plot_composition(adata: ad.AnnData, path, batch_key: str = "sample") -> Path:
    tab = pd.crosstab(adata.obs[batch_key].astype(str), adata.obs["cell_type"].astype(str),
                      normalize="index") * 100
    colours = colour_map(tab.columns)
    cols = sorted([c for c in tab.columns if c != "Unassigned"]) + \
        (["Unassigned"] if "Unassigned" in tab.columns else [])
    fig, ax = plt.subplots(figsize=(7.5, 0.9 + 0.5 * len(tab)))
    left = np.zeros(len(tab))
    for c in cols:
        vals = tab[c].to_numpy()
        ax.barh(tab.index, vals, left=left, color=colours[c], edgecolor=SURFACE, linewidth=2,
                height=0.62, label=c)
        for y, (v, lft) in enumerate(zip(vals, left, strict=True)):
            if v >= 8:
                ax.text(lft + v / 2, y, f"{v:.0f}%", ha="center", va="center", fontsize=8,
                        color=INK if colours[c] in ("#eda100", NEUTRAL) else "#ffffff")
        left += vals
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of cells in the sample")
    ax.invert_yaxis()
    ax.tick_params(axis="y", length=0, colors=SECONDARY)
    ax.spines["left"].set_visible(False)
    ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8, labelcolor=SECONDARY)
    ax.set_title(f"Cell-type composition per {batch_key}")
    return _save(fig, path)


def plot_qc_violins(adata: ad.AnnData, path) -> Path:
    samples = sorted(adata.obs["sample"].astype(str).unique())
    colours = colour_map(samples)
    metrics = [("n_umis", "UMIs per cell"), ("n_genes", "Genes per cell"),
               ("pct_mito", "% mitochondrial")]
    fig, axes = plt.subplots(1, 3, figsize=(4 + 1.3 * len(samples) * 3 * 0.6, 3.6))
    for ax, (key, title) in zip(axes, metrics, strict=True):
        data = [adata.obs.loc[adata.obs["sample"].astype(str) == s, key].to_numpy()
                for s in samples]
        parts = ax.violinplot(data, showextrema=False, widths=0.8)
        for body, s in zip(parts["bodies"], samples, strict=True):
            body.set_facecolor(colours[s])
            body.set_edgecolor("none")
            body.set_alpha(0.85)
        for i, d in enumerate(data, start=1):
            q1, med, q3 = np.percentile(d, [25, 50, 75])
            ax.vlines(i, q1, q3, color=INK, linewidth=2)
            ax.scatter([i], [med], s=22, color=SURFACE, edgecolors=INK, linewidths=1.2, zorder=3)
        ax.set_xticks(range(1, len(samples) + 1), samples, rotation=20, ha="right")
        ax.set_title(title)
        ax.grid(axis="y", color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        ax.set_ylim(bottom=0)
    fig.subplots_adjust(wspace=0.35)
    return _save(fig, path)


def plot_marker_dots(adata: ad.AnnData, markers: pd.DataFrame, path, per_cluster: int = 3) -> Path:
    """Dot plot: dot area = share of the cluster expressing the gene, colour = mean expression."""
    genes: list[str] = []
    for _, grp in markers.sort_values("rank").groupby("cluster", sort=False):
        for g in grp["gene"].head(per_cluster):
            if g in adata.var_names and g not in genes:
                genes.append(g)
    clusters = sorted(adata.obs["cluster"].astype(str).unique(), key=_natural)
    names = (adata.obs.assign(cluster=adata.obs["cluster"].astype(str))
             .groupby("cluster")["cell_type"].agg(lambda s: s.astype(str).mode().iat[0]))
    sub = adata[:, genes].X
    sub = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
    cl = adata.obs["cluster"].astype(str).to_numpy()
    frac = np.array([(sub[cl == c] > 0).mean(axis=0) for c in clusters])
    mean = np.array([sub[cl == c].mean(axis=0) for c in clusters])
    scaled = mean / np.maximum(mean.max(axis=0, keepdims=True), 1e-9)  # per gene, 0 to 1

    fig, ax = plt.subplots(figsize=(2.6 + 0.36 * len(genes), 1.4 + 0.36 * len(clusters)))
    xx, yy = np.meshgrid(range(len(genes)), range(len(clusters)))
    dots = ax.scatter(xx, yy, s=frac * 150 + 2, c=scaled, cmap=SEQUENTIAL, vmin=0, vmax=1,
                      edgecolors=AXIS, linewidths=0.5)
    ax.set_xticks(range(len(genes)), genes, rotation=90)
    ax.set_yticks(range(len(clusters)), [f"{c}  {names.get(c, '')}" for c in clusters])
    ax.set_xlim(-0.7, len(genes) - 0.3)
    ax.set_ylim(len(clusters) - 0.4, -0.6)
    ax.tick_params(length=0, colors=SECONDARY)
    for s in ax.spines.values():
        s.set_visible(False)
    cbar = fig.colorbar(dots, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Mean expression (scaled per gene)", color=SECONDARY)
    cbar.outline.set_visible(False)
    handles = [plt.scatter([], [], s=f * 150 + 2, color=NEUTRAL, edgecolors=AXIS, linewidths=0.5,
                           label=f"{int(f * 100)}%") for f in (0.25, 0.5, 1.0)]
    ax.legend(handles=handles, title="Cells expressing", loc="upper left",
              bbox_to_anchor=(1.18, 1.0), fontsize=8, title_fontsize=8, labelspacing=1.0)
    ax.set_title("Top marker genes per cluster")
    return _save(fig, path)


def plot_sample_qc(raw: ad.AnnData, kept: set[str], summary: dict, path) -> Path:
    """Left: barcode rank plot with the cell-calling threshold. Right: called cells on
    genes against mitochondrial share, showing which ones QC removed."""
    x = raw.X
    totals = np.asarray(x.sum(axis=1)).ravel()
    n_genes = np.asarray((x > 0).sum(axis=1)).ravel()
    mito = raw.var_names.str.upper().str.startswith("MT-")
    pct_mito = 100 * np.asarray(x[:, mito].sum(axis=1)).ravel() / np.maximum(totals, 1)
    threshold = float(summary.get("umi_threshold", 0))
    sample = str(raw.obs["sample"].iloc[0]) if "sample" in raw.obs else ""

    fig, (a, b) = plt.subplots(1, 2, figsize=(10.5, 3.9))
    ranked = np.sort(totals[totals > 0])[::-1]
    ranks = np.arange(1, len(ranked) + 1)
    n_called = int((ranked >= threshold).sum())
    a.loglog(ranks[:n_called], ranked[:n_called], color=CATEGORICAL[0], linewidth=2,
             label=f"called cells ({n_called:,})")
    a.loglog(ranks[n_called:], ranked[n_called:], color=NEUTRAL, linewidth=2,
             label=f"background ({len(ranked) - n_called:,})")
    a.axhline(threshold, color=BAD, linewidth=1, linestyle=(0, (4, 3)))
    a.text(1.2, threshold * 1.25, f"threshold {threshold:,.0f} UMIs", color=BAD, fontsize=8)
    a.set_xlabel("Barcode rank")
    a.set_ylabel("UMIs")
    a.grid(color=GRID, linewidth=0.8)
    a.set_axisbelow(True)
    a.legend(loc="lower left", fontsize=8, labelcolor=SECONDARY)
    a.set_title(f"{sample}: which barcodes are cells")

    called = totals >= threshold
    is_kept = np.array([n in kept for n in raw.obs_names])
    removed = called & ~is_kept
    size = float(np.clip(5000 / max(called.sum(), 1), 2, 12))
    b.scatter(n_genes[is_kept], pct_mito[is_kept], s=size, color=CATEGORICAL[0], alpha=0.6,
              linewidths=0, rasterized=True, label=f"kept ({int(is_kept.sum()):,})")
    b.scatter(n_genes[removed], pct_mito[removed], s=size * 1.6, color=BAD, alpha=0.85,
              linewidths=0, marker="x" if removed.sum() < 2000 else "o", rasterized=True,
              label=f"removed ({int(removed.sum()):,})")
    if "pct_mito_cut" in summary:
        b.axhline(summary["pct_mito_cut"], color=MUTED, linewidth=1, linestyle=(0, (4, 3)))
    b.set_xlabel("Genes detected")
    b.set_ylabel("% mitochondrial")
    b.set_ylim(bottom=0)
    b.grid(color=GRID, linewidth=0.8)
    b.set_axisbelow(True)
    b.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8, labelcolor=SECONDARY)
    b.set_title(f"{sample}: cells removed by QC")
    fig.subplots_adjust(wspace=0.28)
    return _save(fig, path)
