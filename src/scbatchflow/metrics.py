"""Integration metrics: how well batches mix, and how much biology survives.

Batch removal:        iLISI, batch silhouette, kBET acceptance
Biology conservation: cLISI, label silhouette, cluster/label ARI and NMI

All scores are rescaled to 0..1 where 1 is best, following the conventions of
the scIB benchmark (Luecken et al. 2022), so they can be averaged.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import chi2
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_samples
from sklearn.neighbors import NearestNeighbors

BATCH_METRICS = ["ilisi", "batch_asw", "kbet_acceptance"]
BIO_METRICS = ["clisi", "label_asw", "ari", "nmi"]


def lisi_scores(emb: np.ndarray, labels: pd.Series, perplexity: float = 30) -> np.ndarray:
    """Per-cell LISI (Korsunsky et al. 2019): effective number of labels nearby."""
    from harmonypy import compute_lisi

    perplexity = float(min(perplexity, max(5, (emb.shape[0] - 1) // 3)))
    meta = pd.DataFrame({"x": labels.astype(str).to_numpy()})
    return np.asarray(compute_lisi(emb.astype(np.float64), meta, ["x"], perplexity)).ravel()


def kbet_acceptance(emb: np.ndarray, batch: pd.Series, k: int | None = None,
                    alpha: float = 0.05) -> float:
    """Share of neighbourhoods whose batch mix matches the global mix (kBET).

    For each cell, a chi-squared test compares batch counts among its k nearest
    neighbours with the dataset-wide proportions (Buettner et al. 2019).
    """
    codes, uniques = pd.factorize(batch.astype(str))
    n, n_batches = emb.shape[0], len(uniques)
    if n_batches < 2:
        return float("nan")
    k = k or int(np.clip(n // 10, 10, 100))
    k = min(k, n - 1)
    expected = np.bincount(codes, minlength=n_batches) / n * k
    nn = NearestNeighbors(n_neighbors=k + 1).fit(emb)
    idx = nn.kneighbors(emb, return_distance=False)[:, 1:]
    observed = np.stack([(codes[idx] == b).sum(axis=1) for b in range(n_batches)], axis=1)
    stat = ((observed - expected) ** 2 / expected).sum(axis=1)
    pvals = chi2.sf(stat, df=n_batches - 1)
    return float((pvals >= alpha).mean())


def batch_asw(emb: np.ndarray, batch: pd.Series, labels: pd.Series) -> float:
    """1 - |silhouette on batch|, averaged within each cell type then across types."""
    scores = []
    for ct in labels.unique():
        mask = (labels == ct).to_numpy()
        b = batch[mask]
        if mask.sum() < 10 or b.nunique() < 2 or b.value_counts().min() < 2:
            continue
        scores.append(float(np.mean(1 - np.abs(silhouette_samples(emb[mask], b.to_numpy())))))
    return float(np.mean(scores)) if scores else float("nan")


def cluster(emb: np.ndarray, resolution: float = 1.0, n_neighbors: int = 15,
            seed: int = 0) -> np.ndarray:
    tmp = ad.AnnData(X=np.zeros((emb.shape[0], 1), dtype=np.float32))
    tmp.obsm["X_emb"] = emb
    sc.pp.neighbors(tmp, use_rep="X_emb", n_neighbors=min(n_neighbors, emb.shape[0] - 1),
                    random_state=seed)
    sc.tl.leiden(tmp, resolution=resolution, random_state=seed, flavor="igraph",
                 n_iterations=2, directed=False)
    return tmp.obs["leiden"].to_numpy()


def best_clustering(emb: np.ndarray, labels: np.ndarray, seed: int = 0) -> tuple[float, float]:
    """ARI and NMI at the Leiden resolution that best recovers the labels."""
    best = (-1.0, -1.0)
    for res in (0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0):
        cl = cluster(emb, resolution=res, seed=seed)
        nmi = float(normalized_mutual_info_score(labels, cl))
        if nmi > best[1]:
            best = (float(adjusted_rand_score(labels, cl)), nmi)
    return best


def evaluate(emb: np.ndarray, batch: pd.Series, labels: pd.Series, seed: int = 0) -> dict:
    """Score one embedding. Cells labelled 'Unassigned' are left out of bio metrics."""
    batch = batch.astype(str).reset_index(drop=True)
    labels = labels.astype(str).reset_index(drop=True)
    n_batches = batch.nunique()
    out: dict[str, float] = {}

    ilisi = lisi_scores(emb, batch)
    out["ilisi_raw"] = float(np.median(ilisi))
    out["ilisi"] = float((np.median(ilisi) - 1) / (n_batches - 1)) if n_batches > 1 else np.nan
    out["kbet_acceptance"] = kbet_acceptance(emb, batch)

    known = (labels != "Unassigned").to_numpy()
    n_labels = labels[known].nunique()
    if known.sum() >= 20 and n_labels >= 2:
        e, lab, bat = emb[known], labels[known].reset_index(drop=True), batch[known]
        clisi = lisi_scores(e, lab)
        out["clisi_raw"] = float(np.median(clisi))
        out["clisi"] = float((n_labels - np.median(clisi)) / (n_labels - 1))
        out["label_asw"] = float((silhouette_samples(e, lab.to_numpy()).mean() + 1) / 2)
        out["batch_asw"] = batch_asw(e, bat.reset_index(drop=True), lab)
        out["ari"], out["nmi"] = best_clustering(e, lab.to_numpy(), seed=seed)
    else:
        for m in ("clisi_raw", "clisi", "label_asw", "batch_asw", "ari", "nmi"):
            out[m] = float("nan")

    out["batch_score"] = _mean_available([out[m] for m in BATCH_METRICS])
    out["bio_score"] = _mean_available([out[m] for m in BIO_METRICS])
    # scIB weighting: biology counts for more than mixing. Without usable labels the
    # overall score falls back to batch mixing alone, which the summary makes visible
    # through the NaN bio_score.
    if np.isnan(out["bio_score"]):
        out["overall_score"] = out["batch_score"]
    else:
        out["overall_score"] = float(0.4 * out["batch_score"] + 0.6 * out["bio_score"])
    return out


def _mean_available(values: list[float]) -> float:
    vals = [v for v in values if not np.isnan(v)]
    return float(np.mean(vals)) if vals else float("nan")
