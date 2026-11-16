"""The figures and report are built from a small synthetic dataset; the test checks that
every figure is produced and that the report embeds them, not how they look."""

import json

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from scbatchflow import viz
from scbatchflow.annotate import finalise
from scbatchflow.cli import main
from scbatchflow.metrics import evaluate
from scbatchflow.preprocess import label_by_markers, merge_and_normalise
from scbatchflow.qc import cell_qc

MARKERS = {"A": ["A1", "A2", "A3"], "B": ["B1", "B2", "B3"]}


def raw_sample(name, shift, seed):
    rng = np.random.default_rng(seed)
    genes = MARKERS["A"] + MARKERS["B"] + [f"G{i}" for i in range(40)] + ["MT-1", "MT-2"]
    n = 120
    lam = np.full((n, len(genes)), 2.0)
    lam[: n // 2, 0:3] = 25
    lam[n // 2:, 3:6] = 25
    lam[:, 6:26] *= shift
    cells = rng.poisson(lam)
    empties = rng.poisson(0.05, size=(200, len(genes)))
    x = sp.csr_matrix(np.vstack([cells, empties]).astype("float32"))
    obs = pd.DataFrame({"sample": name, "batch": name},
                       index=[f"{name}_{i}" for i in range(x.shape[0])])
    return ad.AnnData(x, obs=obs, var=pd.DataFrame(index=genes))


def test_figures_and_report(tmp_path):
    plots = tmp_path / "in"
    plots.mkdir()
    filtered, summaries, qc_pngs = [], [], []
    for i, (name, shift) in enumerate([("s1", 1.0), ("s2", 2.5)]):
        raw = raw_sample(name, shift, i)
        kept, summary = cell_qc(raw, expected_cells=120, min_genes=10, doublets=False)
        png = viz.plot_sample_qc(raw, set(kept.obs_names), summary, plots / f"{name}.qc.png")
        assert png.stat().st_size > 5000
        filtered.append(kept)
        qc_pngs.append(str(png))
        sfile = plots / f"{name}.summary.json"
        sfile.write_text(json.dumps({"sample": name, "batch": name, "n_reads": 1000,
                                     "pct_pseudoaligned": 95.0, **summary}))
        summaries.append(str(sfile))

    merged = merge_and_normalise(filtered, n_hvg=40, n_pcs=8)
    merged.obs["marker_label"] = pd.Categorical(label_by_markers(merged, MARKERS))
    merged.write_h5ad(tmp_path / "merged.h5ad")
    emb_files, metric_files = [], []
    for method, emb in {"pca": merged.obsm["X_pca"], "other": merged.obsm["X_pca"] * 0.5}.items():
        np.save(tmp_path / f"{method}.embedding.npy", emb)
        scores = evaluate(np.asarray(emb), merged.obs["sample"], merged.obs["marker_label"])
        (tmp_path / f"{method}.metrics.json").write_text(json.dumps({"method": method, **scores}))
        emb_files.append(str(tmp_path / f"{method}.embedding.npy"))
        metric_files.append(str(tmp_path / f"{method}.metrics.json"))

    out = tmp_path / "results"
    main(["finalise", "--input", str(tmp_path / "merged.h5ad"), "--embeddings", *emb_files,
          "--metrics", *metric_files, "--sample-summaries", *summaries, "--qc-plots", *qc_pngs,
          "--run-label", "unit test", "--outdir", str(out)])

    expected = ["umap_sample.png", "umap_cell_type.png", "umap_cluster.png", "umap_overview.png",
                "method_grid.png", "metric_heatmap.png", "integration_tradeoff.png",
                "composition.png", "qc_violins.png", "marker_dots.png", "qc_s1.png", "qc_s2.png"]
    for name in expected:
        assert (out / "plots" / name).stat().st_size > 3000, name
    html = (out / "report.html").read_text()
    assert html.count("data:image/png;base64,") == len(expected) - 3  # overview replaces 3 UMAPs
    for heading in ("QC funnel", "Integration benchmark", "Before and after", "final embedding",
                    "defines each cluster"):
        assert heading in html
    assert "unit test" in html and "<script" not in html


def test_colours_follow_the_entity_not_its_position():
    assert viz.colour_map(["b", "a"])["a"] == viz.colour_map(["a", "c", "b"])["a"]
    assert viz.colour_map(["x", "Unassigned"])["Unassigned"] == viz.NEUTRAL


def test_finalise_names_clusters_from_marker_labels():
    rng = np.random.default_rng(0)
    emb = np.vstack([rng.normal(0, 1, (60, 4)), rng.normal(12, 1, (60, 4))]).astype("float32")
    x = sp.csr_matrix(rng.poisson(3, (120, 30)).astype("float32"))
    adata = ad.AnnData(x, obs=pd.DataFrame({"sample": ["s1", "s2"] * 60,
                                            "marker_label": ["A"] * 60 + ["B"] * 60},
                                           index=[f"c{i}" for i in range(120)]),
                       var=pd.DataFrame(index=[f"g{i}" for i in range(30)]))
    out, _ = finalise(adata, emb, "pca", resolution=0.2)
    assert set(out.obs["cell_type"]) == {"A", "B"}
