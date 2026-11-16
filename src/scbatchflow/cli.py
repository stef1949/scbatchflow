"""Command-line entry point: one sub-command per pipeline step."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import warnings
from pathlib import Path

from . import __version__


def _cmd_simulate(a: argparse.Namespace) -> None:
    from .simulate import SimConfig, simulate

    cfg = SimConfig(cells_per_sample=a.cells, n_samples=a.samples, mean_umis=a.mean_umis,
                    seed=a.seed)
    paths = simulate(a.outdir, cfg)
    print(json.dumps({k: str(v) for k, v in paths.items()}, indent=2))


def _cmd_kb2h5ad(a: argparse.Namespace) -> None:
    from .quant import kb_to_anndata

    adata, stats = kb_to_anndata(a.kb_dir, a.sample, a.batch)
    adata.write_h5ad(a.out)
    Path(a.stats).write_text(json.dumps(stats, indent=2))


def _cmd_qc(a: argparse.Namespace) -> None:
    import anndata as ad

    from .qc import cell_qc

    adata, summary = cell_qc(
        ad.read_h5ad(a.input), expected_cells=a.expected_cells, min_genes=a.min_genes,
        max_pct_mito=a.max_pct_mito, mito_mads=None if a.mito_mads <= 0 else a.mito_mads,
        doublets=not a.no_doublets, seed=a.seed,
    )
    if adata.n_obs == 0:
        sys.exit(f"QC removed every cell in {a.input}; check thresholds and input quality")
    stats = json.loads(Path(a.stats).read_text()) if a.stats else {}
    if a.plot:
        from .viz import plot_sample_qc

        plot_sample_qc(ad.read_h5ad(a.input), set(adata.obs_names), summary, a.plot)
    adata.write_h5ad(a.out)
    Path(a.summary).write_text(json.dumps({**stats, **summary}, indent=2))


def _cmd_merge(a: argparse.Namespace) -> None:
    import anndata as ad

    from .preprocess import label_by_markers, load_markers, merge_and_normalise

    adata = merge_and_normalise([ad.read_h5ad(p) for p in a.inputs], n_hvg=a.n_hvg,
                                n_pcs=a.n_pcs, seed=a.seed)
    adata.obs["marker_label"] = label_by_markers(adata, load_markers(a.markers), seed=a.seed)
    adata.obs["marker_label"] = adata.obs["marker_label"].astype("category")
    adata.write_h5ad(a.out)


def _cmd_integrate(a: argparse.Namespace) -> None:
    import anndata as ad
    import numpy as np

    from .integrate import run_method

    np.save(a.out, run_method(ad.read_h5ad(a.input), a.method, batch_key=a.batch_key,
                              seed=a.seed))


def _cmd_metrics(a: argparse.Namespace) -> None:
    import anndata as ad
    import numpy as np

    from .metrics import evaluate

    adata = ad.read_h5ad(a.input, backed="r")
    scores = evaluate(np.load(a.embedding), adata.obs[a.batch_key], adata.obs["marker_label"],
                      seed=a.seed)
    Path(a.out).write_text(json.dumps({"method": a.method, **scores}, indent=2))


def _cmd_finalise(a: argparse.Namespace) -> None:
    import anndata as ad
    import numpy as np
    import pandas as pd

    from . import viz
    from .annotate import finalise
    from .report import build_report

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame([json.loads(Path(p).read_text()) for p in a.metrics]).set_index("method")
    table = table.sort_values("overall_score", ascending=False)
    method = table.index[0] if a.final_method == "best" else a.final_method
    emb_path = {Path(p).name.split(".")[0]: p for p in a.embeddings}
    if method not in emb_path:
        sys.exit(f"no embedding for final method '{method}'; have {sorted(emb_path)}")

    adata, markers = finalise(ad.read_h5ad(a.input), np.load(emb_path[method]), method,
                              resolution=a.resolution, seed=a.seed)
    adata.write_h5ad(out / "integrated.h5ad")
    cells = adata.obs.copy()
    cells.insert(0, "cell_id", cells.index)
    cells = cells.rename(columns={"sample": "sample_id"})
    cells["umap_1"], cells["umap_2"] = adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1]
    cells.to_csv(out / "cells.csv", index=False)
    markers.to_csv(out / "cluster_markers.csv", index=False)
    table.to_csv(out / "integration_metrics.csv")
    samples = pd.DataFrame([json.loads(Path(p).read_text()) for p in a.sample_summaries])
    samples = samples.rename(columns={"sample": "sample_id"}).sort_values("sample_id")
    samples.to_csv(out / "samples.csv", index=False)
    (out / "final_method.txt").write_text(method + "\n")

    plots = out / "plots"
    plots.mkdir(exist_ok=True)
    for src in a.qc_plots or []:
        shutil.copy(src, plots / ("qc_" + Path(src).name.replace(".qc.png", ".png")))
    viz.plot_final_umaps(adata, plots, batch_key=a.batch_key)
    viz.plot_method_grid(adata, {m: np.load(f) for m, f in emb_path.items()},
                         plots / "method_grid.png", batch_key=a.batch_key, seed=a.seed)
    viz.plot_metric_heatmap(table, plots / "metric_heatmap.png")
    viz.plot_tradeoff(table, plots / "integration_tradeoff.png")
    viz.plot_composition(adata, plots / "composition.png", batch_key=a.batch_key)
    viz.plot_qc_violins(adata, plots / "qc_violins.png")
    if len(markers):
        viz.plot_marker_dots(adata, markers, plots / "marker_dots.png")

    show = ["ilisi", "kbet_acceptance", "batch_asw", "clisi", "label_asw", "ari", "nmi",
            "batch_score", "bio_score", "overall_score"]
    lines = [
        "# scbatchflow run summary", "",
        f"- Cells after QC: {adata.n_obs:,} across {adata.obs['sample'].nunique()} samples",
        f"- Final embedding: **{method}** "
        f"({'best overall score' if a.final_method == 'best' else 'set by --final_method'})",
        f"- Clusters: {adata.obs['cluster'].nunique()}",
        f"- Cells without a marker label: {(adata.obs['marker_label'] == 'Unassigned').mean():.1%}"
        " (left out of the biology metrics)", "",
        "## Integration benchmark (0-1, higher is better)", "",
        "| method | " + " | ".join(show) + " |",
        "|---|" + "---|" * len(show),
    ]
    for m, row in table.iterrows():
        lines.append(f"| {m} | " + " | ".join(f"{row[c]:.3f}" for c in show) + " |")
    lines += ["", "## Cells per sample and type", "",
              pd.crosstab(adata.obs["sample"], adata.obs["cell_type"]).to_markdown()
              if _has_tabulate() else
              pd.crosstab(adata.obs["sample"], adata.obs["cell_type"]).to_string(), ""]
    lines += ["", "Open `report.html` in this folder for the figures.", ""]
    (out / "summary.md").write_text("\n".join(lines))
    build_report(out, run_label=a.run_label or "", batch_key=a.batch_key)


def _cmd_report(a: argparse.Namespace) -> None:
    from .report import build_report

    print(build_report(a.results, run_label=a.run_label or "", batch_key=a.batch_key))


def _cmd_samplesheet(a: argparse.Namespace) -> None:
    from .samplesheet import build_samplesheet

    rows = build_samplesheet(a.library, a.out)
    print(f"wrote {a.out}: {len(rows)} lane(s) across {len({r['sample'] for r in rows})} sample(s)")


def _has_tabulate() -> bool:
    try:
        import tabulate  # noqa: F401
        return True
    except ImportError:
        return False


def _cmd_load_db(a: argparse.Namespace) -> None:
    import pandas as pd

    from .db import load_run

    res = Path(a.results)
    metrics = pd.read_csv(res / "integration_metrics.csv").melt(
        id_vars="method", var_name="metric", value_name="value")
    db_url = a.db_url or os.environ.get("SCBF_DB_URL")
    if not db_url:
        sys.exit("no database: pass --db-url or set SCBF_DB_URL")
    counts = load_run(
        db_url=db_url, run_id=a.run_id, pipeline_version=__version__,
        final_method=(res / "final_method.txt").read_text().strip(),
        params=json.loads(Path(a.params_file).read_text()) if a.params_file else {},
        samples=pd.read_csv(res / "samples.csv"),
        cells=pd.read_csv(res / "cells.csv", dtype={"cluster": str}),
        metrics=metrics,
        markers=pd.read_csv(res / "cluster_markers.csv", dtype={"cluster": str}),
    )
    print(json.dumps({"run_id": a.run_id, "rows_loaded": counts}, indent=2))


def _cmd_db_schema(a: argparse.Namespace) -> None:
    from .db import write_schema

    write_schema(a.out)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="scbf", description=__doc__)
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("simulate", help="write a tiny simulated 10x v3 dataset")
    s.add_argument("--outdir", required=True)
    s.add_argument("--cells", type=int, default=150)
    s.add_argument("--samples", type=int, default=2)
    s.add_argument("--mean-umis", type=int, default=300)
    s.add_argument("--seed", type=int, default=1)
    s.set_defaults(func=_cmd_simulate)

    s = sub.add_parser("kb2h5ad", help="convert kb count output to h5ad")
    s.add_argument("--kb-dir", required=True)
    s.add_argument("--sample", required=True)
    s.add_argument("--batch", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--stats", required=True)
    s.set_defaults(func=_cmd_kb2h5ad)

    s = sub.add_parser("qc", help="cell calling and per-cell QC for one sample")
    s.add_argument("--input", required=True)
    s.add_argument("--stats")
    s.add_argument("--out", required=True)
    s.add_argument("--summary", required=True)
    s.add_argument("--expected-cells", type=int, default=3000)
    s.add_argument("--min-genes", type=int, default=200)
    s.add_argument("--max-pct-mito", type=float, default=20.0)
    s.add_argument("--mito-mads", type=float, default=5.0, help="0 disables the MAD rule")
    s.add_argument("--no-doublets", action="store_true")
    s.add_argument("--plot", help="write a QC figure for this sample to this PNG")
    s.add_argument("--seed", type=int, default=0)
    s.set_defaults(func=_cmd_qc)

    s = sub.add_parser("merge", help="merge samples, normalise, PCA, marker labels")
    s.add_argument("--inputs", nargs="+", required=True)
    s.add_argument("--markers", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--n-hvg", type=int, default=2000)
    s.add_argument("--n-pcs", type=int, default=30)
    s.add_argument("--seed", type=int, default=0)
    s.set_defaults(func=_cmd_merge)

    s = sub.add_parser("integrate", help="run one integration method")
    s.add_argument("--input", required=True)
    s.add_argument("--method", required=True)
    s.add_argument("--batch-key", default="sample")
    s.add_argument("--out", required=True)
    s.add_argument("--seed", type=int, default=0)
    s.set_defaults(func=_cmd_integrate)

    s = sub.add_parser("metrics", help="score one embedding")
    s.add_argument("--input", required=True)
    s.add_argument("--embedding", required=True)
    s.add_argument("--method", required=True)
    s.add_argument("--batch-key", default="sample")
    s.add_argument("--out", required=True)
    s.add_argument("--seed", type=int, default=0)
    s.set_defaults(func=_cmd_metrics)

    s = sub.add_parser("finalise", help="pick the embedding, cluster, annotate, report")
    s.add_argument("--input", required=True)
    s.add_argument("--embeddings", nargs="+", required=True)
    s.add_argument("--metrics", nargs="+", required=True)
    s.add_argument("--sample-summaries", nargs="+", required=True)
    s.add_argument("--qc-plots", nargs="*", help="per-sample QC figures to include")
    s.add_argument("--batch-key", default="sample")
    s.add_argument("--run-label", help="title shown on the report")
    s.add_argument("--final-method", default="best")
    s.add_argument("--resolution", type=float, default=0.8)
    s.add_argument("--outdir", required=True)
    s.add_argument("--seed", type=int, default=0)
    s.set_defaults(func=_cmd_finalise)

    s = sub.add_parser("report", help="rebuild report.html from an existing results directory")
    s.add_argument("--results", required=True)
    s.add_argument("--batch-key", default="sample")
    s.add_argument("--run-label")
    s.set_defaults(func=_cmd_report)

    s = sub.add_parser("samplesheet", help="write a samplesheet from folders of 10x FASTQs")
    s.add_argument("--library", action="append", required=True, metavar="DIR,BATCH[,TECHNOLOGY]",
                   help="one per library folder, e.g. data/pbmc_1k_v3_fastqs,v3,10xv3")
    s.add_argument("--out", default="samplesheet.csv")
    s.set_defaults(func=_cmd_samplesheet)

    s = sub.add_parser("load-db", help="load a results directory into the database")
    s.add_argument("--results", required=True)
    s.add_argument("--db-url", help="SQLAlchemy URL; defaults to $SCBF_DB_URL")
    s.add_argument("--run-id", required=True)
    s.add_argument("--params-file", help="JSON file of run parameters to record")
    s.set_defaults(func=_cmd_load_db)

    s = sub.add_parser("db-schema", help="write the PostgreSQL DDL")
    s.add_argument("--out", required=True)
    s.set_defaults(func=_cmd_db_schema)
    return p


def main(argv: list[str] | None = None) -> None:
    # scanpy 1.13 renamed several keyword arguments; the old names used here still work
    # on every supported version, so the deprecation notices are noise in pipeline logs.
    warnings.filterwarnings("ignore", category=FutureWarning, module="scanpy")
    warnings.filterwarnings("ignore", category=FutureWarning, module="scbatchflow")
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
