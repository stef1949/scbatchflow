"""Assemble one self-contained HTML report from a results directory.

Figures are embedded as base64, so report.html can be copied anywhere (or attached
to an email) and opened without the rest of the output.
"""

from __future__ import annotations

import base64
import html
from pathlib import Path

import pandas as pd

from . import __version__

CSS = """
:root{--bg:#f9f9f7;--surface:#fcfcfb;--ink:#0b0b0b;--muted:#52514e;--line:#e1e0d9;--accent:#2a78d6}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.55 system-ui,-apple-system,
"Segoe UI",sans-serif}
main{max-width:1080px;margin:0 auto;padding:32px 20px 80px;display:flex;flex-direction:column;
gap:28px}
h1{font-size:30px;margin:0 0 4px;line-height:1.15}
h2{font-size:21px;margin:0 0 6px}
p{margin:0;max-width:72ch}
.sub{color:var(--muted)}
section{background:var(--surface);border:1px solid var(--line);border-radius:12px;padding:20px;
display:flex;flex-direction:column;gap:14px}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px}
.tile{background:var(--surface);border:1px solid var(--line);border-radius:12px;padding:14px 16px}
.tile b{display:block;font-size:26px;font-variant-numeric:tabular-nums;line-height:1.2}
.tile span{color:var(--muted);font-size:14px}
img{max-width:100%;height:auto;border-radius:8px}
img.small{max-width:min(100%,540px)}
.left td,.left th{text-align:left}
.scroll{overflow-x:auto}
table{border-collapse:collapse;font-size:14px;font-variant-numeric:tabular-nums}
th,td{padding:6px 12px;border-bottom:1px solid var(--line);text-align:right;white-space:nowrap}
th:first-child,td:first-child{text-align:left}
th{color:var(--muted);font-weight:500}
tr.best td{font-weight:700}
.how{border-left:3px solid var(--accent);padding:2px 0 2px 12px;color:var(--muted);font-size:15px}
code{background:#eeeeea;padding:1px 5px;border-radius:4px;font-size:.9em}
footer{color:var(--muted);font-size:13px}
"""


def _img(path: Path, alt: str, small: bool = False) -> str:
    if not path.exists():
        return ""
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    cls = ' class="small"' if small else ""
    return f'<img{cls} alt="{html.escape(alt)}" src="data:image/png;base64,{data}">'


def _table(df: pd.DataFrame, formats: dict[str, str] | None = None, best: str | None = None,
           index: bool = False, left: bool = False) -> str:
    formats = formats or {}
    if index:
        df = df.reset_index()
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in df.columns)
    rows = []
    for _, r in df.iterrows():
        cells = []
        for c in df.columns:
            v = r[c]
            if isinstance(v, float) and pd.isna(v):
                text = "n/a"
            elif c in formats:
                text = formats[c].format(v)
            else:
                text = str(v)
            cells.append(f"<td>{html.escape(text)}</td>")
        cls = ' class="best"' if best is not None and str(r[df.columns[0]]) == best else ""
        rows.append(f"<tr{cls}>{''.join(cells)}</tr>")
    cls = "scroll left" if left else "scroll"
    return (f'<div class="{cls}"><table><thead><tr>{head}</tr></thead>'
            f'<tbody>{"".join(rows)}</tbody></table></div>')


def build_report(results_dir: str | Path, run_label: str = "", batch_key: str = "sample") -> Path:
    res = Path(results_dir)
    plots = res / "plots"
    samples = pd.read_csv(res / "samples.csv")
    cells = pd.read_csv(res / "cells.csv", dtype={"cluster": str})
    metrics = pd.read_csv(res / "integration_metrics.csv")
    markers = pd.read_csv(res / "cluster_markers.csv", dtype={"cluster": str})
    method = (res / "final_method.txt").read_text().strip()
    unassigned = float((cells["marker_label"] == "Unassigned").mean())

    tiles = [
        (f"{len(cells):,}", "cells after QC"),
        (f"{samples['sample_id'].nunique()}", "samples"),
        (f"{cells['cluster'].nunique()}", "clusters"),
        (method, "embedding chosen"),
        (f"{unassigned:.1%}", "cells without a marker label"),
    ]
    tile_html = "".join(f'<div class="tile"><b>{html.escape(a)}</b><span>{html.escape(b)}</span>'
                        "</div>" for a, b in tiles)

    funnel_cols = {
        "sample_id": "sample", "batch": "batch", "n_reads": "reads",
        "pct_pseudoaligned": "% mapped", "n_barcodes": "barcodes",
        "n_cells_called": "cells called", "n_removed_low_genes": "low genes",
        "n_removed_high_mito": "high mito", "n_removed_doublets": "doublets",
        "n_cells_pass_qc": "cells kept", "median_umis": "median UMIs",
        "median_genes": "median genes", "median_pct_mito": "median % mito",
    }
    funnel = samples[[c for c in funnel_cols if c in samples.columns]].rename(columns=funnel_cols)
    funnel_fmt = {"reads": "{:,.0f}", "% mapped": "{:.1f}", "barcodes": "{:,.0f}",
                  "cells called": "{:,.0f}", "low genes": "{:,.0f}", "high mito": "{:,.0f}",
                  "doublets": "{:,.0f}", "cells kept": "{:,.0f}", "median UMIs": "{:,.0f}",
                  "median genes": "{:,.0f}", "median % mito": "{:.1f}"}

    metric_names = {"method": "method", "ilisi": "iLISI", "kbet_acceptance": "kBET",
                    "batch_asw": "batch ASW", "clisi": "cLISI", "label_asw": "label ASW",
                    "ari": "ARI", "nmi": "NMI", "batch_score": "batch removal",
                    "bio_score": "biology kept", "overall_score": "overall"}
    mt = metrics[[c for c in metric_names if c in metrics.columns]].rename(columns=metric_names)
    mt_fmt = {c: "{:.3f}" for c in mt.columns if c != "method"}

    qc_imgs = "".join(_img(p, f"QC for {p.stem}") for p in sorted(plots.glob("qc_*.png"))
                      if p.name != "qc_violins.png")
    top = (markers[markers["rank"] <= 5].sort_values(["cluster", "rank"])
           .groupby("cluster", sort=False)["gene"].agg(", ".join).reset_index())
    names = cells.groupby("cluster")["cell_type"].agg(lambda s: s.mode().iat[0])
    sizes = cells["cluster"].value_counts()
    top.insert(1, "cell type", top["cluster"].map(names))
    top.insert(2, "cells", top["cluster"].map(sizes))
    top = top.rename(columns={"gene": "top markers"})
    top = top.iloc[sorted(range(len(top)), key=lambda i: (not top["cluster"].iat[i].isdigit(),
                                                         int(top["cluster"].iat[i])
                                                         if top["cluster"].iat[i].isdigit()
                                                         else 0))]

    title = "scbatchflow report" + (f": {run_label}" if run_label else "")
    body = f"""
<header><h1>{html.escape(title)}</h1>
<p class="sub">Raw reads to annotated cells, with batch integration benchmarked on the way.</p>
</header>
<div class="tiles">{tile_html}</div>

<section><h2>1. Read mapping and the QC funnel</h2>
<p class="how">How to read it: follow each sample left to right. A mapping rate far below 90%
usually means the wrong chemistry or reference. A large gap between cells called and cells
kept means QC is removing a lot; look at which rule did it.</p>
{_table(funnel, funnel_fmt)}
{qc_imgs}
{_img(plots / "qc_violins.png", "QC metrics per sample after filtering")}
<p class="sub">Violins show cells that passed QC. The line spans the middle half of cells and
the dot marks the median.</p></section>

<section><h2>2. Integration benchmark</h2>
<p class="how">How to read it: each row is one way of combining the samples, scored 0 to 1.
The first three columns measure how well batches mix, the next four how much cell-type
structure survives. A method that scores high on mixing and low on biology has
over-corrected. <b>{html.escape(method)}</b> was used for everything that follows.</p>
{_table(mt, mt_fmt, best=method)}
{_img(plots / "metric_heatmap.png", "Heatmap of integration scores per method")}
{_img(plots / "integration_tradeoff.png", "Batch removal against biology kept", small=True)}
</section>

<section><h2>3. Before and after correction</h2>
<p class="how">How to read it: the top row is the uncorrected data. In the left column, colours
separated into their own islands mean a batch effect; after a good correction the colours
overlap. In the right column, cell types should stay apart in every row.</p>
{_img(plots / "method_grid.png", "UMAP per method coloured by batch and by marker label")}
</section>

<section><h2>4. The final embedding</h2>
{_img(plots / "umap_overview.png", "Final UMAP by sample, cell type and cluster")}
{_img(plots / "composition.png", "Cell-type composition per sample")}
<p class="sub">Clusters are named after the most common marker label among their cells.
Composition differences between samples can be real biology; they are not corrected.</p>
</section>

<section><h2>5. What defines each cluster</h2>
<p class="how">How to read it: a large, dark dot means most cells in that cluster express the
gene strongly. Check that each cluster's markers match the cell type it was given.</p>
{_img(plots / "marker_dots.png", "Dot plot of top marker genes per cluster")}
{_table(top, left=True)}</section>

<footer>Generated by scbatchflow {__version__}. Integration was run over
<code>{html.escape(batch_key)}</code>. Tables behind every figure are in this folder as CSV,
and in the results database when <code>--db_url</code> was set.</footer>
"""
    page = (f'<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width, initial-scale=1">'
            f"<title>{html.escape(title)}</title><style>{CSS}</style></head>"
            f"<body><main>{body}</main></body></html>")
    out = res / "report.html"
    out.write_text(page)
    return out
