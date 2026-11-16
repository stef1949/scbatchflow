"""Convert kallisto|bustools output into AnnData and collect mapping statistics."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import pandas as pd
import scipy.io
import scipy.sparse as sp


def kb_to_anndata(kb_dir: str | Path, sample: str, batch: str) -> tuple[ad.AnnData, dict]:
    """Read the unfiltered `kb count` matrix; return AnnData and mapping stats."""
    kb_dir = Path(kb_dir)
    counts = kb_dir / "counts_unfiltered"
    mtx = sp.csr_matrix(scipy.io.mmread(counts / "cells_x_genes.mtx"), dtype="float32")
    barcodes = (counts / "cells_x_genes.barcodes.txt").read_text().split()
    gene_ids = (counts / "cells_x_genes.genes.txt").read_text().split()
    names_file = counts / "cells_x_genes.genes.names.txt"
    names = names_file.read_text().split("\n")[: len(gene_ids)] if names_file.exists() else []
    if len(names) != len(gene_ids):
        names = gene_ids
    names = [n if n else g for n, g in zip(names, gene_ids, strict=True)]

    var = pd.DataFrame({"gene_id": gene_ids, "gene_name": names})
    # gene symbols as index, falling back to the stable ID when a symbol repeats
    dup = var["gene_name"].duplicated(keep=False)
    var.index = var["gene_name"].where(~dup, var["gene_name"] + "_" + var["gene_id"])
    var.index.name = None
    obs = pd.DataFrame(
        {"barcode": barcodes, "sample": sample, "batch": batch},
        index=[f"{sample}_{bc}" for bc in barcodes],
    )
    adata = ad.AnnData(X=mtx, obs=obs, var=var)

    stats: dict = {"sample": sample, "batch": batch, "n_barcodes": int(adata.n_obs)}
    run_info = kb_dir / "run_info.json"
    if run_info.exists():
        info = json.loads(run_info.read_text())
        stats["n_reads"] = int(info.get("n_processed", 0))
        stats["n_pseudoaligned"] = int(info.get("n_pseudoaligned", 0))
        stats["pct_pseudoaligned"] = float(info.get("p_pseudoaligned", 0.0))
    return adata, stats
