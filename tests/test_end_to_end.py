"""Reads to labelled cells on simulated data, without Nextflow.

Skipped when kallisto|bustools (kb-python) is not installed.
"""

import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest

from scbatchflow.integrate import METHODS, run_method
from scbatchflow.metrics import evaluate
from scbatchflow.preprocess import label_by_markers, merge_and_normalise
from scbatchflow.qc import cell_qc
from scbatchflow.quant import kb_to_anndata
from scbatchflow.simulate import CELL_TYPE_MARKERS

pytestmark = pytest.mark.skipif(shutil.which("kb") is None, reason="kb-python not installed")


@pytest.fixture(scope="module")
def merged(sim_dir, tmp_path_factory):
    work = tmp_path_factory.mktemp("kb")
    subprocess.run(["kb", "ref", "-i", work / "i.idx", "-g", work / "t2g.txt", "-f1",
                    work / "cdna.fa", sim_dir / "genome.fa", sim_dir / "genes.gtf"],
                   check=True, capture_output=True)
    adatas = []
    for s in ("sim1", "sim2"):
        subprocess.run(["kb", "count", "-i", work / "i.idx", "-g", work / "t2g.txt", "-x", "10xv3",
                        "-w", sim_dir / "whitelist.txt", "-t", "2", "-o", work / s,
                        sim_dir / f"{s}_S1_L001_R1_001.fastq.gz",
                        sim_dir / f"{s}_S1_L001_R2_001.fastq.gz"],
                       check=True, capture_output=True)
        raw, stats = kb_to_anndata(work / s, s, f"batch{s[-1]}")
        assert 90 < stats["pct_pseudoaligned"] < 100, "the 4% junk reads must not align"
        filtered, summary = cell_qc(raw, expected_cells=40, min_genes=15, doublets=False)
        assert summary["n_cells_called"] == 40, "empty droplets must be dropped, cells kept"
        adatas.append(filtered)
    adata = merge_and_normalise(adatas, n_hvg=60, n_pcs=10)
    adata.obs["marker_label"] = label_by_markers(adata, CELL_TYPE_MARKERS)
    return adata


def test_marker_labels_recover_simulated_cell_types(merged, sim_dir):
    truth = pd.read_csv(sim_dir / "truth.csv").set_index("cell_id")["cell_type"]
    labels = merged.obs["marker_label"]
    assert (labels == truth.reindex(labels.index)).mean() > 0.95


@pytest.mark.parametrize("method", sorted(METHODS))
def test_every_method_returns_a_valid_embedding(merged, method):
    emb = run_method(merged, method)
    assert emb.shape[0] == merged.n_obs and np.isfinite(emb).all()


def test_correction_improves_batch_mixing_without_losing_biology(merged):
    scores = {m: evaluate(run_method(merged, m), merged.obs["sample"], merged.obs["marker_label"])
              for m in ("pca", "combat")}
    assert scores["combat"]["batch_score"] > scores["pca"]["batch_score"]
    assert scores["combat"]["bio_score"] > 0.7
