import pandas as pd
import sqlalchemy as sa

from scbatchflow.db import load_run, postgres_ddl


def frames():
    samples = pd.DataFrame([{"sample_id": "s1", "batch": "b1", "n_reads": 1000,
                             "pct_pseudoaligned": 95.5, "n_cells_pass_qc": 2,
                             "median_umis": float("nan"), "not_a_column": 1}])
    cells = pd.DataFrame([
        {"cell_id": "s1_A", "sample_id": "s1", "n_umis": 500, "n_genes": 200, "pct_mito": 2.0,
         "cluster": "0", "cell_type": "T cell", "umap_1": 0.1, "umap_2": 0.2},
        {"cell_id": "s1_B", "sample_id": "s1", "n_umis": 700, "n_genes": 300, "pct_mito": 3.0,
         "cluster": "1", "cell_type": "B cell", "umap_1": 1.1, "umap_2": 1.2},
    ])
    metrics = pd.DataFrame([{"method": "harmony", "metric": "ilisi", "value": 0.8},
                            {"method": "pca", "metric": "ilisi", "value": float("nan")}])
    markers = pd.DataFrame([{"cluster": "0", "gene": "CD3E", "rank": 1, "logfc": 5.0,
                             "pval_adj": 1e-9, "score": 12.0}])
    return samples, cells, metrics, markers


def load(url, run_id="r1"):
    s, c, m, k = frames()
    return load_run(url, run_id, "0.1.0", "harmony", {"n_hvg": 2000}, s, c, m, k)


def test_load_and_query(tmp_path):
    url = f"sqlite:///{tmp_path / 't.db'}"
    assert load(url) == {"sample": 1, "cell": 2, "integration_metric": 2, "cluster_marker": 1}
    with sa.create_engine(url).connect() as conn:
        assert conn.execute(sa.text("select count(*) from cell")).scalar() == 2
        # NaN must arrive as NULL, not as a float
        assert conn.execute(sa.text(
            "select value from integration_metric where method='pca'")).scalar() is None
        assert conn.execute(sa.text("select median_umis from sample")).scalar() is None


def test_reloading_a_run_replaces_it(tmp_path):
    url = f"sqlite:///{tmp_path / 't.db'}"
    load(url)
    load(url)
    load(url, run_id="r2")
    with sa.create_engine(url).connect() as conn:
        assert conn.execute(sa.text("select count(*) from run")).scalar() == 2
        assert conn.execute(sa.text("select count(*) from cell")).scalar() == 4


def test_postgres_ddl_declares_keys():
    ddl = postgres_ddl()
    assert "CREATE TABLE cell" in ddl
    assert "FOREIGN KEY(run_id, sample_id) REFERENCES sample" in ddl
    assert "CREATE INDEX ix_cell_type" in ddl
