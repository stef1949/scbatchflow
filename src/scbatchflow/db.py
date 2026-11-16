"""Results database: schema and loader.

Written against SQLAlchemy Core so the same code targets PostgreSQL in
production and SQLite in unit tests. `db/schema.sql` is the PostgreSQL DDL
generated from these definitions (see `scbf db-schema`).
"""

from __future__ import annotations

import datetime as dt
import json
import math
from pathlib import Path

import pandas as pd
import sqlalchemy as sa
from sqlalchemy.schema import CreateIndex, CreateTable

metadata = sa.MetaData()

run = sa.Table(
    "run", metadata,
    sa.Column("run_id", sa.String(64), primary_key=True),
    sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
    sa.Column("pipeline_version", sa.String(32), nullable=False),
    sa.Column("final_method", sa.String(32)),
    sa.Column("params", sa.Text),
)

sample = sa.Table(
    "sample", metadata,
    sa.Column("run_id", sa.String(64), sa.ForeignKey("run.run_id", ondelete="CASCADE"),
              primary_key=True),
    sa.Column("sample_id", sa.String(128), primary_key=True),
    sa.Column("batch", sa.String(128), nullable=False),
    sa.Column("n_reads", sa.BigInteger),
    sa.Column("pct_pseudoaligned", sa.Float),
    sa.Column("n_barcodes", sa.Integer),
    sa.Column("umi_threshold", sa.Float),
    sa.Column("n_cells_called", sa.Integer),
    sa.Column("n_removed_low_genes", sa.Integer),
    sa.Column("n_removed_high_mito", sa.Integer),
    sa.Column("n_removed_doublets", sa.Integer),
    sa.Column("n_cells_pass_qc", sa.Integer),
    sa.Column("median_umis", sa.Float),
    sa.Column("median_genes", sa.Float),
    sa.Column("median_pct_mito", sa.Float),
)

cell = sa.Table(
    "cell", metadata,
    sa.Column("run_id", sa.String(64), primary_key=True),
    sa.Column("cell_id", sa.String(160), primary_key=True),
    sa.Column("sample_id", sa.String(128), nullable=False),
    sa.Column("n_umis", sa.Integer),
    sa.Column("n_genes", sa.Integer),
    sa.Column("pct_mito", sa.Float),
    sa.Column("doublet_score", sa.Float),
    sa.Column("marker_label", sa.String(128)),
    sa.Column("cluster", sa.String(32)),
    sa.Column("cell_type", sa.String(128)),
    sa.Column("umap_1", sa.Float),
    sa.Column("umap_2", sa.Float),
    sa.ForeignKeyConstraint(["run_id", "sample_id"], ["sample.run_id", "sample.sample_id"],
                            ondelete="CASCADE"),
)
sa.Index("ix_cell_sample", cell.c.run_id, cell.c.sample_id)
sa.Index("ix_cell_type", cell.c.run_id, cell.c.cell_type)

integration_metric = sa.Table(
    "integration_metric", metadata,
    sa.Column("run_id", sa.String(64), sa.ForeignKey("run.run_id", ondelete="CASCADE"),
              primary_key=True),
    sa.Column("method", sa.String(32), primary_key=True),
    sa.Column("metric", sa.String(32), primary_key=True),
    sa.Column("value", sa.Float),
)

cluster_marker = sa.Table(
    "cluster_marker", metadata,
    sa.Column("run_id", sa.String(64), sa.ForeignKey("run.run_id", ondelete="CASCADE"),
              primary_key=True),
    sa.Column("cluster", sa.String(32), primary_key=True),
    sa.Column("gene", sa.String(128), primary_key=True),
    sa.Column("rank", sa.Integer, nullable=False),
    sa.Column("logfc", sa.Float),
    sa.Column("pval_adj", sa.Float),
    sa.Column("score", sa.Float),
)


def postgres_ddl() -> str:
    from sqlalchemy.dialects import postgresql

    dialect = postgresql.dialect()
    parts = []
    for table in metadata.sorted_tables:
        parts.append(str(CreateTable(table).compile(dialect=dialect)).strip() + ";")
        for index in sorted(table.indexes, key=lambda i: i.name):
            parts.append(str(CreateIndex(index).compile(dialect=dialect)).strip() + ";")
    return "\n\n".join(parts) + "\n"


def _clean(value):
    """NaN and numpy scalars are not valid SQL parameters."""
    if value is None:
        return None
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    return value


def _rows(df: pd.DataFrame) -> list[dict]:
    return [{k: _clean(v) for k, v in rec.items()} for rec in df.to_dict("records")]


def load_run(
    db_url: str,
    run_id: str,
    pipeline_version: str,
    final_method: str,
    params: dict,
    samples: pd.DataFrame,
    cells: pd.DataFrame,
    metrics: pd.DataFrame,
    markers: pd.DataFrame,
) -> dict[str, int]:
    """Create tables if needed and load one run in a single transaction.

    Re-loading a run_id replaces it, so a resumed pipeline never duplicates rows.
    """
    engine = sa.create_engine(db_url)
    metadata.create_all(engine)
    sample_cols = [c.name for c in sample.columns if c.name != "run_id"]
    cell_cols = [c.name for c in cell.columns if c.name != "run_id"]
    counts = {}
    with engine.begin() as conn:
        for table in (cluster_marker, integration_metric, cell, sample, run):
            conn.execute(sa.delete(table).where(table.c.run_id == run_id))
        conn.execute(sa.insert(run), [{
            "run_id": run_id, "created_at": dt.datetime.now(dt.timezone.utc),
            "pipeline_version": pipeline_version, "final_method": final_method,
            "params": json.dumps(params, sort_keys=True, default=str),
        }])
        for table, df, cols in (
            (sample, samples.reindex(columns=sample_cols), sample_cols),
            (cell, cells.reindex(columns=cell_cols), cell_cols),
            (integration_metric, metrics, None),
            (cluster_marker, markers, None),
        ):
            rows = _rows(df if cols is None else df[cols])
            for r in rows:
                r["run_id"] = run_id
            for i in range(0, len(rows), 5000):
                conn.execute(sa.insert(table), rows[i : i + 5000])
            counts[table.name] = len(rows)
    engine.dispose()
    return counts


def write_schema(path: str | Path) -> None:
    Path(path).write_text(
        "-- PostgreSQL schema for scbatchflow results.\n"
        "-- Generated by `scbf db-schema`; edit src/scbatchflow/db.py, not this file.\n\n"
        + postgres_ddl()
    )
