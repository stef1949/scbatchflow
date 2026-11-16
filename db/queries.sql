-- Example questions the results database answers. PostgreSQL syntax.

-- 1. Which integration method won each run, and by how much over no correction?
SELECT m.run_id,
       m.method,
       ROUND(m.value::numeric, 3)                    AS overall_score,
       ROUND((m.value - base.value)::numeric, 3)     AS gain_over_uncorrected
FROM integration_metric m
JOIN integration_metric base
  ON base.run_id = m.run_id AND base.method = 'pca' AND base.metric = 'overall_score'
WHERE m.metric = 'overall_score'
ORDER BY m.run_id, m.value DESC;

-- 2. The batch-removal / biology trade-off, one row per method.
SELECT method,
       ROUND(MAX(value) FILTER (WHERE metric = 'batch_score')::numeric, 3) AS batch_removal,
       ROUND(MAX(value) FILTER (WHERE metric = 'bio_score')::numeric, 3)   AS bio_conservation
FROM integration_metric
WHERE run_id = :'run_id'
GROUP BY method
ORDER BY bio_conservation DESC;

-- 3. Per-sample QC funnel: where were barcodes lost?
SELECT sample_id,
       n_reads,
       pct_pseudoaligned,
       n_barcodes,
       n_cells_called,
       n_removed_low_genes + n_removed_high_mito + n_removed_doublets AS n_removed_qc,
       n_cells_pass_qc,
       ROUND(100.0 * n_cells_pass_qc / NULLIF(n_cells_called, 0), 1) AS pct_cells_kept
FROM sample
WHERE run_id = :'run_id'
ORDER BY sample_id;

-- 4. Cell-type composition per sample, as proportions.
SELECT sample_id,
       cell_type,
       COUNT(*)                                                             AS n_cells,
       ROUND(100.0 * COUNT(*) / SUM(COUNT(*)) OVER (PARTITION BY sample_id), 1) AS pct_of_sample
FROM cell
WHERE run_id = :'run_id'
GROUP BY sample_id, cell_type
ORDER BY sample_id, n_cells DESC;

-- 5. Clusters dominated by one sample: candidates for residual batch effect.
WITH per_cluster AS (
    SELECT cluster, sample_id, COUNT(*) AS n,
           SUM(COUNT(*)) OVER (PARTITION BY cluster) AS cluster_total
    FROM cell
    WHERE run_id = :'run_id'
    GROUP BY cluster, sample_id
)
SELECT cluster, sample_id, n, ROUND(100.0 * n / cluster_total, 1) AS pct_of_cluster
FROM per_cluster
WHERE n > 0.9 * cluster_total
ORDER BY cluster;

-- 6. Top five marker genes per cluster.
SELECT cluster, STRING_AGG(gene, ', ' ORDER BY rank) AS top_markers
FROM cluster_marker
WHERE run_id = :'run_id' AND rank <= 5
GROUP BY cluster
ORDER BY cluster;
