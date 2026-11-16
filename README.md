# scbatchflow

[![CI](https://github.com/stef1949/scbatchflow/actions/workflows/ci.yml/badge.svg)](https://github.com/stef1949/scbatchflow/actions/workflows/ci.yml)

A Nextflow pipeline that takes single-cell RNA-seq from **raw 10x Chromium reads** to
**annotated cells**, benchmarks batch-integration methods against each other on the way,
and loads every result into **PostgreSQL**.

Most single-cell pipelines pick one integration method and trust it. This one treats
integration as a measured choice: each method is scored on how well it mixes batches *and*
on how much biology it keeps, and the trade-off is reported before a final embedding is chosen.

```mermaid
flowchart LR
    A[FASTQ + samplesheet] --> B[FastQC / MultiQC]
    A --> C[kallisto / bustools<br/>pseudoalign, correct barcodes, count UMIs]
    R[genome FASTA + GTF] --> C
    C --> D[Cell calling + QC<br/>per sample]
    D --> E[Merge, normalise, HVG, PCA<br/>marker labels per sample]
    E --> F1[PCA baseline]
    E --> F2[Harmony]
    E --> F3[ComBat]
    F1 & F2 & F3 --> G[Metrics<br/>iLISI, kBET, ASW, cLISI, ARI, NMI]
    G --> H[Pick embedding, cluster, UMAP, markers]
    H --> I[(PostgreSQL)]
    H --> J[report.html, plots, h5ad]
```

## What it does

| Stage | Tools | Output |
|---|---|---|
| Read QC | FastQC, MultiQC | Per-file quality report |
| Quantification | kallisto, bustools | Cells x genes UMI counts per sample, mapping statistics |
| Cell calling and QC | Scanpy, Scrublet | Empty droplets, damaged cells and doublets removed, with a per-sample audit |
| Integration benchmark | PCA baseline, Harmony, ComBat | Each method scored on 7 metrics for batch mixing and biology kept |
| Annotation | Leiden, marker scoring, Wilcoxon | Clusters, cell types and marker genes on the best embedding |
| Reporting | Matplotlib, PostgreSQL | A self-contained `report.html`, CSVs, an `.h5ad`, and a results database |

Everything runs from one command, in Docker, locally, on Slurm or in CI.

## Results

### Real data: 10x PBMC, v2 against v3 chemistry

10x Genomics released two libraries of about 1,000 cells from the same PBMC sample, made
with different chemistries. The cells are the same, so any difference between the two
libraries is technical: a clean test of batch correction with a known right answer.

| | v2 library | v3 library |
|---|---:|---:|
| Reads | 76.9 M | 66.6 M |
| Pseudoaligned | 65.0% | 57.4% |
| Cells called | 964 | 1,114 |
| Removed: high mitochondrial / doublets | 7 / 9 | 35 / 7 |
| Cells kept | 948 | 1,072 |
| Median UMIs per cell | 4,038 | 7,407 |
| Median genes per cell | 1,225 | 2,142 |

v3 detected about 1.8 times as many genes per cell as v2, the kind of technical difference
integration has to remove. The mapping rates are in the range usual for a transcriptome-only
index, which does not count reads from introns.

| Method | Batch removal | Biology kept | Overall |
|---|---:|---:|---:|
| **Harmony** (chosen) | 0.807 | 0.840 | **0.827** |
| ComBat | 0.635 | 0.826 | 0.750 |
| PCA, no correction | 0.373 | 0.786 | 0.621 |

Harmony mixed the two libraries best (iLISI 0.87 against 0.58 for ComBat) while keeping cell
types at least as well separated. ComBat, which won on the simulated data below, left the
monocytes partly split by chemistry. The two datasets disagree on the best method, which is
the reason the pipeline benchmarks per dataset instead of fixing one method.

![UMAP per method, coloured by chemistry and by marker label](docs/images/pbmc/method_grid.png)

*Top row: without correction every cell type splits into a v2 half and a v3 half. Bottom row:
after Harmony the two chemistries overlap and the five cell types stay distinct.*

<table>
<tr>
<td><img src="docs/images/pbmc/metric_heatmap.png" alt="Integration scores per method"></td>
<td><img src="docs/images/pbmc/integration_tradeoff.png" alt="Batch removal against biology kept"></td>
</tr>
<tr>
<td><img src="docs/images/pbmc/barcode_rank_v3.png" alt="Barcode rank plot for the v3 library"></td>
<td><img src="docs/images/pbmc/composition.png" alt="Cell-type composition per library"></td>
</tr>
</table>

Composition is almost identical between the libraries (T cells 40% and 39%, monocytes 32%
and 31%, B cells 15% and 17%, NK cells 11% and 10%), as it should be for the same cells.
2,020 cells passed QC, 2.2% had no clear marker label, and Leiden found 14 clusters.

What this run does not settle:

- cLISI was 1.00 for every method, including no correction, so it did not separate the
  methods here; the biology score rests on label ASW, ARI and NMI.
- Two clusters (3 and 6) are defined by ribosomal and mitochondrial genes rather than
  identity markers, so they more likely reflect cell state or quality than distinct types.
- The v3 cells carry more mitochondrial reads (median 10.3% against 3.8%), so the
  per-sample adaptive cut-off differed between libraries (10.8% for v2, 20% for v3).
- Two libraries from one donor is a small test. A multi-donor cohort is the next step.

### Simulated data with known truth

The test profile simulates two 10x v3 samples of 150 cells each, with three cell types,
a gene-wise batch effect, a shift in composition between the samples, empty droplets and
4% unmappable reads. Because the truth is known, the run can be checked rather than
eyeballed:

- 96.2% of reads pseudoaligned, matching the 4% of junk reads put in.
- Cell calling kept exactly the 150 real cells per sample and dropped every empty droplet.
  QC then removed 1 cell, leaving 299.
- Marker-based labels matched the simulated cell type for all 299 cells; the final
  cluster labels matched for 298.

| Method | Batch removal | Biology kept | Overall |
|---|---:|---:|---:|
| **ComBat** (chosen) | 0.751 | 0.904 | **0.843** |
| Harmony | 0.416 | 0.910 | 0.712 |
| PCA, no correction | 0.313 | 0.781 | 0.594 |

Scores run 0 to 1; overall = 0.4 x batch removal + 0.6 x biology kept. Uncorrected data
scores lowest on biology as well as mixing, because the batch effect splits each cell type
in two. ComBat wins here because the simulated effect follows exactly the per-gene shift
and scale that ComBat assumes, so this result validates the scoring, not the method. On the
real PBMC data above, Harmony won instead.

![UMAP per method, coloured by sample and by marker label](docs/images/method_grid.png)

*Top row: without correction the two samples form separate islands within every cell type.
Lower rows: after correction the samples overlap while the three cell types stay apart.*

<table>
<tr>
<td><img src="docs/images/metric_heatmap.png" alt="Integration scores per method"></td>
<td><img src="docs/images/integration_tradeoff.png" alt="Batch removal against biology kept"></td>
</tr>
<tr>
<td><img src="docs/images/qc_sim1.png" alt="Barcode rank plot and QC for one sample"></td>
<td><img src="docs/images/marker_dots.png" alt="Marker genes per cluster"></td>
</tr>
</table>

The full report from this run is in [docs/example_report.html](docs/example_report.html);
download it and open it in a browser.

## Quick start

Requirements: Nextflow 24.04 or newer, and either Docker or Python 3.10+.

```bash
# 1. End-to-end test on simulated reads (about 90 seconds on 2 CPUs)
pip install . multiqc                      # or: docker build -t scbatchflow:0.2.1 .
nextflow run . -profile test               # add ,docker to run inside the container

# 2. With the results database
docker compose up -d
nextflow run . -profile test \
    --db_url postgresql+psycopg2://scbf:scbf@localhost:5432/scbatchflow
psql postgresql://scbf:scbf@localhost:5432/scbatchflow -v run_id=<run name> -f db/queries.sql
```

Then open `results_test/results/report.html` in a browser.

The `test` profile simulates two 10x v3 samples (150 cells each, three cell types, a known
gene-wise batch effect, empty droplets and 4% unmappable reads) and runs every step on them,
pseudoalignment included. Add `--skip_fastqc` if FastQC is not installed locally.

## Running on real data

```bash
# 1. A samplesheet from folders of 10x FASTQs: one --library per folder, as DIR,BATCH,CHEMISTRY
scbf samplesheet --out samplesheet.csv \
    --library data/pbmc_1k_v3_fastqs,v3,10xv3 \
    --library data/pbmc_1k_v2_fastqs,v2,10xv2

# 2. A prebuilt human reference
mkdir -p ref && kb ref -d human -i ref/index.idx -g ref/t2g.txt

# 3. Run
nextflow run . -profile docker \
    --input samplesheet.csv \
    --kb_index ref/index.idx --t2g ref/t2g.txt \
    --batch_key batch --run_id pbmc_v2_v3 \
    --db_url postgresql+psycopg2://scbf:scbf@host.docker.internal:5432/scbatchflow
```

With the `docker` profile, `scbf` and `kb` live inside the image; prefix those two commands
with `docker run --rm -u $(id -u):$(id -g) -v $PWD:$PWD -w $PWD scbatchflow:0.2.1`.

The samplesheet has one row per FASTQ pair, and rows that share a `sample` are treated as
lanes. The `technology` column sets the chemistry per sample (default `--technology 10xv3`),
so libraries made with different kits can be combined in one run. To build the index from a
genome instead, pass `--genome_fasta` and `--gtf` in place of `--kb_index` and `--t2g`.
[docs/real_data.md](docs/real_data.md) walks through a first real run on public 10x PBMC data.

### Main options

| Option | Default | Meaning |
|---|---|---|
| `--methods` | `pca,harmony,combat` | Integration methods to benchmark |
| `--final_method` | `best` | Embedding used for clustering: best overall score, or a method name |
| `--batch_key` | `sample` | Column to integrate over (`sample` or `batch`) |
| `--markers` | `assets/markers_pbmc.json` | Cell type to marker genes; swap for your tissue |
| `--expected_cells` | `3000` | Used by cell calling |
| `--min_genes`, `--max_pct_mito`, `--mito_mads` | `200`, `20`, `5` | Per-cell QC thresholds |
| `--skip_doublets`, `--skip_fastqc` | `false` | Turn off Scrublet or FastQC/MultiQC |
| `--db_url` | none | SQLAlchemy URL; omit to skip the database |

## Outputs

Start with `results/report.html`: one self-contained page with every figure and table,
which can be copied or emailed without the rest of the folder.

| Report section | Figures | Question it answers |
|---|---|---|
| Read mapping and QC funnel | Barcode rank plot and QC scatter per sample, violins | Were the right barcodes called cells, and what did QC remove? |
| Integration benchmark | Score heatmap, trade-off plot | Which method mixes batches without losing cell types? |
| Before and after | UMAP per method, by batch and by marker label | Does the correction look like the scores say? |
| Final embedding | UMAP by sample, cell type and cluster; composition bars | What is in each sample? |
| Cluster markers | Dot plot, top-marker table | Do the labels match the genes? |

```
results/
├── results/
│   ├── report.html                the visual report
│   ├── summary.md                 run summary and the benchmark table as text
│   ├── integrated.h5ad            final AnnData: counts, embedding, UMAP, clusters, labels
│   ├── integration_metrics.csv, samples.csv, cells.csv, cluster_markers.csv
│   └── plots/                     every figure as a PNG
├── qc/, quant/                    per-sample QC and mapping summaries (JSON)
├── fastqc/, multiqc/              read-level quality report
└── pipeline_info/                 Nextflow timeline, report and trace
```

`scbf report --results results/results` rebuilds the page from an existing results folder.
For interactive exploration, open `integrated.h5ad` in cellxgene or a Scanpy notebook.

## Results database

Five tables, defined once in [`src/scbatchflow/db.py`](src/scbatchflow/db.py) and exported
to [`db/schema.sql`](db/schema.sql):

| Table | Grain | Holds |
|---|---|---|
| `run` | one pipeline run | version, parameters, chosen method |
| `sample` | sample per run | reads, mapping rate, QC funnel |
| `cell` | cell per run | QC metrics, cluster, cell type, UMAP coordinates |
| `integration_metric` | method x metric per run | every benchmark score |
| `cluster_marker` | gene per cluster per run | Wilcoxon marker genes |

Loading is one transaction, and reloading a `run_id` replaces it, so `-resume` never
duplicates rows. [`db/queries.sql`](db/queries.sql) has worked examples: the QC funnel per
sample, composition per sample, method ranking against the uncorrected baseline, and
clusters dominated by a single sample.

## How integration is scored

Scores are scaled to 0-1 (higher is better) following the scIB benchmark conventions.

| Group | Metric | Asks |
|---|---|---|
| Batch removal | iLISI | How many batches are represented in each cell's neighbourhood? |
| | kBET acceptance | Does each neighbourhood's batch mix match the global mix? |
| | Batch ASW | Within a cell type, are batches inseparable? |
| Biology | cLISI | Is each neighbourhood a single cell type? |
| | Label ASW | Are cell types compact and separated? |
| | ARI, NMI | Does clustering recover the labels? |

Overall = 0.4 x batch removal + 0.6 x biology.

The biological reference labels come from marker-gene scores computed **inside each sample,
before any correction**, so no integration method can influence the labels it is judged
against. Cells without a clear marker call are labelled `Unassigned` and left out of the
biology metrics.

## Design decisions

- **kallisto|bustools for quantification.** It handles barcode correction and UMI counting
  in a few GB of memory, so a full run fits on a laptop, and it installs from PyPI, which
  keeps the test profile free of large downloads.
- **Cell calling uses the order-of-magnitude rule** rather than EmptyDrops: transparent, fast,
  and adequate for well-loaded libraries. See limitations.
- **One Python package, thin Nextflow processes.** Every step is an `scbf` sub-command with
  unit tests; Nextflow only wires them together. Steps can be run and debugged by hand.
- **Adding a method is one function.** Register it in `METHODS` in
  [`integrate.py`](src/scbatchflow/integrate.py); the workflow, metrics, report and database
  pick it up by name.
- **Simulated reads as the test fixture.** The simulator writes real 10x v3 read layouts
  with known cell types, so the tests check answers (cell calls, labels, mixing), not just
  that the code runs.

## Limitations

- Automated tests use simulated data; real data has been run once, on the 10x PBMC pair
  above.
- Marker-based labels are only as good as the marker file, and the biology metrics inherit
  that. For a new tissue, check the `Unassigned` rate first.
- Order-of-magnitude cell calling misses low-RNA cell types in heterogeneous samples.
- No ambient-RNA correction.
- ComBat runs on dense matrices: fine to tens of thousands of cells, not beyond.
- Batch metrics assume the batches share cell types. Where composition differs strongly,
  read iLISI and kBET alongside batch ASW, which is computed within cell type.

## Roadmap

- [x] Real-data benchmark: 10x PBMC v2 vs v3 chemistry
- [ ] A public multi-donor tumour cohort
- [ ] Add scVI and a HarmonizeNN encoder as integration methods
- [ ] Pseudobulk differential expression per cell type between conditions
- [ ] EmptyDrops-style cell calling and ambient-RNA correction
- [ ] nf-test module tests and a published container on GHCR

## Development

```bash
pip install -e ".[dev]"
ruff check . && pytest
scbf db-schema --out db/schema.sql      # after changing the table definitions
```

## References

- Melsted et al. 2021, *Nat Biotechnol*: kallisto|bustools
- Korsunsky et al. 2019, *Nat Methods*: Harmony and LISI
- Johnson et al. 2007, *Biostatistics*: ComBat
- Büttner et al. 2019, *Nat Methods*: kBET
- Luecken et al. 2022, *Nat Methods*: scIB benchmark and metric scaling
- Wolf et al. 2018, *Genome Biol*: Scanpy
- Wolock et al. 2019, *Cell Syst*: Scrublet

## Licence

MIT
