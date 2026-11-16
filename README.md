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

## Quick start

Requirements: Nextflow 24.04 or newer, and either Docker or Python 3.10+.

```bash
# 1. End-to-end test on simulated reads (about 90 seconds on 2 CPUs)
pip install . multiqc                      # or: docker build -t scbatchflow:0.2.0 .
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
with `docker run --rm -u $(id -u):$(id -g) -v $PWD:$PWD -w $PWD scbatchflow:0.2.0`.

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

- Automated tests use simulated data only. Real-data results belong in
  [docs/real_data.md](docs/real_data.md) once run.
- Marker-based labels are only as good as the marker file, and the biology metrics inherit
  that. For a new tissue, check the `Unassigned` rate first.
- Order-of-magnitude cell calling misses low-RNA cell types in heterogeneous samples.
- No ambient-RNA correction.
- ComBat runs on dense matrices: fine to tens of thousands of cells, not beyond.
- Batch metrics assume the batches share cell types. Where composition differs strongly,
  read iLISI and kBET alongside batch ASW, which is computed within cell type.

## Roadmap

- [ ] Real-data benchmark: 10x PBMC v2 vs v3 chemistry, then a public tumour cohort
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
