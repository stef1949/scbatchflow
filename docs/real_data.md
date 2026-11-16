# First real-data run: 10x PBMC, v2 vs v3 chemistry

10x Genomics publishes two libraries made from the **same PBMC sample** with different
chemistries: `pbmc_1k_v2` (about 1,000 cells, v2) and `pbmc_1k_v3` (about 1,200 cells, v3).
Because the cells are the same, any difference between the two is technical. That makes
the pair a clean batch-correction benchmark with a known right answer: after correction the
two libraries should mix completely within each cell type.

The data are small enough for a laptop (a few GB of FASTQ) and are released under CC BY 4.0.

The commands below assume the `docker` profile and that you are in the repo folder. Define a
shortcut for running tools from the image:

```bash
alias sc='docker run --rm -u $(id -u):$(id -g) -v $HOME:$HOME -w $PWD scbatchflow:0.2.0'
```

## 1. Download the FASTQs

On the 10x Genomics datasets site, open "1k PBMCs from a Healthy Donor (v3 chemistry)" and
the matching "(v2 chemistry)" page and copy the link of each FASTQ archive.

```bash
mkdir -p ~/data ~/ref && cd ~/data
wget "PASTE_V3_FASTQ_LINK"
wget "PASTE_V2_FASTQ_LINK"
for f in *.tar; do tar -xf "$f"; done
cd -
```

## 2. Get a human reference

```bash
sc kb ref -d human -i ~/ref/index.idx -g ~/ref/t2g.txt
```

Create `~/ref` yourself first, as above. If Docker creates the folder it is owned by root
and the download fails with a permission error.

Record which reference you used. It belongs in the write-up.

## 3. Write the samplesheet

```bash
sc scbf samplesheet --out ~/data/samplesheet.csv \
    --library ~/data/pbmc_1k_v3_fastqs,v3,10xv3 \
    --library ~/data/pbmc_1k_v2_fastqs,v2,10xv2
cat ~/data/samplesheet.csv
```

It finds every lane, pairs R1 with R2 and writes full paths. Expect four rows.

## 4. Run

```bash
docker compose up -d
nextflow run . -profile docker \
    --input ~/data/samplesheet.csv \
    --kb_index ~/ref/index.idx --t2g ~/ref/t2g.txt \
    --expected_cells 1000 --batch_key batch \
    --outdir ~/results_pbmc --run_id pbmc_v2_v3 \
    --db_url postgresql+psycopg2://scbf:scbf@host.docker.internal:5432/scbatchflow
```

If a step fails, fix the cause and rerun the same command with `-resume`. Then open
`~/results_pbmc/results/report.html`.

## 5. What to check before believing the numbers

1. **Mapping rate** in `results/quant/*/` and the `sample` table. A low rate on one library
   usually means the wrong `technology`.
2. **Cell calls** against the cell counts 10x reports for each dataset. A large gap means
   `--expected_cells` or the calling rule needs attention.
3. **`Unassigned` rate** at the top of the report. If more than about a fifth of cells have no marker
   label, the biology metrics rest on a biased subset.
4. **The uncorrected row of "Before and after"** in the report: you should see the
   chemistry split. If you do not, there is no batch effect to benchmark.
5. **The trade-off plot.** A method that wins on mixing by losing on biology has
   over-corrected; say so rather than reporting the overall score alone.

## 6. Results

_To be filled in after the run: the benchmark table from `summary.md`, the trade-off plot,
before/after UMAPs, and two or three sentences on which method you would choose and why._

## Next: a tumour dataset

Pick a public multi-patient 10x tumour cohort whose **raw FASTQs are openly downloadable**
(many human cohorts keep raw reads under controlled access, so check before committing),
write a marker file for its expected cell types, and rerun with `--markers`. Patient is then
the batch, and the interesting question changes: malignant cells are expected to separate
by patient, so a method that mixes them completely is removing biology.
