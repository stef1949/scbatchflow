"""Simulate a tiny 10x Chromium v3 experiment from first principles.

The simulator writes a reference (genome FASTA + GTF), a barcode on-list and
paired FASTQ files with the real 10x v3 read layout:

    R1 = 16 bp cell barcode + 12 bp UMI
    R2 = cDNA read from the 3' end of the transcript

It exists so the whole pipeline, including pseudoalignment, can be tested in
seconds with a known ground truth (cell type per barcode, and a known
gene-wise batch effect between samples). It is not a model of real data.
"""

from __future__ import annotations

import csv
import gzip
from dataclasses import dataclass
from pathlib import Path

import numpy as np

BASES = np.array(list("ACGT"))
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

# Real marker names on synthetic sequences, so the default marker file works.
CELL_TYPE_MARKERS: dict[str, list[str]] = {
    "T cell": ["CD3D", "CD3E", "IL7R", "TRAC"],
    "B cell": ["MS4A1", "CD79A", "CD79B", "BANK1"],
    "Monocyte": ["LYZ", "S100A8", "S100A9", "CD14"],
}
MITO_GENES = ["MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP6", "MT-CYB"]


@dataclass
class SimConfig:
    n_background_genes: int = 100
    transcript_length: int = 700
    spacer_length: int = 150
    cells_per_sample: int = 150
    n_samples: int = 2
    mean_umis: int = 300
    read_length: int = 91
    empty_droplets: int = 60
    batch_sd: float = 0.6
    unmappable_fraction: float = 0.04
    seed: int = 1


def _random_seq(rng: np.random.Generator, n: int) -> str:
    return "".join(rng.choice(BASES, size=n))


def _revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def _write_fastq(path: Path, records: list[tuple[str, str]]) -> None:
    with gzip.open(path, "wt", compresslevel=4) as fh:
        for name, seq in records:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def simulate(outdir: str | Path, config: SimConfig | None = None) -> dict[str, Path]:
    """Write reference, FASTQs, samplesheet and truth table. Returns the paths."""
    cfg = config or SimConfig()
    rng = np.random.default_rng(cfg.seed)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    marker_genes = [g for genes in CELL_TYPE_MARKERS.values() for g in genes]
    background = [f"BG{i:03d}" for i in range(1, cfg.n_background_genes + 1)]
    genes = marker_genes + MITO_GENES + background
    n_genes = len(genes)

    # ---- reference: one single-exon transcript per gene on one contig ----
    transcripts: dict[str, str] = {}
    contig_parts: list[str] = []
    gtf_rows: list[str] = []
    pos = 0
    for i, gene in enumerate(genes):
        spacer = _random_seq(rng, cfg.spacer_length)
        contig_parts.append(spacer)
        pos += len(spacer)
        mrna = _random_seq(rng, cfg.transcript_length)
        strand = "-" if i % 3 == 2 else "+"
        genomic = mrna if strand == "+" else _revcomp(mrna)
        start, end = pos + 1, pos + len(genomic)
        contig_parts.append(genomic)
        pos = end
        transcripts[gene] = mrna
        gid, tid = f"SIMG{i:05d}", f"SIMT{i:05d}"
        attrs = f'gene_id "{gid}"; gene_name "{gene}"; gene_biotype "protein_coding";'
        tattrs = f'{attrs} transcript_id "{tid}"; transcript_name "{gene}-201";'
        gtf_rows.append(f"chrT\tsim\tgene\t{start}\t{end}\t.\t{strand}\t.\t{attrs}")
        gtf_rows.append(f"chrT\tsim\ttranscript\t{start}\t{end}\t.\t{strand}\t.\t{tattrs}")
        gtf_rows.append(
            f'chrT\tsim\texon\t{start}\t{end}\t.\t{strand}\t.\t{tattrs} exon_number "1";'
        )
    contig = "".join(contig_parts)
    genome_fa = out / "genome.fa"
    with open(genome_fa, "w") as fh:
        fh.write(">chrT\n")
        for i in range(0, len(contig), 60):
            fh.write(contig[i : i + 60] + "\n")
    gtf = out / "genes.gtf"
    gtf.write_text("\n".join(gtf_rows) + "\n")

    # ---- expression programmes ----
    base = rng.lognormal(mean=0.0, sigma=0.8, size=n_genes)
    base[[genes.index(g) for g in MITO_GENES]] *= 1.5
    cell_types = list(CELL_TYPE_MARKERS)
    profiles = {}
    for ct in cell_types:
        prof = base.copy()
        for other, mk in CELL_TYPE_MARKERS.items():
            idx = [genes.index(g) for g in mk]
            prof[idx] = base.mean() * (12.0 if other == ct else 0.05)
        # a weaker programme across background genes, so types differ beyond markers
        bg_idx = rng.choice(np.arange(len(marker_genes) + len(MITO_GENES), n_genes), 20, False)
        prof[bg_idx] *= rng.lognormal(0.0, 0.9, size=20)
        profiles[ct] = prof

    # ---- barcodes ----
    n_bc = cfg.n_samples * (cfg.cells_per_sample + cfg.empty_droplets) + 500
    barcodes = sorted({_random_seq(rng, 16) for _ in range(n_bc * 2)})
    rng.shuffle(barcodes)
    barcodes = barcodes[:n_bc]
    (out / "whitelist.txt").write_text("\n".join(barcodes) + "\n")

    truth_rows: list[dict[str, str]] = []
    sheet_rows: list[dict[str, str]] = []
    bc_cursor = 0
    for s in range(cfg.n_samples):
        sample = f"sim{s + 1}"
        batch = f"batch{s + 1}"
        # gene-wise multiplicative batch effect; sample 1 is the reference
        effect = np.ones(n_genes) if s == 0 else rng.lognormal(0.0, cfg.batch_sd, size=n_genes)
        depth = 1.0 if s == 0 else 0.7
        # shift composition between batches so batch and biology are partly confounded
        props = np.array([0.5, 0.3, 0.2]) if s % 2 == 0 else np.array([0.25, 0.35, 0.4])
        r1: list[tuple[str, str]] = []
        r2: list[tuple[str, str]] = []

        def emit(bc: str, gene_idx: np.ndarray, sample: str = sample,
                 r1: list = r1, r2: list = r2) -> None:
            for g in gene_idx:
                umi = _random_seq(rng, 12)
                mrna = transcripts[genes[g]]
                for _ in range(int(rng.integers(1, 4))):  # PCR duplicates share a UMI
                    lo = len(mrna) - 400
                    start = int(rng.integers(lo, len(mrna) - cfg.read_length))
                    name = f"{sample}:{len(r1) + 1}"
                    r1.append((name, bc + umi))
                    r2.append((name, mrna[start : start + cfg.read_length]))

        for _ in range(cfg.cells_per_sample):
            bc = barcodes[bc_cursor]
            bc_cursor += 1
            ct = cell_types[int(rng.choice(len(cell_types), p=props))]
            p = profiles[ct] * effect
            n_umi = max(int(rng.lognormal(np.log(cfg.mean_umis * depth), 0.35)), 60)
            emit(bc, rng.choice(n_genes, size=n_umi, p=p / p.sum()))
            truth_rows.append(
                {"cell_id": f"{sample}_{bc}", "sample": sample, "batch": batch, "cell_type": ct}
            )
        ambient = sum(profiles.values()) * effect
        for _ in range(cfg.empty_droplets):  # ambient RNA only: must be removed by cell calling
            bc = barcodes[bc_cursor]
            bc_cursor += 1
            emit(bc, rng.choice(n_genes, size=int(rng.integers(2, 12)), p=ambient / ambient.sum()))

        # reads that match nothing in the reference (adapter dimers, intergenic, contamination)
        cell_bcs = [r["cell_id"].split("_", 1)[1] for r in truth_rows if r["sample"] == sample]
        for _ in range(int(cfg.unmappable_fraction * len(r1))):
            name = f"{sample}:{len(r1) + 1}"
            r1.append((name, cell_bcs[int(rng.integers(len(cell_bcs)))] + _random_seq(rng, 12)))
            r2.append((name, _random_seq(rng, cfg.read_length)))

        order = rng.permutation(len(r1))
        f1 = out / f"{sample}_S1_L001_R1_001.fastq.gz"
        f2 = out / f"{sample}_S1_L001_R2_001.fastq.gz"
        _write_fastq(f1, [r1[i] for i in order])
        _write_fastq(f2, [r2[i] for i in order])
        sheet_rows.append(
            {"sample": sample, "batch": batch, "fastq_1": str(f1.resolve()),
             "fastq_2": str(f2.resolve())}
        )

    truth = out / "truth.csv"
    with open(truth, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["cell_id", "sample", "batch", "cell_type"])
        w.writeheader()
        w.writerows(truth_rows)
    sheet = out / "samplesheet.csv"
    with open(sheet, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "batch", "fastq_1", "fastq_2"])
        w.writeheader()
        w.writerows(sheet_rows)

    return {"genome": genome_fa, "gtf": gtf, "whitelist": out / "whitelist.txt",
            "samplesheet": sheet, "truth": truth}
