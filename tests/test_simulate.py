import csv
import gzip


def _reads(path):
    with gzip.open(path, "rt") as fh:
        lines = fh.read().splitlines()
    return lines[0::4], lines[1::4]


def test_layout_matches_10x_v3(sim_dir):
    names1, r1 = _reads(sim_dir / "sim1_S1_L001_R1_001.fastq.gz")
    names2, r2 = _reads(sim_dir / "sim1_S1_L001_R2_001.fastq.gz")
    assert names1 == names2, "R1 and R2 must stay paired and in the same order"
    assert {len(s) for s in r1} == {28}, "R1 is 16 bp barcode + 12 bp UMI"
    assert {len(s) for s in r2} == {91}


def test_barcodes_are_on_the_whitelist(sim_dir):
    whitelist = set((sim_dir / "whitelist.txt").read_text().split())
    _, r1 = _reads(sim_dir / "sim2_S1_L001_R1_001.fastq.gz")
    assert {s[:16] for s in r1} <= whitelist


def test_truth_and_samplesheet(sim_dir):
    truth = list(csv.DictReader(open(sim_dir / "truth.csv")))
    assert len(truth) == 80
    assert {r["cell_type"] for r in truth} == {"T cell", "B cell", "Monocyte"}
    sheet = list(csv.DictReader(open(sim_dir / "samplesheet.csv")))
    assert [r["sample"] for r in sheet] == ["sim1", "sim2"]


def test_gtf_has_both_strands_and_mito_genes(sim_dir):
    gtf = (sim_dir / "genes.gtf").read_text()
    assert "\t+\t" in gtf and "\t-\t" in gtf
    assert 'gene_name "MT-CO1"' in gtf
