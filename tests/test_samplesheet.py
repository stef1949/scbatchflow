import csv

import pytest

from scbatchflow.samplesheet import build_samplesheet


def touch(folder, *names):
    folder.mkdir(parents=True, exist_ok=True)
    for n in names:
        (folder / n).write_bytes(b"")


def test_lanes_become_rows_with_absolute_paths(tmp_path):
    v3 = tmp_path / "pbmc_1k_v3_fastqs"
    touch(v3, "pbmc_1k_v3_S1_L001_R1_001.fastq.gz", "pbmc_1k_v3_S1_L001_R2_001.fastq.gz",
          "pbmc_1k_v3_S1_L002_R1_001.fastq.gz", "pbmc_1k_v3_S1_L002_R2_001.fastq.gz",
          "pbmc_1k_v3_S1_L001_I1_001.fastq.gz")
    v2 = tmp_path / "pbmc_1k_v2_fastqs"
    touch(v2, "pbmc_1k_v2_S1_L001_R1_001.fastq.gz", "pbmc_1k_v2_S1_L001_R2_001.fastq.gz")
    out = tmp_path / "sheet.csv"
    rows = build_samplesheet([f"{v3},v3,10xv3", f"{v2},v2,10xv2"], out)
    assert [r["sample"] for r in rows] == ["pbmc_1k_v3", "pbmc_1k_v3", "pbmc_1k_v2"]
    assert [r["technology"] for r in rows] == ["10xv3", "10xv3", "10xv2"]
    assert all(r["fastq_1"].startswith("/") and "_R1_" in r["fastq_1"] for r in rows)
    assert all("_R2_" in r["fastq_2"] for r in rows)
    assert list(csv.DictReader(open(out)))[2]["batch"] == "v2"


def test_missing_mate_is_an_error(tmp_path):
    lib = tmp_path / "lib"
    touch(lib, "s_S1_L001_R1_001.fastq.gz")
    with pytest.raises(FileNotFoundError, match="missing mate"):
        build_samplesheet([f"{lib},b1"], tmp_path / "s.csv")


def test_bad_spec_is_an_error(tmp_path):
    with pytest.raises(ValueError):
        build_samplesheet(["just_a_folder"], tmp_path / "s.csv")
