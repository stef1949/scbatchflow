"""Build a samplesheet from folders of 10x FASTQ files."""

from __future__ import annotations

import csv
import re
from pathlib import Path

LANE = re.compile(r"^(?P<sample>.+?)_S\d+_L\d+_R1_\d+\.f(ast)?q\.gz$")


def build_samplesheet(libraries: list[str], out: str | Path) -> list[dict[str, str]]:
    """Each library is 'DIR,BATCH[,TECHNOLOGY]'. Every R1 file in DIR becomes one lane row,
    paired with the R2 file of the same name. Paths are written absolute."""
    rows: list[dict[str, str]] = []
    for spec in libraries:
        parts = [p.strip() for p in spec.split(",")]
        if len(parts) not in (2, 3) or not all(parts):
            raise ValueError(f"expected DIR,BATCH[,TECHNOLOGY], got '{spec}'")
        folder, batch = Path(parts[0]).expanduser().resolve(), parts[1]
        technology = parts[2] if len(parts) == 3 else ""
        if not folder.is_dir():
            raise FileNotFoundError(f"no such folder: {folder}")
        r1_files = sorted(p for p in folder.iterdir() if "_R1_" in p.name
                          and p.name.endswith((".fastq.gz", ".fq.gz")))
        if not r1_files:
            raise FileNotFoundError(f"no *_R1_*.fastq.gz files in {folder}")
        for r1 in r1_files:
            r2 = r1.with_name(r1.name.replace("_R1_", "_R2_"))
            if not r2.exists():
                raise FileNotFoundError(f"missing mate for {r1.name}: expected {r2.name}")
            match = LANE.match(r1.name)
            rows.append({"sample": match["sample"] if match else folder.name, "batch": batch,
                         "technology": technology, "fastq_1": str(r1), "fastq_2": str(r2)})
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample", "batch", "technology", "fastq_1",
                                                "fastq_2"])
        writer.writeheader()
        writer.writerows(rows)
    return rows
