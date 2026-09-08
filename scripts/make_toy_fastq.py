#!/usr/bin/env python3
"""
make_toy_fastq.py

Generates a tiny, fully SYNTHETIC paired-end FASTQ dataset shaped like a
10x Genomics Visium HD 3' library, purely so you can see -- and make --
real, valid FASTQ files by hand.

This is NOT real sequencing data. It won't produce a biologically
meaningful result if you run it through Space Ranger. It exists so you
can open the two files it writes and recognize every line, and so you
have working code to adapt if you ever want to generate test data for
a pipeline.

See docs/fastq_primer.md for the explanation this script is paired with.

Usage:
    python3 make_toy_fastq.py
Produces:
    data/toy_fastq/toy_S1_L001_R1_001.fastq.gz  (barcode + UMI reads)
    data/toy_fastq/toy_S1_L001_R2_001.fastq.gz  (cDNA reads)
"""
import gzip
import os
import random

random.seed(42)  # reproducible output -- same command, same files, every time

BASES = "ACGT"
OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "toy_fastq")

# ---------------------------------------------------------------------------
# Visium HD's 3' chemistry splits each fragment across two reads, the same
# way Chromium 3' single-cell does:
#   R1 = spatial barcode + UMI   (tells you WHERE on the slide and WHICH
#                                  original molecule -- no biology here)
#   R2 = the actual cDNA insert (the biology: a fragment of a transcript)
# A real run also has I1/I2 index reads used to demultiplex samples on the
# flow cell; we skip those here since they don't touch biology or position.
#
# R1_LEN = 43bp is the figure commonly cited for Visium HD (UMI + spatial
# barcode combined). Treat it as illustrative -- for the authoritative
# number, check 10x's own Visium HD library structure documentation before
# relying on it for anything beyond this toy example.
# ---------------------------------------------------------------------------
R1_LEN = 43
R2_LEN = 91
N_READS = 8


def random_seq(length: int) -> str:
    return "".join(random.choice(BASES) for _ in range(length))


def random_qual(length: int, min_q: int = 30, max_q: int = 40) -> str:
    """Phred+33: each quality score Q becomes the ASCII character chr(Q + 33).
    Q30 = ASCII 63 ('?') = a 1-in-1000 chance the base call is wrong.
    Q40 = ASCII 73 ('I') = a 1-in-10000 chance. Real reads vary more than
    this (and usually degrade toward the read's end); we keep it in a
    tight, all-good-quality range just to keep the example readable."""
    return "".join(chr(random.randint(min_q, max_q) + 33) for _ in range(length))


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # A handful of fake "spatial barcodes". In a real Visium HD run these
    # come from a fixed whitelist where every sequence maps to a known
    # physical bin on the slide -- that mapping is what lets Space Ranger
    # turn "which barcode" into "which x,y position." Here we just invent
    # a few so you can see multiple reads sharing one fake bin.
    fake_barcodes = [random_seq(16) for _ in range(4)]

    r1_path = os.path.join(OUT_DIR, "toy_S1_L001_R1_001.fastq.gz")
    r2_path = os.path.join(OUT_DIR, "toy_S1_L001_R2_001.fastq.gz")

    with gzip.open(r1_path, "wt") as r1, gzip.open(r2_path, "wt") as r2:
        for i in range(N_READS):
            # Illumina CASAVA-style header. Same instrument/run/flowcell/
            # lane/tile/x/y between R1 and R2 is what PAIRS them -- it's
            # how you know read i of R1 and read i of R2 came from the
            # same physical molecule.
            base_id = f"SYNTH:1:TOYFLOWCELL:1:1101:{1000+i}:{2000+i}"

            barcode = random.choice(fake_barcodes)
            umi = random_seq(12)
            r1_seq = (barcode + umi)[:R1_LEN].ljust(R1_LEN, "A")
            r2_seq = random_seq(R2_LEN)

            r1.write(f"@{base_id} 1:N:0:AAAAAAAA\n{r1_seq}\n+\n{random_qual(len(r1_seq))}\n")
            r2.write(f"@{base_id} 2:N:0:AAAAAAAA\n{r2_seq}\n+\n{random_qual(len(r2_seq))}\n")

    print(f"Wrote {N_READS} read pairs to:\n  {r1_path}\n  {r2_path}")
    print("Peek at them with: zcat data/toy_fastq/toy_S1_L001_R1_001.fastq.gz | head -8")


if __name__ == "__main__":
    main()
