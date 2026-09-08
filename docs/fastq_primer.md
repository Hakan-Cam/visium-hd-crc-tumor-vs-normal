# FASTQ, and where it fits before this pipeline

Everything else in this repo starts from an *already processed* Visium HD
output (`filtered_feature_bc_matrix.h5`, `tissue_positions.parquet`, ...).
This doc backs up one step and explains what comes before that: the raw
reads a sequencer produces, in FASTQ format, and how they become that
matrix. It also shows you how to build a FASTQ file yourself.

## 1. The format itself

A FASTQ file is plain text (almost always gzipped). Every read is exactly
four lines:

```
@SYNTH:1:TOYFLOWCELL:1:1101:1000:2000 1:N:0:AAAAAAAA
CATGTGCGGCGACCCTGCGACAGTGACGAAAAAAAAAAAAAAA
+
GAHHFBFEB@@IEDEEF?III@?ED@BBBGFAEACFB@FG@?I
```

1. **`@` + header** -- a unique read identifier. Illumina's convention
   packs in `instrument:run:flowcell:lane:tile:x:y`, then (after the
   space) `read-number:filtered-flag:control-number:index-sequence`. You
   rarely need to parse this by hand, but two things matter for 10x data:
   the coordinates before the space are how a downstream tool recognizes
   that an R1 read and an R2 read came from the *same physical cluster*
   on the flow cell (i.e. the same original molecule), and the index
   sequence after `0:` is how samples pooled on one flow cell get sorted
   back into their own files.
2. **sequence line** -- the actual bases called, A/C/G/T (and occasionally
   `N` for "couldn't call this base").
3. **`+`** -- a separator, optionally repeating the header. Always just
   `+` in modern files.
4. **quality line** -- one character per base in the sequence line,
   encoding how confident the base caller was.

### Quality scores (Phred+33)

Each quality character is `chr(Q + 33)`, where `Q` is a Phred score:

| Character | ASCII code | Q  | Meaning                          |
|-----------|-----------|----|-----------------------------------|
| `?`       | 63        | 30 | 1-in-1,000 chance the base is wrong |
| `I`       | 73        | 40 | 1-in-10,000 chance the base is wrong |

Higher Q = more confident. `Q = -10 * log10(P_error)`, so every +10 is a
10x drop in error probability. Real reads typically start high-quality
and drift lower toward the end of the read -- this toy example doesn't
bother simulating that decay, to keep it easy to read by eye.

## 2. Why there are *two* files per sample (R1 / R2)

10x's 3' chemistry (Chromium 3', and Visium HD's 3' assay) sequences each
fragment from both ends into two separate, paired files:

- **R1** -- the spatial barcode + UMI. This is bookkeeping: which bin on
  the slide, and which original molecule (so PCR duplicates of the same
  molecule can be collapsed later). It contains no biology by itself.
- **R2** -- the actual cDNA insert: a fragment of a real transcript. This
  is the biology.

Space Ranger reads R1 and R2 *together*, read-pair by read-pair (matched
by their identical header up to the read-number field), to figure out
"this fragment of this gene came from this position on the slide."

For Visium HD specifically, R1 is commonly cited as 43bp total (UMI +
spatial barcode combined) -- that's the figure this repo's toy generator
uses. Treat it as illustrative rather than gospel: for the authoritative,
current split, check 10x's own [Visium HD documentation](https://www.10xgenomics.com/support/spatial-gene-expression-hd) before relying
on it for anything beyond this example. A real spatial barcode also isn't
a random sequence -- it's drawn from a fixed whitelist where every valid
sequence maps to a known physical (x, y) position on the slide, which is
exactly the lookup table Space Ranger uses to place each read.

## 3. Make your own FASTQ files

[`scripts/make_toy_fastq.py`](../scripts/make_toy_fastq.py) writes a tiny,
fully synthetic paired-end FASTQ dataset shaped like the above -- 8 read
pairs, correctly paired, correctly encoded, gzipped, nothing but the
Python standard library. Run it and look at what it made:

```bash
python3 scripts/make_toy_fastq.py
zcat data/toy_fastq/toy_S1_L001_R1_001.fastq.gz | head -8
zcat data/toy_fastq/toy_S1_L001_R2_001.fastq.gz | head -8
```

That's genuinely it -- a FASTQ file is nothing more than those four lines
repeated per read, written to a `.gz`. The only "trick" is getting the
header pairing and the quality-string length (must exactly match the
sequence length) right, which the script handles.

This toy dataset is **not** real sequencing data and running it through
Space Ranger would not produce a biologically meaningful result (the
sequences and barcodes are random). It exists purely so you have working,
correct code for the mechanics of the format.

### Reading FASTQ back (R side)

Since the rest of this repo is R/Seurat, here's the equivalent read-side,
using Bioconductor's `ShortRead`:

```r
# install.packages("BiocManager"); BiocManager::install("ShortRead")
library(ShortRead)
fq <- readFastq("data/toy_fastq/toy_S1_L001_R1_001.fastq.gz")
sread(fq)              # the sequences, as a DNAStringSet
quality(fq)             # the quality strings
width(fq)                # read lengths -- should all be 43 for R1
```

## 4. Where real Visium HD FASTQs come from, and what's between them and our matrix

The full chain, start to finish:

```
Sequencer's raw signal
        │  base calling
        ▼
BCL files  (Illumina's native per-cycle format)
        │  bcl2fastq / bcl-convert  (or `spaceranger mkfastq`, a thin wrapper)
        ▼
FASTQ files  (R1 + R2 per sample, per lane)  <-- this doc stops here
        │  `spaceranger count`:
        │    1. correct each R1 barcode against the whitelist, correct/collapse UMIs
        │    2. align R2 to the reference transcriptome (STAR under the hood)
        │    3. assign each aligned read to a gene AND a spatial bin
        │    4. sum UMIs per (gene, bin) into a sparse count matrix
        ▼
filtered_feature_bc_matrix.h5 + tissue_positions.parquet + ...  <-- 02_visium_hd_analysis.Rmd starts here
```

Everything in this repo's actual analysis picks up at the last box. That's
deliberate and standard practice -- essentially nobody re-runs the
alignment step for an analysis project; you consume the processed matrix
GEO/10x already published, exactly as we're doing with GSM8594569 /
GSM8594571.

### If you want to see (or run) the real thing

10x publishes a small, real, official Visium HD FASTQ example specifically
for this purpose -- the **"Visium HD Tiny 3' Dataset"** (mouse brain, ~297
MB, downsampled from a real run):
<https://www.10xgenomics.com/support/software/space-ranger/latest/resources/visium-hd-example-data>

Two honest caveats if you want to actually run `spaceranger count` on it
yourself, rather than just open the FASTQs and look:

1. **Space Ranger itself is free but gated** -- 10x requires a (free)
   registration/license click-through before you can download the
   binary, and it isn't on any package registry, so it can't be
   auto-installed as part of this repo's setup.
2. **Compute** -- 10x's own docs recommend well north of what a laptop
   typically has (their full-size Visium HD guidance is in the tens of
   GB of RAM; the "Tiny" dataset above is deliberately downsampled to be
   far more forgiving, but it's still real alignment work, not instant).

Neither of those is a reason not to try it if you're curious -- just
scope it as its own side exercise rather than a dependency of the tumor
vs. normal-adjacent analysis in this repo.
