# scPAISO

scPAISO is an end-to-end APA analysis workflow for 3′ tag-based scRNA-seq. The package covers BAM preprocessing, PAS discovery, isoform quantification, and visualization helpers. This document provides instructions for installing the package and quickly running the included demonstration.

## Requirements

- R ≥ 4.2; using `renv`, Conda, or containerized R is recommended for reproducibility.
- Mandatory packages: `Rsamtools`, `data.table`, `dplyr`, `GenomicAlignments`, `GenomicRanges`, `IRanges`, `Matrix`, `future.apply`, `ggplot2`, `rtracklayer`, `BSgenome`, `fitdistrplus`, etc. See `DESCRIPTION` for the complete list.
- Optional packages:
  - `BSgenome.Hsapiens.UCSC.hg38` (default genome; install different `BSgenome` objects for other species).
  - `nabor` (KNN-based imputation in `pas_impute`).

To install the development version, please use:

```r
remotes::install_github("yongjieliu/scPAISO")
```

Or to install a specific version:

```r
remotes::install_url("https://github.com/yongjieliu/scPAISO/archive/refs/tags/v0.1.0.tar.gz")
```

## Demo and Test Data

Lightweight demo data (chr21/22 subsets, annotations, intermediate results) are hosted on Zenodo (https://zenodo.org/records/17656388)


## Quick Start (chr21/22 subset)

Constructing genome annotation files:

```r
# Make your own gene annotation :
mygtf= "refdata-gex-GRCh38-2024-A/genes/genes.gtf"
gene.anno = MakeAnnoFromGtf(mygtf)
# Or use the file provided in the demo data.:
gene.anno <- readRDS("./scPAISO.test/annotation.gene.hg38.rds")
```

Define the essential parameters:

```r
library(scPAISO)
library(Rsamtools)
library(dplyr)
library(data.table)
library(Matrix)
library(GenomicAlignments)
library(future.apply)
library(rtracklayer)
library(scales)

library(plyranges)

library(fitdistrplus)

library(ggplot2)
library(grid)
library(RColorBrewer)
library(cowplot)
theme_set(theme_cowplot(15))
library(ComplexHeatmap)

library('BSgenome.Hsapiens.UCSC.hg38')

outpath <- "./scPAISO.test"
dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

geneanno <- gene.anno$gene
intronanno <- gene.anno$intron
chrs <- paste0("chr", 21:22)
chrom.size <- read.table("./scPAISO.test/hg38.chrNameLength.txt",
                         col.names = c("seqnames", "seqlengths"))
threads <- 8

bam.r2 <- c("./scPAISO.test/ctl4_R2_chr21_22.bam",
            "./scPAISO.test/itp4_R2_chr21_22.bam")
bam.r1 <- c("./scPAISO.test/ctl4_R1_chr21_22.bam",
            "./scPAISO.test/itp4_R1_chr21_22.bam")
samplenames <- c("ctl4", "itp4")
cell_meta <- readRDS("./scPAISO.test/cell_meta.Rds")
```

Altemative Polyadenylation Analysis:

```r
# splits Read1/Read2 BAM pairs by chromosome to reduce memory use
BamSplit(bam.r1, bam.r2, samplenames, outpath)

# loads split BAMs, annotates reads with gene/intron metadata, and 
# caches intermediate data frames when `cache.bam2df = TRUE`.
ProcessBam(samplenames, chrs, outpath, geneanno, intronanno,
           cell.metadata = cell_meta, chrom.size = chrom.size,
           cache.bam2df = FALSE, threads = threads)

# builds coverage and calls PAS peaks via MACS3; 
# Ensure `macs3` is available in `$PATH` or pass its absolute path via `callpeak`.
# `literal/num_A` control the poly(A) motif filter, while `min.gapwidth` merges nearby peaks.

pas_peak <- findPAS(samplenames, chrs, outpath, threads,
                    callpeak = "macs3", genomeSize = 2.7e9,
                    min.gapwidth = 60, literal = "AAAAAAAA",
                    num_A = 10)

# maps each PAS peak to genes across chromosomes
pas_peak <- pas_peak2gene(pas_peak, samplenames, chrs, outpath, threads)

# labels regions (3′ UTR, exon, intron) using the annotation list generated via `MakeAnnoFromGtf`
pas_peak <- pas_peak_anno(pas_peak, gene.anno, outpath)

# estimates read-length distributions per sample; reducing `sampling.frac` speeds modeling.
pas_fit_model(samplenames, chrs, outpath, gr_pas = pas_peak,
              intronanno = intronanno, intron.remove = TRUE,
              cache.dist = FALSE, sampling.frac = 0.04,
              threads = threads)

# builds the PAS-by-cell matrix using the fitted model.
pas_sc_count <- pas_predict(pas_peak, samplenames, chrs, outpath,
                            geneanno = geneanno, intronanno = intronanno,
                            cell.metadata = cell_meta,
                            intron.remove = TRUE, threads = threads)

# computes per-cell PAS usage scores;
# To obtain the average APA score per cell, set `out.matrix = False` (the default).
# Conversely, setting it to TRUE returns  the APA rank score of each gene in each cell (consumes more computational time and memory).
PasUsageScore(pas_sc_count, gr_pas = pas_peak,
              min.pas.cutoff = 3, by = "utr",out.matrix = F,
              threads = threads)
```

## Usage Notes

- Always supply an explicit `outpath`; the pipeline writes every stage (BAM chunks, PAS calls, model fits, logs) under that directory.
- `pas_predict` dominate runtime. Increasing `threads` shortens runtime but raises memory usage almost linearly.
- Set species-specific genome objects via the `BSg` argument (e.g., `BSg = BSgenome.Mmusculus.UCSC.mm39`).
- Optional features (`pas_impute`, `feature.plot`, `llrplot`) throw explicit errors if `nabor`, `Seurat`, or `Signac` are missing; install them only when needed.

## Troubleshooting

- **Missing gene/transcript columns**: `MakeAnnoFromGtf` tries multiple column names. If it fails, inspect the GTF metadata columns (`colnames(mcols(import(gtf)))`).
- **`macs3` not found**: confirm `which macs3` returns a path or point `callpeak` to the executable.
- **Test data path mismatches**: dump the directory layout (`list.files(outpath, recursive = FALSE)`) and ensure the R session sees the same path (use absolute paths when in doubt).
- **Concurrency issues**: if you hit random failures, rerun with `threads = 1` to verify correctness, then scale up.

When opening issues, include the exact commands, input files, `sessionInfo()`, and relevant logs to keep the workflow reproducible.
