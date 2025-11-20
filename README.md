# Pipeline

## Mapping

``` bash
# for example, you have two scRNA-seq library
## cellranger Read1 + Read2
myindex=S1 S2
for sp in ${myindex[@]}
do
    nice cellranger count --id=${sp}Count \
                 --transcriptome={DIRECTORY_WITH_REFERENCE} \
                 --localcores={NUMBER_OF_CPUS} \
                 --include-introns true \ # both true or false is ok
                 --sample=${sp} \
                 --fastqs={DIRECTORY_WITH_FASTQ_FILES} 
done

## STAR Read1
for sp in ${myindex[@]} 
do
    nice STAR --runThreadN {NUMBER_OF_CPUS} \
        --genomeDir {DIRECTORY_WITH_REFERENCE} \
        --readFilesIn ${sp}_S1_L001_R1_001.fastq.gz \
        --readFilesCommand zcat \
        --clip5pNbases 58 \
        --limitBAMsortRAM 200323256559 \ # depends on your fastq size
        --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.3 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 70 \
        --sjdbGTFfile $mytranscript/genes/genes.gtf \
        --outFileNamePrefix ./$sp.R1 \
        --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate
done
```

## APA analysis

### Setup
``` r
library(Rsamtools)
library(dplyr)
library(data.table)
library(Matrix)
library(GenomicAlignments)
library(future.apply)
library(rtracklayer)
library(scales)
#library(plyr)
library(reshape2)

library(fitdistrplus)

library(ggplot2)
library(grid)
library(RColorBrewer)
library(cowplot)
theme_set(theme_cowplot(15))
library(ComplexHeatmap)

library('BSgenome.Hsapiens.UCSC.hg38')
library('TxDb.Hsapiens.UCSC.hg38.knownGene')
library('org.Hs.eg.db')
'%ni%' <- Negate('%in%')

library(scPAISO)
outpath= "./scPAISO/tests/pe/"
gene.anno  <- readRDS("./scPAISO/tests/hg38/annotation.gene.hg38.rds")
geneanno <- gene.anno$gene
intronanno <-  gene.anno$intron
chrs <- paste0("chr",21:22)
chrom.size <- read.table("./scPAISO/tests/hg38/chrNameLength.txt")
colnames(chrom.size) <- c("seqnames","seqlengths")
rownames(chrom.size) <- chrom.size$seqnames
threads = 4

bam.r2 = c("./scPAISO/tests/pe/ctl4_R2_chr21_22.bam", "./scPAISO/tests/pe/itp4_R2_chr21_22.bam") # subset from cellranger output
bam.r1 = c("./scPAISO/tests/pe/ctl4_R1_chr21_22.bam", "./scPAISO/tests/pe/itp4_R1_chr21_22.bam") # subset from from STAR

samplenames <- c("ctl4","itp4")

# make a data.frame cell_meta
## cb = sample_cellbarcode
## for example, itp4_AAACCCAAGAGCCTGA-1
cell_meta <- readRDS("./scPAISO/tests/pe/cell_meta.Rds")
# cell_meta <- data.frame(cb = colnames(SeuratObj),celltype = SeuratObj$celltype,samplename = sc.SeuratObj$samplename)
```

### Process bams

``` r
# To achieve efficient parallel processing and reduce memory usage, we split the BAM file by chromatin.
BamSplit("./scPAISO/inst/scripts/bam.split.sh",bam.r1,bam.r2,samplenames,outpath)
# bam to data.frame
ProcessBam(samplenames,chrs,outpath,geneanno,intronanno,cell.metadata = cell_meta,chrom.size, cache.bam2df= F,threads = threads)
```

### PAS peaks identification

``` r
# call peak
pas_peak <- findPAS(samplenames,chrs,outpath,threads,cache.covdata = F,
                    callpeak = "~/anaconda3/bin/macs3",genomeSize = 2.7e9, cutOff = 10^-5,additionalParams = "--nomodel",
                    sum.count = 10,cpm.cutoff = 0.5,
                    min.gapwidth = 60,
                    BSg = BSgenome.Hsapiens.UCSC.hg38,literal = "AAAAAAAA",num_A = 10) 

# mapping PAS peaks to genes
pas_peak <- pas_peak2gene(gr_pas = pas_peak,samplenames = samplenames, chrs = chrs,outpath = outpath,threads = threads)

# annotate PAS peaks: 3' UTR, 3' UTR expand, Last Exon, Other Exon, Last Intron, and Other Intron.
pas_peak <- pas_peak_anno(pas_peak,gene.anno,outpath = outpath)

# plot base frequence of surrounding PAS peak
plot_base_freq(gr_pas = pas_peak,extend_base = 100,BSg = BSgenome.Hsapiens.UCSC.hg38,outpath = outpath,outname = "pas.100bp.freq.pdf")
```

### PA isoform quantification

``` r
# fitting model: weibull, gamma, or log-normal distribution
pas_fit_model(samplenames,chrs,outpath,gr_pas = pas_peak,intronanno = intronanno,cache.dist = F,
              intron.remove = T,threads=4,sampling.frac = 4/100)

# PA isoform quantification: PAS peak X cell matrix 
pas_sc_count <- pas_predict(gr_pas = pas_peak,samplenames,chrs,outpath,
                            geneanno = geneanno,intronanno = intronanno,intron.remove = T,
                            cell.metadata = cell_meta,fitmodel = "weibull",threads=threads)
```

### APA rank score

``` r
pas.rank.score <- PasUsageScore(pas_sc_count,gr_pas = pas_peak,min.pas.cutoff = 3,
                                  pseudo_count = 0,threads = threads,by = "utr",out.matrix = F)
```
