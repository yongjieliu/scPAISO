#' scPAISO package
#'
#' Utilities for identifying and quantifying alternative polyadenylation sites
#' using single-cell RNA-seq data.
#'
#' @docType package
#' @name scPAISO
#' @import future.apply
#' @import Matrix
#' @import GenomicAlignments
#' @import GenomicRanges
#' @import IRanges
#' @import Rsamtools
#' @import BSgenome
#' @import ggplot2
#' @import grid
#' @import cowplot
#' @import rtracklayer
#' @import scales
#' @import methods
#' @importFrom data.table ":=" data.table setDT setnames
#' @importFrom data.table as.data.table dcast melt rbindlist
#' @importFrom dplyr arrange filter group_by inner_join mutate summarise
#' @importFrom dplyr row_number
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
#' @importFrom stats aggregate coef dgamma dist dlnorm dweibull predict
#' @importFrom utils write.csv
#' @importFrom grDevices dev.off pdf png
#' @importFrom future plan multisession sequential
#' @importFrom magrittr %>%
#' @importFrom fitdistrplus fitdist denscomp
#' @importFrom Biostrings DNAStringSet letterFrequency
#' @importFrom S4Vectors mcols queryHits subjectHits
#' @keywords internal
"_PACKAGE"
