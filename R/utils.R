## Internal utility operators and helpers

`%ni%` <- function(x, table) {
  !(x %in% table)
}

default_hg38_bsgenome <- function() {
  pkg <- "BSgenome.Hsapiens.UCSC.hg38"
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package ", pkg,
      " must be installed or supply a custom BSg argument.",
      call. = FALSE
    )
  }
  getExportedValue(pkg, pkg)
}

utils::globalVariables(c(
  "..density..", ".", ".N", ":=",
  "AggregateExpression", "AnnotationPlot", "BigwigTrack",
  "CombineTracks", "Count", "DefaultAssay<-", "Description", "FeaturePlot",
  "Freq", "Gene.nCount", "NoLegend", "PAS", "PAS.nCount", "PeakPlot",
  "R2PAS", "Var1", "VlnPlot", "WhichMaxs", "big.path", "celltype",
  "denscomp", "end_r2", "fitdist", "gene_id", "gene_sum", "index_df",
  "intron_width", "map.peak", "mapL_r1", "mapL_r2", "mapq_r1", "mapq_r2",
  "other", "outpath", "overlap.width", "pas_peak", "pas_sc_count", "peak",
  "plotvalue", "pos_r1", "pos_r2", "qname", "queryHits", "rank2", "ratio",
  "readid", "sc.merged", "sc.pas.merged", "seqlevels", "seqlevels<-",
  "strand_r1", "strand_r2", "subjectHits", "tag.CB", "tag.UB", "value",
  "variable", "x", "count"
))
