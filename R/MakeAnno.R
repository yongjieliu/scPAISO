#' Make gene annotation GRanges From gtf
#' 
#' 
#' @param gtf Input gtf file
#' 
#' @export

MakeAnnoFromGtf = function(gtf){

  if (!exists("calc_all_introns", mode = "function")) {
    cpppath <- system.file("calc_introns.cpp", package = "scPAISO")
    if (cpppath == "") stop("calc_introns.cpp not found inside scPAISO installation", call. = FALSE)
    Rcpp::sourceCpp(cpppath)
  }

  
  if (missing(gtf) || length(gtf) != 1) stop("gtf path must be a single string")
  if (!file.exists(gtf)) stop(paste0("gtf file not found: ", gtf))
  
  message("Loading gtf")
  gtf.gr <- rtracklayer::import(gtf)
  gtf.gr <- gtf.gr[,c("type","gene_id","gene_name","transcript_id")]
  
  ######## gene ########
  message("Making gene annotation")
  gene.gene <- gtf.gr[gtf.gr$type %in% c("gene")]
  gene.gene <- gene.gene[,c("gene_id","gene_name")]
  
  ######## exon ########
  message("Making last exon annotation")
  last.exon <-  gtf.gr[gtf.gr$type %in% "exon"]
  last.exon.pos <- sort(last.exon[strand(last.exon) == "+",],decreasing = T)
  last.exon.pos <- last.exon.pos[!duplicated(last.exon.pos$transcript_id),]
  last.exon.neg <- sort(last.exon[strand(last.exon) == "-",],decreasing = F)
  last.exon.neg <- last.exon.neg[!duplicated(last.exon.neg$transcript_id),]
  last.exon <-  c(last.exon.pos,last.exon.neg)
  last.exon <- reduce(split(last.exon, last.exon$gene_id), min.gapwidth = 5)

  last.exon <- stack(last.exon)
  last.exon <- sort(last.exon)
  colnames(last.exon@elementMetadata) = "gene_id"
  last.exon$gene_id = as.character(last.exon$gene_id)

  message("Making other exon annotation")
  other.exon <-  gtf.gr[gtf.gr$type %in% "exon"]
  other.exon <- reduce(split(other.exon, other.exon$gene_id), min.gapwidth = 5)
  other.exon <- stack(other.exon)
  other.exon <- sort(other.exon)
  colnames(other.exon@elementMetadata) = "gene_id"
  other.exon$gene_id = as.character(other.exon$gene_id)

  ######## intron ########
  message("Making intron annotation")
  non.intron <- gtf.gr[gtf.gr$type %ni% c("gene","transcript")]
  non.intron <- reduce(split(non.intron, non.intron$transcript_id), min.gapwidth = 5)
  non.intron <- stack(non.intron)
  colnames(non.intron@elementMetadata) = "transcript_id"

  gene.intron <- calc_all_introns(
    start = start(non.intron),
    end = end(non.intron),
    seqnames = as.character(seqnames(non.intron)),
    strand = as.character(strand(non.intron)),
    transcript_id = as.character(non.intron$transcript_id)
  )

  gene.intron <- as(gene.intron,"GRanges")

  gene.intron$gene_id <- plyr::mapvalues(gene.intron$transcript_id,
                                         from = unique(data.frame(gene_id = gtf.gr$gene_id,transcript_id = gtf.gr$transcript_id))[,2],
                                         to = unique(data.frame(gene_id = gtf.gr$gene_id,transcript_id = gtf.gr$transcript_id))[,1],
                                         warn_missing = F)
  
  last.intro.pos <- sort(gene.intron[strand(gene.intron) == "+",],decreasing = T)
  last.intro.pos <- last.intro.pos[!duplicated(last.intro.pos$transcript_id),]
  last.intro.neg <- sort(gene.intron[strand(gene.intron) == "-",],decreasing = F)
  last.intro.neg <- last.intro.neg[!duplicated(last.intro.neg$transcript_id),]
  last.intro <-  c(last.intro.pos,last.intro.neg)

  last.intro <- reduce(split(last.intro, last.intro$gene_id), min.gapwidth = 5)
  last.intro <- stack(last.intro)
  last.intro <- sort(last.intro)
  colnames(last.intro@elementMetadata) = "gene_id"
  last.intro$gene_id = as.character(last.intro$gene_id)
  
  gene.intron <- disjoin(split(gene.intron, gene.intron$gene_id))
  gene.intron <- stack(gene.intron)
  colnames(gene.intron@elementMetadata) = "gene_id"
  gene.intron$gene_id = as.character(gene.intron$gene_id)
  gene.intron$intron_width = end(gene.intron) - start(gene.intron) + 1
  gene.intron <- sort(gene.intron)
  gene.intron$gene_name = plyr::mapvalues(gene.intron$gene_id,from = gene.gene$gene_id,to = gene.gene$gene_name,warn_missing = F)
  
  ######## utr ########
  message("Making UTR annotation")
  utr = gtf.gr[grepl("UTR|utr", gtf.gr$type)]
  tmp.index = findOverlaps(utr,last.exon, maxgap=5)
  tmp.index = tmp.index[utr$gene_id[queryHits(tmp.index)] == last.exon$gene_id[subjectHits(tmp.index)]]
  
  utr5 <- utr[1:length(utr) %ni% queryHits(tmp.index)]
  utr5 <- reduce(split(utr5, utr5$gene_id), min.gapwidth = 5)
  utr5 <- stack(utr5)
  utr5 <- sort(utr5)
  colnames(utr5@elementMetadata) = "gene_id"
  utr5$gene_id = as.character(utr5$gene_id)
  
  utr3 <- utr[queryHits(tmp.index)]
  utr3 <- reduce(split(utr3, utr3$gene_id), min.gapwidth = 5)
  utr3 <- stack(utr3)
  colnames(utr3@elementMetadata) = "gene_id"
  utr3 <- sort(utr3)
  utr3$gene_id = as.character(utr3$gene_id)

  utr3.label <- c(utr3,last.exon)
  utr3.label <- reduce(split(utr3.label, utr3.label$gene_id), min.gapwidth = 5)
  utr3.label <- stack(utr3.label)
  colnames(utr3.label@elementMetadata) = "gene_id"

  utr3.label$region <- paste0(seqnames(utr3.label),":",start(utr3.label),"-",end(utr3.label))
  tmp.index = findOverlaps(utr3,utr3.label)
  tmp.index = tmp.index[utr3$gene_id[queryHits(tmp.index)] == utr3.label$gene_id[subjectHits(tmp.index)]]
  utr3$region <- paste0(seqnames(utr3),":",start(utr3),"-",end(utr3))
  utr3$region[queryHits(tmp.index)] <- utr3.label$region[subjectHits(tmp.index)]
  utr3 <- sort(utr3)
  
  message("Exporting")
  
  tmp.index = findOverlaps(last.exon,utr3.label)
  tmp.index = tmp.index[last.exon$gene_id[queryHits(tmp.index)] == utr3.label$gene_id[subjectHits(tmp.index)]]
  last.exon$region <- paste0(seqnames(last.exon),":",start(last.exon),"-",end(last.exon))
  last.exon$region[queryHits(tmp.index)] <- utr3.label$region[subjectHits(tmp.index)]
  last.exon <- sort(last.exon)
  
  last.exon[paste0(seqnames(last.exon),":",start(last.exon),"-",end(last.exon)) != last.exon$region]
  
  gtf.anno <- list(gene = gene.gene,
                    utr3 = utr3,
                    last.exon = last.exon,exon = other.exon,
                    last.intron = last.intro,intron = gene.intron,
                    utr5 = utr5)
  return(gtf.anno)
}
