#' Make gene annotation GRanges From gtf
#' 
#' 
#' @param gtf Input gtf file
#' 
#' @export

MakeAnnoFromGtf = function(gtf){
  if (missing(gtf) || length(gtf) != 1) stop("gtf path must be a single string")
  if (!file.exists(gtf)) stop(paste0("gtf file not found: ", gtf))
  
  print("Loading gtf")
  gtf.gr <- rtracklayer::import(gtf)
  gtf.gr <- gtf.gr[,c("type","gene_id","gene_name","transcript_id","exon_number")]
  
  ######## gene ########
  print("Making gene annotation")
  gene.gene <- gtf.gr[gtf.gr$type %in% c("gene")]
  gene.gene <- gene.gene[,c("gene_id","gene_name")]
  
  ######## exon ########
  print("Making last exon annotation")
  last.exon <-  gtf.gr[gtf.gr$type %in% "exon"]
  last.exon.pos <- sort(last.exon[strand(last.exon) == "+",],decreasing = T)
  last.exon.pos <- last.exon.pos[!duplicated(last.exon.pos$transcript_id),]
  last.exon.neg <- sort(last.exon[strand(last.exon) == "-",],decreasing = F)
  last.exon.neg <- last.exon.neg[!duplicated(last.exon.neg$transcript_id),]
  last.exon <-  c(last.exon.pos,last.exon.neg)
  last.exon <- reduce(split(last.exon, last.exon$gene_id), min.gapwidth = 5)
  tmp <- future_lapply(names(last.exon), function(x) {
    gr = last.exon[[x]]
    gr$gene_id <- x
    return(gr)
  })
  last.exon <- do.call("c", tmp)
  last.exon <- sort(last.exon)
  
  print("Making other exon annotation")
  other.exon <-  gtf.gr[gtf.gr$type %in% "exon"]
  # other.exon <- reduce(split(other.exon, other.exon$gene_id), min.gapwidth = 5)
  # tmp <- future_lapply(names(other.exon), function(x) {
  #   gr = other.exon[[x]]
  #   gr$gene_id <- x
  #   return(gr)
  # }) 
  # other.exon <- do.call("c", tmp)
  other.exon <- sort(other.exon)
  
  ######## intron ########
  print("Making intron annotation")
  non.intron <- gtf.gr[gtf.gr$type %ni% c("gene","transcript")]
  non.intron <- reduce(split(non.intron, non.intron$transcript_id), min.gapwidth = 5)
  
  gene.intron <- future_lapply(names(non.intron), function(x) {
    gr = non.intron[[x]]
    if (length(gr) > 1){
      gr = GRanges(seqnames=seqnames(gr)[1], 
                   ranges=IRanges(start= end(gr)[1:(length(gr)-1)] +1,
                                  end=start(gr)[2:length(gr)] -1), 
                   strand=strand(gr)[1])
      gr$transcript_id = x
    } else {
      gr$transcript_id = "exon"
    }
    return(gr)
  })
  rm(non.intron)
  gene.intron <- do.call("c", gene.intron)
  gene.intron <- gene.intron[gene.intron$transcript_id != "exon"]
  gene.intron$gene_id <- plyr::mapvalues(gene.intron$transcript_id,
                                         from = unique(data.frame(gene_id = gtf.gr$gene_id,transcript_id = gtf.gr$transcript_id))[,2],
                                         to = unique(data.frame(gene_id = gtf.gr$gene_id,transcript_id = gtf.gr$transcript_id))[,1],
                                         warn_missing = F)
  
  last.intro.pos <- sort(gene.intron[strand(gene.intron) == "+",],decreasing = T)
  last.intro.pos <- last.intro.pos[!duplicated(last.intro.pos$transcript_id),]
  last.intro.neg <- sort(gene.intron[strand(gene.intron) == "-",],decreasing = F)
  last.intro.neg <- last.intro.neg[!duplicated(last.intro.neg$transcript_id),]
  last.intro <-  c(last.intro.pos,last.intro.neg)
  
  # last.intro <- reduce(split(last.intro, last.intro$gene_id), min.gapwidth = 5)
  # tmp <- future_lapply(names(last.intro), function(x) {
  #   gr = last.intro[[x]]
  #   gr$gene_id <- x
  #   return(gr)
  # }) 
  # last.intro <- do.call("c", tmp)
  last.intro <- sort(last.intro)
  
  gene.intron <- disjoin(split(gene.intron, gene.intron$gene_id))
  tmp <- future_lapply(names(gene.intron), function(x) {
    gr = gene.intron[[x]]
    gr$gene_id <- x
    return(gr)
  }) 
  gene.intron <- do.call("c", tmp)
  gene.intron$intron_width = end(gene.intron) - start(gene.intron) + 1
  gene.intron <- sort(gene.intron)
  
  gene.intron$gene_name = plyr::mapvalues(gene.intron$gene_id,from = gene.gene$gene_id,to = gene.gene$gene_name,warn_missing = F)
  
  ######## utr ########
  print("Making UTR annotation")
  utr = gtf.gr[grepl("UTR|utr", gtf.gr$type)]
  tmp.index = findOverlaps(utr,last.exon, maxgap=5)
  tmp.index = tmp.index[utr$gene_id[queryHits(tmp.index)] == last.exon$gene_id[subjectHits(tmp.index)]]
  
  utr5 <- utr[1:length(utr) %ni% queryHits(tmp.index)]
  # utr5 <- reduce(split(utr5, utr5$gene_id), min.gapwidth = 5)
  # tmp <- future_lapply(names(utr5), function(x) {
  #   gr = utr5[[x]]
  #   gr$gene_id <- x
  #   return(gr)
  # })
  # utr5 <- do.call("c", tmp)
  utr5 <- sort(utr5)
  
  utr3 <- utr[queryHits(tmp.index)]
  utr3 <- reduce(split(utr3, utr3$gene_id), min.gapwidth = 5)
  tmp <- future_lapply(names(utr3), function(x) {
    gr = utr3[[x]]
    gr$gene_id <- x
    return(gr)
  }) 
  utr3 <- do.call("c", tmp)
  utr3 <- sort(utr3)
  
  utr3.label <- c(utr3,last.exon)
  utr3.label <- reduce(split(utr3.label, utr3.label$gene_id), min.gapwidth = 5)
  utr3.label <- future_lapply(names(utr3.label), function(x) {
    gr = utr3.label[[x]]
    gr$gene_id <- x
    return(gr)
  })
  utr3.label <- do.call("c", utr3.label)
  utr3.label$region <- paste0(seqnames(utr3.label),":",start(utr3.label),"-",end(utr3.label))
  tmp.index = findOverlaps(utr3,utr3.label)
  tmp.index = tmp.index[utr3$gene_id[queryHits(tmp.index)] == utr3.label$gene_id[subjectHits(tmp.index)]]
  utr3$region <- paste0(seqnames(utr3),":",start(utr3),"-",end(utr3))
  utr3$region[queryHits(tmp.index)] <- utr3.label$region[subjectHits(tmp.index)]
  utr3 <- sort(utr3)
  
  print("Exporting")
  
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
