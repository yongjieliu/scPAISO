#' Process bam files
#' 
#' 
#' @param bam.split.sh Path to the bam.split.sh file.
#' @param bam.r1 A vector containing bam file path of Read1.
#' @param bam.r2 A vector containing bam file path of Read2.
#' @param samplenames A vector containing sample name of each bam file.
#' @param outpath A character describing where the file will be saved.
#'
#' @param chrs A vector containing chromatin to use. 
#' @param geneanno A ranges describing gene annotation. 
#' Options include "underscore", "dash", "numeric" and "dot". The string portion prior to these will be kept.
#' @param intronanno A ranges describing intron annotation.
#' @param cell.metadata A data.frame containing: cb, cell barcode id; celltype, cell type annotation; samplename, sample name.
#' @param cache.bam2df A boolean value that determines whether cache bam2df file shoulde be used.
#'
#' @param chrom.size A data.frame containing: seqnames, chromatin id; seqlengths, chromatin length. seqnames as rowname of data.frame.
#' @param threads The number of threads to be used for parallel computing.
#' @param isbulk A boolean value that determines whether the RNA-seq library is bulk RNA-seq such as 3'-seq and Poly(A)-seq.
#' 
#' @export

BamSplit <- function(bam.split.sh = NULL,bam.r1 = NULL,bam.r2,samplenames,outpath = "./scPAISO.output"){

  if (is.null(bam.r2) || length(bam.r2) == 0) {
    stop("bam.r2 cannot be empty")
  }
  if (is.null(samplenames) || length(samplenames) == 0) {
    stop("samplenames cannot be empty")
  }
  if (length(unique(c(length(bam.r2), length(samplenames)))) != 1) {
    stop("bam.r2 and samplenames must have the same length")
  }
  if (!is.null(bam.r1) && length(unique(c(length(bam.r1), length(bam.r2), length(samplenames)))) != 1) {
    stop("bam.r1, bam.r2, and samplenames must have the same length")
  }
  if (any(!file.exists(c(bam.r2, bam.r1)))) {
    stop("One or more BAM paths do not exist; please check bam.r1/bam.r2")
  }

  if(is.null(bam.split.sh)) {
    bam.split.sh <- system.file("scripts/bam.split.sh", package = "scPAISO")
    if(!file.exists(bam.split.sh)) {
      stop("Script not found!")
    }
  }
  if (!file.exists(outpath)) {
    dir.create(outpath, recursive = TRUE)
  }
  
  if (is.null(bam.r1)){
    print("Warning: processing by Read 2 only model")
    if (length(unique(length(bam.r2),length(samplenames))) != 1){
      stop("Error: bam.r2 and samplenames should have the same size")
    }
    
    for (i in 1:length(bam.r2)){
      status <- system(paste0(bam.split.sh," ",bam.r2[i]," ",samplenames[i]," ",outpath), wait = TRUE)
      if (status != 0) stop("bam.split.sh failed for sample ", samplenames[i])
      out_dir <- file.path(outpath,"chr.bam.R2",samplenames[i])
      bam_files <- list.files(out_dir, pattern = "\\.bam$", full.names = TRUE)
      if (length(bam_files) == 0) {
        stop("bam.split.sh completed but produced no BAM files for sample ", samplenames[i])
      }
    }
    
  } else {
    print("Warning: processing by Read 1 + Read 2 model")
    if (length(unique(length(bam.r1),length(bam.r2),length(samplenames))) != 1){
      stop("Error: bam.r1, bam.r2, and samplenames should have the same size")
    }
    
    for (i in 1:length(bam.r2)){
      status <- system(paste0(bam.split.sh," ",bam.r1[i]," ",bam.r2[i]," ",samplenames[i]," ",outpath), wait = TRUE)
      if (status != 0) stop("bam.split.sh failed for sample ", samplenames[i])
      out_dir_r2 <- file.path(outpath,"chr.bam.R2",samplenames[i])
      out_dir_r1 <- file.path(outpath,"chr.bam.R1",samplenames[i])
      bam_r2_files <- list.files(out_dir_r2, pattern = "\\.bam$", full.names = TRUE)
      bam_r1_files <- list.files(out_dir_r1, pattern = "\\.bam$", full.names = TRUE)
      if (length(bam_r2_files) == 0 || length(bam_r1_files) == 0) {
        stop("bam.split.sh completed but produced no BAM files for sample ", samplenames[i])
      }
    }
  }
}

bam2df <- function(samplename,chr,outpath,geneanno,intronanno,cell.metadata,chrom.size,cache.bam2df){
  
  if (missing(samplename) || missing(chr) || missing(outpath)) {
    stop("samplename/chr/outpath cannot be missing")
  }
  if (!is.character(samplename) || length(samplename) != 1) stop("samplename must be a single character value")
  if (!is.character(chr) || length(chr) != 1) stop("chr must be a single character value")
  if (!dir.exists(outpath)) stop("outpath directory does not exist")
  if (!all(c("cb","celltype","samplename") %in% colnames(cell.metadata))) {
    stop("cell.metadata must contain columns cb, celltype, and samplename")
  }
  if (!all(c("seqnames","seqlengths") %in% colnames(chrom.size))) {
    stop("chrom.size must contain seqnames and seqlengths columns")
  }
  
  bam_r2 = paste0(outpath,"/chr.bam.R2/",samplename,"/",chr,".bam")
  bam_r1 = paste0(outpath,"/chr.bam.R1/",samplename,"/",chr,".bam")

  if (!file.exists(paste0(outpath,"/stat.qc/"))){
    dir.create(paste0(outpath,"/stat.qc/"))
  }
  
  if (!file.exists(paste0(outpath,"/bam2df/"))){
    dir.create(paste0(outpath,"/bam2df/"))
  }
  
  if (!file.exists(paste0(outpath,"/stat.qc/",samplename))){
    dir.create(paste0(outpath,"/stat.qc/",samplename))
  }
  
  if (!file.exists(paste0(outpath,"/bam2df/",samplename))){
    dir.create(paste0(outpath,"/bam2df/",samplename))
  }
  
  if(!file.exists(paste0(outpath,"/bam2df/",samplename,"/",chr,".bam2df.rds")) | !cache.bam2df){
    
    #### read bam
    df_r2 <- Rsamtools::scanBam(bam_r2, param = Rsamtools::ScanBamParam(what=c('qname',"strand","pos","mapq","cigar"),tag = c("CB","UB"))) #'flag',
    df_r2 <- data.table::as.data.table(do.call(rbind, lapply(df_r2, as.data.frame)))
    if (!"tag.CB" %in% names(df_r2)) stop("tag.CB not found in R2 BAM; ensure Cell Barcodes are tagged")
    df_r2 <- df_r2[!is.na(df_r2$tag.CB),]
    colnames(df_r2)[2:5] <- paste0(colnames(df_r2)[2:5],"_r2")
    df_r2$mapL_r2 <- GenomicAlignments::cigarWidthAlongReferenceSpace(df_r2$cigar_r2)
    
    valid_cb <- gsub(".+_","",cell.metadata$cb[cell.metadata$samplename == samplename])
    df_r2 <- df_r2[tag.CB %in% valid_cb]
  
    df_r2 <- df_r2[mapq_r2 >= 255]
    df_r2[, mapL_r2 := mapL_r2 + pos_r2 - 1]
    df_r2[, c("mapq_r2","cigar_r2") := NULL]
    colnames(df_r2)[6] <- "end_r2"
    
    qc.export <- paste0("R2 reads: ",nrow(df_r2),"\n")
    
    saveRDS(df_r2,file = paste0(outpath,"/bam2df/",samplename,"/",chr,".r2bam2df.rds"))
    
    df_r1 <- Rsamtools::scanBam(bam_r1, param = Rsamtools::ScanBamParam(what=c('qname',"strand","pos","mapq","cigar")))
    df_r1 <- data.table::as.data.table(do.call(rbind, lapply(df_r1, as.data.frame)))
    colnames(df_r1)[2:5] <- paste0(colnames(df_r1)[2:5],"_r1")
    df_r1$mapL_r1 <- GenomicAlignments::cigarWidthAlongReferenceSpace(df_r1$cigar_r1)
    df_r1 <- df_r1[mapq_r1 >= 255]
    df_r1[, mapL_r1 := mapL_r1 + pos_r1 - 1]
    df_r1[, c("mapq_r1","cigar_r1") := NULL]
    colnames(df_r1)[4] <- "end_r1"
    
    qc.export <- paste0(qc.export,
                        "R1 reads: ",nrow(df_r1),"\n")
    
    # df_r1$matched_lengths <- lapply(df_r1$cigar,function(x){
    #   sum(as.numeric(gsub("M","",unlist(regmatches(x, gregexpr("\\d+M", x)))))) 
    # }) %>% Reduce("c",.)
    
    #### merge R1 and R2 & strand filter
    df_result <- merge(df_r1,df_r2, allow.cartesian = FALSE)
    qc.export <- paste0(qc.export,
                        "Pair reads: ",nrow(df_result),"\n")
    
    data.table::setDT(df_result)
    df_result <- df_result[strand_r1 != strand_r2]
    
    qc.export <- paste0(qc.export,
                        "Strand match reads: ",nrow(df_result),"\n")
    
    #### R1/R2 out of range filter
    #       ===========> #R2
    #    <================    #R1 out of range
    df_result <- df_result[!(df_result$pos_r1 < df_result$pos_r2 & df_result$strand_r2 == "+"),]
    #       ===========> #R2 out of range
    #         <======    #R1
    df_result <- df_result[!(df_result$end_r1 < df_result$end_r2 & df_result$strand_r2 == "+"),]
  
    df_result <- df_result[!(df_result$pos_r1 > df_result$pos_r2 & df_result$strand_r2 == "-"),]
    df_result <- df_result[!(df_result$end_r1 > df_result$end_r2 & df_result$strand_r2 == "-"),]
    # dim(df_result) 258763
    
    qc.export <- paste0(qc.export,
                        "Nice mapped reads: ",nrow(df_result),"\n")
    
    #### gene anno
    tmp <- df_result[,.(chr = chr,
                      start = pos_r2,
                      end = end_r2,
                      strand = strand_r2)]
    
    gr.index = findOverlaps(as(tmp,"GRanges"),geneanno)
    
    overlap_length <- pintersect(as(tmp,"GRanges")[queryHits(gr.index)], geneanno[subjectHits(gr.index)])
    
    geneanno_df <- data.frame(geneanno)[,c("seqnames","gene_id","gene_name")]
    df_result <- cbind(df_result[queryHits(gr.index),],
                       geneanno_df[subjectHits(gr.index),])
    
    qc.export <- paste0(qc.export,
                        "Gene mapped reads: ",nrow(df_result),"\n")
    qc.export <- paste0(qc.export,
                        "Gene dup reads: ",table(duplicated(df_result$qname))[2],"\n")
    
    dup_q <- df_result$qname[duplicated(df_result$qname)]
    df_result <- df_result[!qname %in% dup_q]
    
    #### R2 to PAS distance filter
    
    df_result[, start := end_r2]
    df_result[strand_r2 == "-", start := pos_r1]
    df_result$end <- df_result$end_r1
    df_result$end[df_result$strand_r2 == "-"] <- df_result$pos_r2[df_result$strand_r2 == "-"]
    df_result$strand = df_result$strand_r2
    
    df_result <- df_result[df_result$end - df_result$start < 10^6,]
    
    gr.index <- findOverlaps(intronanno, type=c("within"),as(df_result,"GRanges"))
    #gr.index <- findOverlaps(intronanno, type=c("within"),as(df_result[df_result$qname == "A00358:379:HLWFFDSXY:4:1149:3848:22874",],"GRanges"))
    gr.index = gr.index[intronanno[from(gr.index)]$gene_id == df_result$gene_id[to(gr.index)],]
    
    intro_df <- data.frame(index_df = to(gr.index),
                           intron_width = intronanno[from(gr.index)]$intron_width)
    intro_df <- intro_df %>% group_by(index_df) %>% summarise(intron_total = sum(intron_width))
    
    df_result$intron <- 0
    df_result$intron[intro_df$index_df] <- intro_df$intron_total
  
    df_result$R2PAS <- df_result$end - df_result$start - df_result$intron
    
    # table(df_result$gene_name[df_result$R2PAS > 5000])
    # data.frame(df_result[df_result$R2PAS > 5000 & df_result$gene_name == "HMGN1",])[1:4,]
    
    # dim(df_result)  193751
    df_result <- df_result[df_result$R2PAS < 10^5,]
    
    qc.export <- paste0(qc.export,
                        "R2PAS filted reads: ",nrow(df_result),"\n")
    
    #### filter PCR duplicated reads by R1toR2
    # dim(df_result)  193751
    
    df_result <- df_result %>%  group_by(tag.CB,tag.UB,gene_id) %>% filter(R2PAS == min(R2PAS))
    df_result <- df_result[!duplicated(df_result[,c("tag.CB","tag.UB","gene_id")]),]
    #df_result <- data.frame(df_result)
    
    qc.export <- paste0(qc.export,
                        "PCR de-dup reads: ",nrow(df_result),"\n")
  
    writeLines(qc.export, paste0(outpath,"/stat.qc/",samplename,"/stat.",chr,".txt"))
    
    #df_result <- df_result[order(df_result$PAS),]
    
    #### return PAS data
    df_result$PAS <- df_result$end
    df_result$PAS[df_result$strand == "-"] <- df_result$start[df_result$strand == "-"]
    
    saveRDS(df_result,file = paste0(outpath,"/bam2df/",samplename,"/",chr,".bam2df.rds"))
  } else {
    df_result <- readRDS(paste0(outpath,"/bam2df/",samplename,"/",chr,".bam2df.rds"))
  }
  
  #### return merged PAS sites
  df_result$start <- df_result$PAS
  df_result$end <-  df_result$start
  
  df_result_1 <- as(df_result[df_result$strand == "-",],"GRanges")
  df_result_1@seqinfo@seqlengths <- chrom.size[df_result_1@seqinfo@seqnames,]$seqlengths
  df_result_1 <- coverage(df_result_1)
  
  df_result_2 <- as(df_result[df_result$strand == "+",],"GRanges")
  df_result_2@seqinfo@seqlengths <- chrom.size[df_result_2@seqinfo@seqnames,]$seqlengths
  df_result_2 <- coverage(df_result_2)
  
  saveRDS(list(df_result_1 = df_result_1,df_result_2 = df_result_2),
          file = paste0(outpath,"/bam2df/",samplename,"/",chr,".cov.rds"))
  
  tmp <- df_result[df_result$strand == "-",]
  
  if (nrow(tmp) > 1000){
    gr_pas_1 <- as(tmp,"GRanges")
    gr_pas_1 <- reduce(gr_pas_1, min.gapwidth = 0, ignore.strand=FALSE)
  } else {
    gr_pas_1 <- NULL
  }

  tmp <- df_result[df_result$strand == "+",]
  
  if (nrow(tmp) > 1000){
    gr_pas_2 <- as(tmp,"GRanges")
    gr_pas_2 <- reduce(gr_pas_2, min.gapwidth = 0, ignore.strand=FALSE)
  } else {
    gr_pas_2 <- NULL
  }
  
  gr_pas <- c(gr_pas_1,gr_pas_2)
  
  return(gr_pas)
}


bam2df.r2 <- function(samplename,chr,outpath,cell.metadata = NULL,chrom.size,pa.include = F){
  
  bam_r2 = paste0(outpath,"/chr.bam.R2/",samplename,"/",chr,".bam")
  
  if (!file.exists(paste0(outpath,"/stat.qc/"))){
    dir.create(paste0(outpath,"/stat.qc/"))
  }
  
  if (!file.exists(paste0(outpath,"/bam2df/"))){
    dir.create(paste0(outpath,"/bam2df/"))
  }
  
  if (!file.exists(paste0(outpath,"/stat.qc/",samplename))){
    dir.create(paste0(outpath,"/stat.qc/",samplename))
  }
  
  if (!file.exists(paste0(outpath,"/bam2df/",samplename))){
    dir.create(paste0(outpath,"/bam2df/",samplename))
  }
  
  #### read bam
  df_r2 <- Rsamtools::scanBam(bam_r2, 
                   param = Rsamtools::ScanBamParam(what=c('qname',"strand","pos","mapq","cigar","seq","qual"),tag = c("CB","UB"))) #'flag',
  df_r2 <- do.call(rbind, lapply(df_r2, as.data.frame))
  df_r2 <- df_r2[!is.na(df_r2$tag.CB),]
  df_r2 <- df_r2[df_r2$mapq == 255,]
  df_r2$mapL <- GenomicAlignments::cigarWidthAlongReferenceSpace(df_r2$cigar)
  
  if (!is.null(cell.metadata)){
    df_r2 <- df_r2[df_r2$tag.CB %in% gsub(".+_","",cell.metadata$cb[cell.metadata$samplename == samplename]),]
  }
  
  df_r2$mapL <- df_r2$mapL + df_r2$pos - 1
  
  df_r2_pa = df_r2[grepl("[0-9]+S$",df_r2$cigar),]
  
  df_r2$mapq <- NULL
  df_r2$cigar <- NULL
  df_r2$seq <- NULL
  df_r2$qual <- NULL
  colnames(df_r2)[6] <- "end"
  
  qc.export <- paste0("R2 reads: ",nrow(df_r2),"\n")
  
  saveRDS(df_r2,file = paste0(outpath,"/bam2df/",samplename,"/",chr,".r2bam2df.rds"))
  
  rm(df_r2)
  gc()
  
  ### processing PA reads
  if (pa.include){
    len <- regmatches(df_r2_pa$cigar, regexpr("\\d+(?=S$)", df_r2_pa$cigar, perl = TRUE))
    len <- as.integer(len)
    
    df_r2_pa$seq <- substring(as.character(df_r2_pa$seq), nchar(df_r2_pa$seq) - len + 1, nchar(df_r2_pa$seq))
    df_r2_pa$qual <- substring(as.character(df_r2_pa$qual), nchar(df_r2_pa$qual) - len + 1, nchar(df_r2_pa$qual))
    df_r2_pa = df_r2_pa[nchar(df_r2_pa$seq) > 5,]
    
    a.ratio = letterFrequency(DNAStringSet(df_r2_pa$seq),"A")
    df_r2_pa = df_r2_pa[a.ratio >= 5 & a.ratio/nchar(df_r2_pa$seq) > 0.8,]
    
    df_r2_pa$mapq <- NULL
    df_r2_pa$cigar <- NULL
    df_r2_pa$seq <- NULL
    df_r2_pa$qual <- NULL
    colnames(df_r2_pa)[6] <- "end"
    df_r2_pa = df_r2_pa[!duplicated(df_r2_pa[,c("tag.CB","tag.UB")]),]
    
    df_r2_pa$tag.CB <- NULL
    df_r2_pa$tag.UB <- NULL
    df_r2_pa$pos <- NULL
    
    df_r2_pa$chr = chr
    df_r2_pa$start = df_r2_pa$end
    
    qc.export <- paste0(qc.export,
                        "R1 reads: ",nrow(df_r2_pa),"\n")
    
    writeLines(qc.export, paste0(outpath,"/stat.qc/",samplename,"/stat.",chr,".txt"))
    saveRDS(df_r2_pa,file = paste0(outpath,"/bam2df/",samplename,"/",chr,".r2pa.bam2df.rds"))
    
    df_result_1 <- as(df_r2_pa[df_r2_pa$strand == "-",],"GRanges")
    df_result_1@seqinfo@seqlengths <- chrom.size[df_result_1@seqinfo@seqnames,]$seqlengths
    df_result_1 <- coverage(df_result_1)
    
    df_result_2 <- as(df_r2_pa[df_r2_pa$strand == "+",],"GRanges")
    df_result_2@seqinfo@seqlengths <- chrom.size[df_result_2@seqinfo@seqnames,]$seqlengths
    df_result_2 <- coverage(df_result_2)
    
    saveRDS(list(df_result_1 = df_result_1,df_result_2 = df_result_2),
            file = paste0(outpath,"/bam2df/",samplename,"/",chr,".cov.rds"))
  }
  
}



bam2df.bulk <- function(samplename,chr,outpath,geneanno,intronanno,chrom.size,cache.bam2df){
  
  bam_r2 = paste0(outpath,"/chr.bam.R2/",samplename,"/",chr,".bam")
  bam_r1 = paste0(outpath,"/chr.bam.R1/",samplename,"/",chr,".bam")
  
  if (!file.exists(paste0(outpath,"/stat.qc/"))){
    dir.create(paste0(outpath,"/stat.qc/"))
  }
  
  if (!file.exists(paste0(outpath,"/bam2df/"))){
    dir.create(paste0(outpath,"/bam2df/"))
  }
  
  if (!file.exists(paste0(outpath,"/stat.qc/",samplename))){
    dir.create(paste0(outpath,"/stat.qc/",samplename))
  }
  
  if (!file.exists(paste0(outpath,"/bam2df/",samplename))){
    dir.create(paste0(outpath,"/bam2df/",samplename))
  }
  
  if(!file.exists(paste0(outpath,"/bam2df/",samplename,"/",chr,".bam2df.rds")) | !cache.bam2df){
    
    #### read bam
    df_r2 <- Rsamtools::scanBam(file = bam_r2, param = Rsamtools::ScanBamParam(what=c('qname','rname',"strand","pos","mapq","cigar"))) #'flag',
    df_r2 <- do.call(rbind, lapply(df_r2, as.data.frame))
    colnames(df_r2)[2:6] <- paste0(colnames(df_r2)[2:6],"_r2")
    df_r2$mapL_r2 <- GenomicAlignments::cigarWidthAlongReferenceSpace(df_r2$cigar_r2)
    
    df_r2 <- df_r2[df_r2$mapq_r2 >= 255,]
    df_r2$mapL_r2 <- df_r2$mapL_r2 + df_r2$pos_r2 - 1
    df_r2$mapq_r2 <- NULL
    df_r2$cigar_r2 <- NULL
    colnames(df_r2)[5] <- "end_r2"
    df_r2$tag.CB = samplename
    df_r2$tag.UB = paste0(df_r2$rname,"-",df_r2$pos_r2,"-",df_r2$mapL_r2)
    
    qc.export <- paste0("R2 reads: ",nrow(df_r2),"\n")
    
    saveRDS(df_r2,file = paste0(outpath,"/bam2df/",samplename,"/",chr,".r2bam2df.rds"))
    
    df_r1 <- Rsamtools::scanBam(bam_r1, param = Rsamtools::ScanBamParam(what=c('qname','rname',"strand","pos","mapq","cigar")))
    df_r1 <- do.call(rbind, lapply(df_r1, as.data.frame))
    colnames(df_r1)[2:6] <- paste0(colnames(df_r1)[2:6],"_r1")
    df_r1$mapL_r1 <- GenomicAlignments::cigarWidthAlongReferenceSpace(df_r1$cigar_r1)
    df_r1 <- df_r1[df_r1$mapq_r1 >= 255,]
    df_r1$mapL_r1 <- df_r1$mapL_r1 + df_r1$pos_r1 - 1
    df_r1$mapq_r1 <- NULL
    df_r1$cigar_r1 <- NULL
    colnames(df_r1)[5] <- "end_r1"
    
    qc.export <- paste0(qc.export,
                        "R1 reads: ",nrow(df_r1),"\n")
    
    # df_r1$matched_lengths <- lapply(df_r1$cigar,function(x){
    #   sum(as.numeric(gsub("M","",unlist(regmatches(x, gregexpr("\\d+M", x)))))) 
    # }) %>% Reduce("c",.)
    
    #### merge R1 and R2 & strand filter
    df_result <- merge(data.table(df_r1) ,data.table(df_r2))
    qc.export <- paste0(qc.export,
                        "Pair reads: ",nrow(df_result),"\n")
    
    df_result <- df_result[df_result$strand_r1 != df_result$strand_r2,]
    df_result <- df_result[df_result$rname_r1 == df_result$rname_r2,]
    qc.export <- paste0(qc.export,
                        "Strand match reads: ",nrow(df_result),"\n")
    
    #### R1/R2 out of range filter
    #       ===========> #R2
    #    <================    #R1 out of range
    df_result <- df_result[!(df_result$pos_r1 < df_result$pos_r2 & df_result$strand_r2 == "+"),]
    #       ===========> #R2 out of range
    #         <======    #R1
    df_result <- df_result[!(df_result$end_r1 < df_result$end_r2 & df_result$strand_r2 == "+"),]
    
    df_result <- df_result[!(df_result$pos_r1 > df_result$pos_r2 & df_result$strand_r2 == "-"),]
    df_result <- df_result[!(df_result$end_r1 > df_result$end_r2 & df_result$strand_r2 == "-"),]
    # dim(df_result) 258763
    
    qc.export <- paste0(qc.export,
                        "Nice mapped reads: ",nrow(df_result),"\n")
    
    #### gene anno
    tmp <- data.frame(chr = df_result$rname_r2,
                      start = df_result$pos_r2,
                      end = df_result$end_r2,
                      strand = df_result$strand_r2)
    
    gr.index = suppressMessages(findOverlaps(as(tmp,"GRanges"),geneanno))
    
    overlap_length <- pintersect(as(tmp,"GRanges")[queryHits(gr.index)], geneanno[subjectHits(gr.index)])
    
    df_result <- cbind(df_result[queryHits(gr.index),],
                       data.frame(geneanno)[subjectHits(gr.index),c("seqnames","gene_id","gene_name")])
    
    # df_result$overlap.width = width(overlap_length)/(df_result$end_r2 - df_result$pos_r2 + 1)
    
    # dim(df_result) 243314
    # quantile(df_result$overlap.width ,seq(0,1,0.05))
    
    df_result <- data.table(df_result)

    df_result$rname_r1 <- NULL
    df_result$rname_r2 <- NULL
    # dim(df_result) 235926
    # df_result2$gene_name[duplicated(df_result2$qname)]
    
    # df_result[df_result$qname == "A00511:334:HHVY3DSXY:2:2676:11433:10927",]
    # MRPS6 / SLC5A3
    
    qc.export <- paste0(qc.export,
                        "Gene mapped reads: ",nrow(df_result),"\n")
    qc.export <- paste0(qc.export,
                        "Gene dup reads: ",table(duplicated(df_result$qname))[2],"\n")
    
    qname.dup = df_result$qname[duplicated(df_result$qname)]
    df_result = df_result[df_result$qname %ni% qname.dup,]
    
    #### R2 to PAS distance filter
    
    df_result$start <- df_result$end_r2
    df_result$start[df_result$strand_r2 == "-"] <- df_result$pos_r1[df_result$strand_r2 == "-"]
    df_result$end <- df_result$end_r1
    df_result$end[df_result$strand_r2 == "-"] <- df_result$pos_r2[df_result$strand_r2 == "-"]
    df_result$strand = df_result$strand_r2
    df_result$strand_r2 <- NULL
    df_result$strand_r1 <- NULL
    
    df_result <- df_result[df_result$end - df_result$start < 10^6,]
    
    gr.index <- suppressMessages(findOverlaps(intronanno, type=c("within"),as(df_result,"GRanges")))
    gr.index = gr.index[intronanno[from(gr.index)]$gene_id == df_result$gene_id[to(gr.index)],]
    
    intro_df <- data.frame(index_df = to(gr.index),
                           intron_width = intronanno[from(gr.index)]$intron_width)
    intro_df <- intro_df %>% group_by(index_df) %>% summarise(intron_total = sum(intron_width))
    
    df_result$intron <- 0
    df_result$intron[intro_df$index_df] <- intro_df$intron_total
    
    df_result$R2PAS <- df_result$end - df_result$start - df_result$intron
    
    # table(df_result$gene_name[df_result$R2PAS > 5000])
    # data.frame(df_result[df_result$R2PAS > 5000 & df_result$gene_name == "HMGN1",])[1:4,]
    
    # dim(df_result)  193751
    df_result <- df_result[df_result$R2PAS < 10^5,]
    
    qc.export <- paste0(qc.export,
                        "R2PAS filted reads: ",nrow(df_result),"\n")
    
    writeLines(qc.export, paste0(outpath,"/stat.qc/",samplename,"/",chr,".stat.txt"))
    
    #df_result <- df_result[order(df_result$PAS),]
    
    #### return PAS data
    df_result$PAS <- df_result$end
    df_result$PAS[df_result$strand == "-"] <- df_result$start[df_result$strand == "-"]
    
    saveRDS(df_result,file = paste0(outpath,"/bam2df/",samplename,"/",chr,".bam2df.rds"))
  } else {
    df_result <- readRDS(paste0(outpath,"/bam2df/",samplename,"/",chr,".bam2df.rds"))
  }
  
  #### return merged PAS sites
  df_result$start <- df_result$PAS
  df_result$end <-  df_result$start
  
  df_result_1 <- as(df_result[df_result$strand == "-",],"GRanges")
  df_result_1@seqinfo@seqlengths <- chrom.size[df_result_1@seqinfo@seqnames,]$seqlengths
  df_result_1 <- coverage(df_result_1)
  
  df_result_2 <- as(df_result[df_result$strand == "+",],"GRanges")
  df_result_2@seqinfo@seqlengths <- chrom.size[df_result_2@seqinfo@seqnames,]$seqlengths
  df_result_2 <- coverage(df_result_2)
  
  saveRDS(list(df_result_1 = df_result_1,df_result_2 = df_result_2),
          file = paste0(outpath,"/bam2df/",samplename,"/",chr,".cov.rds"))
  
  tmp <- df_result[df_result$strand == "-",]
  
  if (nrow(tmp) > 1000){
    gr_pas_1 <- as(tmp,"GRanges")
    gr_pas_1 <- reduce(gr_pas_1, min.gapwidth = 0, ignore.strand=FALSE)
  } else {
    gr_pas_1 <- NULL
  }
  
  tmp <- df_result[df_result$strand == "+",]
  
  if (nrow(tmp) > 1000){
    gr_pas_2 <- as(tmp,"GRanges")
    gr_pas_2 <- reduce(gr_pas_2, min.gapwidth = 0, ignore.strand=FALSE)
  } else {
    gr_pas_2 <- NULL
  }
  
  gr_pas <- c(gr_pas_1,gr_pas_2)
  
  return(gr_pas)
}

PAScalling  <- function(samplenames,chrs,outpath,cell.metadata = NULL,chrom.size,threads = 24,cache.covdata = F,
                        callpeak = "macs3",genomeSize = 2.7e9, cutOff = 10^-5,additionalParams = "--nomodel",
                        min.gapwidth = 50,sum.count = 5,min.count = 5,cpm.cutoff = 0.5,
                        BSg = BSgenome.Hsapiens.UCSC.hg38,literal = "AAAAAAAA",num_A = 10){
  
  if (!cache.covdata){
    combinations <- expand.grid(samplenames, chrs)
    plan(sequential)
    plan(multisession, workers= threads)
    print("Reading bam")
    future_lapply(1:nrow(combinations),function(i){
      
      samplename <- combinations$Var1[i]
      chr <- combinations$Var2[i]
      
      bam2df.r2(samplename,chr,outpath,cell.metadata,chrom.size,pa.include =T)
    })
    plan(sequential)
  }
  print("Peak calling")
  pas_peak_gr <- findPAS(samplenames,chrs,outpath,threads,cache.covdata = F,
                      callpeak,genomeSize, cutOff,additionalParams,
                      min.gapwidth,sum.count,min.count,cpm.cutoff,
                      BSg = BSgenome.Hsapiens.UCSC.hg38,literal,num_A)
  return(pas_peak_gr)
}

ProcessBam <- function(samplenames,chrs,outpath,geneanno,intronanno,cell.metadata = NULL,chrom.size,cache.bam2df = F,threads = 24,isbulk = F){
  gr_pas.list <- list()
  combinations <- expand.grid(samplenames, chrs)
  
  check_split_bam <- function(sample, chr, read = c("R1","R2")) {
    read <- match.arg(read)
    bam_path <- file.path(outpath, paste0("chr.bam.",read), sample, paste0(chr,".bam"))
    if (!file.exists(bam_path)) {
      stop("Split BAM not found: ", bam_path, ". Please run BamSplit first.")
    }
    info <- file.info(bam_path)
    if (is.na(info$size) || info$size == 0) {
      stop("Split BAM has zero size: ", bam_path)
    }
    if (difftime(Sys.time(), info$mtime, units = "secs") < 5) {
      stop("Split BAM appears to be still being modified: ", bam_path)
    }
  }
  
  if (!isbulk) {
    apply(combinations, 1, function(row) {
      check_split_bam(row[1], row[2], "R1")
      check_split_bam(row[1], row[2], "R2")
    })
  } else {
    apply(combinations, 1, function(row) {
      check_split_bam(row[1], row[2], "R2")
    })
  }
  
  plan(sequential)
  if (threads > 1) {
    plan(multisession, workers= threads)
  }
  gr_pas_data <- future_lapply(1:nrow(combinations),function(i){
    
    samplename <- as.character(combinations$Var1[i]) 
    chr <- as.character(combinations$Var2[i]) 
    
    if (isbulk){
      bam2df.bulk(samplename,chr,outpath,geneanno,intronanno,chrom.size,cache.bam2df = cache.bam2df)
    } else {
      bam2df(samplename,chr,outpath,geneanno,intronanno,cell.metadata,chrom.size,cache.bam2df = cache.bam2df)
    }
  })
  
  plan(sequential)
  
  gr_pas.list <- do.call("c",gr_pas_data)
  gr_pas.list <- reduce(gr_pas.list, min.gapwidth = 0, ignore.strand=FALSE)
  
  saveRDS(gr_pas.list,
          file = paste0(outpath,"/bam2df/gr_pas.rds"))
}
