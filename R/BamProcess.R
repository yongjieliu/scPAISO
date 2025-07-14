#' Process bam files
#' 
#' 
#' @param bam.split.sh Path to the bam.split.sh file.
#' @param bam.r1 A vector containing bam file path of Read1.
#' @param bam.r1 A vector containing bam file path of Read2.
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

BamSplit <- function(bam.split.sh,bam.r1,bam.r2,samplenames,outpath){
  for (i in 1:length(bam.r2)){
    system(paste0(bam.split.sh," ",bam.r1[i]," ",bam.r2[i]," ",samplenames[i]," ",outpath)) 
  }
}

bam2df <- function(samplename,chr,outpath,geneanno,intronanno,cell.metadata,chrom.size,cache.bam2df){
  
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
    df_r2 <- scanBam(bam_r2, param = ScanBamParam(what=c('qname',"strand","pos","mapq","cigar"),tag = c("CB","UB"))) #'flag',
    df_r2 <- do.call(rbind, lapply(df_r2, as.data.frame))
    df_r2 <- df_r2[!is.na(df_r2$tag.CB),]
    colnames(df_r2)[2:5] <- paste0(colnames(df_r2)[2:5],"_r2")
    df_r2$mapL_r2 <- GenomicAlignments::cigarWidthAlongReferenceSpace(df_r2$cigar_r2)
    
    df_r2 <- df_r2[df_r2$tag.CB %in% gsub(".+_","",cell.metadata$cb[cell.metadata$samplename == samplename]),]
  
    df_r2 <- df_r2[df_r2$mapq_r2 >= 255,]
    df_r2$mapL_r2 <- df_r2$mapL_r2 + df_r2$pos_r2 - 1
    df_r2$mapq_r2 <- NULL
    df_r2$cigar_r2 <- NULL
    colnames(df_r2)[6] <- "end_r2"
    
    qc.export <- paste0("R2 reads: ",nrow(df_r2),"\n")
    
    saveRDS(df_r2,file = paste0(outpath,"/bam2df/",samplename,"/",chr,".r2bam2df.rds"))
    
    df_r1 <- scanBam(bam_r1, param = ScanBamParam(what=c('qname',"strand","pos","mapq","cigar")))
    df_r1 <- do.call(rbind, lapply(df_r1, as.data.frame))
    colnames(df_r1)[2:5] <- paste0(colnames(df_r1)[2:5],"_r1")
    df_r1$mapL_r1 <- GenomicAlignments::cigarWidthAlongReferenceSpace(df_r1$cigar_r1)
    df_r1 <- df_r1[df_r1$mapq_r1 >= 255,]
    df_r1$mapL_r1 <- df_r1$mapL_r1 + df_r1$pos_r1 - 1
    df_r1$mapq_r1 <- NULL
    df_r1$cigar_r1 <- NULL
    colnames(df_r1)[4] <- "end_r1"
    
    qc.export <- paste0(qc.export,
                        "R1 reads: ",nrow(df_r1),"\n")
    
    # df_r1$matched_lengths <- lapply(df_r1$cigar,function(x){
    #   sum(as.numeric(gsub("M","",unlist(regmatches(x, gregexpr("\\d+M", x)))))) 
    # }) %>% Reduce("c",.)
    
    #### merge R1 and R2 & strand filter
    df_result <- merge(df_r1,df_r2)
    qc.export <- paste0(qc.export,
                        "Pair reads: ",nrow(df_result),"\n")
    
    df_result <- df_result[df_result$strand_r1 != df_result$strand_r2,]
    
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
    tmp <- data.frame(chr = chr,
                      start = df_result$pos_r2,
                      end = df_result$end_r2,
                      strand = df_result$strand_r2)
    
    gr.index = findOverlaps(as(tmp,"GRanges"),geneanno)
    
    overlap_length <- pintersect(as(tmp,"GRanges")[queryHits(gr.index)], geneanno[subjectHits(gr.index)])
    
    df_result <- cbind(df_result[queryHits(gr.index),],
                       data.frame(geneanno)[subjectHits(gr.index),c("seqnames","gene_id","gene_name")])
    
    #df_result$overlap.width = width(overlap_length)/(df_result$end_r2 - df_result$pos_r2 + 1)
  
    # dim(df_result) 243314
    # quantile(df_result$overlap.width ,seq(0,1,0.05))
    
    df_result <- data.table(df_result)
    
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
    df_r2 <- scanBam(file = bam_r2, param = ScanBamParam(what=c('qname','rname',"strand","pos","mapq","cigar"))) #'flag',
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
    
    df_r1 <- scanBam(bam_r1, param = ScanBamParam(what=c('qname','rname',"strand","pos","mapq","cigar")))
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

ProcessBam <- function(samplenames,chrs,outpath,geneanno,intronanno,cell.metadata = NULL,chrom.size,cache.bam2df = F,threads = 24,isbulk = F){
  gr_pas.list <- list()
  combinations <- expand.grid(samplenames, chrs)
  
  plan(sequential)
  plan(multisession, workers= threads)
  gr_pas_data <- future_lapply(1:nrow(combinations),function(i){
    
    samplename <- combinations$Var1[i]
    chr <- combinations$Var2[i]
    
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
