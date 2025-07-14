#' PAS peaks identification
#' 
#' @param samplenames A vector containing sample names.
#' @param chrs A vector containing chromatin to use. 
#' @param outpath A character describing where the file will be saved.
#' @param threads The number of threads to be used for parallel computing.
#' @param cache.covdata A boolean value that determines whether cache coverage file shoulde be used.
#' @param callpeak Method to identify PAS peak. Options include "macs2","macs3", or "trim". 
#' @param genomeSize The number of genome size.
#' @param cutOff The number of peak threshold for macs-based mathod.
#' @param additionalParams Additional params for macs-based mathod.
#' @param min.gapwidth The minimum distance between PAS peaks. Peaks gaps smaller than this distance will be merged.
#' @param sum.count The minimum count of PAS peak. Peaks count smaller than this will be removed.
#' @param cpm.cutoff The minimum CPM normalized count of PAS peak. PAS peaks with a count smaller than this will be removed.
#'
#' @param min.count The minimum count of PAS peak identification for "trim" method.
#' 
#' @param literal If the literal sequence is detected within 5bp downstream and 20bp upstream of PAS 3' end, the PAS peak will be deleted.
#' @param num_A If more than the number of A base is detected within 5bp downstream and 20bp upstream of PAS 3' end, the PAS peak will be deleted.
#' 
#' @export

findPAS <- function(samplenames,chrs,outpath,threads,cache.covdata = F,
                    callpeak = "macs3",genomeSize = 2.7e9, cutOff = 10^-5,additionalParams = "--nomodel",
                    min.gapwidth = 50,sum.count = 10,min.count = 10,cpm.cutoff = 0.5,
                    BSg = BSgenome.Hsapiens.UCSC.hg38,literal = "AAAAAAAA",num_A = 10){
  
  # bed for macs3
  # covrage for trim and macs3
  if (!file.exists(paste0(outpath,"/pas.gr"))){
    dir.create(paste0(outpath,"/pas.gr"))
  }
  
  if (!cache.covdata){
    print("PAS Coverage")
    bam_df_1 <- NULL
    bam_df_2 <- NULL
    for (samplename in samplenames){
      for (i in chrs) {
        cov.data <- readRDS(paste0(outpath,"/bam2df/",samplename,"/",i,".cov.rds"))
        if (is.null(bam_df_1)){
          bam_df_1 <- cov.data$df_result_1
          bam_df_2 <- cov.data$df_result_2
        } else {
          bam_df_1 <- bam_df_1 + cov.data$df_result_1
          bam_df_2 <- bam_df_2 + cov.data$df_result_2
        }
      }
    }
    
    print("Bed writing")
    
    rtracklayer::export.bw(object = bam_df_1, con = paste0(outpath,"/pas.gr/cov1.bigwig"))
    bam_bed_1 <- as(bam_df_1,"GRanges")
    bam_bed_1 <- bam_bed_1[bam_bed_1$score > 0,]
    bam_bed_1 <- data.frame(bam_bed_1)
    data.table::fwrite(bam_bed_1[,c("seqnames","start","end","score")], paste0(outpath,"/pas.gr/cov1.bed"), 
                       sep = "\t", col.names = FALSE, append = F)
    
    rtracklayer::export.bw(object = bam_df_2, con = paste0(outpath,"pas.gr/cov2.bigwig"))
    bam_bed_2 <- as(bam_df_2,"GRanges")
    bam_bed_2 <- bam_bed_2[bam_bed_2$score > 0,]
    bam_bed_2 <- data.frame(bam_bed_2)
    data.table::fwrite(bam_bed_2[,c("seqnames","start","end","score")], paste0(outpath,"/pas.gr/cov2.bed"), 
                       sep = "\t", col.names = FALSE, append = F)
    
    bam_df = list(bam_df_1 = bam_df_1,bam_df_2 = bam_df_2)
    saveRDS(bam_df,
            file = paste0(outpath,"/pas.gr/cov.rds"))
  }
  
  bam_df <- readRDS(paste0(outpath,"/pas.gr/cov.rds"))
  bam_df_1 <- bam_df$bam_df_1
  bam_df_2 <- bam_df$bam_df_2
  
  if (grepl("macs3",callpeak)){
    print("Peak calling by macs3")
    cmd <- sprintf("callpeak -g %s --name merged1 --treatment %s --outdir %s --format BED --max-gap 10 --min-length 5 --extsize 3 --shift -1 --call-summits --keep-dup all -p %s %s", #  
                   genomeSize, paste0(outpath,"/pas.gr/cov1.bed"), paste0(outpath,"/pas.gr"), cutOff, additionalParams) # 
    
    run <- system2(callpeak, cmd, wait=TRUE) #, stdout=NULL, stderr=NULL
    
    cmd <- sprintf("callpeak -g %s --name merged2 --treatment %s --outdir %s --format BED --max-gap 10 --min-length 5 --extsize 3 --shift -1 --call-summits --keep-dup all -p %s %s", #  
                   genomeSize, paste0(outpath,"/pas.gr/cov2.bed"), paste0(outpath,"/pas.gr"), cutOff, additionalParams) # 
    
    run <- system2(callpeak, cmd, wait=TRUE)
    
    tmp1 <- read.table(paste0(outpath,"pas.gr/merged1_peaks.narrowPeak"))
    colnames(tmp1) <- c("seqnames","start","end","peakname","score","strand","signalValue","pValue","qValue","peak")
    tmp1$strand = "-"
    tmp2 <- read.table(paste0(outpath,"pas.gr/merged2_peaks.narrowPeak"))
    colnames(tmp2) <- c("seqnames","start","end","peakname","score","strand","signalValue","pValue","qValue","peak")
    tmp2$strand = "+"
    
    merge.peak <- as(rbind(tmp1,tmp2),"GRanges")
    ### filter PAS by counts
    gr_pas_1 <- merge.peak[strand(merge.peak) == "-"]
    #gr_pas_1 <- reduce(gr_pas_1, min.gapwidth = 0, ignore.strand=FALSE)
    gr_pas_1$sum <-  grangesSums(gr_pas_1, bam_df_1)
    gr_pas_1 <- gr_pas_1[gr_pas_1$sum >= sum.count,]
    gr_pas_1 <- reduce(gr_pas_1, min.gapwidth = 30, ignore.strand=FALSE)
    gr_pas_1$sum <-  grangesSums(gr_pas_1, bam_df_1)
    gr_pas_1$max <-  grangesMaxs(gr_pas_1, bam_df_1)
    gr_pas_1$WhichMaxs <- grangesWhichMaxs(gr_pas_1, bam_df_1)
    
    gr_pas_2 <- merge.peak[strand(merge.peak) == "+"]
    #gr_pas_2 <- reduce(gr_pas_2, min.gapwidth = 0, ignore.strand=FALSE)
    gr_pas_2$sum <-  grangesSums(gr_pas_2, bam_df_2)
    gr_pas_2 <- gr_pas_2[gr_pas_2$sum >= sum.count,]
    gr_pas_2 <- reduce(gr_pas_2, min.gapwidth = 30, ignore.strand=FALSE)
    gr_pas_2$sum <-  grangesSums(gr_pas_2, bam_df_2)
    gr_pas_2$max <-  grangesMaxs(gr_pas_2, bam_df_2)
    gr_pas_2$WhichMaxs <- grangesWhichMaxs(gr_pas_2, bam_df_2)
    
    gr_pas <- c(gr_pas_1,gr_pas_2)
    callpeak.label = callpeak
    saveRDS(gr_pas,
            file = paste0(outpath,"/pas.gr/gr_pas_peak_macs.rds"))
  } else if (grepl("trim",callpeak)){
    print("Peak calling by count")
    
    gr_pas <- readRDS(paste0(outpath,"/bam2df/gr_pas.rds"))
    
    gr_pas_1 <- gr_pas[strand(gr_pas) == "-"]
    gr_pas_1 <- reduce(gr_pas_1, min.gapwidth = 0, ignore.strand=FALSE)
    gr_pas_1$sum <-  grangesSums(gr_pas_1, bam_df_1)
    gr_pas_1 <- gr_pas_1[gr_pas_1$sum >= min.count,]
    gr_pas_1 <- reduce(gr_pas_1, min.gapwidth = 30, ignore.strand=FALSE)
    gr_pas_1$sum <-  grangesSums(gr_pas_1, bam_df_1)
    gr_pas_1 <- gr_pas_1[gr_pas_1$sum >= sum.count,]
    gr_pas_1$max <-  grangesMaxs(gr_pas_1, bam_df_1)
    #gr_pas_1 <- gr_pas_1[gr_pas_1$max > top.count,]
    gr_pas_1$WhichMaxs <- grangesWhichMaxs(gr_pas_1, bam_df_1)
    
    gr_pas_2 <- gr_pas[strand(gr_pas) == "+"]
    gr_pas_2 <- reduce(gr_pas_2, min.gapwidth = 0, ignore.strand=FALSE)
    gr_pas_2$sum <-  grangesSums(gr_pas_2, bam_df_2)
    gr_pas_2 <- gr_pas_2[gr_pas_2$sum >= min.count,]
    gr_pas_2 <- reduce(gr_pas_2, min.gapwidth = 30, ignore.strand=FALSE)
    gr_pas_2$sum <-  grangesSums(gr_pas_2, bam_df_2)
    gr_pas_2 <- gr_pas_2[gr_pas_2$sum >= sum.count,]
    gr_pas_2$max <-  grangesMaxs(gr_pas_2, bam_df_2)
    #gr_pas_2 <- gr_pas_2[gr_pas_2$max > top.count,]
    gr_pas_2$WhichMaxs <- grangesWhichMaxs(gr_pas_2, bam_df_2)
    
    gr_pas <- c(gr_pas_1,gr_pas_2)
    callpeak.label = callpeak
    saveRDS(gr_pas,
            file = paste0(outpath,"/pas.gr/gr_pas_peak_trim.rds"))
  }
  
  total_count <- sum(gr_pas$sum)
  gr_pas <- gr_pas[gr_pas$sum > total_count/10^6*cpm.cutoff]
  
  ### filter PAS by polyA in DNA
  print("filting PAS by polyA")
  
  # by peak end
  pas.index_1 <- gr_pas[strand(gr_pas) == "-"]
  end(pas.index_1) <- start(pas.index_1) + 5 #
  start(pas.index_1) <- start(pas.index_1) - 20

  pas.index_2 <- gr_pas[strand(gr_pas) == "+"]
  start(pas.index_2) <- end(pas.index_2) - 5 #
  end(pas.index_2) <- end(pas.index_2) + 20

  pas.index <- getSeq(BSg, c(pas.index_1,pas.index_2))
  
  gr_pas$filter <- "keep"
  gr_pas$filter[grepl(literal,data.frame(pas.index)$pas.index)] <- "remove"
  gr_pas$filter[letterFrequency(pas.index,"A") >= num_A & grepl("AAAAA",data.frame(pas.index)$pas.index)] <- "remove"
  
  # by peak which max
  pas.index <- gr_pas
  
  pas.index_1 <- gr_pas[strand(gr_pas) == "-"]
  end(pas.index_1) <- pas.index$WhichMaxs + 12
  start(pas.index_1) <- pas.index$WhichMaxs - 12
  
  pas.index_2 <- gr_pas[strand(gr_pas) == "+"]
  start(pas.index) <- pas.index$WhichMaxs - 12
  end(pas.index) <- pas.index$WhichMaxs + 12
  
  pas.index <- getSeq(BSg, pas.index)
  
  gr_pas$filter[grepl(literal,data.frame(pas.index)$pas.index)] <- "remove"
  gr_pas$filter[letterFrequency(pas.index,"A") >= num_A & 
                  grepl("AAAAA",data.frame(pas.index)$pas.index)] <- "remove"
  
  gr_pas <- sort(gr_pas)
  
  saveRDS(gr_pas[gr_pas$filter == "remove"],
          file = paste0(outpath,"/pas.gr/gr_pas_peak_remove.",callpeak,".rds"))
  
  gr_pas <- gr_pas[gr_pas$filter == "keep"]
  
  # reduce 
  gr_pas_1 <- gr_pas[strand(gr_pas) == "-"]
  gr_pas_1 <- reduce(gr_pas_1, min.gapwidth = min.gapwidth, ignore.strand=FALSE)
  gr_pas_1$sum <-  grangesSums(gr_pas_1, bam_df_1)
  gr_pas_1$max <-  grangesMaxs(gr_pas_1, bam_df_1)
  gr_pas_1$WhichMaxs <- grangesWhichMaxs(gr_pas_1, bam_df_1)
  
  gr_pas_2 <- gr_pas[strand(gr_pas) == "+"]
  gr_pas_2 <- reduce(gr_pas_2, min.gapwidth = min.gapwidth, ignore.strand=FALSE)
  gr_pas_2$sum <-  grangesSums(gr_pas_2, bam_df_2)
  gr_pas_2$max <-  grangesMaxs(gr_pas_2, bam_df_2)
  gr_pas_2$WhichMaxs <- grangesWhichMaxs(gr_pas_2, bam_df_2)
  
  gr_pas <- c(gr_pas_1,gr_pas_2)
  
  gr_pas$filter <- NULL
  
  gr_pas$peak <- paste0(seqnames(gr_pas),":",start(gr_pas),"-",end(gr_pas))
  
  gr_pas$callpeak = callpeak.label
  
  saveRDS(gr_pas,
          file = paste0(outpath,"/pas.gr/gr_pas_peak_keep.",callpeak,".rds"))
  return(gr_pas)
}

#' link pas peak to gene
#' @param gr_pas A range of PAS peak.
#' @param samplenames A vector containing sample names.
#' @param chrs A vector containing chromatin to use. 
#' @param threads The number of threads to be used for parallel computing.
#' 
#' @export

pas_peak2gene <- function(gr_pas = pas_peak,samplenames = samplenames, 
                          chrs = chrs,threads= 24){
  
  if (!file.exists(paste0(outpath,"/pas.df"))){
    dir.create(paste0(outpath,"/pas.df"))
  }
  
  combinations <- expand.grid(samplenames, chrs)
  
  plan(sequential)
  plan(multisession, workers= threads)
  
  start(gr_pas) = start(gr_pas) - 16
  end(gr_pas) = end(gr_pas) + 16
  
  # filter not mapped reads and anno
  gr_pas.list <- future_lapply(1:nrow(combinations),function(i){
    
    samplename <- combinations$Var1[i]
    chr <- combinations$Var2[i]
    
    if (!file.exists(paste0(outpath,"/pas.df/",samplename))){
      dir.create(paste0(outpath,"/pas.df/",samplename))
    }
    
    bam_df <- readRDS(paste0(outpath,"/bam2df/",samplename,"/",chr,".bam2df.rds"))
    
    bam_df$start <- bam_df$PAS
    bam_df$end <-  bam_df$start
    
    bam_df <- data.frame(bam_df)
    
    gr.index = findOverlaps(as(bam_df,"GRanges"),gr_pas, ignore.strand=FALSE)
    gr_pas <- gr_pas[subjectHits(gr.index)]
    gr_pas$gene_id <- bam_df$gene_id[queryHits(gr.index)]
    gr_pas <- unique(gr_pas)
    
    gr.index = findOverlaps(as(bam_df,"GRanges"),gr_pas, ignore.strand=FALSE)
    bam_df <- bam_df[queryHits(gr.index),]
    bam_df$WhichMaxs <- gr_pas$WhichMaxs[subjectHits(gr.index)]
    bam_df$peak <- gr_pas$peak[subjectHits(gr.index)]
    bam_df$start <- start(gr_pas)[subjectHits(gr.index)] + 16
    bam_df$end <- end(gr_pas)[subjectHits(gr.index)] - 16
    bam_df <- bam_df[paste0(bam_df$gene_id,bam_df$peak) %in% paste0(gr_pas$gene_id,gr_pas$peak),]
    
    saveRDS(bam_df,file = paste0(outpath,"/pas.df/",samplename,"/",chr,".pas.df.rds"))
    
    plot_dist_data <- bam_df[,c("seqnames","start","end","strand","tag.CB","gene_id","R2PAS","PAS","pos_r2","end_r2")]
    saveRDS(plot_dist_data,file = paste0(outpath,"/pas.df/",samplename,"/",chr,".plot_dist_data.rds"))
    
    return(gr_pas)
  })
  
  gr_pas <- do.call("c",gr_pas.list )
  gr_pas <- unique(gr_pas)
  
  plan(sequential)

  start(gr_pas) = start(gr_pas) + 16
  end(gr_pas) = end(gr_pas) - 16
  
  gr_pas_pos <- data.frame(gr_pas[strand(gr_pas) == "+",])
  gr_pas_pos <- gr_pas_pos %>% 
    group_by(gene_id) %>% 
    arrange(WhichMaxs) %>% mutate(rank = row_number())
  
  gr_pas_neg <- data.frame(gr_pas[strand(gr_pas) == "-",])
  gr_pas_neg <- gr_pas_neg %>% 
    group_by(gene_id) %>% 
    arrange(desc(WhichMaxs)) %>% mutate(rank = row_number())
  
  gr_pas <- as(rbind(gr_pas_pos,gr_pas_neg),"GRanges") 
  gr_pas <- sort(gr_pas)
  
  saveRDS(gr_pas,file = paste0(outpath,"/pas.gr/gr_pas_peak_anno.rds"))
  
  return(gr_pas)
}


#' pas peak annotation
#' @param gr_pas A range of PAS peak.
#' @param gene.anno A ranges list describing 3' UTR, exon, and intron. 
#' @param chrs A vector containing chromatin to use. 
#' @param threads The number of threads to be used for parallel computing.
#' 
#' @export

pas_peak_anno <- function(gr_pas = pas_peak,gene.anno,
                          priority = c("utr3","last.exon","last.intron","exon","intron","utr5")){
  gr_pas$anno <- "utr3"
  gr_pas$anno.region <- ""
  
  for (i in rev(priority)){
    tmp = gene.anno[[i]]
    tmp.index = findOverlaps(gr_pas,tmp)
    tmp.index = tmp.index[gr_pas$gene_id[queryHits(tmp.index)] == tmp$gene_id[subjectHits(tmp.index)]]
    gr_pas$anno[queryHits(tmp.index)] = i
    if (i %in% c("utr3","last.exon")){
      gr_pas$anno.region[queryHits(tmp.index)] = tmp$region[subjectHits(tmp.index)]
    } else {
      gr_pas$anno.region[queryHits(tmp.index)] = paste0(seqnames(tmp),":",start(tmp),"-",end(tmp))[subjectHits(tmp.index)]
    }
  }
  
  if (length(gr_pas[gr_pas$anno.region == ""]) > 0){
    gr_pas$anno[gr_pas$anno.region == ""] = "utr3.new"
    tmp = gene.anno$utr3
    tmp.pos <- sort(tmp[strand(tmp) == "+",],decreasing = T)
    tmp.pos <-tmp.pos[!duplicated(tmp.pos$gene_id),]
    tmp.neg <- sort(tmp[strand(tmp) == "-",],decreasing = F)
    tmp.neg <-tmp.neg[!duplicated(tmp.neg$gene_id),]
    tmp <-  c(tmp.pos,tmp.neg)
    tmp$anno.region = paste0(seqnames(tmp),":",start(tmp),"-",end(tmp))
    gr_pas$anno.region[gr_pas$anno.region == ""] <- plyr::mapvalues(gr_pas$gene_id[gr_pas$anno.region == ""],
                                                                    from = tmp$gene_id,
                                                                    to = tmp$anno.region,warn_missing = F)
    
  }
  
  gr_pas$anno <-  plyr::mapvalues(gr_pas$anno,
                                  from = c("utr3","utr3.new","last.exon","last.intron","exon","intron","utr5"),
                                  to = c("3' UTR","3' UTR expand","Last Exon","Last Intron","Other Exon","Other Intron","5' UTR"),
                                  warn_missing = F)  
  
  gr_pas$gene_name <- plyr::mapvalues(gr_pas$gene_id, from = gene.anno$gene$gene_id,to = gene.anno$gene$gene_name,warn_missing = F) 
  gr_pas$peak = paste0(seqnames(gr_pas),":",start(gr_pas),"-",end(gr_pas))
  
  plotdata = data.frame(round(prop.table(table(gr_pas$anno))*100,digits = 2))

  plot.level = c("3' UTR","3' UTR expand","Last Exon","Last Intron","Other Exon","Other Intron","5' UTR")
  plotdata$Var1 <- factor(plotdata$Var1,levels = plot.level[plot.level %in% plotdata$Var1])
  
  plotdata <- plotdata[order(plotdata$Var1),]
  pdf(paste0(outpath,"/peakanno.pdf"),height = 4,width = 6)
  print(ggplot(plotdata, aes(x="", y=Freq, fill=Var1)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    scale_fill_discrete(name = "Type",labels=paste0(plotdata$Var1, " (",plotdata$Freq,"%)")) +
    theme_void())
  dev.off()

  return(gr_pas)
}


