#' Model fitting
#' 
#' @param samplenames A vector containing sample names.
#' @param chrs A vector containing chromatin to use. 
#' @param outpath A character describing where the file will be saved.
#' @param gr_pas A range of PAS peak.
#' @param cache.dist A boolean value that determines whether cache distance file shoulde be used.
#' @param intronanno A ranges describing intron annotation.
#' @param threads The number of threads to be used for parallel computing.
#' @param intron.remove A boolean value that determines whether remove intron when caculate the distance between Read2 and PAS. 
#' @param sampling.frac The proportion of reads used for distance fitting.
#' 
#' @export

pas_fit_model <- function(samplenames,chrs,outpath,gr_pas = pas_peak,cache.dist = F,
                          intronanno = intronanno,threads=24,
                          intron.remove = T,sampling.frac = 4/100){
 
  if (!file.exists(paste0(outpath,"/model"))){
    dir.create(paste0(outpath,"/model"))
  }
  
  combinations <- expand.grid(samplenames, chrs)
  
  print("Reading pas df")
  
  if (!cache.dist){
    plan(sequential)
    plan(multisession, workers= threads)
    plot_dist_data <- future_lapply(1:nrow(combinations),function(i){
      
      samplename <- combinations$Var1[i]
      chr <- combinations$Var2[i]
      
      tmp <- readRDS(paste0(outpath,"/pas.df/",samplename,"/",chr,".plot_dist_data.rds"))
      if (!is.null(sampling.frac)){
        tmp <- tmp[sample(1:nrow(tmp),floor(nrow(tmp)*sampling.frac)),]
      }
      return(tmp)
    })
    
    plot_dist_data <- rbindlist(plot_dist_data)
    
    saveRDS(plot_dist_data,file = paste0(outpath,"/model/all.data.rds"))
    
    train.index = sample(1:nrow(plot_dist_data),100000)
    train_data <- plot_dist_data[train.index,]
    
    saveRDS(train_data,file = paste0(outpath,"/model/train_data.rds"))
    
    evaluate_data <- plot_dist_data[1:nrow(plot_dist_data) %ni% train.index,]
    evaluate_data$readid = paste0("readid",1:nrow(evaluate_data))
    saveRDS(evaluate_data,file = paste0(outpath,"/model/evaluate_data.rds"))
    
    plan(sequential)
  } else {
    plot_dist_data = readRDS(paste0(outpath,"/model/all.data.rds"))
    train_data = readRDS(paste0(outpath,"/model/train_data.rds"))
    evaluate_data = readRDS(paste0(outpath,"/model/evaluate_data.rds"))
  }
  
  print("Remove intron")
  train_data$start = train_data$start - 16
  train_data$end = train_data$end + 16
  
  train_data$start[train_data$strand == "+"] = train_data$end_r2[train_data$strand == "+"]
  train_data$end[train_data$strand == "-"] = train_data$pos_r2[train_data$strand == "-"]
  
  gr.index <- findOverlaps(intronanno, type=c("within"),as(train_data,"GRanges"))
  gr.index = gr.index[intronanno[from(gr.index)]$gene_id == train_data$gene_id[to(gr.index)],]
  
  intro_df <- data.table(index_df = to(gr.index),
                         intron_width = intronanno[from(gr.index)]$intron_width)
  intro_df <- intro_df %>% group_by(index_df) %>% summarise(intron_total = sum(intron_width))
  
  train_data$intron <- 0
  
  if (intron.remove){
    train_data$intron[intro_df$index_df] <- intro_df$intron_total
  }
  
  train_data$dist = train_data$end - train_data$start - train_data$intron
  
  train_data <- log2(train_data$dist + 1) + 1
  
  print("Model fitting")
  
  # fitdistrplus::fitdist
  fw <- fitdist(train_data, "weibull")
  fg <- fitdist(train_data, "gamma")
  fln <- fitdist(train_data, "lnorm")
  
  if (!file.exists(paste0(outpath,"/model/"))){
    dir.create(paste0(outpath,"/model/"))
  }
  
  saveRDS(fg,file = paste0(outpath,"/model/fit.gamma.rds"))
  saveRDS(fw,file = paste0(outpath,"/model/fit.weibull.rds"))
  saveRDS(fln,file = paste0(outpath,"/model/fit.fln.rds"))
  
  dcomp <- denscomp(list(fw, fln, fg), legendtext = c("weibull", "lognormal", "gamma"),
                    xlab = "-log2(dist + 16)", 
                    plotstyle = "ggplot", breaks = 80, addlegend = T)
  
  pdf(paste0(outpath,"fitdist.fit.pdf"))
  print(dcomp + ggplot2::ggtitle("Fits"))
  dev.off()
  
  print("Model evaluating")
  
  start(gr_pas) = start(gr_pas) - 16
  end(gr_pas) = end(gr_pas) + 16
  
  pas_peak_index <- setnames(data.frame(gr_pas)[,c("start","end","peak","gene_id","rank")],
                             c("map.start","map.end","map.peak","gene_id","rank"))  
  
  evaluate_data$peak = paste0(evaluate_data$seqnames,":",evaluate_data$start ,"-",evaluate_data$end)
  
  evaluate_data <- inner_join(evaluate_data,pas_peak_index, by = 'gene_id')
  
  evaluate_data$start <- evaluate_data$end_r2
  evaluate_data$end <- evaluate_data$map.end
  
  evaluate_data$start[evaluate_data$strand == "-"] <- evaluate_data$map.start[evaluate_data$strand == "-"]
  evaluate_data$end[evaluate_data$strand == "-"] <- evaluate_data$pos_r2[evaluate_data$strand == "-"]
  
  evaluate_data <- evaluate_data[evaluate_data$start <= evaluate_data$end,]
  
  gr.index <- findOverlaps(intronanno, type=c("within"),as(evaluate_data,"GRanges"))
  gr.index = gr.index[intronanno[from(gr.index)]$gene_id == evaluate_data$gene_id[to(gr.index)],]
  
  intro_df <- data.table(index_df = to(gr.index),
                         intron_width = intronanno[from(gr.index)]$intron_width)
  
  intro_df <- intro_df %>% group_by(index_df) %>% summarise(intron_total = sum(intron_width))
  
  evaluate_data$intron <- 0
  if (intron.remove){
    evaluate_data$intron[intro_df$index_df] <- intro_df$intron_total
  }

  evaluate_data$dist = evaluate_data$end - evaluate_data$start - evaluate_data$intron
  
  for (fitmodel in c("weibull", "lognormal", "gamma")) {
    if (fitmodel == "weibull"){
      evaluate_data$predict <- dweibull(log2(evaluate_data$dist + 1) + 1, coef(fw)[1], coef(fw)[2])
    } else if (fitmodel == "gamma") {
      evaluate_data$predict <- dgamma(log2(evaluate_data$dist + 1) + 1, coef(fg)[1], coef(fg)[2]) 
    } else if (fitmodel == "lognormal"){
      evaluate_data$predict <- dlnorm(log2(evaluate_data$dist + 1) +1, coef(fln)[1], coef(fln)[2]) 
    }

    saveRDS(evaluate_data,file = paste0(outpath,"/model/evaluate_data.",fitmodel,".rds"))
    
    plot_data <- evaluate_data %>%
      group_by(readid,gene_id) %>%
      filter(predict == max(predict))
    
    # cutoff.data <- lapply(seq(0, 0.4, by=0.01),function(i){
    #   plot_data_max <- plot_data[plot_data$predict >= i,]
    #   fdr <- table(plot_data_max$peak == plot_data_max$map.peak)/nrow(plot_data_max)
    #   return(c(nrow(plot_data_max)/nrow(plot_data),fdr[2]))
    # }) %>% Reduce("rbind",.)
    
    # plot_data <- plot_data %>%
    #   arrange(qname,gene_id, dist) %>%
    #   group_by(qname,gene_id) %>%
    #   mutate(dist.order = row_number())
    # plot_data$dist.order[plot_data$dist.order > 4] <- 4
    # ord.dis <- table(plot_dist_data$pas_order)/nrow(plot_dist_data) 
    # plot_data$predict2 <- plot_data$predict*ord.dis[plot_data$dist.order]
    # plot_data$predict2 <- as.numeric(round(plot_data$predict2,digits = 6))
    # plot_data_max <- plot_data %>%
    #   group_by(qname,gene_id) %>%
    #   filter(predict2 == max(predict2))
    
    plot_data_max <- plot_data[plot_data$peak == plot_data$map.peak,]
    
    fdr <- table(plot_data$peak == plot_data$map.peak)/nrow(plot_data)
    
    gene.stat = data.frame(table(plot_data$gene_id))
    
    gene.stat2 <- data.frame(table(plot_data_max$gene_id))
    
    gene.stat <- merge(gene.stat,gene.stat2,by = "Var1",all = T)
    gene.stat$Freq.y[is.na(gene.stat$Freq.y)] = 0
    gene.stat$stat = gene.stat$Freq.y/gene.stat$Freq.x
    
    grob <- grobTree(textGrob(paste0("Overall accuracy ",round(fdr[2],digits = 4)), x=0.1,  y=0.95, hjust=0,
                              gp=gpar(col="red", fontsize=13, fontface="italic")))
    
    
    pdf(paste0(outpath,"fdr.gene.dist.",fitmodel,".pdf"),height = 5,width = 5)
    print(ggplot(gene.stat,aes(x = stat)) + 
            geom_histogram(aes(y=..density..),
                           binwidth=0.01,
                           colour="black", fill="white") + 
            geom_density(alpha=.2, fill="#FF6666")+ xlab("Accuracy") +
            ggtitle("Distribution of gene accuracy level") + theme(plot.title = element_text(hjust = 0.5)) +
            annotation_custom(grob)
    )
    dev.off()
  } 
  
  ## plot distance dis
  start(gr_pas) = start(gr_pas) + 16
  end(gr_pas) = end(gr_pas) - 16
  
  saveRDS(gr_pas,file = paste0(outpath,"/model/plot.pas.gr.rds"))
  plot_pas_stat(plotdata = plot_dist_data,gr_pas = gr_pas)
}

pas_bulk <- function(samplenames,chrs,outpath,cell.metadata = cell.metadata,threads=24){
  
  combinations <- expand.grid(samplenames, chrs)
  
  plan(sequential)
  plan(multisession, workers= threads)
  
  print("Counting")
  bam_df.list <- future_lapply(1:nrow(combinations),function(i){
    
    samplename <- combinations$Var1[i]
    chr <- combinations$Var2[i]
    
    bam_df <- readRDS(paste0(outpath,"/pas.df/",samplename,"/",chr,".plot_dist_data.rds"))

    bam_df$peak = paste0(bam_df$seqnames,":",bam_df$start,"-",bam_df$end)
    bam_df$tag.CB = paste0(samplename,"_",bam_df$tag.CB)
    bam_df <- bam_df[bam_df$tag.CB %in% cell.metadata$cb,]
    bam_df$celltype <- plyr::mapvalues(bam_df$tag.CB,from = cell.metadata$cb,to = cell.metadata$celltype,warn_missing = F)
    
    bam_df <- suppressMessages(
      bam_df[,c("celltype","peak")] %>%
        group_by(celltype,peak) %>%
        summarise(count=n()))
  })
  plan(sequential)
  
  bam_df.list <- rbindlist(bam_df.list)
  bam_df.list = bam_df.list %>%
    group_by(celltype,peak) %>%
    summarise(count = sum(count))
  
  bam_df.list <- dcast(bam_df.list ,peak ~ celltype)
  rownames(bam_df.list) <- bam_df.list$peak
  bam_df.list$peak <- NULL
  bam_df.list[is.na(bam_df.list)] <- 0
  
  saveRDS(bam_df.list,file = paste0(outpath,"/pas_count.bulk.rds"))
  return(bam_df.list)
}

#' PA isoform quantification
#' @param gr_pas A range of PAS peak.
#' @param samplenames A vector containing sample names.
#' @param chrs A vector containing chromatin to use. 
#' @param outpath A character describing where the file will be saved.
#' @param geneanno A ranges describing gene annotation.
#' @param intronanno A ranges describing intron annotation.
#' @param intron.remove A boolean value that determines whether remove intron when caculate the distance between Read2 and PAS. 
#' @param cell.metadata A data.frame containing: cb, cell barcode id; celltype, cell type annotation; samplename, sample name.
#' @param fitmodel The distribution used for distance fitting. Options include "weibull", "gamma", or "lognormal". 
#' @param threads The number of threads to be used for parallel computing.
#' 
#' @export

pas_predict <- function(gr_pas = pas_peak,samplenames,chrs,outpath,
                        geneanno = geneanno,intronanno = intronanno,
                        intron.remove = T,remove.bad = "dist",
                        cell.metadata = cell.metadata,fitmodel = "gamma",threads=24){
  
  start(gr_pas) = start(gr_pas) - 16
  end(gr_pas) = end(gr_pas) + 16
  
  cell.index <- 1:nrow(cell.metadata)
  names(cell.index) <- cell.metadata$cb
  peak.index <- 1:length(gr_pas)  
  names(peak.index) <- gr_pas$peak
  
  pas_peak_index <- setnames(data.frame(gr_pas)[,c("start","end","peak","gene_id","rank")],c("map.start","map.end","map.peak","gene_id","rank"))
  
  if (fitmodel == "weibull"){
    fit_para <- readRDS(paste0(outpath,"/model/fit.weibull.rds"))
  } else if (fitmodel == "gamma") {
    fit_para <- readRDS(paste0(outpath,"/model/fit.gamma.rds"))
  } else if (fitmodel == "lognormal"){
    fit_para <- readRDS(paste0(outpath,"/model/fit.fln.rds"))
  }
  
  combinations <- expand.grid(samplenames, chrs)
  
  plan(sequential)
  plan(multisession, workers= threads)
  
  index = unique(gr_pas$callpeak)
  
  if (remove.bad == "dist"){
    bad.reads <- lapply(index, function(x){
      readRDS(paste0(outpath,"pas.gr/gr_pas_peak_remove.",x,".rds"))
    })  %>% Reduce("c",.)
    start(bad.reads)[as.logical(strand(bad.reads) == "+")] <- start(bad.reads)[as.logical(strand(bad.reads) == "+")] - 400
    end(bad.reads)[as.logical(strand(bad.reads) == "-")] <- end(bad.reads)[as.logical(strand(bad.reads) == "-")] + 400
  }
  
  print("Predicting")
  bam_df.list <- future_lapply(1:nrow(combinations),function(i){
    
    samplename <- combinations$Var1[i]
    chr <- combinations$Var2[i]

    bam_df <- readRDS(paste0(outpath,"/bam2df/",samplename,"/",chr,".r2bam2df.rds"))
    
    qc.export <- paste0("R2 reads: ",nrow(bam_df),"\n")
    
    bam_df$seqnames = chr
    bam_df$start <- bam_df$pos_r2
    bam_df$end <- bam_df$end_r2
    bam_df$strand = bam_df$strand_r2
    
    gr.index = findOverlaps(as(bam_df,"GRanges"),geneanno)
    
    overlap_length <- pintersect(as(bam_df,"GRanges")[queryHits(gr.index)], geneanno[subjectHits(gr.index)])
    
    bam_df <- cbind(bam_df[queryHits(gr.index),],
                    data.frame(geneanno)[subjectHits(gr.index),c("gene_id","gene_name")])
    
    bam_df$overlap.width = width(overlap_length)/(bam_df$end_r2 - bam_df$pos_r2 + 1)
    
    bam_df <- data.table(bam_df)
    bam_df <- bam_df %>%  group_by(qname) %>% filter(overlap.width == max(overlap.width))
    
    qc.export <- paste0(qc.export,
                        "Gene mapped reads: ",nrow(bam_df),"\n")
    
    if ("tag.UB" %in% colnames(bam_df)){
      #bam_df$tag.UB <- paste0(bam_df$rname_r2,"-",bam_df$pos_r2,"-",bam_df$end_r2)
      bam_df <- bam_df[!duplicated(bam_df[,c("tag.CB","tag.UB","gene_id")]),]
    }
    
    qc.export <- paste0(qc.export,
                        "Dedup reads: ",nrow(bam_df),"\n")
    
    if (remove.bad == "dist"){
      bad.reads.index <- findOverlaps(as(bam_df,"GRanges"),bad.reads,type=c("within"))
      bam_df <- bam_df[1:nrow(bam_df) %ni% queryHits(bad.reads.index),]
      qc.export <- paste0(qc.export,"MissPriming removed R2 reads: ",nrow(bam_df),"\n")
    }
    
    bam_df <- inner_join(bam_df,pas_peak_index, by = 'gene_id')
    
    bam_df$start <- bam_df$end_r2
    bam_df$end <- bam_df$map.end
    
    bam_df$start[bam_df$strand == "-"] <- bam_df$map.start[bam_df$strand == "-"]
    bam_df$end[bam_df$strand == "-"] <- bam_df$pos_r2[bam_df$strand == "-"]
    
    bam_df <- bam_df[bam_df$end >= bam_df$start,]
    
    if (nrow(bam_df) > 3e6){
      bam_df_split <- split(bam_df, sample(1:ceiling(nrow(bam_df)/2e6),nrow(bam_df),replace = T))
      
      bam_df <- lapply(1:length(bam_df_split),function(x){
        tmp = bam_df_split[[x]]
        gr.index <- findOverlaps(intronanno, as(tmp,"GRanges"),type="within")
        gr.index = gr.index[intronanno[from(gr.index)]$gene_id == tmp$gene_id[to(gr.index)],]
        intro_df <- data.table(index_df = to(gr.index),
                               intron_width = intronanno[from(gr.index)]$intron_width)
        intro_df <- intro_df %>% group_by(index_df) %>% summarise(intron_total = sum(intron_width))
        tmp$intron <- 0
        if (intron.remove){
          tmp$intron[intro_df$index_df] <- intro_df$intron_total
        }
        tmp$dist = tmp$end - tmp$start - tmp$intron
        tmp <- tmp[tmp$dist < 10^5,]
        tmp <- tmp[is.numeric(tmp$dist),]
        return(tmp)
      }) %>% Reduce("rbind",.)
      
      rm(bam_df_split)
    } else {
     
        gr.index <- findOverlaps(intronanno, as(bam_df,"GRanges"),type="within")
        gr.index = gr.index[intronanno[from(gr.index)]$gene_id == bam_df$gene_id[to(gr.index)],]
        intro_df <- data.table(index_df = to(gr.index),
                               intron_width = intronanno[from(gr.index)]$intron_width)
        intro_df <- intro_df %>% group_by(index_df) %>% summarise(intron_total = sum(intron_width))
        bam_df$intron <- 0
        if (intron.remove){
          bam_df$intron[intro_df$index_df] <- intro_df$intron_total
        }
        bam_df$dist = bam_df$end - bam_df$start - bam_df$intron
        bam_df <- bam_df[bam_df$dist < 10^5,]
        bam_df <- bam_df[is.numeric(bam_df$dist),]
    }
    
    if (fitmodel == "weibull"){
      bam_df$predict <- dweibull(log2(bam_df$dist + 1) + 1, fit_para$estimate[1],fit_para$estimate[2])
    } else if (fitmodel == "gamma") {
      bam_df$predict <- dgamma(log2(bam_df$dist + 1) + 1, fit_para$estimate[1],fit_para$estimate[2])
    } else if (fitmodel == "lognormal"){
      bam_df$predict <- dlnorm(log2(bam_df$dist + 1) + 1, fit_para$estimate[1],fit_para$estimate[2])
    }
    
    bam_df <- bam_df %>%
      group_by(qname,gene_id) %>%
      filter(predict == max(predict))
    
    
    
    qc.export <- paste0(qc.export,
                        "PAS mapped: ",nrow(bam_df))
    
    writeLines(qc.export, paste0(outpath,"stat.qc/",samplename,"/stat.predict.",chr,".txt"))
    
    bam_df <-  suppressMessages(
      bam_df[,c("tag.CB","map.peak")] %>%
      group_by(tag.CB,map.peak) %>%
      summarise(count=n())
      )
    
    bam_df$tag.CB <- paste0(samplename,"_",bam_df$tag.CB)
    bam_df <- bam_df[bam_df$tag.CB %in% cell.metadata$cb,]
    
    bam_df$tag.CB <- cell.index[bam_df$tag.CB]
    bam_df$map.peak <- peak.index[bam_df$map.peak]
    saveRDS(bam_df,file = paste0(outpath,"/bam2df/",samplename,"/predict.",chr,".rds"))
    
    return(bam_df)
  })
  plan(sequential)
  
  print("To dgcMatrix")
  
  bam_df.list <- rbindlist(bam_df.list)
  
  bam_df.list <- sparseMatrix(i = bam_df.list$map.peak, 
                              j = bam_df.list$tag.CB, 
                              x = bam_df.list$count)
  
  rownames(bam_df.list) <- gr_pas$peak[1:nrow(bam_df.list)]
  colnames(bam_df.list) <- cell.metadata$cb[1:ncol(bam_df.list)]
  # rownames(bam_df.list) = plyr::mapvalues(rownames(bam_df.list),
  #                                         from = peak.index,to = names(peak.index))
  # colnames(bam_df.list) = plyr::mapvalues(rownames(bam_df.list),
  #                                         from = cell.index,to = names(cell.index))
  
  saveRDS(bam_df.list,file = paste0(outpath,"/pas_count.rds"))
  return(bam_df.list)
}

########## imputation ##########

pas_impute <- function(matDR = sc.pas.merged@reductions$pca@cell.embeddings,
                       mat = sc.pas.merged@assays$RNA@counts,
                       samplenames,threads = length(samplenames),
                       k = 15,td = 10,ka = 4,epsilon = 1){
  
  blocks <- lapply(samplenames, function(x){
    matDR[grepl(paste0(x,"_"),rownames(matDR)),]
  })
  names(blocks) = samplenames
  
  plan(sequential)
  plan(multisession, workers= threads)
  
  blockList <- future_lapply(seq_along(blocks), function(x){
    
    if (!file.exists(paste0(outpath,"/imputation/"))){
      dir.create(paste0(outpath,"/imputation/"))
    }
    
    mat.sub = blocks[[x]]
    ix <- rownames(mat.sub)
    Nx <- length(ix)
    
    #Compute KNN
    knnObj <- nabor::knn(data = mat.sub, query = mat.sub, k = k)
    knnIdx <- knnObj$nn.idx
    knnDist <- knnObj$nn.dists
    rm(knnObj)
    
    if(ka > 0){
      knnDist <- knnDist / knnDist[,ka]
    }
    
    if(epsilon > 0){
      W <- Matrix::sparseMatrix(rep(seq_len(Nx), k), c(knnIdx), x=c(knnDist), dims = c(Nx, Nx))
    } else {
      W <- Matrix::sparseMatrix(rep(seq_len(Nx), k), c(knnIdx), x=1, dims = c(Nx, Nx)) # unweighted kNN graph
    }
    W <- W + Matrix::t(W)
    
    #Compute Kernel
    if(epsilon > 0){
      W@x <- exp(-(W@x / epsilon^2))
    }
    
    #Markov normalization
    W <- W / Matrix::rowSums(W)
    
    #Initialize Matrix
    Wt <- W
    
    #Computing Diffusion Matrix
    for(i in seq_len(td)){
      Wt <- Wt %*% W
    }
    rownames(Wt) <- ix
    colnames(Wt) <- ix
    
    rm(knnIdx)
    rm(knnDist)
    rm(W)
    gc()
    
    saveRDS(Wt,file = paste0(outpath,"/imputation/block", x))
    return(Wt)
  })
  
  mat.impute <- future_lapply(seq_along(blockList), function(x){
    tmp <- Matrix::t(as.matrix(blockList[[x]]) %*% Matrix::t(mat[, colnames(blockList[[x]]), drop = FALSE]))
    saveRDS(tmp,file = paste0(outpath,"/imputation/block",x,".impute.rds"))
    return(tmp)
  })
  names(mat.impute) = samplenames
  return(mat.impute)
}


#' Model fitting
#' 
#' @param sc.pas PAS count matrix: could be PAS peak X cell matrix or PAS peak X celltype.
#' @param gr_pas A range of PAS peak.
#' @param min.pas.cutoff If the rowSums(PAS count) of a PAS peak is less than this value, it will not be used for calculation.
#' @param pseudo_count The number of count that will be add to PAS count.
#' @param threads The number of threads to be used for parallel computing.
#' @param out.matrix A boolean value that determines whether  to export the APA rank score of each 3' UTR region in each cell or celltype. 
#' 
#' @export

PasUsageScore <- function(sc.pas = pas_sc_count,gr_pas,min.pas.cutoff = 5,pseudo_count = 0.1,threads = 24,by = "utr",out.matrix = F){
  
  plan(sequential)
  plan(multisession, workers= threads)
  peak.index = data.frame(gr_pas)[,c("peak","gene_id","anno", "anno.region","rank")]
  peak.index = peak.index[peak.index$peak %in% rownames(sc.pas),]
  
  if (by == "all"){
    gene.index <- table(peak.index$gene_id)
    gene.index <- names(gene.index)[gene.index > 1]
    
    peak.index <- peak.index[peak.index$gene_id %in% gene.index,]
  } else if (by == "utr"){
    
    peak.index = peak.index[peak.index$anno %in% c("3' UTR", "3' UTR expand","Last Exon"),]
    peak.index$gene_id = NULL
    gene.index <- table(peak.index$anno.region)
    gene.index <- names(gene.index)[gene.index > 1]
    
    peak.index <- peak.index[peak.index$anno.region %in% gene.index,]
    colnames(peak.index)[colnames(peak.index) == "anno.region"] = "gene_id"
  }
  
  
  if (out.matrix){
      pas_sc_score <-  future_lapply(colnames(sc.pas),function(x){
        
        tmp1 <- data.frame(peak = rownames(sc.pas), value = sc.pas[,x])
        tmp1 <- merge(tmp1,peak.index)
        
        tmp1 <- tmp1[order(tmp1$gene_id),]
        
        ## filter low identified genes
        tmp2 <- tmp1 %>%
          group_by(gene_id) %>%
          mutate(score = sum(value))  %>%
          filter(score < min.pas.cutoff)  %>%
          summarise(score = sum(value))
        tmp2$score <- NA
        
        ## apa score
        tmp3 <- tmp1 %>%
          group_by(gene_id) %>%
          mutate(score = sum(value))  %>%
          filter(score  >= min.pas.cutoff)  %>%
          group_by(gene_id) %>%
          arrange(rank) %>% mutate(rank2 = row_number()) %>% #desc(rank)
          #mutate(score = (max(rank)-rank)*(value + pseudo_count)/sum(value+ pseudo_count)/max(rank))  %>%
          mutate(score = (rank2 -1)/(max(rank2) -1)*(value + pseudo_count)/sum(value + pseudo_count))  %>%
          summarise(score = sum(score))
        
        # order and export
        tmp1 <- rbind(data.frame(tmp2),data.frame(tmp3))
        tmp1 <- tmp1[order(tmp1$gene_id),]
        rownames(tmp1) <- tmp1$gene_id
        tmp1$gene_id <- NULL
        return(tmp1)
      })
      
      pas_sc_score <- do.call("cbind",pas_sc_score)
      colnames(pas_sc_score) = colnames(sc.pas)
      pas_sc_score_mean <- colMeans(pas_sc_score,na.rm = T)
    
      pas_sc_score <- list(score.mt = pas_sc_score,
                           score.mean = pas_sc_score_mean)
      
      plan(sequential)

  }
    
  if (!out.matrix){
    pas_sc_score <-  future_lapply(colnames(sc.pas),function(x){
      
      tmp1 <- data.frame(peak = rownames(sc.pas), value = sc.pas[,x])
      tmp1 <- merge(tmp1,peak.index)
      
      ## apa score
      tmp1 <- tmp1 %>%
        group_by(gene_id) %>%
        mutate(score = sum(value))  %>%
        filter(score  >= min.pas.cutoff)  %>%
        group_by(gene_id) %>%
        arrange(rank) %>% mutate(rank2 = row_number()) %>% #desc(rank)
        #mutate(score = (max(rank)-rank)*(value + pseudo_count)/sum(value+ pseudo_count)/max(rank))  %>%
        mutate(score = (rank2 -1)/(max(rank2) -1)*(value + pseudo_count)/sum(value + pseudo_count))  %>%
        summarise(score = sum(score))
      
      gc()
      return(mean(tmp1$score))
    })
    
    pas_sc_score <- do.call("c",pas_sc_score)
    
    plan(sequential)
    names(pas_sc_score) <- colnames(sc.pas)
  }
  
  return(pas_sc_score)
}


