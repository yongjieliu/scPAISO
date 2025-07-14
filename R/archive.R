##########  APA score ##########  
PasUsageScore.arv <- function(sc.pas.list = pas_sc_count,gr_pas,min.pas.cutoff = 2,pseudo_count = 0.1,threads = 24,return.gene = F){
  
  plan(sequential)
  plan(multisession, workers= threads)
  peak.index = data.frame(gr_pas)[,c("peak","gene_id", "rank")]
  gene.index <- table(peak.index$gene_id)
  gene.index <- names(gene.index)[gene.index > 1]
  peak.index <- peak.index[peak.index$gene_id %in% gene.index,]
  
  out.list <- list()
  pas_sc_score_mean <- c()
  name.index = c()
  
  for (i in seq_along(sc.pas.list)){
    sc.pas <- sc.pas.list[[i]]
    name.index = c(name.index,colnames(sc.pas))
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
      tmp3$score <- 1 - tmp3$score
      
      # order and export
      tmp1 <- rbind(data.frame(tmp2),data.frame(tmp3))
      tmp1 <- tmp1[order(tmp1$gene_id),]
      rownames(tmp1) <- tmp1$gene_id
      tmp1$gene_id <- NULL
      return(tmp1)
    })
    pas_sc_score <- do.call("cbind",pas_sc_score)
    pas_sc_score_mean <- c(pas_sc_score_mean,colMeans(pas_sc_score,na.rm = T))
    if (return.gene){
      out.list <- c(out.list,list(pas_sc_score))  
    }
  }
  
  names(pas_sc_score_mean) <- name.index
  #pas_sc_score <- do.call("cbind",out.list)
  
  plan(sequential)

  # colmean.scale
  pas_sc_score_mean.scale <- scale(pas_sc_score_mean) 
  # gene.div.mean
  # pas_sc_score_div <- sweep(pas_sc_score,1,rowMeans(pas_sc_score,na.rm = T),FUN = "/")
  # pas_sc_score_div <- colMeans(pas_sc_score_div,na.rm = T)
  
  # scale by each gene
  # pas_sc_score_scale <- t(scale(t(pas_sc_score)))
  # pas_sc_score_scale <- colMeans(pas_sc_score_scale,na.rm = T)
  
  pas.rank.score <- list()
  
  pas.rank.score[["matrix"]] <- out.list
  pas.rank.score[["apa.mean.scale"]] <- pas_sc_score_mean.scale
  pas.rank.score[["apa.mean"]] <- pas_sc_score_mean
  # pas.rank.score[["apa.div"]] <- pas_sc_score_div[colnames(sc.pas)]
  # pas.rank.score[["apa.scale"]] <- pas_sc_score_scale[colnames(sc.pas)]
  
  return(pas.rank.score)
}
