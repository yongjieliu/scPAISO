########### plot ###################
enrichplot <- function(enrichout = enrichout,showsize = 10,showtitle = "Enrichments",xgap = 0.1,xvalue = "p.adjust",xline = 2){
  func.enrichout <- enrichout@result[1:showsize,]
  func.enrichout$plotvalue <- -log10(func.enrichout[,xvalue])
  func.enrichout$x <- paste0("level_",1:nrow(func.enrichout)) 
  func.enrichout$x <- factor(func.enrichout$x,levels = func.enrichout$x)
  
  func.enrichplot <- ggplot(func.enrichout, aes(x=plotvalue,y=x,fill = Count)) +
    geom_vline(xintercept = xline, linetype="dotted", 
               color = "red", size=1.5) +
    geom_bar(stat = "identity",width = 0.8) +
    scale_fill_gradient(low = "green", high = "red") +
    geom_text(aes(x=xgap, label=Description),hjust=0) + 
    scale_y_discrete(limits=rev) + 
    ggtitle(showtitle) + theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(expand = c(0,0)) +
    ylab("") + xlab(paste0("-log10(",xvalue,")")) + 
    theme(axis.text.y =element_blank(),
          #axis.text.x =element_blank(),
          axis.ticks.y = element_blank(),
          #axis.ticks.x = element_blank(),
          #strip.text.x = element_blank(),
          panel.spacing = unit(1.5, "lines"))
  return(func.enrichplot)
}

plot_base_freq <- function(gr_pas = pas_peak,extend_base = 100,BSg = BSgenome.Hsapiens.UCSC.hg38,
                           outpath, outname = "pas.100bp.freq.pdf",ylim = 0.8){
  if (!methods::is(gr_pas, "GRanges")) stop("gr_pas  GRanges ")
  if (!methods::is(BSg, "BSgenome")) stop("BSg  BSgenome ")
  if (!is.numeric(extend_base) || extend_base <= 0) stop("extend_base ")
  if (!is.numeric(ylim) || ylim <= 0) stop("ylim ")
  if (missing(outpath) || !dir.exists(outpath)) stop("outpath ")
  plot_pas <- gr_pas
  start(plot_pas) <- gr_pas$WhichMaxs - extend_base
  end(plot_pas) <- gr_pas$WhichMaxs + extend_base
  
  plot_pas <- getSeq(BSg, plot_pas)

  counts <- consensusMatrix(plot_pas, as.prob = T,baseOnly=T)
  counts <- data.frame(t(counts))
  counts$other <- 1:nrow(counts)
  counts$other <- counts$other - extend_base - 1
  counts <- reshape2::melt(counts,id = "other")

  myplot <- ggplot(counts,aes(x = other,y = value, group=variable, color=variable)) +
    geom_line() + xlab("Position") + ylab("Base fraction") + ylim(0,ylim)
    guides(color=guide_legend(title="Base"))
  
  pdf(file.path(outpath,outname),height = 5,width = 4)
  print(myplot)
  dev.off()
}

plot_pas_stat <- function(plotdata,gr_pas = gr_pas,outpath){ 
  if (missing(outpath) || !dir.exists(outpath)) stop("outpath ")

  print("Ploting")
  pdf(file.path(outpath,"pas.dendity.R2PAS.pdf"),height = 5,width = 8)
  print(ggplot(plotdata, aes(x=log2(R2PAS+1))) + 
          geom_histogram(aes(y=..density..),
                         binwidth=.2,
                         colour="black", fill="white") +
          geom_density(alpha=.2, fill="#FF6666"))
  dev.off()
  
  #quantile(plotdata$R2PAS,seq(0,1,0.05))
  
  plotdata$dist  <- plotdata$R2PAS + plotdata$end -  plotdata$PAS
  plotdata$dist[plotdata$strand == "-"]  <- 
    c(plotdata$R2PAS + plotdata$PAS -  plotdata$start)[plotdata$strand == "-"]
  
  pdf(file.path(outpath,"pas.dendity.R2Peak.max.pdf"),height = 5,width = 8)
  print(ggplot(plotdata, aes(x=log2(dist + 1))) + 
          geom_histogram(aes(y=..density..),
                         binwidth=.2,
                         colour="black", fill="white") +
          geom_density(alpha=.2, fill="#FF6666") + xlab("log2(R2PASmax + 1)"))
  dev.off()
  
  #quantile(plotdata$dist,seq(0,1,0.05))
  
  plotdata$start[plotdata$strand == "+"] <- plotdata$end_r2[plotdata$strand == "+"]
  plotdata$end[plotdata$strand == "-"] <- plotdata$pos_r2[plotdata$strand == "-"]
  
  plotdata = plotdata[plotdata$end >= plotdata$start,]
  myindex = findOverlaps(as(plotdata,"GRanges"),gr_pas)
  myindex = myindex[plotdata$gene_id[from(myindex)] == gr_pas$gene_id[to(myindex)],]
  plotdata$pas_order <- as.numeric(table(queryHits(myindex))[1:nrow(plotdata)]) 
  
  write.csv(table(plotdata$pas_order),file =  file.path(outpath,"model","sta.closest.csv"))
  pdf(file.path(outpath,"pas.dendity.pas.order.pdf"),height = 5,width = 4)
  print(ggplot(plotdata, aes(x=plotdata$pas_order -1)) + 
          geom_histogram(binwidth = 1,colour="black", fill="#FF6666") +
          xlim(NA,7) +
          xlab("Skipped PASs"))
  dev.off()
}

plot_aggr_cor <- function(pas_count_predict = pas_sc_count,
                          sc.rna = sc.rna@assays$RNA@counts,
                          cell.metadata,geneanno,threads =24){
  
  cell.id <- intersect(colnames(sc.rna),colnames(pas_count_predict))
  plan(sequential)
  plan(multisession, workers= threads)
  cor.aggr <- future_lapply(cell.id,function(i){
    cor.data <- pas_count_predict[,i]
    cor.data <- aggregate(cor.data,by = list(pas_peak$gene_name),sum)
    rownames(cor.data) <- cor.data$Group.1
    cor.gene <- sc.rna[,i]
    cor.gene <- cor.gene[names(cor.gene) %in% cor.data$Group.1]
    cor.data <- cor.data[cor.data$Group.1 %in% names(cor.gene),]
    
    cor(cor.gene[cor.data$Group.1],cor.data$x)
  })
  plan(sequential)
  cor.aggr <- do.call("c",cor.aggr)
  cor.aggr <- data.frame(cor = cor.aggr,cb = cell.id)
  cor.aggr$Gene.nCount <- colSums(sc.rna)[cell.id]
  cor.aggr$PAS.nCount <- colSums(pas_count_predict)[cell.id]
  
  pdf(paste0(outpath,"cor.aggr.pdf"),height = 10,width = 10)
  print(ggplot(cor.aggr,aes(x = PAS.nCount,y = Gene.nCount)) + geom_point() + 
          xlab("Total PAS counts") + ylab("Total mRNA counts"))
  print(ggplot(cor.aggr,aes(x = "Cells",y = cor)) + 
          geom_boxplot(outlier.colour="black", outlier.shape=16,notch=FALSE))
  dev.off()
  
  return(cor.aggr)
}

feature.plot <- function(object = sc.merged,feature = x,plot.type="vln"){
  DefaultAssay(object) = "RNA"
  plot1 <- FeaturePlot(object, features = feature) + coord_fixed() 
  
  plot.peak = pas_peak$peak[pas_peak$gene_name == feature]
  
  DefaultAssay(object) = "polyA"
  plot.peak <- plot.peak[plot.peak %in% rownames(object)]
  plot2 <- lapply(plot.peak, function(x){
    if (plot.type =="vln"){
      myplot = VlnPlot(object, features = x,pt.size = 0) + NoLegend()
    }
    if (plot.type =="umap"){
      myplot = FeaturePlot(object, features = x) + coord_fixed()
    }
    myplot
  })
  plot1 <- c(list(plot1),plot2)
  png(file = paste0(outpath,feature,".pas.png")  ,height =  ceiling(length(plot1)/3)*1200 ,width = 3600,res = 200)
  print(plot_grid(plotlist = plot1,ncol = 3))
  dev.off()
}

llrplot <- function(object,gene,bulk.df = NULL,region = NULL,bw.path = big.path,
                    smooth = 100,hei_rel = c(4,1),height = 4,width = 10,
                    assays = "ISO",group.by = "ident",bp.expand  = 2000){
  if (is.null(bulk.df)){
    bulk.df = AggregateExpression(object,assays = assays,group.by = group.by)[[assays]]
    bulk.df <- data.frame(bulk.df)
  }
  
  group.id = colnames(bulk.df)
  bulk.df$peak <- rownames(bulk.df)
  
  llr.df = pas_peak[pas_peak$gene_name == gene,c("peak","anno","anno.region")]
  llr.df <- merge(llr.df,bulk.df)
  
  llr.df <- llr.df[order(llr.df$start,decreasing = F),]
  
  if (!is.null(region)){
    llr.df <- llr.df[llr.df$start > as.numeric(strsplit(region,"-")[[1]][2]),]
    llr.df <- llr.df[llr.df$start < as.numeric(strsplit(region,"-")[[1]][3]),]
  }
  
  plot.region = paste0(llr.df$seqnames[1],"-",
                       min(llr.df$start) - bp.expand,"-",
                       max(llr.df$start) + bp.expand)
  
  gene.str = llr.df$strand[1]
  
  llr.df <- reshape2::melt(llr.df[,c("start",group.id)],id = c("start"))
  if (gene.str == "+"){
    llr.df <- llr.df[order(llr.df$start,decreasing = T),]
    plot.direction = "vh"
  } else {
    llr.df <- llr.df[order(llr.df$start,decreasing = F),]
    plot.direction = "hv"
  }
  
  llr.df <- llr.df %>%
    group_by(variable) %>% 
    mutate(ratio = value/sum(value))  %>% 
    mutate(cumsum = cumsum(ratio))
  
  covplot = BigwigTrack(
    region = plot.region,type = "coverage",
    bigwig = bw.path,bigwig.scale = "separate",
    extend.upstream = 0, extend.downstream = 0,
    smooth = smooth,downsample.rate = 1 #,ymax = "q99.9"
  )  + scale_fill_manual(values = hue_pal()(length(big.path)))
  
  peak_plot <- PeakPlot(
    object = object,
    region = plot.region
  )
  
  gene_plot <- AnnotationPlot(
    object = object,
    region = plot.region )
  
  llr_plot <- ggplot(llr.df, aes(x = start, y = cumsum,color = variable)) +
    geom_step(direction = plot.direction) + #xlim(min(llr.df$x) - expand,max(llr.df$x) + expand) + 
    geom_point() + guides(color=guide_legend(nrow=10, byrow=F)) + NoLegend() +
    xlim(as.numeric(strsplit(plot.region,"-")[[1]][2]),
         as.numeric(strsplit(plot.region,"-")[[1]][3])) + 
    xlab("") + 
    theme_classic() +
    theme(
      legend.title = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    ) 
  # scale_y_reverse()
  pdf(paste0(outpath,gene,".llr.pdf"),height = height,width = width)
  print(CombineTracks(
    plotlist = list(covplot,llr_plot,peak_plot,gene_plot),
    heights = hei_rel))
  dev.off()
}

########### export ###################

# normalize by celltype R1 profile or merged R1 profile
export_pas_bigwig <- function(samplenames,chrs,predict.bulk,outpath,cell.metadata,track.nor = "all", # "celltype"
                               chrom.size,threads=24){
  
  if (!file.exists(paste0(outpath,"celltype.bw"))){
    dir.create(paste0(outpath,"celltype.bw"))
  }
  
  pseudo.count = 0.2
  combinations <- expand.grid(samplenames, chrs)
  
  plan(sequential)
  plan(multisession, workers= threads)
  
  colnames(predict.bulk) = gsub("-|_",".",colnames(predict.bulk))
  cell.metadata$celltype = gsub("-|_",".",cell.metadata$celltype)
  
  print("Read raw peak counts")
  bam_df.list <- future_lapply(1:nrow(combinations),function(i){
    
    samplename <- combinations$Var1[i]
    chr <- combinations$Var2[i]
    
    bam_df <- readRDS(paste0(outpath,"/pas.df/",samplename,"/",chr,".plot_dist_data.rds"))
    
    bam_df$peak = paste0(bam_df$seqnames,":",bam_df$start,"-",bam_df$end)
    bam_df$tag.CB = paste0(samplename,"_",bam_df$tag.CB)
    bam_df <- bam_df[bam_df$tag.CB %in% cell.metadata$cb,]
    bam_df$celltype <- plyr::mapvalues(bam_df$tag.CB,from = cell.metadata$cb,to = as.character(cell.metadata$celltype) ,warn_missing = F)
    
    bam_df <- suppressMessages(
      bam_df %>%
        group_by(seqnames,start,end,strand,celltype,peak) %>%
        summarise(count=n())
      )
    
    return(bam_df)
  })
  
  bam_df.list <- rbindlist(bam_df.list)
  bam_df.list = suppressMessages(
    bam_df.list %>%
    group_by(seqnames,start,end,strand,celltype,peak) %>%
    summarise(count = sum(count)))
  
  if (track.nor == "celltype"){
    bam_df.list <- data.table::dcast(bam_df.list ,seqnames + start+end+strand+peak ~ celltype)
    bam_df.list[is.na(bam_df.list)] <- 0
    bam_df.list = data.frame(bam_df.list)
    rownames(bam_df.list) <- bam_df.list$peak
    bam_df.list$peak <- NULL
    
    peaks.used = rownames(predict.bulk)[rownames(predict.bulk) %in% rownames(bam_df.list)]
    cellused = colnames(predict.bulk)[colnames(predict.bulk) %in% colnames(bam_df.list)]
    
    bam_df.list <- bam_df.list[peaks.used,]
    
    peak.width = bam_df.list$end - bam_df.list$start + 1
    
    sc.plus = bam_df.list
    sc.plus[,cellused] = round((data.frame(predict.bulk[peaks.used,cellused]))/(bam_df.list[,cellused] + peak.width*pseudo.count),digits = 5)
    
  } else if (track.nor == "all"){

    bam_df.list = bam_df.list %>%
      group_by(seqnames,start,end,strand,peak) %>%
      summarise(count = sum(count))
    
    bam_df.list = data.frame(bam_df.list)
    rownames(bam_df.list) <- bam_df.list$peak
    
    peaks.used = rownames(predict.bulk)[rownames(predict.bulk) %in% rownames(bam_df.list)]
    
    bam_df.list <- bam_df.list[peaks.used,] # raw.count.R1
    
    sc.plus = bam_df.list
    # FC to multiply, celltype fc = cell.predict.count/raw.count.R1
    for (x in colnames(predict.bulk)) {
      sc.plus[,x] = round((predict.bulk[peaks.used,x] + 1)/(bam_df.list$count + 1),digits = 5)
    }
  }
  
  saveRDS(sc.plus,file = paste0(outpath,"/pas.df/fc2multiply.rds"))
  
  print("Peak counts to covrage")
  
  cov.data.list <- future_lapply(1:nrow(combinations),function(x){
    
    cov.list <- list() # samplename.chr.celltype
    samplename <- combinations$Var1[x]
    chr <- combinations$Var2[x]
    
    # R1 and R2
    tmp <- readRDS(paste0(outpath,"/pas.df/",samplename,"/",chr,".plot_dist_data.rds"))
    tmp <- data.frame(tmp)
    tmp$tag.CB <- paste0(samplename,"_",tmp$tag.CB)
    tmp <- tmp[tmp$tag.CB %in% cell.metadata$cb,]
    tmp$celltype <- plyr::mapvalues(tmp$tag.CB,from = cell.metadata$cb,
                                    to = as.character(cell.metadata$celltype),warn_missing = F)
    tmp$peak = paste0(tmp$seqnames,":",tmp$start,"-",tmp$end)
    tmp <- tmp[tmp$peak %in% rownames(sc.plus),]
    
    tmp.exp <- tmp[tmp$strand == "-",]
    tmp.exp$start <- tmp.exp$PAS
    tmp.exp$end <- tmp.exp$PAS
    tmp.exp <-  as(tmp.exp,"GRanges")
    tmp.exp@seqinfo@seqlengths <- chrom.size[tmp.exp@seqinfo@seqnames,]$seqlengths
    
    cov.list[[paste0(samplename,".",chr,".R1.all.neg")]] = coverage(tmp.exp)
    
    tmp.exp <- tmp[tmp$strand == "+",]
    tmp.exp$start <- tmp.exp$PAS
    tmp.exp$end <- tmp.exp$PAS
    tmp.exp <-  as(tmp.exp,"GRanges")
    tmp.exp@seqinfo@seqlengths <- chrom.size[tmp.exp@seqinfo@seqnames,]$seqlengths
    
    cov.list[[paste0(samplename,".",chr,".R1.all.pos")]] = coverage(tmp.exp)
    
    tmp.exp <- tmp[tmp$strand == "-",]
    tmp.exp$start <- tmp.exp$pos_r2
    tmp.exp$end <- tmp.exp$pos_r2
    tmp.exp <-  as(tmp.exp,"GRanges")
    tmp.exp@seqinfo@seqlengths <- chrom.size[tmp.exp@seqinfo@seqnames,]$seqlengths
    
    cov.list[[paste0(samplename,".",chr,".R2.all.neg")]] = coverage(tmp.exp)
    
    tmp.exp <- tmp[tmp$strand == "+",]
    tmp.exp$start <- tmp.exp$end_r2
    tmp.exp$end <- tmp.exp$end_r2
    tmp.exp <-  as(tmp.exp,"GRanges")
    tmp.exp@seqinfo@seqlengths <- chrom.size[tmp.exp@seqinfo@seqnames,]$seqlengths
    
    cov.list[[paste0(samplename,".",chr,".R2.all.pos")]] = coverage(tmp.exp)
    
    # celltype R1
    tmp <- split(tmp,tmp$celltype)
    
    for (i in names(tmp)) {
      
        tmp.exp <- tmp[[i]] # raw R1 count
        tmp.exp <- tmp.exp[tmp.exp$strand == "-",]
        tmp.exp$start <- tmp.exp$PAS
        tmp.exp$end <- tmp.exp$PAS
        tmp.exp <-  as(tmp.exp,"GRanges")
        tmp.exp@seqinfo@seqlengths <- chrom.size[tmp.exp@seqinfo@seqnames,]$seqlengths
        
        cov.list[[paste0(samplename,".",chr,".",i,".raw.neg")]] = coverage(tmp.exp)

        tmp.exp <- tmp[[i]]
        tmp.exp <- tmp.exp[tmp.exp$strand == "+",]
        tmp.exp$start <- tmp.exp$PAS
        tmp.exp$end <- tmp.exp$PAS
        tmp.exp <-  as(tmp.exp,"GRanges")
        tmp.exp@seqinfo@seqlengths <- chrom.size[tmp.exp@seqinfo@seqnames,]$seqlengths
        
        cov.list[[paste0(samplename,".",chr,".",i,".raw.pos")]] = coverage(tmp.exp)
    }
    return(cov.list)
  })
  
  cov.data.list <- cov.data.list %>% Reduce("c",.)
  saveRDS(cov.data.list,file = paste0(outpath,"/pas.df/cov.data.list.rds"))
  
  print("Writing raw profile")
  # cov.data.list[R.raw.strand // samplename.chr.celltype.raw.strand ]
  write.index = apply(expand.grid(unique(cell.metadata$celltype), 
                                  c("raw.pos","raw.neg")), 
                      1, paste, collapse=".")
  
  cov.data.merge <- future_lapply(c("R1.all.neg","R1.all.pos","R2.all.pos","R2.all.neg",write.index),function(i){
    sample.index =  apply(expand.grid(samplenames,chrs,i), 1, paste, collapse=".")
    tmp <- cov.data.list[names(cov.data.list) %in% sample.index]
    tmp <- Reduce("+",tmp)
    
    cov.data.list[names(cov.data.list) %in% sample.index] <- NULL
    cov.data.list[[i]] <- tmp
    return(tmp)
    
    tmp <- tmp*1000/nrow(cell.metadata)
    rtracklayer::export.bw(object = tmp, con = paste0(outpath,"celltype.bw/",i,".bigwig"))
  })
  
  names(cov.data.merge) = c("R1.all.neg","R1.all.pos","R2.all.pos","R2.all.neg",write.index)
  saveRDS(cov.data.merge,file = paste0(outpath,"/pas.df/cov.data.merge.rds"))
  
  # cov.data.merge[R.raw.strand // celltype.raw.strand ]
  
  write.index = table(cell.metadata$celltype)
  
  future_lapply(names(cov.data.merge),function(i){
    tmp <- cov.data.merge[[i]]
    rtracklayer::export.bw(object = tmp, con = paste0(outpath,"celltype.bw/",i,".bigwig"))
  }) 
  
  print("Writing celltype profile")
  
  if (track.nor == "all"){
    future_lapply(names(write.index),function(i){
      tmp <- cov.data.merge[["R1.all.pos"]]
      fc2multi = coverage(as(sc.plus[sc.plus$strand == "+",],"GRanges"),
                          weight = sc.plus[sc.plus$strand == "+",i])
      tmp <- tmp*fc2multi*1000/write.index[i]
      rtracklayer::export.bw(object = tmp, con = paste0(outpath,"celltype.bw/",i,".predict.pos.bigwig"))
      
      tmp <- cov.data.merge[["R1.all.neg"]]
      fc2multi = coverage(as(sc.plus[sc.plus$strand == "-",],"GRanges"),
                          weight = sc.plus[sc.plus$strand == "-",i])
      tmp <- tmp*fc2multi*1000/write.index[i]
      rtracklayer::export.bw(object = tmp, con = paste0(outpath,"celltype.bw/",i,".predict.neg.bigwig"))
    }) 
  } else if (track.nor == "celltype"){
    future_lapply(names(write.index),function(i){
      tmp <- cov.data.merge[[paste0(i,".raw.pos")]]
      tmp = tmp + coverage(as(sc.plus[sc.plus$strand == "+",],"GRanges"),weight = pseudo.count)
      fc2multi = coverage(as(sc.plus[sc.plus$strand == "+",],"GRanges"),
                          weight = sc.plus[sc.plus$strand == "+",i])
      tmp = tmp*fc2multi*1000/write.index[i]
      rtracklayer::export.bw(object = tmp, con = paste0(outpath,"celltype.bw/",i,".predict.pos.bigwig"))
      
      tmp <- cov.data.merge[[paste0(i,".raw.neg")]]
      tmp = tmp + coverage(as(sc.plus[sc.plus$strand == "-",],"GRanges"),weight = pseudo.count)
      fc2multi = coverage(as(sc.plus[sc.plus$strand == "-",],"GRanges"),
                          weight = sc.plus[sc.plus$strand == "-",i])
      tmp = tmp*fc2multi*1000/write.index[i]
      rtracklayer::export.bw(object = tmp, con = paste0(outpath,"celltype.bw/",i,".predict.neg.bigwig"))
    })
  }
  
  plan(sequential)
}


write.igv.xml <- function(igv.pre = "ipst",outpath,cell.metadata){

  xml.export <- paste0('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n',
                       '<Session genome="hg38" locus="chr2:231452062-231458313" nextAutoscaleGroup="5" version="8">\n',
                       '    <Resources>\n')

  cell.metadata$celltype = gsub("-",".",cell.metadata$celltype)
  write.index = apply(expand.grid(unique(cell.metadata$celltype), c("predict.pos","predict.neg")), 1, paste, collapse=".")

  for (i in c("R1.all.neg","R1.all.pos","R2.all.pos","R2.all.neg",write.index)) {
    xml.export <- paste0(xml.export,
                         '        <Resource path="view\\',igv.pre,'\\','celltype.bw','\\',i,'.bigwig" type="bigwig"/>\n')
  }
  xml.export <- paste0(xml.export,
                       '    </Resources>\n')
  xml.export <- paste0(xml.export,
                       '    <Panel height="1008" name="DataPanel" width="2541">\n')

  for (i in c("R1.all","R2.all")){

    xml.export <- paste0(xml.export,
                         '    <Track attributeKey="Overlay" autoScale="true" clazz="org.broad.igv.track.MergedTracks" fontSize="15" height="60" id="',i,'.bigwig" ',
                         'name="',i,'" renderer="BAR_CHART" visible="true">\n')

    xml.export <- paste0(xml.export,
                         '            <Track altColor="0,0,178" attributeKey="',i,'.pos.bigwig" ',
                         'autoScale="true" color="255,51,51" fontSize="10" height="40" id="',
                         'E:\\projects\\apa\\view\\',igv.pre,'\\celltype.bw\\',i,'.pos.bigwig" ','name="',i,'.pos.bigwig" ',' renderer="BAR_CHART" visible="true" windowFunction="mean"/>\n')
    xml.export <- paste0(xml.export,
                         '            <Track altColor="0,0,178" attributeKey="',i,'.neg.bigwig" ',
                         'autoScale="true" color="0,0,178" fontSize="10" height="40" id="',
                         'E:\\projects\\apa\\view\\',igv.pre,'\\celltype.bw\\',i,'.neg.bigwig" ','name="',i,'.neg.bigwig" ',' renderer="BAR_CHART" visible="true" windowFunction="mean"/>\n')

    xml.export <- paste0(xml.export,
                         '            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="2.4320135" minimum="0.0" type="LINEAR"/>\n')
    xml.export <- paste0(xml.export,
                         '        </Track>\n')
  }

  for (i in sort(unique(cell.metadata$celltype))){
    
    xml.export <- paste0(xml.export,
                         '    <Track attributeKey="Overlay" autoScale="true" clazz="org.broad.igv.track.MergedTracks" fontSize="15" height="60" id="',i,'.bigwig" ',
                         'name="',i,'" renderer="BAR_CHART" visible="true">\n')
    
    xml.export <- paste0(xml.export,
                         '            <Track altColor="0,0,178" attributeKey="',i,'.pos.bigwig" ',
                         'autoScale="true" color="255,51,51" fontSize="10" height="40" id="',
                         'E:\\projects\\apa\\view\\',igv.pre,'\\celltype.bw\\',i,'.predict.pos.bigwig" ','name="',i,'.predict.pos.bigwig" ',' renderer="BAR_CHART" visible="true" windowFunction="mean"/>\n')
    xml.export <- paste0(xml.export,
                         '            <Track altColor="0,0,178" attributeKey="',i,'.neg.bigwig" ',
                         'autoScale="true" color="0,0,178" fontSize="10" height="40" id="',
                         'E:\\projects\\apa\\view\\',igv.pre,'\\celltype.bw\\',i,'.predict.neg.bigwig" ','name="',i,'.predict.neg.bigwig" ',' renderer="BAR_CHART" visible="true" windowFunction="mean"/>\n')
    
    xml.export <- paste0(xml.export,
                         '            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="2.4320135" minimum="0.0" type="LINEAR"/>\n')
    xml.export <- paste0(xml.export,
                         '        </Track>\n')
  }
  
  xml.export <- paste0(xml.export,
                       '    </Panel>\n',
                       '    <Panel height="4350" name="FeaturePanel" width="2541">\n',
                       '        <Track attributeKey="Reference sequence" clazz="org.broad.igv.track.SequenceTrack" fontSize="10" id="Reference sequence" name="Reference sequence" sequenceTranslationStrandValue="POSITIVE" shouldShowTranslation="false" visible="true"/>\n',
                       '        <Track attributeKey="MANE Transcripts" clazz="org.broad.igv.track.FeatureTrack" colorScale="ContinuousColorScale;0.0;1.2590000629425049;255,255,255;0,0,178" fontSize="10" groupByStrand="false" id="https://hgdownload.soe.ucsc.edu/gbdb/hg38/mane/mane.bb" name="MANE Transcripts" visible="true"/>\n',
                       '        <Track attributeKey="Refseq Genes" clazz="org.broad.igv.track.FeatureTrack" colorScale="ContinuousColorScale;0.0;836.0;255,255,255;0,0,178" displayMode="EXPANDED" fontSize="10" groupByStrand="false" id="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz" name="Refseq Genes" visible="true"/>\n',
                       '    </Panel>\n',
                       '    <PanelLayout dividerFractions="0.8684436801375752"/>\n',
                       '    <HiddenAttributes>\n',
                       '        <Attribute name="DATA FILE"/>\n',
                       '        <Attribute name="DATA TYPE"/>\n',
                       '        <Attribute name="NAME"/>\n',
                       '    </HiddenAttributes>\n',
                       '</Session>\n'
  )

  writeLines(xml.export,paste0(outpath,igv.pre,".igv_session.cluster.xml"))
}

