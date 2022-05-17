options(future.globals.maxSize = 10737418240)  # 10240*1024^2 = 10737418240  for 10GB

library(Seurat)
library(dplyr)
library(Matrix)
library(Matrix.utils)
library(ggplot2)
library(reshape2)
library(pryr)
library(RColorBrewer)
library(tidyr)
library(cowplot)


YLINKED_GENES=c("USP9Y","DDX3Y","EIF1AY","RPS4Y1")
ReadBowen <- function(bdir,bid,btype){
  bdata <- Read10X(data.dir = bdir)
  bsample <- CreateSeuratObject(counts = bdata, min.cells = 3, min.features = 100, project = "10X")
  #bsample$individual <- individual
  bsample$sampid <- bid
  bsample$samptype <- btype
  return(bsample)
}

add.souporcell <- function(sfile, sobject){
  pool1.info <- read.table(file=sfile, header=T, sep="\t", stringsAsFactors=F)
  #pool1.info$barcode <- gsub("-1", "", pool1.info$barcode)
  row.names(pool1.info) <- pool1.info$barcode
  pool1.info <- pool1.info[,c(2,3)]
  clono_seurat <- AddMetaData(object=sobject, metadata=pool1.info)
  return(clono_seurat)
}


run.pool <- function(i, p){
  sfile1 <- paste("/vol/projects/CIIM/Covid_50MHH/scRNA/souporcell/out_",i,"/clusters.tsv",sep="")
  bdir1 = paste("/vol/projects/CIIM/Covid_50MHH/scRNA/align_out/align_",i,"/outs/filtered_feature_bc_matrix",sep="")
  fig1 = paste("raw_stat_",p,".pdf",sep="")
  fig2 = paste("assign_stat_",p,".pdf",sep="")
  fig3 = paste("Y_link_",p,".pdf",sep="")
  pbmc.pool <- ReadBowen(bdir = bdir1, bid=p, btype=p)
  pbmc.pool <- add.souporcell(sfile1,pbmc.pool)
  pbmc.pool[["percent.mt"]] <- PercentageFeatureSet(object = pbmc.pool, pattern = "^MT-")
  
  Idents(pbmc.pool)<-"status"
  try(dev.off())
  pdf(fig1,height=5,width=12)
  a<-VlnPlot(pbmc.pool,features = c("nFeature_RNA","nCount_RNA","percent.mt"),pt.size=0,ncol=3)
  plot(a)
  dev.off()
  pbmc.pool <- subset(x = pbmc.pool, subset = nFeature_RNA > 100 & nFeature_RNA < 3000 & nCount_RNA > 500 & percent.mt < 15 & status == "singlet")
  Idents(pbmc.pool)<-"assignment"
  pdf(fig2,height=5,width=12)
  a<-VlnPlot(pbmc.pool,features = c("nFeature_RNA","nCount_RNA","percent.mt"),pt.size=0,ncol=3)
  plot(a)
  dev.off()
  
  pbmc.pool$pool = p
  pbmc.pool$sample_assign <- paste(pbmc.pool$pool,pbmc.pool$assignment, sep="_")
  
  pbmc.pool <- NormalizeData(pbmc.pool)
  
  
YLINKED_GENES=c("FAM197Y5","TSPY6P","TSPY10","CDY1B","BPY2B","FAM197Y1","SLC9B1P1","DAZ3","USP9Y","DAZ4","DDX3Y","UTY","BPY2C","CDY1","TMSB4Y","VCY","VCY1B",
            "NLGN4Y","CDY2B","CDY2A","HSFY1","HSFY2","KDM5D","EIF1AY","CYorf17","RBMY1B","RBMY1A1","RBMY1E","RBMY1D","SRY","RPS4Y1","PRY2","RBMY1F","ZFY","TGIF2LY",
                "PCDH11Y","AC012067.1","TSPY2","RBMY1J","PRY","BPY2","AMELY","DAZ1","TBL1Y","TSPY4","TSPY8","TSPY3","TSPY1","DAZ2")

   pdf(fig3,height=6,width=12)
  a<- DotPlot(pbmc.pool,features = YLINKED_GENES,group.by="sample_assign")+RotatedAxis()
	plot(a)
  dev.off()
  
  saveRDS(pbmc.pool, paste("pbmc.pool.",p,".rds",sep=""))
  return(pbmc.pool)
}

#pbmc.pool2a <- run.pool("2a","P2a")

#blood samples
#pbmc.pool1a <- run.pool("1a","P1a")
#pbmc.pool2a <- run.pool("2a","P2a")
#pbmc.pool3a <- run.pool("3a","P3a")
#pbmc.pool4a <- run.pool("4a","P4a")

pbmc.pool1b <- run.pool("1b","P1b")
pbmc.pool2b <- run.pool("2b","P2b")
pbmc.pool3b <- run.pool("3b","P3b")
pbmc.pool4b <- run.pool("4b","P4b")

pbmc.pool1c <- run.pool("1c","P1c")
pbmc.pool2c <- run.pool("2c","P2c")
pbmc.pool3c <- run.pool("3c","P3c")
pbmc.pool4c <- run.pool("4c","P4c")

pbmc.pool1d <- run.pool("1d","P1d")
pbmc.pool2d <- run.pool("2d","P2d")
pbmc.pool3d <- run.pool("3d","P3d")
pbmc.pool4d <- run.pool("4d","P4d")

 # i="4d"
  

pbmc.poolb <- merge(x= pbmc.pool1b, y=c(pbmc.pool2b,pbmc.pool3b,pbmc.pool4b),project = "10X", 
                    add.cell.ids = c("p1b","p2b","p3b","p4b"))
pbmc.poolb$batch <- "b"

pbmc.poolc <- merge(x= pbmc.pool1c, y=c(pbmc.pool2c,pbmc.pool3c,pbmc.pool4c),project = "10X", 
                    add.cell.ids = c("p1c","p2c","p3c","p4c"))
pbmc.poolc$batch <- "c"
pbmc.poold <- merge(x= pbmc.pool1d, y=c(pbmc.pool2d,pbmc.pool3d,pbmc.pool4d),project = "10X", 
                    add.cell.ids = c("p1d","p2d","p3d","p4d"))
pbmc.poold$batch <- "d"


pbmc.pools <- merge(x= pbmc.poolb, y=c(pbmc.poold, pbmc.poolc), project = "10X")
saveRDS(pbmc.pools, "pbmc.pools.raw.rds")
pbmc.pools.test <- pbmc.pools

pbmc.pools.test <- NormalizeData(object = pbmc.pools.test, normalization.method = "LogNormalize")      
pbmc.pools.test <- FindVariableFeatures(object = pbmc.pools.test, selection.method = "vst", nfeatures = 2000)
pbmc.pools.test <- ScaleData(object = pbmc.pools.test, features = VariableFeatures(pbmc.pools.test))
pbmc.pools.test <- RunPCA(pbmc.pools.test, verbose = FALSE)
pbmc.pools.test <- RunUMAP(pbmc.pools.test, dims = 1:20)
pbmc.pools.test <- FindNeighbors(object = pbmc.pools.test, reduction = "pca", dims = 1:20) #used to be 18
pbmc.pools.test <- FindClusters(pbmc.pools.test, resolution = 0.4)    
saveRDS(pbmc.pools.test, "pbmc.pools.test.rds")
pdf("UMAP_merge.pdf",height=6,width=6) 
p1 <- DimPlot(pbmc.pools.test,label=T)
                   plot(p1)
                   dev.off()
                   
pdf("UMAP_merge.split.pdf",height=12,width=12) 
p2 <- DimPlot(pbmc.pools.test,label=T,split.by="sample_assign",ncol=4)   
                   plot(p2)
                   dev.off()
                   
                   
pdf("UMAP_merge.pool.pdf",height=6,width=12) 
p3 <- DimPlot(pbmc.pools.test,label=T,split.by="samptype",ncol=4)   
plot(p3)
dev.off()
                   
                   
mks <- c("CCR7","TCF7","IL7R","GZMK",
         "CD8A","NKG7","GZMB","CD14","LYZ","FCGR3A","CST3",
         "CD86","CD163","TNF","IL1B","CD79A","CD27","SDC1","PPBP")


get.violin.data1 <- function(seurat, genes) {
  #  output = data.frame(gene = character(0), value= numeric(0), ident = character(0))
  #####################################################################################################
  output = data.frame(gene = character(0), value= numeric(0), ident = character(0), tech=character(0))
  #####################################################################################################
  for (gene in genes) {
    if(any(gene == seurat@assays$RNA@data@Dimnames[[1]])){
      data.use = data.frame(FetchData(seurat,gene))
      data.use = t(data.use)
      data.melt=data.frame(rep(gene, length(seurat@active.ident)))
      colnames(data.melt)[1]="gene"
      data.melt$value=as.numeric(data.use[1,1:length(seurat@active.ident)])
      data.melt$id=names(data.use)[1:length(seurat@active.ident)]
      data.melt$ident=seurat@active.ident
      ############################################
      #data.melt$tech=seurat$stim # ???What was used to stimulate the cells?
      ############################################
     # if(any(data.melt$value != 0)) noise = rnorm(length(data.melt$value))/100000 else noise = 0
      data.melt$value=as.numeric(as.character(data.melt$value))#+noise
      output = rbind(output, data.melt)
    } else {
      data.melt=data.frame(
        gene = rep(gene, seurat@assays$RNA@data@Dim[2]),
        value = rep(0, seurat@assays$RNA@data@Dim[2]), 
        ident = seurat@active.ident
      )
      output = rbind(output, data.melt)
    }
  }
  return(output)
}

genes <- intersect(mks,row.names(pbmc.pools.test@assays$RNA@data))

violin.plot.data <- get.violin.data1(pbmc.pools.test, rev(genes))
pv <- ggplot(violin.plot.data, aes(factor(ident),value)) +
  geom_violin(scale="width",adjust=1,trim=TRUE,aes(fill=factor(ident)),show.legend = F) + 
  ylab("") + xlab("") +
  coord_flip() +facet_wrap(~ gene,scales = "free_x", ncol = length(levels(violin.plot.data$gene))) +
  theme(strip.text.x = element_text(size=12, angle=-50),
        axis.text.y = element_text(size=15),
        strip.background = element_blank(),
        panel.spacing.x = unit(c(-0.2), "lines"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0.5,0.4,0.5,1.3), "cm")) +
  scale_x_discrete(limits = rev(levels(violin.plot.data$ident)), position = "top")

pdf("VlnMarkers_merge.pdf",height=6,width=12)   
plot(pv)
dev.off()


pdf("UMAPMks_merge.pdf",height=12,width=12)          
mks <- c(  "CD3E", "CD14", "FCGR3A", "CD79A"  )
FeaturePlot(pbmc.pools.test, features=mks, order=T, ncol=2)      
dev.off()
                   
                   
                   
                   
library(harmony)
pbmc.pools <- NormalizeData(object = pbmc.pools, normalization.method = "LogNormalize")
pbmc.pools <- FindVariableFeatures(object = pbmc.pools,
                                   assay="RNA",
                                   selection.method = 'vst')
pbmc.pools <- ScaleData(object = pbmc.pools,
                        vars.to.regress = c("nCount_RNA"))
pbmc.pools <- RunPCA(object = pbmc.pools, features = VariableFeatures(object = pbmc.pools), verbose = FALSE)
pbmc.pools <- RunHarmony(pbmc.pools, group.by.vars = "batch", reduction = "pca", dims.use = 1:20)
pbmc.pools <- RunUMAP(pbmc.pools, reduction = "harmony", dims = 1:20, seed.use = 42)
pbmc.pools <- FindNeighbors(object = pbmc.pools, dims = 1:20, reduction="harmony", force.recalc = TRUE)
pbmc.pools <- FindClusters(object = pbmc.pools, resolution = 0.4, algorithm = 1)

saveRDS(pbmc.pools, "pbmc.pools.harmony.rds")

pdf("UMAP_harmony.pdf",height=6,width=6) 
p1 <- DimPlot(pbmc.pools,label=T)
                   plot(p1)
                   dev.off()
                   
pdf("UMAP_harmony.split.pdf",height=12,width=12) 
p2 <- DimPlot(pbmc.pools,label=T,split.by="sample_assign",ncol=4)   
                   plot(p2)
                   dev.off()
                   
                   
pdf("UMAP_harmony.pool.pdf",height=6,width=12) 
p3 <- DimPlot(pbmc.pools,label=T,split.by="samptype",ncol=4)   
plot(p3)
dev.off()


mks <- c("CCR7","TCF7","IL7R","GZMK",
         "CD8A","NKG7","GZMB","CD14","LYZ","FCGR3A","CST3",
         "CD86","CD163","TNF","IL1B","CD79A","CD27","SDC1","PPBP")
genes <- intersect(mks,row.names(pbmc.pools@assays$RNA@data))

violin.plot.data <- get.violin.data1(pbmc.pools, genes)
pv <- ggplot(violin.plot.data, aes(factor(ident),value)) +
  geom_violin(scale="width",adjust=1,trim=TRUE,aes(fill=factor(ident)),show.legend = F) + 
  ylab("") + xlab("") +
  coord_flip() +facet_wrap(~ gene,scales = "free_x", ncol = length(levels(violin.plot.data$gene))) +
  theme(strip.text.x = element_text(size=12, angle=-90),
        axis.text.y = element_text(size=15),
        strip.background = element_blank(),
        panel.spacing.x = unit(c(-0.2), "lines"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0.5,0.4,0.5,1.3), "cm")) +
  scale_x_discrete(limits = rev(levels(violin.plot.data$ident)), position = "top")

pdf("VlnMarkers_harmony.pdf",height=6,width=12)   
plot(pv)
dev.off()


pdf("UMAPMks_harmony.pdf",height=12,width=12)          
mks <- c(  "CD3E", "CD14", "FCGR3A", "CD79A"  )
FeaturePlot(pbmc.pools, features=mks, order=T, ncol=2)      
dev.off()



