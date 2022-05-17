library(Seurat)
library(ggplot2)

library(dplyr)
library(Matrix)
library(Matrix.utils)
library(reshape2)
library(pryr)
library(RColorBrewer)
library(cowplot)
options(future.globals.maxSize=20368709120)



#mono.pools[["percent.mt"]] <- PercentageFeatureSet(object = mono.pools, pattern = "^MT-")
VlnPlot(mono.pools,features = c("nFeature_RNA","nCount_RNA","percent.mt"),pt.size=0,group.by="patient")



mono.pools <- FindVariableFeatures(mono.pools, nfeatures = 1000, selection.method = "vst", assay="RNA")
mono.pools <- ScaleData(mono.pools, features=VariableFeatures(mono.pools))#, vars.to.regress = c("percent.mt"))

######
#top10 <- head(VariableFeatures(mono.pools), 10)

# plot variable features with and without labels
#plot1 <- VariableFeaturePlot(mono.pools)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1 + plot2


###### PCA
mono.pools <- RunPCA(mono.pools, features=VariableFeatures(mono.pools), verbose = FALSE)
#DimHeatmap(mono.pools, dims = 11:20, cells = 500, balanced = TRUE)

#mono.pools <- ScaleData(mono.pools, verbose = FALSE)

# t-SNE and Clustering
mono.pools <- RunUMAP(mono.pools, reduction = "pca", dims = 1:10)
mono.pools <- FindNeighbors(mono.pools, reduction = "pca", dims = 1:10)
mono.pools <- FindClusters(mono.pools,resolution = 0.3)

DimPlot(mono.pools,split.by="Severity",label=T)#p2 <- DimPlot(mono.pools,group.by="SampleID",split.by="Severity2")+NoLegend()



#####
tmp<-data.frame(cell_type= mono.pools@active.ident,sample= mono.pools@meta.data$SampleID,donor= mono.pools@meta.data$Severity2)
res<-matrix(NA,ncol=4,nrow=length(unique(tmp$cell_type))*length(unique(tmp$sample)))
res<-as.data.frame(res)
names(res)<-c("cell_type","freq","sample","Condition")
n=1;
for (d in unique(tmp$sample)){
  t1<-subset(tmp,sample==d)
  t2<-as.data.frame(table(t1$cell_type))
  res[n:(n+length(unique(tmp$cell_type))-1),1]<-as.character(t2$Var1)
  res[n:(n+length(unique(tmp$cell_type))-1),2]<-t2$Freq
  res[n:(n+length(unique(tmp$cell_type))-1),3]<-rep(as.character(d),length(unique(tmp$cell_type)))
  res[n:(n+length(unique(tmp$cell_type))-1),4]<-rep(as.character(unique(t1$donor)),length(unique(tmp$cell_type)))
  n=n+length(unique(tmp$cell_type));
}
x<-split(res,res$sample)
y<-lapply(x, function(w) {w=transform(w,f=freq/sum(freq))})
y2<-do.call(rbind,y)

#p <- ggplot(data=y2, aes(x=Condition,y=f))+geom_boxplot(aes(fill=Condition),outlier.shape=NA)
#p2 <- p + facet_grid(.~cell_type)+geom_jitter(width = 0.2)+ theme_bw()+
#theme(axis.text.x=element_text(angle = 70,hjust = 1),
#strip.background = element_blank()
#)

y2$Condition <- factor(y2$Condition,levels=c("Severe","Mild","Post"))
p <- ggplot(data=y2, aes(x=Condition,y=f))+geom_boxplot(aes(fill=Condition),outlier.shape=NA)
p2 <- p + facet_grid(.~cell_type)+geom_jitter(width = 0.2)+ theme_classic()+theme(axis.text.x=element_text(angle = 70,hjust = 1))+scale_fill_manual(values=c("#9e0142","#fee08b","#4393c3"))+xlab("")+ylab("percent")

pdf("scRNA_cpp.pdf",width=10,height=3)
plot(p2)
dev.off()

jpeg("testUMAP_1.jpeg",width=1200,height=500)
DimPlot(mono.pools,split.by="Severity2",label=T)
dev.off()

jpeg("testUMAP_noreg_nFeature1000_dim10_res0.3_Fig1.jpeg",width=1200,height=700)
plot_grid(p1,p2,ncol=1)
dev.off()

jpeg("testUMAP_noreg_nFeature1000_dim10_res0.3_Fig2.jpeg",width=1200,height=380)
FeaturePlot(mono.pools,features=c("CD14","FCGR3A","IFITM3","CD163"),ncol=4,order = T)
dev.off()




mono.mks <- FindAllMarkers(object = mono.pools, only.pos = TRUE, min.pct = 0.2, min.diff.pct = 0.05, logfc.threshold = 0.25)

top10.mono <- mono.mks %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DotPlot(mono.pools,features = unique(top10.mono$gene),col.min= -1.5, col.max= 1.5)+
scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n =100, name = "RdBu")))+coord_flip()+xlab("Monocytes clusters")+ylab("")


Idents(mono.pools) <- "seurat_clusters"
mono.pools <- RenameIdents(mono.pools,
`0` = "R1",
`1` = "R2",
`3` = "R3",
`7` = "R4",
`6` = "R5",
`4` = "R6",
`2` = "R7",
`5` = "R8"
)
mono.pools$new.id <- Idents(mono.pools)


Idents(mono.pools) <- "Severity2"
mono.pools <- RenameIdents(mono.pools,
`Severe` = "Severe",
`Mild` = "Mild",
`Post` = "PostCovid")
mono.pools$Severity <- Idents(mono.pools)
mcols <- brewer.pal(8,"Set1")
names(mcols)<-c("R1","R2","R3","R4","R5","R6","R7","R8")
DimPlot(mono.pools,group.by="new.id",split.by="Severity",label=F,cols= mcols)


saveRDS(mono.pools, "MHH50_monocytes.rds")



####################

mono.pools <- subset(mono.pools, seurat_clusters != "5")
mono.pools <- FindVariableFeatures(mono.pools, nfeatures = 1000, selection.method = "vst", assay="RNA")
mono.pools <- ScaleData(mono.pools, features=VariableFeatures(mono.pools))#, vars.to.regress = c("percent.mt"))

###### PCA
mono.pools <- RunPCA(mono.pools, features=VariableFeatures(mono.pools), verbose = FALSE)

# t-SNE and Clustering
mono.pools <- RunUMAP(mono.pools, reduction = "pca", dims = 1:10)
mono.pools <- FindNeighbors(mono.pools, reduction = "pca", dims = 1:10)
mono.pools <- FindClusters(mono.pools,resolution = 0.3)

DimPlot(mono.pools,split.by="Severity",label=T)



Idents(mono.pools) <- "seurat_clusters"
mono.pools <- RenameIdents(mono.pools,
`2` = "R1",
`0` = "R2",
`1` = "R3",
`3` = "R4",
`4` = "R5",
`5` = "R6",
`6` = "R7",
`7` = "R8"
)
mono.pools$new.id <- Idents(mono.pools)


Idents(mono.pools) <- "Severity"
mono.pools <- RenameIdents(mono.pools,
`severe` = "Active",
`mild` = "Active",
`post` = "Convalescent")
mono.pools$Disease <- Idents(mono.pools)

mcols <- brewer.pal(8,"Set1")
names(mcols)<-c("R1","R2","R3","R4","R5","R6","R7","R8")
DimPlot(mono.pools,group.by="new.id",split.by="Disease",label=F,cols= mcols)

Idents(mono.pools) <- "new.id"
mono.mks <- FindAllMarkers(object = mono.pools, only.pos = TRUE, min.pct = 0.2, min.diff.pct = 0.05, logfc.threshold = 0.25)

top10.mono <- mono.mks %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DotPlot(mono.pools,features = unique(top10.mono$gene),col.min= -1.5, col.max= 1.5)+
scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n =100, name = "RdBu")))+coord_flip()+xlab("Monocytes clusters")+ylab("")
