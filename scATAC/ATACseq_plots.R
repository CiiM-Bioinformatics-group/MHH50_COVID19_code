### covid MHH50 plots

library(ArchR)

library(Seurat)
library(dplyr)
library(Matrix)
library(Matrix.utils)
library(ggplot2)
library(reshape2)
library(pryr)
library(RColorBrewer)
library(cowplot)
source("/vol/projects/CIIM/Covid_50MHH/scATAC/scripts/getMarkers.raw.R")

##### set up ####
addArchRGenome("hg38")
addArchRThreads(threads=4)

####
setwd("/vol/projects/CIIM/Covid_50MHH/scATAC/analysis_new")

projMono <- loadArchRProject("./projMono7")
projDm <- loadArchRProject("./projDm5/")

pbmc.pools <- readRDS("../../scRNA/analysis/pbmc.harmony.annotev2.rds")
mono.pools <- readRDS("../../scRNA/analysis/mono.pools.regMT.nF2k.dim10.reg0.3.rds")
mono.mks <- read.table("mono.markers.txt",sep="\t",header = T,stringsAsFactors = F)

pdf("Feature_monos.pdf",width=5.5,height=6.3)
FeaturePlot(mono.pools,features=c("CD14","FCGR3A","CD1C","HLA-DRA","CD163","IFI6"),ncol=2)
dev.off()


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

top10 <- mono.mks %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10$cluster <- factor(top10$cluster, levels=c(0,1,3,7,6,4,2,5))
top10 <- top10[order(top10$cluster),]

pdf("mono_marker_plots.newid.pdf",width=6.5,height=12)
DotPlot(mono.pools,features = rev(unique(top10$gene)),col.min= -1.5, col.max= 1.5,group.by="new.id")+
scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n =100, name = "RdBu")))+coord_flip()+xlab("Monocytes clusters")+ylab("")
dev.off()

Idents(mono.pools) <- "Severity2"
mono.pools <- RenameIdents(mono.pools,
`Severe` = "Severe",
`Mild` = "Mild",
`Post` = "PostCovid")
mono.pools$Severity <- Idents(mono.pools)

mcols <- brewer.pal(8,"Set1")
names(mcols)<-c("R1","R2","R3","R4","R5","R6","R7","R8")
DimPlot(mono.pools,group.by="new.id",split.by="Severity",label=F,cols= mcols)




#####
tmp<-data.frame(cell_type= mono.pools@active.ident,sample= mono.pools@meta.data$SampleID,donor= mono.pools@meta.data$Severity)
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
y2$Condition <- factor(y2$Condition,levels=c("Severe","Mild","PostCovid"))
#####
library(DirichletReg)
library(ggsignif)
y2$other <- 1-y2$f
dat.test <- NULL
for(i in levels(as.factor(y2$cell_type))){
  print(i)
  y2.i <- y2[y2$cell_type==i,]
  y2.i$Smp <- DR_data(y2.i[,5:6])
  res <- DirichReg(Smp~Condition,y2.i,model="alternative",base=2)
  x <- summary(res)
  sc <- x$coef.mat["ConditionMild",3:4]
  hc <- x$coef.mat["ConditionPostCovid",3:4]
  dat.test <- rbind(dat.test, c(i,sc,"Severe","Mild"))
  dat.test <- rbind(dat.test, c(i,hc,"Severe","PostCovid"))
  
  y2.i <- y2.i[y2.i$Condition != "Severe",]
  y2.i$Condition <- factor(y2.i$Condition, levels=c("Mild","PostCovid"))
   y2.i$Smp <- DR_data(y2.i[,5:6])
  res <- DirichReg(Smp~Condition,y2.i,model="alternative",base=2)
  x <- summary(res)
  mc <- x$coef.mat["ConditionPostCovid",3:4]
  dat.test <- rbind(dat.test, c(i,mc,"Mild","PostCovid"))
}

dat.test<- as.data.frame(dat.test)
names(dat.test) <- c("cell_type","Z","Pr","Comp1","Comp2")
dat.test$Pr <- as.numeric(as.character(dat.test$Pr))
dat.test$Z <- as.numeric(as.character(dat.test$Z))

dat.test$sig <- "ns"
dat.test$sig[dat.test$Pr < 0.05]="*"
dat.test$sig[dat.test$Pr < 0.01]="**"
dat.test$sig[dat.test$Pr < 0.001]="***"

dat.sig <- dat.test[dat.test$sig != "ns",]
maxs <- NULL
for(i in 1:nrow(dat.sig)){
	c <- as.character(dat.sig[i,"cell_type"])
	if(dat.sig[i, "Comp1"]=="Severe" & dat.sig[i,"Comp2"]=="Mild"){
	        maxs <- c(maxs,max(y2$f[y2$cell_type == c])+0.03)
	}
	if(dat.sig[i, "Comp1"]=="Mild" & dat.sig[i,"Comp2"]=="PostCovid"){
	        maxs <- c(maxs,max(y2$f[y2$cell_type == c])+0.09)
	}
	if(dat.sig[i, "Comp1"]=="Severe" & dat.sig[i,"Comp2"]=="PostCovid"){
	        maxs <- c(maxs,max(y2$f[y2$cell_type == c])+0.06)
	}
	}
dat.sig$maxs <- maxs

#pdf("cpp_mono3.pdf",width=15,height=4.5)
p <- ggplot(data=y2, aes(x=Condition,y=f))+geom_boxplot(aes(fill=Condition),outlier.shape=NA)
p2 <- p + facet_grid(.~cell_type)+geom_jitter(width = 0.2)+ theme_classic()+scale_fill_manual(values=c("#9e0142","#fee08b","#4393c3"))+xlab("")+ylab("percent")
plot(p2)
p2 + geom_signif(data=dat.sig, aes(xmin=Comp1,xmax=Comp2, annotations=sig, y_position=maxs), manual = TRUE, vjust=0.8, tip_length=0.005)+theme(axis.text.x=element_text(angle = 70,hjust = 1))
#dev.off()



###### ATAC figures ####

projMono$predictedGroup_Monoflt3 <- projMono$predictedGroup_Mono3
projMono$predictedGroup_Monoflt3[projMono$predictedScore_Mono3 < 0.6] = "undefine"

mcols <- brewer.pal(9,"Set1")
names(mcols)<-c("R1","R2","R3","R4","R5","R6","R7","R8","undefine")

p1 <- plotEmbedding(ArchRProj = projMono, colorBy = "cellColData", name = "Clusters", embedding="UMAP",labelAsFactors=F,size = 0.5)
p2 <- plotEmbedding(ArchRProj = projMono, colorBy = "cellColData", name = "predictedGroup_Monoflt3", embedding="UMAP", labelAsFactors=F,pal= mcols,size = 0.5)

plotPDF(p1,p2, name = "Plots-UMAP_linkMono4.pdf", ArchRProj = projMono,addDoc=F, width=5, height=5)


idx <- BiocGenerics::which(projMono$Severity %in% "Mild")
mild <- projMono$cellNames[idx]

idx <- BiocGenerics::which(projMono$Severity %in% "Severe")
severe <- projMono$cellNames[idx]

idx <- BiocGenerics::which(projMono$Severity %in% "Post")
post <- projMono$cellNames[idx]

p1 <- plotEmbedding(ArchRProj = projMono, colorBy = "cellColData", highlightCells =mild, name = "Clusters", embedding="UMAP",labelAsFactors=F,size = 0.8,labelSize =0)
p2 <- plotEmbedding(ArchRProj = projMono, colorBy = "cellColData", highlightCells = severe, name = "Clusters", embedding="UMAP",labelAsFactors=F,size = 0.8,labelSize =0)
p3 <- plotEmbedding(ArchRProj = projMono, colorBy = "cellColData", highlightCells = post, name = "Clusters", embedding="UMAP",labelAsFactors=F,size = 0.8,labelSize =0)
plotPDF(p1,p2,p3, name = "Plots-UMAP_condition.pdf", ArchRProj = projMono,addDoc=F, width=5, height=5)


#whitePurple

#projMono <- addImputeWeights(projMono)
markerGenes  <- c(
    "CD14","FCGR3A","CD1C","HLA-DRA","IFI6","CD163",
    "TNF","IL1B","IL6","CXCL9","CXCL10","CXCL11"
  )
p <- plotEmbedding(
    ArchRProj = projMono, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    continuousSet= "whitePurple",
    embedding = "UMAP",
    imputeWeights = NULL #getImputeWeights(projMono)
)

p2 <- lapply(p, function(x){
x + guides(color = FALSE, fill = FALSE) + theme_ArchR(baseSize = 6.5) +
theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + theme(
axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()
) })
p3 <- do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

plotPDF(p3,
name = "Plot-UMAP-Mono1-W-Imputation.pdf", ArchRProj = projMono,
addDOC = FALSE, width = 10, height = 12)



markersGS <- getMarkerFeatures(
	ArchRProj=projMono,
	useMatrix="GeneScoreMatrix",
	groupBy="Clusters",
	bias=c("TSSEnrichment","log10(nFrags)"),
	testMethod="wilcoxon")
	
markerList <- getMarkers(markersGS,cutOff="FDR <= 0.05 ",n=20)
c1 <- markerList$C1$name
c2 <- markerList$C2$name
c3 <- markerList$C3$name
c4 <- markerList$C4$name
c5 <- markerList$C5$name
c6 <- markerList$C6$name
cc <- unique(c(c1,c2,c3,c4,c5,c6))


heatmapGS <- markerHeatmap(
	seMarker = markersGS,
	cutOff = "FDR<0.1",
	nPrint=20,
	nLabel=20,
	transpose=T,
	labelRows=T,
	clusterCols =F,
	returnMatrix=T
)
hetmapGS.top20 <- heatmapGS[,cc]
intersectGS <- t(heatmapGS[,intersect(colnames(heatmapGS),mono.mks$gene)])
g1 <- names(sort(intersectGS[,1],decreasing = T)[1:20])
g2 <- names(sort(intersectGS[,2],decreasing = T)[1:20])
g3 <- names(sort(intersectGS[,3],decreasing = T)[1:20])
g4 <- names(sort(intersectGS[,4],decreasing = T)[1:20])
g5 <- names(sort(intersectGS[,5],decreasing = T)[1:20])
g6 <- names(sort(intersectGS[,6],decreasing = T)[1:20])
gg <- unique(c(g1,g2,g3,g4,g5,g6))
hetmapGS.gg <- heatmapGS[,gg]
pheatmap(t(hetmapGS.gg),cluster_rows = F,cluster_cols = F)


DotPlot(mono.pools,features = rev(gg),col.min= -1.5, col.max= 1.5,group.by="new.id")+
scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n =100, name = "RdBu")))+coord_flip()+xlab("Monocytes clusters")+ylab("")

flt <- unique(c(g1,g2,g4,g6))
hetmapGS.flt <- heatmapGS[c(1,2,4,6),flt]
mono.pools.flt <- subset(mono.pools, new.id %in% c("R7","R8","R2","R3","R5"))
mono.pools.flt$sort.id <- factor(mono.pools.flt$new.id,levels=c("R7","R8","R2","R3","R5"))
DotPlot(mono.pools.flt,features = rev(flt),col.min= -1.5, col.max= 1.5,group.by="sort.id")+
scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n =100, name = "RdBu")))+coord_flip()+xlab("")+ylab("")+theme(axis.text.y=element_text(size = 8))
pheatmap(t(hetmapGS.gg),cluster_rows = F,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name =
  "RdBu")))(100),breaks=-2:2)




####

pdf("testplot.pdf", height=20, width=10)
DotPlot(pbmc.sub,features = unique(gene.spi),col.min= -1.5, col.max= 1.5,group.by="cell.cond")+
scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n =100, name = "RdBu")))+coord_flip()+xlab("")+ylab("")+theme(axis.text.x = element_text(angle = 45,hjust=1), axis.text.y=element_text(size = 2))
dev.off()

gene.cebpd <- read.table("./gene.CEBPD.txt",stringsAsFactors=F)[,1]
pdf("testplot2.pdf", height=20, width=10)
DotPlot(pbmc.sub,features = unique(gene.cebpd),col.min= -1.5, col.max= 1.5,group.by="cell.cond")+
scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n =100, name = "RdBu")))+coord_flip()+xlab("")+ylab("")+theme(axis.text.x = element_text(angle = 45,hjust=1), axis.text.y=element_text(size = 2))
dev.off()

##############

	
	
markersPeaks <- getMarkerFeatures(
    ArchRProj = projDm2, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)	
	
	
	
markerList <- getMarkers(markersPeaks,cutOff="FDR <= 0.1 ")

a <- data.frame(markerList$cMono ,cell="cMono")
b <- data.frame(markerList$ncMono,cell="ncMono")


c <- data.frame(markerList$`CD4 T`,cell="CD4 T")
d <- data.frame(markerList$`CD8 T`,cell="CD8 T")
e <- data.frame(markerList$NK,cell="NK")

f <- data.frame(markerList$B,cell="B")
g <- data.frame(markerList$Plasmablast,cell="Plasmablast")


cc <- rbind(a,b,c,d,e,f,g)
write.table(cc,"Cluster.markers.projDm.txt",sep="\t",quote=F)

##############
markersPeaks <- getMarkerFeatures(
    ArchRProj = projMono, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)	
	
	
	
markerList <- getMarkers(markersPeaks,cutOff="FDR <= 0.1 ")

a <- data.frame(markerList$C1,cell="C1")
b <- data.frame(markerList$C2,cell="C2")
c <- data.frame(markerList$C3,cell="C3")
d <- data.frame(markerList$C4,cell="C4")

cc <- rbind(a,b,c,d)
write.table(cc,"Cluster.markers.projMono.txt",sep="\t",quote=F)

