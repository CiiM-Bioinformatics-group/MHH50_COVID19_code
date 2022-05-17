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

#projCovid2 <- loadArchRProject("./Save_projCovid")
#projMonoFlt <- loadArchRProject("projMonoFlt2")
#projMono <- loadArchRProject("./projMono3")
##### set up ####
addArchRGenome("hg38")
addArchRThreads(threads=4)


#idx <- BiocGenerics::which(projMono$ClustersHm %in% c("C14","C12","C13"))
#mono.cells <- projMono$cellNames[idx]
#projMono <- subsetArchRProject(ArchRProj= projMono, cells= mono.cells, outputDirectory="projMono2")

#idx <- BiocGenerics::which(projMono$ClustersHm %in% c("C12","C13"))
#mono.cells <- projMono$cellNames[idx]
#projMono <- subsetArchRProject(ArchRProj= projMono, cells= mono.cells, outputDirectory="projMono3")


projMono <- loadArchRProject("./projMono4")


projMono <- addIterativeLSI(
	ArchRProj= projMono,
	useMatrix="TileMatrix",
	name="IterativeLSI",
	iterations=2,
	clusterParams=list(
		resolution=c(0.3),
		sampleCells=10000,
		n.start=10
	),
	varFeatures=25000,
	dimsToUse=1:30,
	force=TRUE
)

projMono <- addClusters(
	input= projMono,
	reducedDims="IterativeLSI",
	method="Seurat",
	name="Clusters",
	resolution=0.8,
	force=TRUE
)

projMono <- addUMAP(
	ArchRProj= projMono,
	reducedDims="IterativeLSI",
	name="UMAP",
	nNeighbors=30,
	minDist=0.5,
	metric="cosine",
	force=TRUE
)

p1 <- plotEmbedding(ArchRProj = projMono, colorBy = "cellColData", name = "Severity", embedding="UMAP")
p2 <- plotEmbedding(ArchRProj = projMono, colorBy = "cellColData", name = "Clusters", embedding="UMAP")

plotPDF(p1,p2, name = "Plots-UMAP_Mono-severe-Clusters.pdf", ArchRProj = projMono,addDoc=F, width=5, height=5)


library(DirichletReg)
library(ggsignif)
tmp<-data.frame(cell_type= projMono$Clusters,sample= projMono$patient,donor= projMono$Severity)
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

y2$Condition <- factor(y2$Condition,levels=c("Severe","Mild","Post"))

y2$other <- 1-y2$f


dat.test <- NULL
for(i in levels(as.factor(y2$cell_type))){
#for(i in c("PLA1","PLA2","PLA3")){  
  print(i)
  y2.i <- y2[y2$cell_type==i,]
  y2.i$Smp <- DR_data(y2.i[,5:6])
  res <- DirichReg(Smp~Condition,y2.i,model="alternative",base=2)
  x <- summary(res)
  sc <- x$coef.mat["ConditionMild",3:4]
  hc <- x$coef.mat["ConditionPost",3:4]
  dat.test <- rbind(dat.test, c(i,sc,"Severe","Mild"))
  dat.test <- rbind(dat.test, c(i,hc,"Severe","Post"))
  
  y2.i <- y2.i[y2.i$Condition != "Severe",]
  y2.i$Condition <- factor(y2.i$Condition, levels=c("Mild","Post"))
   y2.i$Smp <- DR_data(y2.i[,5:6])
  res <- DirichReg(Smp~Condition,y2.i,model="alternative",base=2)
  x <- summary(res)
  mc <- x$coef.mat["ConditionPost",3:4]
  dat.test <- rbind(dat.test, c(i,mc,"Mild","Post"))
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
	if(dat.sig[i, "Comp1"]=="Mild" & dat.sig[i,"Comp2"]=="Post"){
	        maxs <- c(maxs,max(y2$f[y2$cell_type == c])+0.09)
	}
	if(dat.sig[i, "Comp1"]=="Severe" & dat.sig[i,"Comp2"]=="Post"){
	        maxs <- c(maxs,max(y2$f[y2$cell_type == c])+0.06)
	}
	}
dat.sig$maxs <- maxs


y2$Cond <- as.character(y2$Condition)
y2$Cond[y2$Cond=="Post"]="PostCovid"
dat.sig$Comp2 <- as.character(dat.sig$Comp2)
dat.sig$Comp2[dat.sig$Comp2=="Post"]="PostCovid"
dat.sig$Comp2 <- as.factor(dat.sig$Comp2)
y2$Cond <- factor(y2$Cond,levels=c("Severe","Mild","PostCovid"))


pdf("cpp_mono3.pdf",width=10,height=3)
p <- ggplot(data=y2, aes(x= Cond,y=f))+geom_boxplot(aes(fill= Cond),outlier.shape=NA)
p2 <- p + facet_grid(.~cell_type)+geom_jitter(width = 0.2)+ theme_classic()+scale_fill_manual(values=c("#9e0142","#fee08b","#4393c3"))+xlab("")+ylab("percent")
plot(p2)
p2 + geom_signif(data=dat.sig, aes(xmin=Comp1,xmax=Comp2, annotations=sig, y_position=maxs), manual = TRUE, vjust=0.8, tip_length=0.005)+theme(axis.text.x=element_text(angle = 70,hjust = 1))
dev.off()



markersGS <- getMarkerFeatures(
	ArchRProj=projMono,
	useMatrix="GeneScoreMatrix",
	groupBy="Clusters",
	bias=c("TSSEnrichment","log10(nFrags)"),
	testMethod="wilcoxon"

)

markerList <- getMarkers(markersGS,cutOff="FDR <= 0.05 ",n=20)


c1 <- markerList$C1$name
c2 <- markerList$C2$name
c3 <- markerList$C3$name
c4 <- markerList$C4$name
c5 <- markerList$C5$name
c6 <- markerList$C6$name
cc <- unique(c(c1,c2,c3,c4,c5,c6))

deg <- NULL
for(i in names(markerList)){
	aaa <- as.data.frame(markerList[[i]])
	aaa$Cluster <- i
	deg <- rbind(deg,aaa)
}



marker.gene <- c("KCNQ1","LYN","BRAT1","IFITM3","HLA-DQA1","HLA-DRB1","HLA-DRA","LMNA","CLU","CSMD3","FAM187B")

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
#ComplexHeatmap::draw(heatmapGS, heatmap_legend_side="bot", annotation_legend_side="bot")

plotPDF(heatmapGS, name = "Plot-geneScores_test.pdf", ArchRProj = projMono,addDoc=F, width=8, height=6)

saveArchRProject(ArchRProj= projMono, outputDirectory="projMono3",load=F)


#whitePurple

projMono <- addImputeWeights(projMono)
markerGenes  <- c(
    "CD14","FCGR3A","CD68","CD163","IFI6","IFITM3",
    "TNF","IL1B","IL6","CXCL9","CXCL10","CXCL11"
  )
  
markerGenes  <- c(
    "CD14","FCGR3A","OAS3","CD163",
    "BATF3" ,"CDKN1C", "LINC02432", "MEG3", "NECAB3"
  )  

p <- plotEmbedding(
    ArchRProj = projMono, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projMono)
)

p2 <- lapply(p, function(x){
x + guides(color = FALSE, fill = FALSE) + theme_ArchR(baseSize = 6.5) +
theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + theme(
axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()
) })
p3 <- do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

plotPDF(p3,
name = "Plot-UMAP-Mono2-W-Imputation.pdf", ArchRProj = projMono,
addDOC = FALSE, width = 10, height = 12)






projMono$cluster.cond <- paste(projMono$Clusters, projMono$Severity,sep=".")
projMono <- addGroupCoverages(ArchRProj = projMono, groupBy = "cluster.cond")

pathToMacs2 <- "/vol/projects/CIIM/resources/tools/encode-atac-seq-pipeline/bin/macs2"
projMono <- addReproduciblePeakSet(
	ArchRProj = projMono,
	groupBy="cluster.cond",
	pathToMacs2= pathToMacs2
	)
getPeakSet(projMono)
projMono <- addPeakMatrix(projMono)
saveArchRProject(ArchRProj= projMono, outputDirectory="projMono4",load=F)

table(projMono$cluster.cond)

markerTest.cluster <- getMarkerFeatures(
  ArchRProj = projMono, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
)

markerList <- getMarkers(markerTest.cluster, cutOff = "FDR <= 0.1", returnGR=F)
markerList <- getMarkers(markerTest.cluster, cutOff = "FDR <= 0.1 & Log2FC >= 0.5", returnGR=F)
markerList$C4

markerTest.cluster.cond <- getMarkerFeatures(
  ArchRProj = projMono, 
  useMatrix = "PeakMatrix",
  groupBy = "cluster.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
)
markerList <- getMarkers(markerTest.cluster.cond, cutOff = "FDR <= 0.1 & Log2FC >= 0.5", returnGR=F)


aa <- getMarkers.raw(markerTest.cluster.cond, cutOff="Pval < 0.05 & Log2FC >=0.5", returnGR=F)


markerList$C4.Severe


peakset = getPeakSet(projMono)
promoter_peak = peakset[which(peakset$peakType == "Promoter"), ]
projMono = addFeatureMatrix(
  input = projMono,
  features = promoter_peak,
  matrixName = "PromoterMatrix"
)
promoter_cm <- getMarkerFeatures(
  ArchRProj = projMono,
  useMatrix = "PromoterMatrix",
  groupBy = "cell_dis2",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Classical Monocytes_post",
  bgdGroups = "Classical Monocytes_healthy"
)



heatmapPeaks <- markerHeatmap(
  seMarker = markerTest.cluster, 
  cutOff = "Pval <= 0.05 & Log2FC >= 1",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = projMono, addDOC = FALSE)





markerGenes <- c("ABO","CXCR6")
p <- plotBrowserTrack(
	ArchRProj= projMono,
	groupBy="cluster.cond",
	geneSymbol=markerGenes,
	upstream=500000,
	downstream=500000
)
plotPDF(plotList=p, name="plot_track_Mono_subs.pdf",ArchRProj= projMono, addDOC=F, width=10,height=5)

plotPDF(plotList=p, name="plot_track_Mono_subs2.pdf",ArchRProj= projMono, addDOC=F, width=10,height=8)

library(chromVARmotifs)
projMono <- addMotifAnnotations(ArchRProj = projMono, motifSet = "cisbp", name = "Motif", force=TRUE)
projMono <- addArchRAnnotations(ArchRProj = projMono, collection = "EncodeTFBS",force=TRUE)
saveArchRProject(ArchRProj= projMono, outputDirectory="projMonoPeak2",load=F)



enrichMotifs <- peakAnnoEnrichment(
    seMarker = markerTest.cluster.cond,
    ArchRProj = projMono,
    peakAnnotation = "Motif",
    cutOff = "Pval <= 0.05 & Log2FC >= 1"
  )



enrichMotifs
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 5, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap2", width = 8, height = 6, ArchRProj = projMono, addDOC = FALSE)


enrichMotifs <- peakAnnoEnrichment(
    seMarker = markerTest.cluster,
    ArchRProj = projMono,
    peakAnnotation = "Motif",
    cutOff = "Pval <= 0.05 & Log2FC >= 1"
  )

enrichMotifs
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 5, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap3", width = 8, height = 6, ArchRProj = projMono, addDOC = FALSE)





markerList <- getMarkers(markerTest.cmono, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR=F)
head(markerList$cMono.Mild)

markerList <- getMarkers.raw(markerTest.cmono, cutOff = "Pval <= 0.05 & abs(Log2FC) >= 1", returnGR=F)
head(markerList$cMono.Mild)
write.table(markerList$cMono.Mild, file="mk.cMono.MvS.txt",sep="\t",quote=F,)

pma <- plotMarkers(seMarker = markerTest.cmono, name = "cMono.Mild", cutOff = "Pval <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
plotPDF(pma, name = "Mild-vs-Seve-cMono-MA", width = 5, height = 5, ArchRProj = projMono, addDOC = FALSE)
















########## Try other peaks.  #######
pathToMacs2 <- "/vol/projects/CIIM/resources/tools/encode-atac-seq-pipeline/bin/macs2"
projMono <- addGroupCoverages(ArchRProj = projMono, groupBy = "Clusters")
projMono <- addReproduciblePeakSet(
	ArchRProj = projMono,
	groupBy="Clusters",
	pathToMacs2= pathToMacs2
	)
getPeakSet(projMono)
projMono <- addPeakMatrix(projMono)
saveArchRProject(ArchRProj= projMono, outputDirectory="projMono5",load=F)

markerTest.cluster2 <- getMarkerFeatures(
  ArchRProj = projMono, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
)
markerList2 <- getMarkers(markerTest.cluster2, cutOff = "FDR <= 0.1 & Log2FC >= 0.5", returnGR=F)
markerList2$C4

library(chromVARmotifs)
projMono <- addMotifAnnotations(ArchRProj = projMono, motifSet = "cisbp", name = "Motif",force = TRUE)
projMono <- addArchRAnnotations(ArchRProj = projMono, collection = "EncodeTFBS",force = TRUE)
saveArchRProject(ArchRProj= projMono, outputDirectory="projMono5",load=F)

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markerTest.cluster2,
    ArchRProj = projMono,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )


enrichMotifs
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, cutOff = 2, pMax=5, pal = paletteContinuous(set = "comet", n = 10), transpose = TRUE, clusterCols=F)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap5", width = 8, height = 6, ArchRProj = projMono, addDOC = FALSE)


enrichTFBS <- peakAnnoEnrichment(
    seMarker = markerTest.cluster2,
    ArchRProj = projMono,
    peakAnnotation = "EncodeTFBS",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

heatmapEncode <- plotEnrichHeatmap(enrichTFBS, n = 7, cutOff = 1, pal = paletteContinuous(set = "comet", n = 10), pMax=10, transpose = TRUE, clusterCols=F)
plotPDF(heatmapEncode, name = "EncodeTFBS-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = projMono, addDOC = FALSE)
#################

markerList2 <- getMarkers(markerTest.cluster2, cutOff = "FDR <= 0.1 & Log2FC >= 0.5", returnGR=F)
aa <- getMarkers.raw(markerTest.cluster2, cutOff = "Pval <= 0.05 & Log2FC >= 0.5", returnGR=F)

markerList2 <- aa
df1 <- NULL
df2 <- NULL
df3 <- NULL
df1 <- data.frame(markerList2$C1.Post, "cond"="Post")
df2 <- data.frame(markerList2$C1.Mild, "cond"="Mild")
df3 <- data.frame(markerList2$C1.Severe, "cond"="Severe")
df <- rbind(df1,df2,df3)
write.table(df,"peakcond.C1.txt",sep="\t", quote=F, row.names=F)
df$cell = "C1"
dfc1 <- df


df1 <- NULL
df2 <- NULL
df3 <- NULL

df1 <- data.frame(markerList2$C2.Post, "cond"="Post")
df2 <- data.frame(markerList2$C2.Mild, "cond"="Mild")
df3 <- data.frame(markerList2$C2.Severe, "cond"="Severe")
df <- rbind(df1,df2,df3)
write.table(df,"peakcond.C2.txt",sep="\t", quote=F, row.names=F)
df$cell = "C2"
dfc2 <- df

df1 <- NULL
df2 <- NULL
df3 <- NULL

df1 <- data.frame(markerList2$C3.Post, "cond"="Post")
df2 <- data.frame(markerList2$C3.Mild, "cond"="Mild")
df3 <- data.frame(markerList2$C3.Severe, "cond"="Severe")
df <- rbind(df1,df2,df3)
write.table(df,"peakcond.C3.txt",sep="\t", quote=F, row.names=F)
df$cell = "C3"
dfc3 <- df

df1 <- NULL
df2 <- NULL
df3 <- NULL

df1 <- data.frame(markerList2$C4.Post, "cond"="Post")
df2 <- data.frame(markerList2$C4.Mild, "cond"="Mild")
df3 <- data.frame(markerList2$C4.Severe, "cond"="Severe")
df <- rbind(df1,df2,df3)
write.table(df,"peakcond.C4.txt",sep="\t", quote=F, row.names=F)
df$cell = "C4"
dfc4 <- df


df1 <- NULL
df2 <- NULL
df3 <- NULL

df1 <- data.frame(markerList2$C5.Post, "cond"="Post")
df2 <- data.frame(markerList2$C5.Mild, "cond"="Mild")
df3 <- data.frame(markerList2$C5.Severe, "cond"="Severe")
df <- rbind(df1,df2,df3)
write.table(df,"peakcond.C5.txt",sep="\t", quote=F, row.names=F)
df$cell = "C5"
dfc5 <- df

df1 <- NULL
df2 <- NULL
df3 <- NULL
df1 <- data.frame(markerList2$C6.Post, "cond"="Post")
df2 <- data.frame(markerList2$C6.Mild, "cond"="Mild")
df3 <- data.frame(markerList2$C6.Severe, "cond"="Severe")
df <- rbind(df1,df2,df3)
write.table(df,"peakcond.C6.txt",sep="\t", quote=F, row.names=F)
df$cell = "C6"
dfc6 <- df



####################   var motif

projMono <- addBgdPeaks(projMono)
projMono <- addDeviationsMatrix(
  ArchRProj = projMono, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(projMono, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = projMono, addDOC = FALSE)
motifs <- c("ELF2_326", "IRF4_632", "SMARCC1_651","FOSL2_105", "FOSL2_105","FOSL1_142", "NR2F2_694", "NR2E1_665")
motifs <- c("NR2C1_670", "SPIB_336", "NR4A1_671","NR2F1_692", "RARA_675","BCL11A_194", "BCL11B_825", "SPI1_322","RUNX1_733")
markerMotifs <- getFeatures(projMono, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs


projMono <- addImputeWeights(projMono)
p <- plotGroups(ArchRProj = projMono, 
  groupBy = "cluster.cond", 
  colorBy = "MotifMatrix", 
  name = markerMotifs,
  imputeWeights = getImputeWeights(projMono)
)

p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    ) + ylab("")
  }
})
do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))

plotPDF(p, name = "Monosubs-Groups2-Deviations-w-Imputation", width = 5, height = 5, ArchRProj = projMono, addDOC = FALSE)

########## TF footprint
########## TF footprint
########## TF footprint




#markerMotifs <- getFeatures(projMono, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
#markerMotifs

motifPositions <- getPositions(projMono)
motifs <- c("ELF2_326", "IRF4_632", "SMARCC1_651","FOSL2_105", "FOSL2_105","FOSL1_142", "NR2F2_694", "NR2E1_665")
motifs <- c("NR2C1_670", "SPIB_336", "NR4A1_671","NR2F1_692", "RARA_675","BCL11A_194", "BCL11B_825", "SPI1_322","RUNX1_733")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

encode <- c("RXRA-HepG","102.PAX5_N1","83.BCL11A","278.CEBPB")
markerEncode <- unlist(lapply(encode, function(x) grep(x, names(motifPositions), value = TRUE)))

markerMotifs

seFoot <- getFootprints(
  ArchRProj = projMono, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Clusters"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = projMono, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)


saveArchRProject(ArchRProj= projMono, outputDirectory="projMono7",load=F)

motifs <- c("ELF2_326","NR2F2_694","SPIB_336","BCL11A_194","SPI1_322", "SMARCC1_651","FOSL2_105", "FOSL1_142","JUND_124","JUNB_139","FOS_137","JUN_143")
#motifs <- c("NR2C1_670", "SPIB_336", "NR4A1_671","NR2F1_692", "RARA_675","BCL11A_194", "BCL11B_825", "SPI1_322","RUNX1_733")

markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs

projMono <- addGroupCoverages(ArchRProj = projMono, groupBy = "Severity")
seFoot <- getFootprints(
  ArchRProj = projMono, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Severity"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = projMono, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)



#saveArchRProject(ArchRProj= projMono, outputDirectory="projMono5",load=F)

########################
seGroupMotif <- getGroupSE(ArchRProj = projDm, useMatrix = "MotifMatrix", groupBy = "Clusters2")
seGroupMotif
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

corGSM_MM <- correlateMatrices(
    ArchRProj = projHeme5,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "IterativeLSI"
)
corGIM_MM <- correlateMatrices(
    ArchRProj = projHeme5,
    useMatrix1 = "GeneIntegrationMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "IterativeLSI"
)

corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]


corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])
p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )

p



corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])
p <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGIM_MM$maxDelta)*1.05)
  )

p

#############

projMono <- addCoAccessibility(
    ArchRProj = projMono,
    reducedDims = "IterativeLSI"
)
markerGenes <- c("ABO","CXCR6")
markerGenes <- c("LUCAT1","LINC00211")
p <- plotBrowserTrack(
	ArchRProj= projMono,
	groupBy="cluster.cond",
	geneSymbol=markerGenes,
	upstream=100000,
	downstream=100000,
	loops = getPeak2GeneLinks(projMono)
)
plotPDF(plotList=p, name="plot_track_Mono_lincRNA_cond.pdf",ArchRProj= projMono, addDOC=F, width=10,height=8)
plotPDF(plotList=p, name="plot_track_Mono_subs_coacc.pdf",ArchRProj= projMono, addDOC=F, width=10,height=5)

plotPDF(plotList=p, name="plot_CoAccess_track_Mono_OAS3.pdf",ArchRProj= projMono, addDOC=F, width=10,height=10)

projMono <- addPeak2GeneLinks(
    ArchRProj = projMono,
    reducedDims = "IterativeLSI"
)

p2g <- getPeak2GeneLinks(
    ArchRProj = projMono,
    returnLoops = FALSE
)

p2g <- getPeak2GeneLinks(
    ArchRProj = projMono,
    corCutOff = 0.45,
    resolution = 1,
    returnLoops = FALSE
)

markerGenes <- c("S100A10","IL10RA","IL15RA")
p <- plotBrowserTrack(
	ArchRProj= projMono,
	groupBy="Clusters",
	geneSymbol=markerGenes,
	upstream=30000,
	downstream=30000,
	loops = getPeak2GeneLinks(projMono)
)
plotPDF(plotList=p, name="plot_track_Mono_subs_links.pdf",ArchRProj= projMono, addDOC=F, width=10,height=5)

saveArchRProject(ArchRProj= projMono, outputDirectory="projMono6",load=F)



markerGenes <- c("FOSL1","IFI6","IFITM3")
p <- plotBrowserTrack(
	ArchRProj= projMono,
	groupBy="severity.patient",
	geneSymbol=markerGenes,
	upstream=50000,
	downstream=50000,
	loops = getPeak2GeneLinks(projMono)
)
plotPDF(plotList=p, name="test_track_links.pdf",ArchRProj= projMono, addDOC=F, width=10,height=10)



mono.pools <- readRDS("/vol/projects/CIIM/Covid_50MHH/scRNA/analysis/mono.pools.regMT.nF2k.dim10.reg0.3.rds")
projMono <- addGeneIntegrationMatrix(
    ArchRProj = projMono, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = mono.pools,
    addToArrow = FALSE,
    groupRNA = "seurat_clusters",
    nameCell = "predictedCell_Mono",
    nameGroup = "predictedGroup_Mono",
    nameScore = "predictedScore_Mono"
)

projMono <- addGeneIntegrationMatrix(
    ArchRProj = projMono, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    sampleCellsATAC=11548,
    sampleCellsRNA=30000,
    seRNA = mono.pools,
    addToArrow = FALSE,
    groupRNA = "seurat_clusters",
    nameCell = "predictedCell_Mono2",
    nameGroup = "predictedGroup_Mono2",
    nameScore = "predictedScore_Mono2"
)

projMono <- addGeneIntegrationMatrix(
    ArchRProj = projMono, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    sampleCellsATAC=11548,
    sampleCellsRNA=64003,
    seRNA = mono.pools,
    addToArrow = TRUE,
    groupRNA = "new.id",
    nameCell = "predictedCell_Mono3",
    nameGroup = "predictedGroup_Mono3",
    nameScore = "predictedScore_Mono3"
)

projMono <- addGeneIntegrationMatrix(
    ArchRProj = projMono, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    sampleCellsATAC=11548,
    sampleCellsRNA=64003,
    seRNA = mono.pools,
    addToArrow = TRUE,
 	force =TRUE,
    groupRNA = "new.id",
    nameCell = "predictedCell_Mono4",
    nameGroup = "predictedGroup_Mono4",
    nameScore = "predictedScore_Mono4"
)

projMono$predictedGroup_Monoflt <- projMono$predictedGroup_Mono
projMono$predictedGroup_Monoflt[projMono$predictedScore_Mono < 0.6] = "undefine"
table(projMono$predictedGroup_Mono)

projMono$predictedGroup_Monoflt2 <- projMono$predictedGroup_Mono2
projMono$predictedGroup_Monoflt2[projMono$predictedScore_Mono2 < 0.6] = "undefine"

projMono$predictedGroup_Monoflt3 <- projMono$predictedGroup_Mono3
projMono$predictedGroup_Monoflt3[projMono$predictedScore_Mono3 < 0.6] = "undefine"



projMono <- addClusters(
	input= projMono,
	reducedDims="IterativeLSI",
	method="Seurat",
	name="Clusters",
	resolution=0.8,
	force=TRUE
)

p1 <- plotEmbedding(ArchRProj = projMono, colorBy = "cellColData", name = "Clusters", embedding="UMAP")
p2 <- plotEmbedding(ArchRProj = projMono, colorBy = "cellColData", name = "predictedGroup_Monoflt3", embedding="UMAP")

plotPDF(p1,p2, name = "Plots-UMAP_linkMono3.pdf", ArchRProj = projMono,addDoc=F, width=5, height=5)


cM <- as.matrix(confusionMatrix(projMono$Clusters, projMono$predictedGroup_Mono))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments


test <- read.table("~/Downloads/Book2.txt",row.names = 1, header=T,sep="\t",stringsAsFactors = F)
test2 <- cM/rowSums(cM)
library(pheatmap)
#heatmap(test2)
pdf("Match_heatmap1.pdf",width=5,height=4)

dev.off()


#####################

###################

markerGenes  <- c(
    "CD14","FCGR3A","OAS3","CD163",
    "BATF3" ,"CDKN1C", "LINC02432", "MEG3", "NECAB3"
  )  


p1 <- plotEmbedding(
    ArchRProj = projMono, 
    colorBy = "GeneIntegrationMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projMono)
)

p2 <- lapply(p1, function(x){
x + guides(color = FALSE, fill = FALSE) + theme_ArchR(baseSize = 6.5) +
theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + theme(
axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()
) })
p3 <- do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

plotPDF(p3,
name = "Plot-UMAP-Mono2-W-Integration.pdf", ArchRProj = projMono,
addDOC = FALSE, width = 10, height = 12)


p <- plotEmbedding(
    ArchRProj = projMono, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projMono)
)


p2 <- lapply(p, function(x){
x + guides(color = FALSE, fill = FALSE) + theme_ArchR(baseSize = 6.5) +
theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + theme(
axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()
) })
p3 <- do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

plotPDF(p3,
name = "Plot-UMAP-Mono2-W-Imputation.pdf", ArchRProj = projMono,
addDOC = FALSE, width = 10, height = 12)

####################

seGroupMotif <- getGroupSE(ArchRProj = projMono, useMatrix = "MotifMatrix", groupBy = "Clusters")

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs


corGSM_MM <- correlateMatrices(
    ArchRProj = projMono,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "IterativeLSI"
)

corGIM_MM <- correlateMatrices(
    ArchRProj = projMono,
    useMatrix1 = "GeneIntegrationMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "IterativeLSI"
)

corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.1 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.5))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])


p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_bw() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )
p



corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])

library(ggrepel)
p <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  geom_label_repel(aes(label = ifelse(TFRegulator  == "YES", as.character(GeneIntegrationMatrix_name), '')), hjust = 1.25, vjust = 0, size = 3) +
  theme_bw() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGIM_MM$maxDelta)*1.05)
  )

p



###########################
###########################
peaks <- readRDS("./peaks.signac.rds")
projMono <- addPeakSet(projMono, peakSet=peaks,force=TRUE)
projMono <- addPeakMatrix(projMono)

saveArchRProject(ArchRProj= projMono, outputDirectory="projMonoPeak2",load=F)


markerTest.cluster <- getMarkerFeatures(
  ArchRProj = projMono, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
)

markerList <- getMarkers(markerTest.cluster, cutOff = "FDR <= 0.1", returnGR=F)
markerList <- getMarkers(markerTest.cluster, cutOff = "Pval <= 0.05 & Log2FC >= 1", returnGR=F)
markerList <- getMarkers(markerTest.cluster, cutOff = "FDR <= 0.1 & Log2FC >= 0.5", returnGR=F)
markerList$C4

markerTest.cluster.cond <- getMarkerFeatures(
  ArchRProj = projMono, 
  useMatrix = "PeakMatrix",
  groupBy = "cluster.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
)
markerList <- getMarkers(markerTest.cluster.cond, cutOff = "Pval <= 0.05 & Log2FC >= 1", returnGR=F)
markerList$C4.Severe


peakset = getPeakSet(projMono)
promoter_peak = peakset[which(peakset$peakType == "Promoter"), ]
projMono = addFeatureMatrix(
  input = projMono,
  features = promoter_peak,
  matrixName = "PromoterMatrix"
)
promoter_cm <- getMarkerFeatures(
  ArchRProj = projMono,
  useMatrix = "PromoterMatrix",
  groupBy = "cell_dis2",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Classical Monocytes_post",
  bgdGroups = "Classical Monocytes_healthy"
)




library(chromVARmotifs)
projMono <- addMotifAnnotations(ArchRProj = projMono, motifSet = "cisbp", name = "Motif", force=TRUE)
projMono <- addArchRAnnotations(ArchRProj = projMono, collection = "EncodeTFBS",force=TRUE)
saveArchRProject(ArchRProj= projMono, outputDirectory="projMonoPeak2",load=F)



enrichMotifs <- peakAnnoEnrichment(
    seMarker = markerTest.cluster.cond,
    ArchRProj = projMono,
    peakAnnotation = "Motif",
    cutOff = "Pval <= 0.05 & Log2FC >= 1"
  )

enrichMotifs
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 5, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap2", width = 8, height = 6, ArchRProj = projMono, addDOC = FALSE)

heatmapEM <- plotEnrichHeatmap(enrichMotifs_bak, n = 7, transpose = TRUE, 
           pal = paletteContinuous(set = "whitePurple", n = 20)[c(1,5,12:20)], 
          clusterCols=F)
#ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-cluster_cond", width = 8, height = 6, ArchRProj = projMono, addDOC = FALSE)
#returnMatrix
heatmapEM <- plotEnrichHeatmap(enrichMotifs_bak, n = 20, transpose = TRUE, returnMatrix=T,
                   pal = paletteContinuous(set = "whitePurple", n = 20)[c(1,5,12:20)],clusterCols=F)

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markerTest.cluster,
    ArchRProj = projMono,
    peakAnnotation = "Motif",
    cutOff = "Pval <= 0.05 & Log2FC >= 1"
  )

enrichMotifs
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE, 
           pal = paletteContinuous(set = "whitePurple", n = 20)[c(1,5,12:20)], 
          clusterCols=F)
#ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap4", width = 8, height = 6, ArchRProj = projMono, addDOC = FALSE)


heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 20, transpose = TRUE,  returnMatrix=T,
           pal = paletteContinuous(set = "whitePurple", n = 20)[c(1,5,12:20)], 
          clusterCols=F)




#############



###############


motifPositions <- getPositions(projMono)
#motifs <- c("SPIB_336", "SPIC_344","POU3F4_619", "BCL11A_194", "RUNX1_733", "STAT2_778",
# "FOS_137", "JUNB_139","JUN_143","FOSL2_105","FOSL1_142","JUND_124","CEBPA_155","CEBPB_140")
#motifs <- c("POU2F2","SPI1","SP3","ZBTB7A","FOSL2","JUNB","JUND","JDP2","BACH1")
motifs <- c("CEBPA_155","CEBPB_140", "CEBPD")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
#markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs


projMono <- addGroupCoverages(ArchRProj = projMono, groupBy = "cluster.cond",force=T)

seFoot <- getFootprints(
  ArchRProj = projMono, 
  positions = motifPositions[markerMotifs], 
  groupBy = "cluster.cond"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = projMono, 
  normMethod = "Subtract",
  plotName = "subFootprints-Subtract-Bias5",
  addDOC = FALSE,
  smoothWindow = 5
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = projMono, 
  normMethod = "Subtract",
  plotName = "subFootprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 10
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = projMono, 
  normMethod = "Divide",
  plotName = "subFootprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 10
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = projMono, 
  normMethod = "None",
  plotName = "subFootprints-No-Normalization",
  addDOC = FALSE,
  smoothWindow = 10
)



#######
DotPlot(pbmc.sub, features=uniq(gene.spi), group.by="cell.cond", cols="RdBu")+ coord_flip()