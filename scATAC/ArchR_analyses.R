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
#projDmFlt <- loadArchRProject("projDmFlt2")
#projDm <- loadArchRProject("./projDm3")
##### set up ####
addArchRGenome("hg38")
addArchRThreads(threads=4)


meta <- read.table("/vol/projects/CIIM/Covid_50MHH/scATAC/souporcell/extract.new.txt",header=T,sep="\t",stringsAsFactors=F,comment.char = "")
tmp <- projCovid2$cellNames

meta2 <- meta[meta$cellname %in% tmp,]
projDm <- subsetArchRProject(ArchRProj=projCovid2, cells= meta2$cellname, outputDirectory="projDm")


projDm <- addCellColData(
ArchRProj = projDm,
  data = meta2$patient,
  name = "patient",
  cells = meta2$cellname
)


projDm <- addCellColData(
ArchRProj = projDm,
  data = meta2$id,
  name = "IDs",
  cells = meta2$cellname
)

projDm <- addCellColData(
ArchRProj = projDm,
  data = meta2$gender,
  name = "Gender",
  cells = meta2$cellname
)

projDm <- addCellColData(
ArchRProj = projDm,
  data = meta2$age,
  name = "Age",
  cells = meta2$cellname
)


projDm <- addCellColData(
ArchRProj = projDm,
  data = meta2$Severity2,
  name = "Severity",
  cells = meta2$cellname
)

#saveArchRProject(ArchRProj= projDm, outputDirectory="projDm",load=F)
#### clustering

projDm <- addIterativeLSI(
	ArchRProj= projDm,
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

projDm <- addHarmony(
	ArchRProj= projDm,
	reducedDims="IterativeLSI",
	name="Harmony",
	groupBy="Sample",
	force=TRUE
)

projDm <- addClusters(
	input= projDm,
	reducedDims="IterativeLSI",
	method="Seurat",
	name="Clusters",
	resolution=0.8,
	force=TRUE
)

projDm <- addUMAP(
	ArchRProj= projDm,
	reducedDims="IterativeLSI",
	name="UMAP",
	nNeighbors=30,
	minDist=0.5,
	metric="cosine",
	force=TRUE
)



projDm <- addClusters(
	input= projDm,
	reducedDims="Harmony",
	method="Seurat",
	name="ClustersHm",
	resolution=0.8,
	force=TRUE
)


projDm <- addUMAP(
	ArchRProj= projDm,
	reducedDims="Harmony",
	name="UMAPHarmony",
	nNeighbors=30,
	minDist=0.5,
	metric="cosine",
	force=TRUE
)

saveArchRProject(ArchRProj= projDm, outputDirectory="projDm",load=F)

p1 <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", name = "Sample", embedding="UMAP")
p2 <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", name = "Clusters", embedding="UMAP")
p3 <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", name = "Sample", embedding="UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", name = "ClustersHm", embedding="UMAPHarmony")
p5 <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", name = "Clusters2", embedding="UMAPHarmony")

plotPDF(p1,p2,p3,p4, name = "Plots-UMAP2Harmony-Sample-Clusters.pdf", ArchRProj = projDm,addDoc=F, width=5, height=5)

plotPDF(p4,p5, name = "Plots-UMAP2Harmony-Clusters-celltype.pdf", ArchRProj = projDm,addDoc=F, width=5, height=5)

p2 <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", name = "ClustersHm", embedding="UMAPHarmony")

plotPDF(p2, name = "Plot-UMAPHm-Clusters.pdf", ArchRProj = projDm, addDoc=F, width=4, height=4)


#######################
idx <- BiocGenerics::which(projDm$Severity %in% "mild")
mild <- projDm$cellNames[idx]

idx <- BiocGenerics::which(projDm$Severity %in% "severe")
severe <- projDm$cellNames[idx]

idx <- BiocGenerics::which(projDm$Severity %in% "post")
post <- projDm$cellNames[idx]

p1 <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =mild, name = "Clusters2", embedding="UMAPHarmony",labelAsFactors=F,size = 0.8,labelSize =0)
p2 <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells = severe, name = "Clusters2", embedding="UMAPHarmony",labelAsFactors=F,size = 0.8,labelSize =0)
p3 <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells = post, name = "Clusters2", embedding="UMAPHarmony",labelAsFactors=F,size = 0.8,labelSize =0)
plotPDF(p1,p2,p3, name = "Plots-UMAP_Harmony_BYcondition.pdf", ArchRProj = projDm,addDoc=F, width=5, height=5)





#######################

idx <- BiocGenerics::which(projDm$Sample %in% "pool1a")
pool1a <- projDm$cellNames[idx]

idx <- BiocGenerics::which(projDm$Sample %in% "pool1b")
pool1b <- projDm$cellNames[idx]

idx <- BiocGenerics::which(projDm$Sample %in% "pool1d")
pool1d <- projDm$cellNames[idx]

idx <- BiocGenerics::which(projDm$Sample %in% "pool2a")
pool2a <- projDm$cellNames[idx]

idx <- BiocGenerics::which(projDm$Sample %in% "pool2d")
pool2d <- projDm$cellNames[idx]

idx <- BiocGenerics::which(projDm$Sample %in% "pool3a")
pool3a <- projDm$cellNames[idx]

idx <- BiocGenerics::which(projDm$Sample %in% "pool3b")
pool3b <- projDm$cellNames[idx]

idx <- BiocGenerics::which(projDm$Sample %in% "pool3d")
pool3d <- projDm$cellNames[idx]

idx <- BiocGenerics::which(projDm$Sample %in% "pool4a")
pool4a <- projDm$cellNames[idx]

idx <- BiocGenerics::which(projDm$Sample %in% "pool4b")
pool4b <- projDm$cellNames[idx]

idx <- BiocGenerics::which(projDm$Sample %in% "pool4c")
pool4c <- projDm$cellNames[idx]

idx <- BiocGenerics::which(projDm$Sample %in% "pool4d")
pool4d <- projDm$cellNames[idx]

p1a <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool1a, name = "ClustersHm",embedding="UMAPHarmony",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool1a")
p1b <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool1b, name = "ClustersHm",embedding="UMAPHarmony",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool1b")
p1d <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool1d, name = "ClustersHm",embedding="UMAPHarmony",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool1d")
p2a <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool2a, name = "ClustersHm",embedding="UMAPHarmony",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool2a")
p2d <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool2d, name = "ClustersHm",embedding="UMAPHarmony",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool2d")
p3a <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool3a, name = "ClustersHm",embedding="UMAPHarmony",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool3a")
p3b <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool3b, name = "ClustersHm",embedding="UMAPHarmony",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool3b")
p3d <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool3d, name = "ClustersHm",embedding="UMAPHarmony",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool3d")
p4a <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool4a, name = "ClustersHm",embedding="UMAPHarmony",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool4a")
p4b <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool4b, name = "ClustersHm",embedding="UMAPHarmony",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool4b")
p4c <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool4c, name = "ClustersHm",embedding="UMAPHarmony",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool4c")
p4d <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool4d, name = "ClustersHm",embedding="UMAPHarmony",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool4d")

plotPDF(
p1a,
p1b,
p1d,
p2a,
p2d,
p3a,
p3b,
p3d,
p4a,
p4b,
p4c,
p4d,name = "Plots-UMAP_Harmony_BYpool.pdf", ArchRProj = projDm,addDoc=F, width=5, height=5)

ggAlignPlots(p1a,p1b,p1d,p2a, type="h")
ggAlignPlots(p2d, p3a, p3b, p3d, type="h")
ggAlignPlots(p4a, p4b, p4c, p4d, type="h")



p1a <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool1a, name = "Clusters",embedding="UMAP",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool1a")
p1b <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool1b, name = "Clusters",embedding="UMAP",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool1b")
p1d <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool1d, name = "Clusters",embedding="UMAP",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool1d")
p2a <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool2a, name = "Clusters",embedding="UMAP",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool2a")
p2d <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool2d, name = "Clusters",embedding="UMAP",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool2d")
p3a <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool3a, name = "Clusters",embedding="UMAP",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool3a")
p3b <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool3b, name = "Clusters",embedding="UMAP",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool3b")
p3d <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool3d, name = "Clusters",embedding="UMAP",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool3d")
p4a <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool4a, name = "Clusters",embedding="UMAP",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool4a")
p4b <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool4b, name = "Clusters",embedding="UMAP",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool4b")
p4c <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool4c, name = "Clusters",embedding="UMAP",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool4c")
p4d <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", highlightCells =pool4d, name = "Clusters",embedding="UMAP",labelAsFactors=F,size = 0.8,labelSize =0)+ggtitle("pool4d")

ggAlignPlots(p1a,p1b,p1d,p2a, type="h")
ggAlignPlots(p2d, p3a, p3b, p3d, type="h")
ggAlignPlots(p4a, p4b, p4c, p4d, type="h")


plotPDF(
p1a,
p1b,
p1d,
p2a,
p2d,
p3a,
p3b,
p3d,
p4a,
p4b,
p4c,
p4d,name = "Plots-UMAP_raw_BYpool.pdf", ArchRProj = projDm,addDoc=F, width=5, height=5)


#######################
#######################



markersGS <- getMarkerFeatures(
	ArchRProj=projDm,
	useMatrix="GeneScoreMatrix",
	groupBy="ClustersHm",
	bias=c("TSSEnrichment","log10(nFrags)"),
	testMethod="wilcoxon"

)

markerList <- getMarkers(markersGS,cutOff="FDR <= 0.05 & Log2FC >= 0.7")




marker.gene <- c(
"IL7R","CCR7","TCF7",  #IL7R,CCR7-,TCF7-             Memory CD4 T   
"CD8A","CD8B","GZMK", #CD8A+,CD8B+,GZMK+             CD8 T    
"NKG7","NCAM1", #NKG7+,NCAM1+,CD8A-             NK    
"CD14","LYZ","FCGR3A", #CD14+, LYZ+, FCGR3A-            CD14+ monocytes   
#FCGR3A+, CD14-              CD16+ monocytes   
"CST3","CD86","HLA-DRA",  #CST3+, CD86+, HLA-DRA+, CD14-, FCGR3A-   mDC, pDC (dendritic cells)    
"IL1B","TNF","IFI6","IFITM3",  #IL1B,TNF,IL6               proinflammatory markers  # macrophages
"IFNG","GZMB","CCL5","CST7","CLEC9A","CD1D","CXCL9","CXCL10","CXCL11",
"CD79A","CD27","SDC1", #CD79A+, CD27-, SDC1-           B cell              
"PPBP", #PPBP+                                 megakaryocyte/platelet   
)

heatmapGS <- markerHeatmap(
	seMarker = markersGS,
	cutOff = "FDR<0.01 & Log2FC >= 1",
	labelMarkers=marker.gene,
	transpose=T
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side="bot", annotation_legend_side="bot")


plotPDF(heatmapGS, name = "Plot-geneScores.pdf", ArchRProj = projDm,addDoc=F, width=8, height=6)



projDm <- addImputeWeights(projDm)


#plotGene <- c("IL1B","IFITM3","CD3E","CD14","FCGR3A","IFNG","HLA-DRA","CCL5","GZMB")
plotGene <- c("IL7R","CCR7","CD8A","GZMK","NKG7","NCAM1","IFNG","HLA-DRA","CD86")

plotGene <- c("CD14","FCGR3A","CD68","CD163","HLA-DRA",
             "CD3E","CD3D","IL7R","CD8A","IFNG",
             "GZMB","MKI67","CD79A","MS4A1","CD74")


p <- plotEmbedding( 
	ArchRProj = projDm,
	colorBy="GeneScoreMatrix",
	name=plotGene,
	embedding="UMAPHarmony",
	imputeWeights = getImputeWeights(projDm)
)


p2 <- lapply(p, function(x){
x + guides(color = FALSE, fill = FALSE) + theme_ArchR(baseSize = 6.5) +
theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + theme(
axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()
) })
p3 <- do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

plotPDF(p3,
name = "Plot-UMAP-Marker-Genes2-W-Imputation.pdf", ArchRProj = projDm,
addDOC = FALSE, width = 10, height = 10)
plotPDF(plotList = p,
name = "Plot-UMAP-Marker-Genes-W-Imputation_each.pdf", ArchRProj = projDm,
addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj= projDm, outputDirectory="projDm2",load=F)


plotGene <- c("CD14","FCGR3A",  # monocytes markers 
             "IFITM3","IFI6",   # Interferon markers
             "HLA-DRA","CD163","CD68",    #HLA,  CD163, CD68
             "CCR7","IL7R","CD3E","CD3D",  # T cells
			"MKI67",  ## 
			"CCL5","GZMB","CD8A","NKG7",  #CD8 NK
			"MS4A1","CD74","CD79A","CD27","SDC1")  #B cells ...
			
			
p <- plotEmbedding( 
	ArchRProj = projDm,
	colorBy="GeneScoreMatrix",
	name=plotGene,
	embedding="UMAPHarmony",
	imputeWeights = getImputeWeights(projDm)
)
plotPDF(plotList = p,
name = "Plot-UMAP-select-Genes-W-Imputation_each.pdf", ArchRProj = projDm,
addDOC = FALSE, width = 5, height = 5)


p <- plotEmbedding( 
	ArchRProj = projDm,
	colorBy="GeneScoreMatrix",
	name=plotGene,
	embedding="UMAP",
	imputeWeights = getImputeWeights(projDm)
)
plotPDF(plotList = p,
name = "Plot-UMAPraw-select-Genes-W-Imputation_each.pdf", ArchRProj = projDm,
addDOC = FALSE, width = 5, height = 5)


scRNA <- readRDS("../../scRNA/analysis/pbmc.harmony.annotev2.rds")



projDm <- addGeneIntegrationMatrix(
    ArchRProj = projDm, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = scRNA,
    addToArrow = FALSE,
    groupRNA = "celltypeL0",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)


cM <- as.matrix(confusionMatrix(projDm$ClustersHm, projDm$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments


df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment"))
#df2 <- df[df$TSSEnrichment >3,]
p <- ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

pdf("test1.pdf")
plot(p)
dev.off()





saveArchRProject(ArchRProj= projDm, outputDirectory="projDm2",load=F)

tmp <- names(grep("cMono",scRNA$celltypeL0,value=T))
rna.ncmono <- names(grep("ncMono",scRNA$celltypeL0,value=T))
rna.cmono <- setdiff(tmp,rna.ncmono)

rna.tm <- names(grep("Tm",scRNA$celltypeL0,value=T))
rna.b <- names(grep("B",scRNA$celltypeL0,value=T))

rna.nk <- names(grep("NK",scRNA$celltypeL0,value=T))
rna.pb <- names(grep("Plasmablast",scRNA$celltypeL0,value=T))
rna.mdc <- names(grep("mDC",scRNA$celltypeL0,value=T))
rna.pdc <- names(grep("pDC",scRNA$celltypeL0,value=T))

rna.un <- names(grep("undefined|Megakaryocyte|Erythrocytes|mDC",scRNA$celltypeL0,value=T))
#Erythrocytes Megakaryocyte     undefined

groupList <- SimpleList(
    cMono = SimpleList(
        ATAC = projDm$cellNames[projDm$ClustersHm %in% c("C12")],
        RNA = rna.cmono
    ),
    ncMono = SimpleList(
        ATAC = projDm$cellNames[projDm$ClustersHm %in% c("C13")],
        RNA = rna.ncmono
    ),
    Ts = SimpleList(
        ATAC = projDm$cellNames[projDm$ClustersHm %in% c("C6","C10","C8","C9")],
        RNA = rna.tm
    ),
    Plasmablast = SimpleList(
        ATAC = projDm$cellNames[projDm$ClustersHm %in% c("C4")],
        RNA = rna.pb
    ),
    B = SimpleList(
        ATAC = projDm$cellNames[projDm$ClustersHm %in% c("C3")],
        RNA = rna.b
    ),
    NK = SimpleList(
        ATAC = projDm$cellNames[projDm$ClustersHm %in% c("C7")],
        RNA = rna.nk
    ),
    pDC = SimpleList(
        ATAC = projDm$cellNames[projDm$ClustersHm %in% c("C5")],
        RNA = rna.pdc
    ),
    undefined = SimpleList(
        ATAC = projDm$cellNames[projDm$ClustersHm %in% c("C11","C15","C14","C2","C1")],
        RNA = rna.un
    )
)


projDm <- addGeneIntegrationMatrix(
    ArchRProj = projDm, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = scRNA,
    addToArrow = TRUE,
    force=TRUE,
    groupList = groupList,
    groupRNA = "celltypeL0",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore",
    k.filter = 30
)




cM <- confusionMatrix(projDm$ClustersHm, projDm$predictedGroup)
labelOld <- rownames(cM)
labelOld

labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew


preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments

labelNew[c(5,6,13,14,15)]="undefine"
labelNew[labelNew == "Tm"]="CD4 T"
labelNew[4]="CD8 T"

#idx <- BiocGenerics::which(projDm$ClustersHm %in% c("C2","C3","C5","C6","C7","C8","C9","C10","C15"))
#flt.cells <- projDm$cellNames[idx]
#projDmFlt <- subsetArchRProject(ArchRProj=projDm, cells= flt.cells, outputDirectory="projDmFlt")
#projDmFlt$Clusters2 <- mapLabels(projDmFlt$ClustersHm, newLabels = labelNew, oldLabels = labelOld)
#projDmFlt <- addUMAP(
#	ArchRProj= projDmFlt,
#	reducedDims="Harmony",
#	name="UMAPHarmony2",
#	nNeighbors=30,
#	minDist=0.5,
#	metric="cosine"
#)
#p2 <- plotEmbedding(ArchRProj = projDmFlt, colorBy = "cellColData", name = "Clusters2", embedding="UMAPHarmony2")
#p3 <- plotEmbedding(ArchRProj = projDmFlt, colorBy = "cellColData", name = "Severity", embedding="UMAPHarmony2")
#plotPDF(p2,p3, name = "Plot-UMAPHm-ClustersNew.pdf", ArchRProj = projDmFlt, addDoc=F, width=4, height=4)






#### raw??######
projDm$Clusters2 <- mapLabels(projDm$ClustersHm, newLabels = labelNew, oldLabels = labelOld)
p2 <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", name = "Clusters2", embedding="UMAPHarmony")
p3 <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", name = "Severity", embedding="UMAPHarmony")
plotPDF(p2,p3, name = "Plot-UMAPHm-ClustersNew.pdf", ArchRProj = projDm, addDoc=F, width=4, height=4)
#######

saveArchRProject(ArchRProj= projDm, outputDirectory="projDm3",load=F)

aaa <- confusionMatrix(projDmFlt$Clusters2, projDmFlt$Severity)


#tmp<-data.frame(cell_type= projDmFlt$Clusters2,sample= projDmFlt$patient,donor= projDmFlt$Severity)
tmp<-data.frame(cell_type= as.factor(projDm$Clusters2),sample= projDm$patient,donor= projDm$Severity)
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

y2$Condition <- factor(y2$Condition,levels=c("severe","mild","post"))


###
#test <- read.table("~/Downloads/Book2.txt",row.names = 1, header=T,sep="\t",stringsAsFactors = F)
#test2 <- test/rowSums(test)
#library(pheatmap)
#heatmap(test2)
#pheatmap(test2,cluster_rows = F,cluster_cols = F)

#
#
#


pdf("cpp_merge2.pdf",width=15,height=4.5)
p <- ggplot(data=y2, aes(x=Condition,y=f))+geom_boxplot(aes(fill=Condition),outlier.shape=NA)
p2 <- p + facet_grid(.~cell_type)+geom_jitter(width = 0.2)+ theme_classic()+scale_fill_manual(values=c("#9e0142","#fee08b","#4393c3"))+xlab("")+ylab("percent")
plot(p2)
#p2 + geom_signif(data=dat.sig, aes(xmin=Comp1,xmax=Comp2, annotations=sig, y_position=maxs), manual = TRUE, vjust=0.8, tip_length=0.005)+theme(axis.text.x=element_text(angle = 70,hjust = 1))
#filename <- paste("BoxPlot_cohort2_cluster",i,".pdf",sep="")
dev.off()


y2$other <- 1-y2$f

i<-'B'
  print(i)
  y2.i <- y2[y2$cell_type==i,]
  y2.i$Smp <- DR_data(y2.i[,5:6])
  res <- DirichReg(Smp~Condition,y2.i,model="alternative",base=2)
  x <- summary(res)
  sc <- x$coef.mat["Conditionmild",3:4]
  hc <- x$coef.mat["Conditionpost",3:4]
  dat.test <- rbind(dat.test, c(i,sc,"severe.vs.mild"))
  dat.test <- rbind(dat.test, c(i,hc,"severe.vs.post"))

##> head(y2)

dat.test <- NULL
for(i in levels(as.factor(y2$cell_type))){
#for(i in c("PLA1","PLA2","PLA3")){  
  print(i)
  y2.i <- y2[y2$cell_type==i,]
  y2.i$Smp <- DR_data(y2.i[,5:6])
  res <- DirichReg(Smp~Condition,y2.i,model="alternative",base=2)
  x <- summary(res)
  sc <- x$coef.mat["Conditionmild",3:4]
  hc <- x$coef.mat["Conditionpost",3:4]
  dat.test <- rbind(dat.test, c(i,sc,"severe","mild"))
  dat.test <- rbind(dat.test, c(i,hc,"severe","post"))
  
  y2.i <- y2.i[y2.i$Condition != "severe",]
  y2.i$Condition <- factor(y2.i$Condition, levels=c("mild","post"))
   y2.i$Smp <- DR_data(y2.i[,5:6])
  res <- DirichReg(Smp~Condition,y2.i,model="alternative",base=2)
  x <- summary(res)
  mc <- x$coef.mat["Conditionpost",3:4]
  dat.test <- rbind(dat.test, c(i,mc,"mild","post"))
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
	if(dat.sig[i, "Comp1"]=="severe" & dat.sig[i,"Comp2"]=="mild"){
	        maxs <- c(maxs,max(y2$f[y2$cell_type == c])+0.03)
	}
	if(dat.sig[i, "Comp1"]=="mild" & dat.sig[i,"Comp2"]=="post"){
	        maxs <- c(maxs,max(y2$f[y2$cell_type == c])+0.09)
	}
	if(dat.sig[i, "Comp1"]=="severe" & dat.sig[i,"Comp2"]=="post"){
	        maxs <- c(maxs,max(y2$f[y2$cell_type == c])+0.06)
	}
	}
dat.sig$maxs <- maxs

pdf("cpp_merge3.pdf",width=15,height=4.5)
p <- ggplot(data=y2, aes(x=Condition,y=f))+geom_boxplot(aes(fill=Condition),outlier.shape=NA)
p2 <- p + facet_grid(.~cell_type)+geom_jitter(width = 0.2)+ theme_classic()+scale_fill_manual(values=c("#9e0142","#fee08b","#4393c3"))+xlab("")+ylab("percent")
plot(p2)
p2 + geom_signif(data=dat.sig, aes(xmin=Comp1,xmax=Comp2, annotations=sig, y_position=maxs), manual = TRUE, vjust=0.8, tip_length=0.005)+theme(axis.text.x=element_text(angle = 70,hjust = 1))
dev.off()




projDm <- addGroupCoverages(ArchRProj=projDm, groupBy="Clusters2")
pathToMacs2 <- "/vol/projects/CIIM/resources/tools/encode-atac-seq-pipeline/bin/macs2"
projDm <- addReproduciblePeakSet(
	ArchRProj = projDm,
	groupBy="Clusters2",
	pathToMacs2= pathToMacs2
	)
getPeakSet(projDm)
projDm <- addPeakMatrix(projDm)
saveArchRProject(ArchRProj= projDm, outputDirectory="projDm3",load=F)


mks <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)")
)

markerList <- getMarkers(mks, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
markerList


heatmapPeaks <- markerHeatmap(
  seMarker = mks, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.5",
  transpose = TRUE
)

plotPDF(heatmapPeaks, name = "Peak-cell.cond-Marker-Heatmap", width = 14, height = 10, ArchRProj = projDm, addDOC = FALSE)


mks <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)")
)

markerList <- getMarkers(mks, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
markerList

peak1 <- as.data.frame(markerList$cMono)
peak1$celltypes = "cMono"
peak2 <- as.data.frame(markerList$ncMono)
peak2$celltypes = "ncMono"
peak3 <- as.data.frame(markerList$NK)
peak3$celltypes = "NK"
peak4 <- as.data.frame(markerList$B)
peak4$celltypes = "B"
peak5 <- as.data.frame(markerList$`CD4 T`)
peak5$celltypes = "CD4T"
peak6 <- as.data.frame(markerList$`CD8 T`)
peak6$celltypes = "CD8T"

combinePeaks <- rbind(peak1,peak2,peak3,peak4,peak5,peak6)
write.table(combinePeaks, "Cluster.markers.projDm_main.txt",sep="\t",quote=F,row.names=F)

##

heatmapPeaks <- markerHeatmap(
  seMarker = mks, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.1",
  transpose = TRUE
)

plotPDF(heatmapPeaks, name = "Peak-cell.cond-Marker-Heatmap2", width = 14, height = 10, ArchRProj = projDm, addDOC = FALSE)

heatmapPeaks2 <- markerHeatmap(
  seMarker = mks, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.1",
  transpose = TRUE,
  returnMatrix = TRUE
)


#########

peaks <- getPeakSet(projDm)
peaks2 <- peaks[peaks$peakType == "Promoter",]
projDmSub <- addPeakSet(
  ArchRProj = projDm,
  peakSet = peaks2,
  force = TRUE
)
projDmSub <- addPeakMatrix(projDmSub,force = TRUE)




mk.cmono.mp <- getMarkerFeatures(
  ArchRProj = projDmSub, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "cMono.Mild",
  bgdGroups = "cMono.Post"
)

mk.cmono.sp <- getMarkerFeatures(
  ArchRProj = projDmSub, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "cMono.Severe",
  bgdGroups = "cMono.Post"
)

markerList <- getMarkers.raw(mk.cmono.mp, cutOff = "Pval <= 1", returnGR=F)
#markerList <- getMarkers(mk.cmono.mp, cutOff = "FDR <= 0.1", returnGR=F)
markerList$cMono.Mild
write.table(markerList$cMono.Mild, file="mkPromoter.cMono.MvP.txt",sep="\t",quote=F)

markerList <- getMarkers.raw(mk.cmono.sp, cutOff = "Pval <= 1", returnGR=F)
#markerList <- getMarkers(mk.cmono.sp, cutOff = "FDR <= 0.1", returnGR=F)
markerList$cMono.Severe
write.table(markerList$cMono.Severe, file="mkPromoter.cMono.SvP.txt",sep="\t",quote=F)

mk.cmono.ms <- getMarkerFeatures(
  ArchRProj = projDmSub, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "cMono.Mild",
  bgdGroups = "cMono.Severe"
)

markerList <- getMarkers.raw(mk.cmono.ms, cutOff = "Pval <= 1", returnGR=F)
#markerList <- getMarkers(mk.cmono.mp, cutOff = "FDR <= 0.1", returnGR=F)
markerList$cMono.Mild
write.table(markerList$cMono.Mild, file="mkPromoter.cMono.MvS.txt",sep="\t",quote=F)



###########################



mk.ncmono.mp <- getMarkerFeatures(
  ArchRProj = projDmSub, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "ncMono.Mild",
  bgdGroups = "ncMono.Post"
)

mk.ncmono.sp <- getMarkerFeatures(
  ArchRProj = projDmSub, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "ncMono.Severe",
  bgdGroups = "ncMono.Post"
)

markerList <- getMarkers.raw(mk.ncmono.mp, cutOff = "Pval <= 1", returnGR=F)
#markerList <- getMarkers(mk.cmono.mp, cutOff = "FDR <= 0.1", returnGR=F)
markerList$ncMono.Mild
write.table(markerList$ncMono.Mild, file="mkPromoter.ncMono.MvP.txt",sep="\t",quote=F)

markerList <- getMarkers.raw(mk.ncmono.sp, cutOff = "Pval <= 1", returnGR=F)
#markerList <- getMarkers(mk.cmono.sp, cutOff = "FDR <= 0.1", returnGR=F)
markerList$ncMono.Severe
write.table(markerList$ncMono.Severe, file="mkPromoter.ncMono.SvP.txt",sep="\t",quote=F)

mk.ncmono.ms <- getMarkerFeatures(
  ArchRProj = projDmSub, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "ncMono.Mild",
  bgdGroups = "ncMono.Severe"
)

markerList <- getMarkers.raw(mk.ncmono.ms, cutOff = "Pval <= 1", returnGR=F)
#markerList <- getMarkers(mk.cmono.mp, cutOff = "FDR <= 0.1", returnGR=F)
markerList$ncMono.Mild
write.table(markerList$ncMono.Mild, file="mkPromoter.ncMono.MvS.txt",sep="\t",quote=F)


####################################


mk.mp <- getMarkerFeatures(
  ArchRProj = projDmSub, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CD4 T.Mild",
  bgdGroups = "CD4 T.Post"
)

mk.sp <- getMarkerFeatures(
  ArchRProj = projDmSub, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CD4 T.Severe",
  bgdGroups = "CD4 T.Post"
)

markerList <- getMarkers.raw(mk.mp, cutOff = "Pval <= 1", returnGR=F)
write.table(markerList$`CD4 T.Mild`, file="mkPromoter.cd4t.MvP.txt",sep="\t",quote=F)

markerList <- getMarkers.raw(mk.sp, cutOff = "Pval <= 1", returnGR=F)
write.table(markerList$`CD4 T.Severe`, file="mkPromoter.cd4t.SvP.txt",sep="\t",quote=F)


mk.ms <- getMarkerFeatures(
  ArchRProj = projDmSub, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CD4 T.Mild",
  bgdGroups = "CD4 T.Severe"
)

markerList <- getMarkers.raw(mk.ms, cutOff = "Pval <= 1", returnGR=F)
write.table(markerList$`CD4 T.Mild`, file="mkPromoter.cd4t.MvS.txt",sep="\t",quote=F)

#########
mk.mp <- getMarkerFeatures(
  ArchRProj = projDmSub, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CD8 T.Mild",
  bgdGroups = "CD8 T.Post"
)

mk.sp <- getMarkerFeatures(
  ArchRProj = projDmSub, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CD8 T.Severe",
  bgdGroups = "CD8 T.Post"
)

markerList <- getMarkers.raw(mk.mp, cutOff = "Pval <= 1", returnGR=F)
write.table(markerList$`CD8 T.Mild`, file="mkPromoter.cd8t.MvP.txt",sep="\t",quote=F)

markerList <- getMarkers.raw(mk.sp, cutOff = "Pval <= 1", returnGR=F)
write.table(markerList$`CD8 T.Severe`, file="mkPromoter.cd8t.SvP.txt",sep="\t",quote=F)


mk.ms <- getMarkerFeatures(
  ArchRProj = projDmSub, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CD8 T.Mild",
  bgdGroups = "CD8 T.Severe"
)

markerList <- getMarkers.raw(mk.ms, cutOff = "Pval <= 1", returnGR=F)
write.table(markerList$`CD8 T.Mild`, file="mkPromoter.cd8t.MvS.txt",sep="\t",quote=F)

######################

#########
mk.mp <- getMarkerFeatures(
  ArchRProj = projDmSub, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Plasmablast.Mild",
  bgdGroups = "Plasmablast.Post"
)

mk.sp <- getMarkerFeatures(
  ArchRProj = projDmSub, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Plasmablast.Severe",
  bgdGroups = "Plasmablast.Post"
)

markerList <- getMarkers.raw(mk.mp, cutOff = "Pval <= 1", returnGR=F)
write.table(markerList$`Plasmablast.Mild`, file="mkPromoter.Plasmablast.MvP.txt",sep="\t",quote=F)

markerList <- getMarkers.raw(mk.sp, cutOff = "Pval <= 1", returnGR=F)
write.table(markerList$`Plasmablast.Severe`, file="mkPromoter.Plasmablast.SvP.txt",sep="\t",quote=F)


mk.ms <- getMarkerFeatures(
  ArchRProj = projDmSub, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Plasmablast.Mild",
  bgdGroups = "Plasmablast.Severe"
)

markerList <- getMarkers.raw(mk.ms, cutOff = "Pval <= 1", returnGR=F)
write.table(markerList$`Plasmablast.Mild`, file="mkPromoter.Plasmablast.MvS.txt",sep="\t",quote=F)

#########
mk.mp <- getMarkerFeatures(
  ArchRProj = projDmSub, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "NK.Mild",
  bgdGroups = "NK.Post"
)

mk.sp <- getMarkerFeatures(
  ArchRProj = projDmSub, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "NK.Severe",
  bgdGroups = "NK.Post"
)

markerList <- getMarkers.raw(mk.mp, cutOff = "Pval <= 1", returnGR=F)
write.table(markerList$`NK.Mild`, file="mkPromoter.NK.MvP.txt",sep="\t",quote=F)

markerList <- getMarkers.raw(mk.sp, cutOff = "Pval <= 1", returnGR=F)
write.table(markerList$`NK.Severe`, file="mkPromoter.NK.SvP.txt",sep="\t",quote=F)


mk.ms <- getMarkerFeatures(
  ArchRProj = projDmSub, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "NK.Mild",
  bgdGroups = "NK.Severe"
)

markerList <- getMarkers.raw(mk.ms, cutOff = "Pval <= 1", returnGR=F)
write.table(markerList$`NK.Mild`, file="mkPromoter.NK.MvS.txt",sep="\t",quote=F)

#####################
saveArchRProject(ArchRProj= projDmSub, outputDirectory="projDmSub",load=F)


#idx <- BiocGenerics::which(projDmFlt$ClustersHm %in% c("C2","C3"))
#cmono.cells <- projDmFlt$cellNames[idx]
#projcMono <- subsetArchRProject(ArchRProj= projDmFlt, cells= mono.cells, outputDirectory="projMono")

projDm$cell.cond <- paste(projDm$Clusters2, projDm$Severity,sep=".")

table(projDm$cell.cond)
 #           B.Mild             B.Post           B.Severe         cMono.Mild 
 #              240                239                146                352 
 #       cMono.Post       cMono.Severe           mDC.Mild           mDC.Post 
 #              402                602                 14                 13 
 #       mDC.Severe           mix.Mild           mix.Post         mix.Severe 
 #                8                434                472                319 
 #      ncMono.Mild        ncMono.Post      ncMono.Severe            NK.Mild 
 #               44                 80                 49                359 
 #          NK.Post          NK.Severe           pDC.Mild           pDC.Post 
 #              383                386                  4                 10 
 ##       pDC.Severe   Plasmablast.Mild   Plasmablast.Post Plasmablast.Severe 
#                 5                 13                 12                 33 
#            T.Mild             T.Post           T.Severe 
#              1255               1353                732 


markerTest.cmono <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "cMono.Mild",
  bgdGroups = "cMono.Severe"
)

markerList <- getMarkers(markerTest.cmono, cutOff = "FDR <= 0.1", returnGR=F)
head(markerList$cMono.Mild)

markerList <- getMarkers.raw(markerTest.cmono, cutOff = "Pval <= 0.05 & abs(Log2FC) >= 1", returnGR=F)
head(markerList$cMono.Mild)
write.table(markerList$cMono.Mild, file="mk.cMono.MvS.txt",sep="\t",quote=F,)

pma <- plotMarkers(seMarker = markerTest.cmono, name = "cMono.Mild", cutOff = "Pval <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
plotPDF(pma, name = "Mild-vs-Seve-cMono-MA", width = 5, height = 5, ArchRProj = projDm, addDOC = FALSE)





markerTest.ncmono <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "ncMono.Mild",
  bgdGroups = "ncMono.Severe"
)

markerList <- getMarkers(markerTest.ncmono, cutOff = "FDR <= 0.1", returnGR=F)
head(markerList$ncMono.Mild)


markerList <- getMarkers.raw(markerTest.ncmono, cutOff = "Pval <= 0.05 & abs(Log2FC) >= 1", returnGR=F)
head(markerList$ncMono.Mild)
write.table(markerList$ncMono.Mild, file="mk.ncMono.MvS.txt",sep="\t",quote=F,)

pma <- plotMarkers(seMarker = markerTest.ncmono, name = "ncMono.Mild", cutOff = "Pval <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
plotPDF(pma, name = "Mild-vs-Seve-ncMono-MA", width = 5, height = 5, ArchRProj = projDm, addDOC = FALSE)

markerTest.nk <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "NK.Mild",
  bgdGroups = "NK.Severe"
)

markerList <- getMarkers(markerTest.nk, cutOff = "FDR <= 0.1", returnGR=F)
head(markerList$NK.Mild)


markerList <- getMarkers.raw(markerTest.nk, cutOff = "Pval <= 0.05 & abs(Log2FC) >= 1", returnGR=F)
head(markerList$NK.Mild)

write.table(markerList$NK.Mild, file="mk.NK.MvS.txt",sep="\t",quote=F,)

pma <- plotMarkers(seMarker = markerTest.nk, name = "NK.Mild", cutOff = "Pval <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
plotPDF(pma, name = "Mild-vs-Seve-NK-MA", width = 5, height = 5, ArchRProj = projDm, addDOC = FALSE)

markerTest.b <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "B.Mild",
  bgdGroups = "B.Severe"
)

markerList <- getMarkers(markerTest.b, cutOff = "FDR <= 0.1", returnGR=F)
head(markerList$B.Mild)


markerList <- getMarkers.raw(markerTest.b, cutOff = "Pval <= 0.05 & abs(Log2FC) >= 1", returnGR=F)
head(markerList$B.Mild)

write.table(markerList$B.Mild, file="mk.B.MvS.txt",sep="\t",quote=F,)

pma <- plotMarkers(seMarker = markerTest.b, name = "B.Mild", cutOff = "Pval <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
plotPDF(pma, name = "Mild-vs-Seve-B-MA", width = 5, height = 5, ArchRProj = projDm, addDOC = FALSE)

############
markerTest.cd4t <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CD4 T.Mild",
  bgdGroups = "CD4 T.Severe"
)

markerList <- getMarkers(markerTest.cd4t, cutOff = "FDR <= 0.1", returnGR=F)
head(markerList$`CD4 T.Mild`)


markerList <- getMarkers.raw(markerTest.cd4t, cutOff = "Pval <= 0.05 & abs(Log2FC) >= 1", returnGR=F)
head(markerList$`CD4 T.Mild`)

write.table(markerList$`CD4 T.Mild`, file="mk.cd4t.MvS.txt",sep="\t",quote=F,)

pma <- plotMarkers(seMarker = markerTest.cd4t, name = "CD4 T.Mild", cutOff = "Pval <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
plotPDF(pma, name = "Mild-vs-Seve-cd4T-MA", width = 5, height = 5, ArchRProj = projDm, addDOC = FALSE)

###########

############
markerTest.cd8t <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CD8 T.Mild",
  bgdGroups = "CD8 T.Severe"
)

markerList <- getMarkers(markerTest.cd8t, cutOff = "FDR <= 0.1", returnGR=F)
head(markerList$`CD8 T.Mild`)


markerList <- getMarkers.raw(markerTest.cd8t, cutOff = "Pval <= 0.05 & abs(Log2FC) >= 1", returnGR=F)
head(markerList$`CD8 T.Mild`)

write.table(markerList$`CD8 T.Mild`, file="mk.cd8t.MvS.txt",sep="\t",quote=F,)

pma <- plotMarkers(seMarker = markerTest.cd8t, name = "CD8 T.Mild", cutOff = "Pval <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
plotPDF(pma, name = "Mild-vs-Seve-cd8T-MA", width = 5, height = 5, ArchRProj = projDm, addDOC = FALSE)
###########


############
markerTest.pla <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Plasmablast.Mild",
  bgdGroups = "Plasmablast.Severe"
)

markerList <- getMarkers(markerTest.pla, cutOff = "FDR <= 0.1", returnGR=F)
head(markerList$Plasmablast.Mild)


markerList <- getMarkers.raw(markerTest.pla, cutOff = "Pval <= 0.05 & abs(Log2FC) >= 1", returnGR=F)
head(markerList$Plasmablast.Mild)

write.table(markerList$Plasmablast.Mild, file="mk.plasma.MvS.txt",sep="\t",quote=F,)

pma <- plotMarkers(seMarker = markerTest.pla, name = "Plasmablast.Mild", cutOff = "Pval <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
plotPDF(pma, name = "Mild-vs-Seve-Plasmablasts-MA", width = 5, height = 5, ArchRProj = projDm, addDOC = FALSE)
###########







p <- plotBrowserTrack(
    ArchRProj = projDm, 
    groupBy = "cell.cond", 
    geneSymbol = c("ABO","DDX1","DPP9","IFNAR2","OAS1"),
    features =  getMarkers(markerTest.cmono, cutOff = "Pval <= 0.05 & abs(Log2FC) >= 1", returnGR = TRUE)["cMono.Mild"],
    upstream = 250000,
    downstream = 250000
)
plotPDF(p, name = "Plot-cmono-With-GWASregion", width = 8, height = 6, ArchRProj = projDm, addDOC = FALSE)
## NULL
## [1] 0

p <- plotBrowserTrack(
    ArchRProj = projDm, 
    groupBy = "cell.cond", 
    geneSymbol = c("ABO","DDX1","DPP9","IFNAR2","OAS1"),
    features =  getMarkers(markerTest.ncmono, cutOff = "Pval <= 0.05 & abs(Log2FC) >= 1", returnGR = TRUE)["ncMono.Mild"],
    upstream = 250000,
    downstream = 250000
)
plotPDF(p, name = "Plot-ncmono-With-GWASregion", width = 8, height = 6, ArchRProj = projDm, addDOC = FALSE)



p <- plotBrowserTrack(
    ArchRProj = projDm, 
    groupBy = "pat.gender", 
    geneSymbol = c("USP9Y","DDX3Y","EIF1AY","RPS4Y1"),
    upstream = 50000,
    downstream = 50000
)
plotPDF(p, name = "Plot-Y-link", width = 5, height = 15, ArchRProj = projDm, addDOC = FALSE)


pma <- plotMarkers(seMarker = markerTest2, name = "cMono.Mild", cutOff = "Pval <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
plotPDF(pma, name = "Mild-vs-Seve-Markers-MA", width = 5, height = 5, ArchRProj = projDm, addDOC = FALSE)


library(chromVARmotifs)
projDm <- addMotifAnnotations(ArchRProj = projDm, motifSet = "cisbp", name = "Motif")
projDm <- addArchRAnnotations(ArchRProj = projDm, collection = "EncodeTFBS")
saveArchRProject(ArchRProj= projDm, outputDirectory="projDm4",load=F)



##############################################
############################



markerTest.cmono <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "cMono.Mild",
  bgdGroups = "cMono.Severe"
)

markerList <- getMarkers.raw(markerTest.cmono, cutOff = "Pval <= 0.05 & abs(Log2FC) >= 1", returnGR=F)
head(markerList$cMono.Mild)
#####
markerTest2 <- markerTest.cmono

motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest2,
    ArchRProj = projDm,
    peakAnnotation = "Motif",  #EncodeTFBS
    cutOff = "Pval <= 0.05 & Log2FC >= 1"
  )
  
  


df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

motifsDo <- peakAnnoEnrichment(
    seMarker = markerTest2,
    ArchRProj = projDm,
    peakAnnotation = "Motif", #EncodeTFBS
    cutOff = "Pval <= 0.05 & Log2FC <= -1"
  )
  
df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))
  
plotPDF(ggUp, ggDo, name = "Mild-vs-Seve-cMono-Motifs-Enriched", width = 5, height = 5, ArchRProj = projDm, addDOC = FALSE)

############### TFBS
motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest2,
    ArchRProj = projDm,
    peakAnnotation = "EncodeTFBS",  #EncodeTFBS
    cutOff = "Pval <= 0.05 & Log2FC >= 1"
  )

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) EncodeTFBS Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

motifsDo <- peakAnnoEnrichment(
    seMarker = markerTest2,
    ArchRProj = projDm,
    peakAnnotation = "EncodeTFBS", #EncodeTFBS
    cutOff = "Pval <= 0.05 & Log2FC <= -1"
  )
  
df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) EncodeTFBS Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))
  
plotPDF(ggUp, ggDo, name = "Mild-vs-Seve-cMono-EncodeTFBS-Enriched", width = 5, height = 5, ArchRProj = projDm, addDOC = FALSE)






projDm <- addBgdPeaks(projDm)
projDm <- addDeviationsMatrix(
  ArchRProj = projDm, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(projDm, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = projDm, addDOC = FALSE)
motifs <- c("KLF5", "CTCF", "SP1","CTCFL", "SMARCC1","FOSL1", "FOSL2", "JUN")
markerMotifs <- getFeatures(projDm, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs <- markerMotifs[markerMotifs %ni% "z:MESP1_69"]
markerMotifs <- markerMotifs[markerMotifs %ni% "z:JUND_124"]
markerMotifs <- markerMotifs[markerMotifs %ni% "z:JUNB_139"]
markerMotifs

projDm <- addImputeWeights(projDm)
p <- plotGroups(ArchRProj = projDm, 
  groupBy = "cell.cond", 
  colorBy = "MotifMatrix", 
  name = markerMotifs,
  imputeWeights = getImputeWeights(projDm)
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

plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation", width = 5, height = 5, ArchRProj = projDm, addDOC = FALSE)


markerMotifs <- getFeatures(projDm, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

motifPositions <- getPositions(projDm)
motifs <- c("KLF5", "CTCF", "SP1", "SMARCC1","FOSL1", "FOSL2", "JUN")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% "MESP1_69"]
markerMotifs <- markerMotifs[markerMotifs %ni% "JUND_124"]
markerMotifs <- markerMotifs[markerMotifs %ni% "JUNB_139"]
markerMotifs



###
motifs <- c("PRPF8","EFTUD2","HNRNPC","SF3B4","U2AF1","U2AF2","SRSF1")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
###
projDm <- loadArchRProject("./projDm_main2")
motifPositions <- getPositions(projDm)
all_motif <- ""
for(i in names(motifPositions)){
	aa <- as.data.frame(motifPositions[[i]])
	aa$motif <- i
	all_motif <- rbind(all_motif,aa)
}

saveRDS(all_motif, "all_motif.rds")
write.table(all_motif,"all_motif.txt",sep="\t",quote=F,row.names=F)

####

projDm <- addGroupCoverages(ArchRProj = projDm, groupBy = "cell.cond")
seFoot <- getFootprints(
  ArchRProj = projDm, 
  positions = motifPositions[markerMotifs], 
  groupBy = "cell.cond"
)


plotFootprints(
  seFoot = seFoot,
  ArchRProj = projDm, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)
saveArchRProject(ArchRProj= projDm, outputDirectory="projDm5",load=F)

projDm <- addCoAccessibility(
    ArchRProj = projDm,
    reducedDims = "IterativeLSI"
)

projDm <- addPeak2GeneLinks(
    ArchRProj = projDm,
    reducedDims = "IterativeLSI"
)

p2g <- getPeak2GeneLinks(
    ArchRProj = projDm,
    corCutOff = 0.45,
    resolution = 10000,
    returnLoops = TRUE
)

p2g[[1]]

markerGenes  <- c(
    "ISG15", "IFI27","IFI6", "IFITM3",#Early Progenitor
    "MPO","PLAC8", #Erythroid
    "DUSP6", "FCER1G"
  )
  
  
markerGenes  <- c("CCR1","CCR2","CCR3","DPP9")

p <- plotBrowserTrack(
    ArchRProj = projDm, 
    groupBy = "cell.cond", 
    geneSymbol = markerGenes, 
    upstream = 100000,
    downstream = 100000,
    loops = getPeak2GeneLinks(projDm)
)

plotPDF(plotList = p, 
    name = "Plot-Tracks2-Peak2GeneLinks.pdf", 
    ArchRProj = projDm, 
    addDOC = FALSE, width = 10, height = 7)


###
markerGenes  <- c("LUCAT1")

p <- plotBrowserTrack(
    ArchRProj = projDm2, 
    groupBy = "cell.cond", 
    geneSymbol = markerGenes, 
    upstream = 200000,
    downstream = 200000,
    loops = getPeak2GeneLinks(projDm2)
)


p <- plotBrowserTrack(
    ArchRProj = projDm2, 
    groupBy = "Clusters2", 
    geneSymbol = markerGenes, 
    #plotSummary = c("scTrack", "loopTrack", "geneTrack"),
    sizes = c(12, 1.5, 3, 2),
    upstream = 200000,
    downstream = 200000,
    loops = getPeak2GeneLinks(projDm2)
)

plotPDF(plotList = p, 
    name = "Plot-Tracks_LUCAT-Peak2GeneLinks_t4.pdf", 
    ArchRProj = projDm2, 
    addDOC = FALSE, width = 10, height = 5)





##############################

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks2,
    ArchRProj = projDm,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )




idx <- BiocGenerics::which(projDmFlt$ClustersHm %in% c("C2","C3"))
mono.cells <- projDmFlt$cellNames[idx]
projMono <- subsetArchRProject(ArchRProj= projDmFlt, cells= mono.cells, outputDirectory="projMono")



markersPeaks <- getMarkerFeatures(
    ArchRProj = projDm, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markersPeaks2 <- getMarkerFeatures(
    ArchRProj = projDm, 
    useMatrix = "PeakMatrix", 
    groupBy = "cell.cond",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


markersPeaks3 <- getMarkerFeatures(
    ArchRProj = projDm, 
    useMatrix = "PeakMatrix", 
    groupBy = "cell.disease",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList$cMono

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
markerList


enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks2,
    ArchRProj = projDm,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.05 & Log2FC >= 1"
  )
  
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 20, transpose = TRUE, clusterCols = F, pMax = 100)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap_main", width = 8, height = 6, ArchRProj = projDm, addDOC = FALSE)


enrichEncode <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = projDm,
    peakAnnotation = "EncodeTFBS",
    cutOff = "FDR <= 0.1 & Log2FC >= 1"
  )
heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEncode, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEncode, name = "EncodeTFBS-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = projDm, addDOC = FALSE)



saveArchRProject(ArchRProj= projDm, outputDirectory="projDm4",load=F)



projDm <- addCoAccessibility(
    ArchRProj = projDm,
    reducedDims = "IterativeLSI"
)

projDm <- addPeak2GeneLinks(
    ArchRProj = projDm,
    reducedDims = "IterativeLSI"
)



idx <- BiocGenerics::which(projDmFlt$ClustersHm %in% c("C9"))
nk.cells <- projDmFlt$cellNames[idx]
projNK <- subsetArchRProject(ArchRProj= projDmFlt, cells= nk.cells, outputDirectory="projNK")

projMono <- addGroupCoverages(ArchRProj=projMono, groupBy="patient")

pathToMacs2 <- "/vol/projects/CIIM/resources/tools/encode-atac-seq-pipeline/bin/macs2"
projMono <- addReproduciblePeakSet(
	ArchRProj = projMono,
	groupBy="patient",
	pathToMacs2= pathToMacs2
	)
getPeakSet(projMono)
projMono <- addPeakMatrix(projMono)
saveArchRProject(ArchRProj= projMono, outputDirectory="projMono2",load=F)


markerGenes <- c("CLEC12A","SNHG5","DUSP6","IFI27","FCER1G","AIF1","IFITM1","MPO")
p <- plotBrowserTrack(
	ArchRProj= projMono,
	groupBy="Severity",
	geneSymbol=markerGenes,
	upstream=50000,
	downstream=50000
)
plotPDF(plotList=p, name="plot_track_monocytes_sample.pdf",ArchRProj= projDmFlt, addDOC=F, width=5,height=4)

markerGenes <- c("IGKC","IGHA1","IGHM","IGLC2","S100A8","FAU","GNLY","HLA-C","LAIR2")
p <- plotBrowserTrack(
	ArchRProj= projNK,
	groupBy="Severity",
	geneSymbol=markerGenes,
	upstream=50000,
	downstream=50000
)
plotPDF(plotList=p, name="plot_track_NKs_sample.pdf",ArchRProj= projDmFlt, addDOC=F, width=5,height=4)














getMarkers.raw <- function (seMarker = NULL, cutOff = "Pval <= 0.05 & Log2FC >= 0.5", 
    n = NULL, returnGR = FALSE) 
{

    assayNames <- names(SummarizedExperiment::assays(seMarker))
    for (an in assayNames) {
        eval(parse(text = paste0(an, " <- ", "SummarizedExperiment::assays(seMarker)[['", 
            an, "']]")))
    }
    passMat <- eval(parse(text = cutOff))
    for (an in assayNames) {
        eval(parse(text = paste0("rm(", an, ")")))
    }
    if (returnGR) {
        if (metadata(seMarker)$Params$useMatrix != "PeakMatrix") {
            stop("Only markers can be returned as GRanges when PeakMatrix!")
        }
        rr <- GRanges(rowData(seMarker)$seqnames, IRanges(rowData(seMarker)$start, 
            rowData(seMarker)$end))
        grL <- lapply(seq_len(ncol(passMat)), function(x) {
            idx <- which(passMat[, x])
            rrx <- rr[idx]
            rrx$Log2FC <- SummarizedExperiment::assays(seMarker[idx, 
                ])[["Log2FC"]][, x]
            rrx$Pval <- SummarizedExperiment::assays(seMarker[idx, 
                ])[["Pval"]][, x]
            if ("MeanDiff" %in% assayNames) {
                rrx$MeanDiff <- SummarizedExperiment::assays(seMarker[idx, 
                  ])[["MeanDiff"]][, x]
            }
            rrx <- rrx[order(rrx$Pval), , drop = FALSE]
            if (!is.null(n)) {
                if (n < nrow(rrx)) {
                  rrx <- rrx[seq_len(n), , drop = FALSE]
                }
            }
            rrx
        }) %>% SimpleList
        names(grL) <- colnames(seMarker)
        grL <- grL[gtools::mixedsort(names(grL))]
        return(grL)
    }
    else {
        markerList <- lapply(seq_len(ncol(passMat)), function(x) {
            idx <- which(passMat[, x])
            rrx <- SummarizedExperiment::rowData(seMarker[idx, 
                ])
            rrx$Log2FC <- SummarizedExperiment::assays(seMarker[idx, 
                ])[["Log2FC"]][, x]
            rrx$Pval <- SummarizedExperiment::assays(seMarker[idx, 
                ])[["Pval"]][, x]
            if ("MeanDiff" %in% assayNames) {
                rrx$MeanDiff <- SummarizedExperiment::assays(seMarker[idx, 
                  ])[["MeanDiff"]][, x]
            }
            rrx <- rrx[order(rrx$Pval), , drop = FALSE]
            if (!is.null(n)) {
                if (n < nrow(rrx)) {
                  rrx <- rrx[seq_len(n), , drop = FALSE]
                }
            }
            rrx
        }) %>% SimpleList
        names(markerList) <- colnames(seMarker)
        return(markerList)
    }
}







##########


#projMono <- addGroupCoverages(ArchRProj=projMono, groupBy="patient")#

#pathToMacs2 <- "/vol/projects/CIIM/resources/tools/encode-atac-seq-pipeline/bin/macs2"
#projMono <- addReproduciblePeakSet(
#	ArchRProj = projMono,
#	groupBy="patient",
#	pathToMacs2= pathToMacs2
#	)
#getPeakSet(projMono)
#projMono <- addPeakMatrix(projMono)
#saveArchRProject(ArchRProj= projMono, outputDirectory="projMono2",load=F)



markerTest1 <- getMarkerFeatures(
  ArchRProj = projMono3, 
  useMatrix = "PeakMatrix",
  groupBy = "Severity",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Severe",
  bgdGroups = "Post"
)
pv1 <- plotMarkers(seMarker = markerTest1, name = "Severe", cutOff = "Pval <= 0.05", plotAs = "Volcano")


markerTest2 <- getMarkerFeatures(
  ArchRProj = projMono3, 
  useMatrix = "PeakMatrix",
  groupBy = "Severity",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Severe",
  bgdGroups = "Mild"
)
pv2 <- plotMarkers(seMarker = markerTest2, name = "Severe", cutOff = "Pval <= 0.05", plotAs = "Volcano")

markerTest3 <- getMarkerFeatures(
  ArchRProj = projMono, 
  useMatrix = "PeakMatrix",
  groupBy = "Severity",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Mild",
  bgdGroups = "Post"
)
pv3 <- plotMarkers(seMarker = markerTest3, name = "Mild", cutOff = "Pval <= 0.05", plotAs = "Volcano")

markerList <- getMarkers(markerTest1, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", returnGR=T)
markerList$Severe


plotPDF(pv1,pv2,pv3, name="plot_DP.pdf",ArchRProj= projDmFlt, addDOC=F, width=5,height=4)





######################
#> table(projDm$Clusters2)
#
#          B       CD4 T       CD8 T       cMono      ncMono          NK 
#       3728       11145        7224       10324        1224        5838 
#        pDC Plasmablast    undefine 
#        189         350        6668 


idx <- BiocGenerics::which(projDm$Clusters2 %in% c("cMono") & projDm$Severity %in% "Severe")
mono.cells <- projDm$cellNames[idx]
projSub <- subsetArchRProject(ArchRProj= projDm, cells= mono.cells, outputDirectory="proj_cMonoSev")

###
idx <- BiocGenerics::which(projDm$Clusters2 %in% c("cMono") & projDm$Severity %in% "Mild")
mono.cells <- projDm$cellNames[idx]
projSub <- subsetArchRProject(ArchRProj= projDm, cells= mono.cells, outputDirectory="proj_cMonoMild")

seGroupMotif <- getGroupSE(ArchRProj = projSub, useMatrix = "MotifMatrix", groupBy = "patient")
seGroupMotif

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

saveRDS(seZ,"seZ.cmono.Mild.rds")
write.table(rowData(seZ),"seZ.cmono.mild.txt",sep="\t",quote=F,row.names=F)

###
idx <- BiocGenerics::which(projDm$Clusters2 %in% c("cMono") & projDm$Severity %in% "Post")
mono.cells <- projDm$cellNames[idx]
projSub <- subsetArchRProject(ArchRProj= projDm, cells= mono.cells, outputDirectory="proj_cMonoPost")

seGroupMotif <- getGroupSE(ArchRProj = projSub, useMatrix = "MotifMatrix", groupBy = "patient")
seGroupMotif

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

saveRDS(seZ,"seZ.cmono.Post.rds")
write.table(rowData(seZ),"seZ.cmono.post.txt",sep="\t",quote=F,row.names=F)

########################

seZ.cmono.Severe <- readRDS("seZ.mono.rds")
seZ.cmono.Post <- readRDS("seZ.cmono.Post.rds")
seZ.cmono.Mild <- readRDS("seZ.cmono.Mild.rds")

seZ.1 <- read.table("seZ.cmono.Sev.txt",header = T,stringsAsFactors = F)
seZ.2 <- read.table("seZ.cmono.Post.txt",header = T,stringsAsFactors = F)
seZ.3 <- read.table("seZ.cmono.Mild.txt",header = T,stringsAsFactors = F)


seZ.1$cond <- "Severe"
seZ.3$cond <- "Mild"
seZ.2$cond <- "Post"

seZ.1$dev <- "low"
seZ.1$dev[seZ.1$maxDelta > quantile(seZ.1$maxDelta, 0.75)]= "high"
seZ.2$dev <- "low"
seZ.2$dev[seZ.2$maxDelta > quantile(seZ.2$maxDelta, 0.75)]="high"
seZ.3$dev <- "low"
seZ.3$dev[seZ.3$maxDelta > quantile(seZ.3$maxDelta, 0.75)]="high"


hi.dev.severe <- seZ.1$name[seZ.1$maxDelta > quantile(seZ.1$maxDelta, 0.75)]
hi.dev.post <- seZ.2$name[seZ.2$maxDelta > quantile(seZ.2$maxDelta, 0.75)]
hi.dev.mild <- seZ.3$name[seZ.3$maxDelta > quantile(seZ.3$maxDelta, 0.75)]

write.table(hi.dev.post,"hi.dev.post.txt",row.names=F, col.names=F, quote=F)
write.table(hi.dev.mild,"hi.dev.mild.txt",row.names=F, col.names=F, quote=F)
write.table(hi.dev.severe,"hi.dev.sev.txt",row.names=F, col.names=F, quote=F)

seZ <- rbind(seZ.3,seZ.1,seZ.2)
seZ$cond <- factor(seZ$cond, levels=c("Severe","Mild","Post"))
#p <- ggplot(data.frame(seZ), aes(x=cond, y=maxDelta)) +
#  geom_violin()
  
  ggplot(data.frame(seZ), aes(x=cond, y=maxDelta)) +
  geom_violin(scale="width",adjust=1,aes(fill=factor(cond)),show.legend = F,draw_quantiles =  c(0.25,0.5,0.75)) + 
  ylab("Motif Deviation") + xlab("") + theme_bw()+scale_fill_manual(values=c("#9e0142","#fee08b","#4393c3"))

######
######
######


###

idx <- BiocGenerics::which(projDm$Clusters2 %in% c("ncMono") & projDm$Severity %in% "Severe")
mono.cells <- projDm$cellNames[idx]
projSub <- subsetArchRProject(ArchRProj= projDm, cells= mono.cells, outputDirectory="proj_ncMonoSevere")

seGroupMotif <- getGroupSE(ArchRProj = projSub, useMatrix = "MotifMatrix", groupBy = "patient")
seGroupMotif

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

saveRDS(seZ,"seZ.ncmono.Sev.rds")
write.table(rowData(seZ),"seZ.ncmono.sev.txt",sep="\t",quote=F,row.names=F)

idx <- BiocGenerics::which(projDm$Clusters2 %in% c("ncMono") & projDm$Severity %in% "Mild")
mono.cells <- projDm$cellNames[idx]
projSub <- subsetArchRProject(ArchRProj= projDm, cells= mono.cells, outputDirectory="proj_ncMonoMild")

seGroupMotif <- getGroupSE(ArchRProj = projSub, useMatrix = "MotifMatrix", groupBy = "patient")
seGroupMotif

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

saveRDS(seZ,"seZ.ncmono.Mild.rds")
write.table(rowData(seZ),"seZ.ncmono.mild.txt",sep="\t",quote=F,row.names=F)

###
idx <- BiocGenerics::which(projDm$Clusters2 %in% c("ncMono") & projDm$Severity %in% "Post")
mono.cells <- projDm$cellNames[idx]
projSub <- subsetArchRProject(ArchRProj= projDm, cells= mono.cells, outputDirectory="proj_ncMonoPost")

seGroupMotif <- getGroupSE(ArchRProj = projSub, useMatrix = "MotifMatrix", groupBy = "patient")
seGroupMotif

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

saveRDS(seZ,"seZ.ncmono.Post.rds")
write.table(rowData(seZ),"seZ.ncmono.post.txt",sep="\t",quote=F,row.names=F)
########################

seZ.1 <- read.table("seZ.ncmono.Sev.txt",header = T,stringsAsFactors = F)
seZ.2 <- read.table("seZ.ncmono.Post.txt",header = T,stringsAsFactors = F)
seZ.3 <- read.table("seZ.ncmono.Mild.txt",header = T,stringsAsFactors = F)


seZ.1$cond <- "Severe"
seZ.3$cond <- "Mild"
seZ.2$cond <- "Post"

seZ.1$dev <- "low"
seZ.1$dev[seZ.1$maxDelta > quantile(seZ.1$maxDelta, 0.75)]= "high"
seZ.2$dev <- "low"
seZ.2$dev[seZ.2$maxDelta > quantile(seZ.2$maxDelta, 0.75)]="high"
seZ.3$dev <- "low"
seZ.3$dev[seZ.3$maxDelta > quantile(seZ.3$maxDelta, 0.75)]="high"


nc.dev.severe <- seZ.1$name[seZ.1$maxDelta > quantile(seZ.1$maxDelta, 0.75)]
nc.dev.post <- seZ.2$name[seZ.2$maxDelta > quantile(seZ.2$maxDelta, 0.75)]
nc.dev.mild <- seZ.3$name[seZ.3$maxDelta > quantile(seZ.3$maxDelta, 0.75)]

write.table(nc.dev.post,"nc.dev.post.txt",row.names=F, col.names=F, quote=F)
write.table(nc.dev.mild,"nc.dev.mild.txt",row.names=F, col.names=F, quote=F)
write.table(nc.dev.severe,"nc.dev.sev.txt",row.names=F, col.names=F, quote=F)

seZ <- rbind(seZ.3,seZ.1,seZ.2)
seZ$cond <- factor(seZ$cond, levels=c("Severe","Mild","Post"))
#p <- ggplot(data.frame(seZ), aes(x=cond, y=maxDelta)) +
#  geom_violin()
  
  ggplot(data.frame(seZ), aes(x=cond, y=maxDelta)) +
  geom_violin(scale="width",adjust=1,aes(fill=factor(cond)),show.legend = F,draw_quantiles =  c(0.25,0.5,0.75)) + 
  ylab("Motif Deviation") + xlab("") + theme_bw()+scale_fill_manual(values=c("#9e0142","#fee08b","#4393c3"))

c.mvp <- unique(gsub("_[0-9]+","",setdiff(hi.dev.mild,hi.dev.post)))
c.svp <- unique(gsub("_[0-9]+","",setdiff(hi.dev.severe,hi.dev.post)))
nc.mvp <- unique(gsub("_[0-9]+","",setdiff(nc.dev.mild,nc.dev.post)))
nc.svp <- unique(gsub("_[0-9]+","",setdiff(nc.dev.severe,nc.dev.post)))
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
de1 <- bitr(c.mvp, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
de2 <- bitr(c.svp, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
de3 <- bitr(nc.mvp, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
de4 <- bitr(nc.svp, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")


nt1_kk <- enrichKEGG(gene = de1$ENTREZID,
                     organism     = "hsa",
                     pvalueCutoff = 0.1)
nt2_kk <- enrichKEGG(gene = de2$ENTREZID,
                     organism     = "hsa",
                     pvalueCutoff = 0.1)
nt3_kk <- enrichKEGG(gene = de3$ENTREZID,
                     organism     = "hsa",
                     pvalueCutoff = 0.1)                     
nt4_kk <- enrichKEGG(gene = de4$ENTREZID,
                     organism     = "hsa",
                     pvalueCutoff = 0.1)
                     
nt1_kk<- setReadable(nt1_kk, "org.Hs.eg.db", "ENTREZID")
nt2_kk<- setReadable(nt2_kk, "org.Hs.eg.db", "ENTREZID")
nt3_kk<- setReadable(nt3_kk, "org.Hs.eg.db", "ENTREZID")
nt4_kk<- setReadable(nt4_kk, "org.Hs.eg.db", "ENTREZID")

k = merge_result(list('cMono.Mild_vs_Post' = nt1_kk,'cMono.Severe_vs_Post' = nt2_kk ,'ncMono.Mild_vs_Post' = nt3_kk , 'ncMono.Severe_vs_Post' = nt4_kk)) %>%
clusterProfiler::dotplot(.,showCategory=10) + theme(axis.text.x=element_text(angle=60, hjust=1))

##############
corGSM_MM <- correlateMatrices(
    ArchRProj = projDm,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "IterativeLSI"
)
corGIM_MM <- correlateMatrices(
    ArchRProj = projDm,
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





###############
idx <- BiocGenerics::which(projDm$Clusters2 %in% c("B", "CD4 T","CD8 T"  , "cMono" ,  "ncMono" ,"NK"))
mono.cells <- projDm$cellNames[idx]
projMain <- subsetArchRProject(ArchRProj= projDm, cells= mono.cells, outputDirectory="projDm_main")



projMain$Clusters3 <- projMain$Clusters2
projMain$Clusters3[projMain$Clusters3 == "cMono"]="Mono"
projMain$Clusters3[projMain$Clusters3 == "ncMono"]="Mono"
projMain$Clusters3[projMain$Clusters3 == "CD4 T"]="T"
projMain$Clusters3[projMain$Clusters3 == "CD8 T"]="T"
projMain$Clusters3.cond <- paste(projMain$Clusters3, projMain$Severity, sep=".")
projMain <- addGroupCoverages(ArchRProj = projMain, groupBy = "Clusters3.cond",force=T)


motifPositions <- getPositions(projMain)
#motifs <- c("SPIB_336", "SPIC_344","POU3F4_619", "BCL11A_194", "RUNX1_733", "STAT2_778",
# "FOS_137", "JUNB_139","JUN_143","FOSL2_105","FOSL1_142","JUND_124","CEBPA_155","CEBPB_140")
#motifs <- c("POU2F2","SPI1","SP3","ZBTB7A","FOSL2","JUNB","JUND","JDP2","BACH1")
motifs <- c("BCL11A","BCL11B","CEBPA_155","CEBPB_140", "CEBPD")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
#markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs

seFoot <- getFootprints(
  ArchRProj = projMain, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Clusters3.cond"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = projMain, 
  normMethod = "Subtract",
  plotName = "main4Footprints-Subtract-Bias3",
  addDOC = FALSE,
  smoothWindow = 5
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = projMain, 
  normMethod = "Divide",
  plotName = "main4Footprints-Divide-Bias3",
  addDOC = FALSE,
  smoothWindow = 5
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = projMain, 
  normMethod = "None",
  plotName = "main4Footprints-No-Normalization3",
  addDOC = FALSE,
  smoothWindow = 5
)

saveArchRProject(ArchRProj= projMain, outputDirectory="projDm_main",load=F)


#######

projMain <- addGroupCoverages(ArchRProj = projMain, groupBy = "cell.disease")

markersPeaks3 <- getMarkerFeatures(
    ArchRProj = projMain, 
    useMatrix = "PeakMatrix", 
    groupBy = "cell.disease",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)



enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks3,
    ArchRProj = projMain,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.05 & Log2FC >= 1"
  )
  
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 20, transpose = TRUE, clusterCols = F, pMax = 100)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap_main", width = 12, height = 6, ArchRProj = projMain, addDOC = FALSE)




#######
seGroupMotif <- getGroupSE(ArchRProj = projSub, useMatrix = "PeakMatrix", groupBy = "patient")
seGroupMotif

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

saveRDS(seZ,"seZ.ncmono.Post.rds")
write.table(rowData(seZ),"seZ.ncmono.post.txt",sep="\t",quote=F,row.names=F)




#########

uniq <- read.table("/vol/projects/CIIM/Covid_50MHH/scATAC/analysis_new/stats/uniq_atac.txt")[,1]

idx <- BiocGenerics::which(projDm$patient %in% uniq)
mono.cells <- projDm$cellNames[idx]
projMain <- subsetArchRProject(ArchRProj= projDm, cells= mono.cells, outputDirectory="projDm_uniq")




peak_m = getMatrixFromProject(
  ArchRProj = projMain,
  useMatrix = "PeakMatrix",
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)
ids = unique(projMain$cell.cond)
int_m = peak_m@assays@data$PeakMatrix
rownames(int_m) = paste(as.vector(projMain@peakSet@seqnames), as.vector(projMain@peakSet@ranges@start), sep = "_")
#colnames(int_m) = rownames(peak_m@colData)
colnames(int_m) = peak_m@colData$cell.cond
bulk_m = matrix(0, dim(int_m)[1], length(ids))
for(i in 1:length(ids))
{
  bulk_m[,i] = apply(int_m[,which(colnames(int_m) == ids[i])],1,sum)
}
rownames(bulk_m) = rownames(int_m)
colnames(bulk_m) = ids



