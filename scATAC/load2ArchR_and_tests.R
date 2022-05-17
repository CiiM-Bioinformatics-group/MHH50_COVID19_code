

#### use ArchR

#devtools::install_github("GreenleafLab/ArchR",ref="master",repos=BiocManager::repositories())
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


###load  (if saved before)
projMono3 <- loadArchRProject("./projMono3")


##### set up ####
addArchRGenome("hg38")
addArchRThreads(threads=4)


########Creating Arrow Files #######
inputfiles <- c(
"pool1a"="/vol/projects/CIIM/Covid_50MHH/scATAC/align_out/atac_1a/outs/fragments.tsv.gz",
"pool1b"="/vol/projects/CIIM/Covid_50MHH/scATAC/align_out/atac_1b/outs/fragments.tsv.gz",
"pool1d"="/vol/projects/CIIM/Covid_50MHH/scATAC/align_out/atac_1d/outs/fragments.tsv.gz",
"pool2a"="/vol/projects/CIIM/Covid_50MHH/scATAC/align_out/atac_2a/outs/fragments.tsv.gz",
"pool2d"="/vol/projects/CIIM/Covid_50MHH/scATAC/align_out/atac_2d/outs/fragments.tsv.gz",
"pool3a"="/vol/projects/CIIM/Covid_50MHH/scATAC/align_out/atac_3a/outs/fragments.tsv.gz",
"pool3b"="/vol/projects/CIIM/Covid_50MHH/scATAC/align_out/atac_3b/outs/fragments.tsv.gz",
"pool3d"="/vol/projects/CIIM/Covid_50MHH/scATAC/align_out/atac_3d/outs/fragments.tsv.gz",
"pool4a"="/vol/projects/CIIM/Covid_50MHH/scATAC/align_out/atac_4a/outs/fragments.tsv.gz",
"pool4c"="/vol/projects/CIIM/Covid_50MHH/scATAC/align_out/atac_4c/outs/fragments.tsv.gz"
)




ArrowFiles <- createArrowFiles(
	inputFiles = inputfiles,
	sampleNames = names(inputfiles),
	filterTSS = 4,
	filterFrags=1000,
	addTileMat=T,
	addGeneScoreMat=T
)

#ArrowFiles
########Done create Arrow Files #######

## Doublet
doubScores <- addDoubletScores(
	input=ArrowFiles,
	k=10,
	knnMethod="UMAP",
	LSIMethod=1
)

###create ArchRProject
projCovid <- ArchRProject(
	ArrowFiles=ArrowFiles,
	outputDirectory="postCovid",
	copyArrows=T
)

p1 <- plotFragmentSizes(ArchRProj = projCovid)
p2 <- plotTSSEnrichment(ArchRProj = projCovid)
plotPDF(p1,p2,name="QC_sample_fragSizes_TSS.pdf", ArchRProj=projCovid, addDoc=F, width=5, height=5)
saveArchRProject(ArchRProj=projCovid, outputDirectory="Save_projCovid",load=F)


projCovid2 <- filterDoublets(projCovid)
projCovid2 <- addIterativeLSI(
	ArchRProj=projCovid2,
	useMatrix="TileMatrix",
	name="IterativeLSI",
	iterations=2,
	clusterParams=list(
		resolution=c(0.3),
		sampleCells=10000,
		n.start=10
	),
	varFeatures=25000,
	dimsToUse=1:30
)

projCovid2 <- addHarmony(
	ArchRProj=projCovid2,
	reducedDims="IterativeLSI",
	name="Harmony",
	groupBy="Sample"
)

projCovid2 <- addClusters(
	input=projCovid2,
	reducedDims="IterativeLSI",
	method="Seurat",
	name="Clusters",
	resolution=0.8
)

projCovid2 <- addUMAP(
	ArchRProj=projCovid2,
	reducedDims="IterativeLSI",
	name="UMAP",
	nNeighbors=30,
	minDist=0.5,
	metric="cosine"
)

projCovid2 <- addUMAP(
	ArchRProj=projCovid2,
	reducedDims="Harmony",
	name="UMAPHarmony",
	nNeighbors=30,
	minDist=0.5,
	metric="cosine"
)

p1 <- plotEmbedding(ArchRProj = projCovid2, colorBy = "cellColData", name = "Sample", embedding="UMAP")
p2 <- plotEmbedding(ArchRProj = projCovid2, colorBy = "cellColData", name = "Clusters", embedding="UMAP")
p3 <- plotEmbedding(ArchRProj = projCovid2, colorBy = "cellColData", name = "Sample", embedding="UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = projCovid2, colorBy = "cellColData", name = "Clusters", embedding="UMAPHarmony")

#plotPDF(p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = projCovid2,addDoc=F, width=4, height=4)

plotPDF(p1,p2,p3,p4, name = "Plots-UMAP2Harmony-Sample-Clusters.pdf", ArchRProj = projCovid2,addDoc=F, width=5, height=5)

saveArchRProject(ArchRProj=projCovid2, outputDirectory="Save_projCovid2",load=F)

#devtools::install_github("immunogenomics/presto")

markersGS <- getMarkerFeatures(
	ArchRProj=projCovid2,
	useMatrix="GeneScoreMatrix",
	groupBy="Clusters",
	bias=c("TSSEnrichment","log10(nFrags)"),
	testMethod="wilcoxon"

)

markerList <- getMarkers(markersGS,cutoff="FDR <= 0.01 & Log2FC >= 1.25")




marker.gene <- c(
"IL7R","CCR7","TCF7",  #IL7R,CCR7-,TCF7-             Memory CD4 T   
"CD8A","CD8B","GZMK", #CD8A+,CD8B+,GZMK+             CD8 T    
"NKG7","NCAM1", #NKG7+,NCAM1+,CD8A-             NK    
"CD14","LYZ","FCGR3A", #CD14+, LYZ+, FCGR3A-            CD14+ monocytes   
#FCGR3A+, CD14-              CD16+ monocytes   
"CST3","CD86","HLA-DRA",  #CST3+, CD86+, HLA-DRA+, CD14-, FCGR3A-   mDC, pDC (dendritic cells)    
"IL1B","TNF","IFI6","IFITM3",  #IL1B,TNF,IL6               proinflammatory markers  # macrophages
"IFNG","GZMB","CCL5","CST7","CLEC9A","CD1D","CXCL9","CXCL10","CXCL11",
"CD79A", #CD79A+, CD27-, SDC1-           B cell        
#CD79A+, CD27+, SDC1+           plasmablast        
"PPBP", #PPBP+                                 megakaryocyte/platelet   
"KIT","TPSAB1" #KIT+, TPSAB1+                 Mast (edited)
)

heatmapGS <- markerHeatmap(
	seMarker = markersGS,
	cutOff = "FDR<0.01 & Log2FC >= 1",
	labelMarkers=marker.gene,
	transpose=T
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side="bot", annotation_legend_side="bot")


plotPDF(heatmapGS, name = "Plot-geneScores.pdf", ArchRProj = projCovid2,addDoc=F, width=8, height=6)



meta <- read.table("/vol/projects/CIIM/Covid_50MHH/scATAC/souporcell/extract.impute.txt",header=T,sep="\t",stringsAsFactors=F,comment.char = "&")
tmp <- projCovid2$cellNames

meta2 <- meta[meta$cellname %in% tmp,]
projDm <- subsetArchRProject(ArchRProj=projCovid2, cells= meta2$cellname, outputDirectory="projDm")



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

projDm <- addUMAP(
	ArchRProj= projDm,
	reducedDims="Harmony",
	name="UMAPHarmony",
	nNeighbors=30,
	minDist=0.5,
	metric="cosine",
	force=TRUE
)

p1 <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", name = "Sample", embedding="UMAP")
p2 <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", name = "Clusters", embedding="UMAP")
p3 <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", name = "Sample", embedding="UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", name = "Clusters", embedding="UMAPHarmony")

#plotPDF(p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = projCovid2,addDoc=F, width=4, height=4)

plotPDF(p1,p2,p3,p4, name = "Plots-UMAP2Harmony-Sample-Clusters.pdf", ArchRProj = projDm,addDoc=F, width=5, height=5)



markersGS <- getMarkerFeatures(
	ArchRProj=projDm,
	useMatrix="GeneScoreMatrix",
	groupBy="Clusters",
	bias=c("TSSEnrichment","log10(nFrags)"),
	testMethod="wilcoxon"

)

markerList <- getMarkers(markersGS,cutoff="FDR <= 0.01 & Log2FC >= 1.25")




marker.gene <- c(
"IL7R","CCR7","TCF7",  #IL7R,CCR7-,TCF7-             Memory CD4 T   
"CD8A","CD8B","GZMK", #CD8A+,CD8B+,GZMK+             CD8 T    
"NKG7","NCAM1", #NKG7+,NCAM1+,CD8A-             NK    
"CD14","LYZ","FCGR3A", #CD14+, LYZ+, FCGR3A-            CD14+ monocytes   
#FCGR3A+, CD14-              CD16+ monocytes   
"CST3","CD86","HLA-DRA",  #CST3+, CD86+, HLA-DRA+, CD14-, FCGR3A-   mDC, pDC (dendritic cells)    
"IL1B","TNF","IFI6","IFITM3",  #IL1B,TNF,IL6               proinflammatory markers  # macrophages
"IFNG","GZMB","CCL5","CST7","CLEC9A","CD1D","CXCL9","CXCL10","CXCL11",
"CD79A", #CD79A+, CD27-, SDC1-           B cell        
#CD79A+, CD27+, SDC1+           plasmablast        
"PPBP", #PPBP+                                 megakaryocyte/platelet   
"KIT","TPSAB1" #KIT+, TPSAB1+                 Mast (edited)
)

heatmapGS <- markerHeatmap(
	seMarker = markersGS,
	cutOff = "FDR<0.01 & Log2FC >= 1",
	labelMarkers=marker.gene,
	transpose=T
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side="bot", annotation_legend_side="bot")


plotPDF(heatmapGS, name = "Plot-geneScores.pdf", ArchRProj = projDm,addDoc=F, width=8, height=6)



tmp <- meta2$id
names(tmp)=meta2$cellname
#                   cellname   id patient gender age Severity2
#1 pool3b#AAACGAACAACGGGTA-1 3b_0      31      M  46      Mild
#2 pool3b#AAACGAAGTAACGGCA-1 3b_0      31      M  46      Mild
#3 pool3b#AAACGAAGTTGAAGCG-1 3b_0      31      M  46      Mild

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

projDm@cellColData


saveArchRProject(ArchRProj= projDm, outputDirectory="projDm2",load=F)



idx <- BiocGenerics::which(projDm$Clusters %in% c("C11","C12"))
mono.cells <- projDm$cellNames[idx]
projMono <- subsetArchRProject(ArchRProj=projDm, cells= mono.cells, outputDirectory="projMono")


p1 <- plotEmbedding(ArchRProj = projDm, colorBy = "cellColData", name = "Severity", embedding="UMAPHarmony")
p2 <- plotEmbedding(ArchRProj = projMono, colorBy = "cellColData", name = "Severity", embedding="UMAPHarmony")

plotPDF(p1,p2, name = "Plot-UMAP-condition-Clusters.pdf", ArchRProj = projDm,addDoc=F, width=4, height=4)


idx <- BiocGenerics::which(projCovid2$Clusters %in% c("C8","C9"))
NK.cells <- projCovid2$cellNames[idx]
projNK <- subsetArchRProject(ArchRProj=projCovid2, cells= NK.cells, outputDirectory="projNK")



markerGenes <- c("IL1B","TNF","CXCL10","HLA-DRA","IL6","PTSG2","GBP1","IFI6","IFITM3")
p <- plotBrowserTrack(
	ArchRProj= projMono,
	groupBy="Severity",
	geneSymbol=markerGenes,
	upstream=50000,
	downstream=50000
)
plotPDF(plotList=p, name="plot_track_monocytes_sample.pdf",ArchRProj= projDm, addDOC=F, width=5,height=4)


markerGenes <- c("IFNG","TGFB1","CCL5")
p <- plotBrowserTrack(
	ArchRProj= projNK,
	groupBy="Sample",
	geneSymbol=markerGenes,
	upstream=50000,
	downstream=50000
)
plotPDF(plotList=p, name="plot_track_NK_sample.pdf",ArchRProj=projCovid2, addDOC=F, width=5,height=4)


#projMono <- addGroupCoverages(ArchRProj=projMono, groupBy="Clusters")
#projMono <- addGroupCoverages(ArchRProj=projMono, groupBy="patient")



projMono <- addGroupCoverages(ArchRProj=projMono, groupBy="Severity")

pathToMacs2 <- "/vol/projects/CIIM/resources/tools/encode-atac-seq-pipeline/bin/macs2"
projMono <- addReproduciblePeakSet(
	ArchRProj = projMono,
	groupBy="Severity",
	pathToMacs2= pathToMacs2
	)
getPeakSet(projMono)
projMono <- addPeakMatrix(projMono)
saveArchRProject(ArchRProj= projMono, outputDirectory="projMono2",load=F)



markerTest1 <- getMarkerFeatures(
  ArchRProj = projMono3, 
  useMatrix = "PeakMatrix",
  groupBy = "Severity",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Severe",
  bgdGroups = "Post"
)
pv1 <- plotMarkers(seMarker = markerTest1, name = "Severe", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")


markerTest2 <- getMarkerFeatures(
  ArchRProj = projMono3, 
  useMatrix = "PeakMatrix",
  groupBy = "Severity",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Severe",
  bgdGroups = "Mild"
)
pv2 <- plotMarkers(seMarker = markerTest2, name = "Severe", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")

markerTest3 <- getMarkerFeatures(
  ArchRProj = projMono, 
  useMatrix = "PeakMatrix",
  groupBy = "Severity",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Mild",
  bgdGroups = "Post"
)
pv3 <- plotMarkers(seMarker = markerTest3, name = "Mild", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")

markerList <- getMarkers(markerTest1, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", returnGR=T)
markerList$Severe


plotPDF(pv1,pv2,pv3, name="plot_DP.pdf",ArchRProj= projDm, addDOC=F, width=5,height=4)


markerGenes <- c("IL1B","TNF","CXCL10","HLA-DRA","IL6","PTSG2","GBP1","IFI6","IFITM3")

gene1 <- c("WNT2B","CRISPLD2","LPP","FAM19A3","FTH1","SH3BP2","FAM85B","TMEM53")
p <- plotBrowserTrack(
	ArchRProj= projMono3,
	groupBy="Severity",
	geneSymbol= gene1,
	upstream=100000,
	downstream=100000,
	features =  markerList$Severe
)
plotPDF(plotList=p, name="plot_tracks.pdf",ArchRProj= projMono3, addDOC=F, width=7,height=4)


gene1 <- c("LPP")
p <- plotBrowserTrack(
	ArchRProj= projMono3,
	groupBy="Severity",
	geneSymbol= gene1,
	upstream=400000,
	downstream=400000,
	features =  markerList$Severe
)
plotPDF(plotList=p, name="plot_track3.pdf",ArchRProj= projMono3, addDOC=F, width=7,height=4)






####
markerGenes <- c("LUCAT1","LINC00211","S100A8","S100A9","CHI3L1")
p <- plotBrowserTrack(
	ArchRProj= projDm,
	groupBy="cell.cond",
	geneSymbol=markerGenes,
	upstream=100000,
	downstream=100000,
	loops = getPeak2GeneLinks(projDm)
)
plotPDF(plotList=p, name="plot_track_pbmc_lincRNA_cond.pdf",ArchRProj= projMono, addDOC=F, width=10,height=8)





