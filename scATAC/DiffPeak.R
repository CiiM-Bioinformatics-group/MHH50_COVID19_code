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


#projDm <- loadArchRProject("./projDm_main2")
projDm <- loadArchRProject("./projDm_uniq")
projDm <- addGroupCoverages(ArchRProj=projDm, groupBy="cell.cond",force=TRUE)


mk.cmono.mp <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "cMono.mild",
  bgdGroups = "cMono.post"
)

mk.cmono.sp <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "cMono.severe",
  bgdGroups = "cMono.post"
)

markerList <- getMarkers.raw(mk.cmono.mp, cutOff = "Pval <= 0.05", returnGR=F)
markerList$cMono.mild
write.table(markerList$cMono.mild, file="DiffPeak.cMono.MvP.txt",sep="\t",quote=F)

markerList <- getMarkers.raw(mk.cmono.sp, cutOff = "Pval <= 0.05", returnGR=F)
markerList$cMono.severe
write.table(markerList$cMono.severe, file="DiffPeak.cMono.SvP.txt",sep="\t",quote=F)

mk.cmono.ms <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "cMono.mild",
  bgdGroups = "cMono.severe"
)

markerList <- getMarkers.raw(mk.cmono.ms, cutOff = "Pval <= 0.05", returnGR=F)
markerList$cMono.mild
write.table(markerList$cMono.mild, file="DiffPeak.cMono.MvS.txt",sep="\t",quote=F)


###########################



mk.ncmono.mp <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "ncMono.mild",
  bgdGroups = "ncMono.post"
)

mk.ncmono.sp <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "ncMono.severe",
  bgdGroups = "ncMono.post"
)

markerList <- getMarkers.raw(mk.ncmono.mp, cutOff = "Pval <= 0.05", returnGR=F)
markerList$ncMono.mild
write.table(markerList$ncMono.mild, file="DiffPeak.ncMono.MvP.txt",sep="\t",quote=F)

markerList <- getMarkers.raw(mk.ncmono.sp, cutOff = "Pval <= 0.05", returnGR=F)
markerList$ncMono.severe
write.table(markerList$ncMono.severe, file="DiffPeak.ncMono.SvP.txt",sep="\t",quote=F)

mk.ncmono.ms <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "ncMono.mild",
  bgdGroups = "ncMono.severe"
)

markerList <- getMarkers.raw(mk.ncmono.ms, cutOff = "Pval <= 0.05", returnGR=F)
markerList$ncMono.mild
write.table(markerList$ncMono.mild, file="DiffPeak.ncMono.MvS.txt",sep="\t",quote=F)


####################################


mk.mp <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CD4 T.mild",
  bgdGroups = "CD4 T.post"
)

mk.sp <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CD4 T.severe",
  bgdGroups = "CD4 T.post"
)

markerList <- getMarkers.raw(mk.mp, cutOff = "Pval <= 0.05", returnGR=F)
write.table(markerList$`CD4 T.mild`, file="DiffPeak.cd4t.MvP.txt",sep="\t",quote=F)

markerList <- getMarkers.raw(mk.sp, cutOff = "Pval <= 0.05", returnGR=F)
write.table(markerList$`CD4 T.severe`, file="DiffPeak.cd4t.SvP.txt",sep="\t",quote=F)


mk.ms <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CD4 T.mild",
  bgdGroups = "CD4 T.severe"
)

markerList <- getMarkers.raw(mk.ms, cutOff = "Pval <= 0.05", returnGR=F)
write.table(markerList$`CD4 T.mild`, file="DiffPeak.cd4t.MvS.txt",sep="\t",quote=F)

#########
mk.mp <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CD8 T.mild",
  bgdGroups = "CD8 T.post"
)

mk.sp <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CD8 T.severe",
  bgdGroups = "CD8 T.post"
)

markerList <- getMarkers.raw(mk.mp, cutOff = "Pval <= 0.05", returnGR=F)
write.table(markerList$`CD8 T.mild`, file="DiffPeak.cd8t.MvP.txt",sep="\t",quote=F)

markerList <- getMarkers.raw(mk.sp, cutOff = "Pval <= 0.05", returnGR=F)
write.table(markerList$`CD8 T.severe`, file="DiffPeak.cd8t.SvP.txt",sep="\t",quote=F)


mk.ms <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CD8 T.mild",
  bgdGroups = "CD8 T.severe"
)

markerList <- getMarkers.raw(mk.ms, cutOff = "Pval <= 0.05", returnGR=F)
write.table(markerList$`CD8 T.mild`, file="DiffPeak.cd8t.MvS.txt",sep="\t",quote=F)

######################

#########
mk.mp <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "B.mild",
  bgdGroups = "B.post"
)

mk.sp <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "B.severe",
  bgdGroups = "B.post"
)

markerList <- getMarkers.raw(mk.mp, cutOff = "Pval <= 0.05", returnGR=F)
write.table(markerList$`B.mild`, file="DiffPeak.B.MvP.txt",sep="\t",quote=F)

markerList <- getMarkers.raw(mk.sp, cutOff = "Pval <= 0.05", returnGR=F)
write.table(markerList$`B.severe`, file="DiffPeak.B.SvP.txt",sep="\t",quote=F)


mk.ms <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "B.mild",
  bgdGroups = "B.severe"
)

markerList <- getMarkers.raw(mk.ms, cutOff = "Pval <= 0.05", returnGR=F)
write.table(markerList$`B.mild`, file="DiffPeak.B.MvS.txt",sep="\t",quote=F)

#########
mk.mp <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "NK.mild",
  bgdGroups = "NK.post"
)

mk.sp <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "NK.severe",
  bgdGroups = "NK.post"
)

markerList <- getMarkers.raw(mk.mp, cutOff = "Pval <= 0.05", returnGR=F)
write.table(markerList$`NK.mild`, file="DiffPeak.NK.MvP.txt",sep="\t",quote=F)

markerList <- getMarkers.raw(mk.sp, cutOff = "Pval <= 0.05", returnGR=F)
write.table(markerList$`NK.severe`, file="DiffPeak.NK.SvP.txt",sep="\t",quote=F)


mk.ms <- getMarkerFeatures(
  ArchRProj = projDm, 
  useMatrix = "PeakMatrix",
  groupBy = "cell.cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "NK.mild",
  bgdGroups = "NK.severe"
)

markerList <- getMarkers.raw(mk.ms, cutOff = "Pval <= 0.05", returnGR=F)
write.table(markerList$`NK.mild`, file="DiffPeak.NK.MvS.txt",sep="\t",quote=F)


