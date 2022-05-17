##### calculate co-expression
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

library(Hmisc)
library(corrplot)

pbmc <- readRDS("pbmc.filter.annotev3.rds")


#features = c("STAT1","STAT2","SPI1","LUCAT1","LINC00211","S100A8","S100A9","CHI3L1")

ISG <- grep("^IF",x=VariableFeatures(pbmc), value=T)[1:10]


features = c("STAT1","STAT2","SPI1","LUCAT1","CHI3L1","CEBPA","CEBPB","CEBPD","CEBPE","CEBPG", "IFI27","IFI30","IFITM3")


pbmc$cell.cond <- paste(pbmc$celltypeL1, pbmc$Severity, sep=".")

pdf("plot_ncRNA.pdf",height=12,width=10)
DotPlot(pbmc,features=rev(features),cols="RdBu",group.by="cell.cond")+RotatedAxis()
dev.off()

mono <- readRDS("mono.annoteV5.rds")
mono <- subset(mono, celltypeL1 == "merged.cMonos")

features = c("LUCAT1","CEBPE","CEBPD","SPI1","JUN","FOS","IFI30","IFI27","CCR2")  ### or other features
DotPlot(mono,features=rev(features),cols="RdBu",group.by="Severity")+RotatedAxis()

mono.sub <- mono[features,]
mono1 <- subset(mono.sub, Severity %in% c("mild","severe") )   #### or post 
exp.hla <- as.data.frame(mono1@assays$RNA@data)
t.exp.hla <- t(exp.hla)
#t.exp.hla[t.exp.hla==0]=NA
r.exp.hla <- rcorr(t.exp.hla, type="pearson")
r.exp.hla$logP <- -log10(r.exp.hla$P)
r.exp.hla$logP[r.exp.hla$logP > 10] = 10

#pdf("~/plot_rcorr_2.pdf",width=5,height=5)
corrplot(r.exp.hla$r, p.mat = r.exp.hla$P, sig.level=0.05, method = "square",order="original",pch.col = "grey",pch.cex = 0.4, col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(10),tl.col = "black",tl.cex = 0.7)
#dev.off()




