#### DE genes between conditions 
options(future.globals.maxSize = 10737418240)  # 10240*1024^2 = 10737418240  for 10GB

library(Seurat)
library(dplyr)
library(Matrix)
library(Matrix.utils)
library(ggplot2)
library(reshape2)
library(pryr)
library(RColorBrewer)
library(cowplot)


pbmc.pools <- readRDS("path to seurat obj")

Idents(pbmc.pools)="celltypeL0"
#cells <- levels(Idents(pbmc.pools))
cells <- c(
"cMono",  #1
"CD163.cMono",  #4
"ncMono",
"mDC",  #5
"pDC",  #6
"CD4.T",  #7
"NK",  #9
"CD8.T",  #8
"B",  #13
"Plasmablast"  #14
)




all.de <- NULL
for(cell in cells){
  print(cell)
  sub <- subset(pbmc.pools, idents=cell)

  de.uhA <- ""
  try(de.uhA <- FindMarkers(subA, min.pct = 0.1, logfc.threshold = 0.05, ident.1 = "Mild", ident.2 = "post", test.use="wilcox"))
  de.uhA$cell= cell
  de.uhA$comp="mild.vs.post"
  de.uhA$gene= rownames(de.uhA)
  
  de.uhP <- ""
  try(de.uhP <- FindMarkers(subP, min.pct = 0.1, logfc.threshold = 0.05, ident.1 = "Severe", ident.2 = "post", test.use="wilcox"))
  de.uhP$cell= cell
  de.uhP$comp="severe.vs.post"
  de.uhP$gene= rownames(de.uhP)
     
  all.de <- rbind(all.de,de.uhA,de.uhP)
}
pbmc.pools.de <- all.de
pbmc.pools.de$method="wilcox"
saveRDS(pbmc.pools.de,file="DE.wilcox.rds")

all.de <- NULL
for(cell in cells){
  print(cell)
  sub <- subset(pbmc.pools, idents=cell)

  de.uhA <- ""
  try(de.uhA <- FindMarkers(subA, min.pct = 0.1, logfc.threshold = 0.05, ident.1 = "Mild", ident.2 = "post", min.cells.feature = 1, min.cells.group = 1, test.use = "MAST", latent.vars = c("Age")))
  de.uhA$cell= cell
  de.uhA$comp="Mild.vs.light"
  de.uhA$gene= rownames(de.uhA)
  
  de.uhP <- ""
  try(de.uhP <- FindMarkers(subP, min.pct = 0.1, logfc.threshold = 0.05, ident.1 = "Severe", ident.2 = "post", min.cells.feature = 1, min.cells.group = 1, test.use = "MAST", latent.vars = c("Age")))
  de.uhP$cell= cell
  de.uhP$comp="Severe.vs.light"
  de.uhP$gene= rownames(de.uhP)
     
  all.de <- rbind(all.de,de.uhA,de.uhP)
}
pbmc.pools.ma <- all.de
pbmc.pools.ma$method="MAST.age"
saveRDS(pbmc.pools.ma,file="DE.mast.age.rds")


de.all <- rbind(pbmc.pools.de, pbmc.pools.ma)
de.all <- de.all[grep("^RPS|^RPL|^MT-", de.all$gene,invert = T),]
de.all$direction = "Down"
de.all$direction[de.all$avg_logFC>0] <- "Up"
de.sig <- de.all[abs(de.all$avg_logFC) > 0.7 & de.all$p_val_adj < 0.05, ]
de <- de.all[de.all$p_val_adj < 0.05, ]
write.table(de.sig, file="table2.SevsMild.txt",sep="\t",row.names=F,quote=F)
write.table(de, file="table2.SevsMild.all.txt",sep="\t",row.names=F,quote=F)
saveRDS(de.all,file="table2.all.rds")

