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


pbmc.pools <- readRDS("../pbmc.uniq.annoteV4.rds")

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




############# active vs conv. ######
Idents(pbmc.pools)<-"Severity"
pbmc.pools <- RenameIdents(pbmc.pools,
`severe`="active",
`mild`="active",
`post`="convalescent")
pbmc.pools$Disease <- Idents(pbmc.pools)

table(pbmc.pools$Severity)
table(pbmc.pools$Disease)


Idents(pbmc.pools)="celltypeL0"

all.de <- NULL
for(cell in cells){
  print(cell)
  sub <- subset(pbmc.pools, idents=cell)
  Idents(sub)="Disease"
  de.uhA <- ""
  try(de.uhA <- FindMarkers(sub, min.pct = 0.1, logfc.threshold = 0.05, ident.1 = "active", ident.2 = "convalescent", test.use="wilcox"))
  de.uhA$cell= cell
  de.uhA$comp="active.vs.convalescent"
  de.uhA$gene= rownames(de.uhA)

  all.de <- rbind(all.de,de.uhA)
}

#### test among all cMonos #####
  print("cMono + CD163.cMono = merged.cMonos")
  cell = "merged.cMonos"
  sub <- subset(pbmc.pools, celltypeL0 %in% c("cMono", "CD163.cMono"))
  Idents(sub)="Disease"
  de.uhA <- ""
  try(de.uhA <- FindMarkers(sub, min.pct = 0.1, logfc.threshold = 0.05, ident.1 = "active", ident.2 = "convalescent", test.use="wilcox"))
  de.uhA$cell= cell
  de.uhA$comp="active.vs.convalescent"
  de.uhA$gene= rownames(de.uhA)


  all.de <- rbind(all.de,de.uhA)

#########

pbmc.pools.de <- all.de
pbmc.pools.de$method="wilcox"
saveRDS(pbmc.pools.de,file="DE2cond.wilcox.rds")

de.all <- pbmc.pools.de
de.all <- de.all[grep("^RPS|^RPL|^MT-", de.all$gene,invert = T),]
de.all$direction = "Down"
de.all$direction[de.all$avg_log2FC>0] <- "Up"
de.sig <- de.all[abs(de.all$avg_log2FC) > 1 & de.all$p_val_adj < 0.05, ]
de <- de.all[de.all$p_val_adj < 0.05, ]
write.table(de.sig, file="DEtable.2cond.log2fc1.txt",sep="\t",row.names=F,quote=F)
write.table(de, file="DEtable.2cond.txt",sep="\t",row.names=F,quote=F)





#############  compare between severity #####

all.de <- NULL
for(cell in cells){
  print(cell)
  sub <- subset(pbmc.pools, idents=cell)
  Idents(sub)="Severity"
  de.uhA <- ""
  try(de.uhA <- FindMarkers(sub, min.pct = 0.1, logfc.threshold = 0.05, ident.1 = "mild", ident.2 = "post", test.use="wilcox"))
  de.uhA$cell= cell
  de.uhA$comp="Mild.vs.Post"
  de.uhA$gene= rownames(de.uhA)

  de.uhP <- ""
  try(de.uhP <- FindMarkers(sub, min.pct = 0.1, logfc.threshold = 0.05, ident.1 = "severe", ident.2 = "post", test.use="wilcox"))
  de.uhP$cell= cell
  de.uhP$comp="Severe.vs.Post"
  de.uhP$gene= rownames(de.uhP)

  de.uhQ <- ""
  try(de.uhQ <- FindMarkers(sub, min.pct = 0.1, logfc.threshold = 0.05, ident.1 = "severe", ident.2 = "mild", test.use="wilcox"))
  de.uhQ$cell= cell
  de.uhQ$comp="Severe.vs.Mild"
  de.uhQ$gene= rownames(de.uhQ)

  all.de <- rbind(all.de,de.uhA,de.uhP,de.uhQ)
}

#### test among all cMonos #####
  print("cMono + CD163.cMono = merged.cMonos")
  cell = "merged.cMonos"
  sub <- subset(pbmc.pools, celltypeL0 %in% c("cMono", "CD163.cMono"))
  Idents(sub)="Severity"
  de.uhA <- ""
  try(de.uhA <- FindMarkers(sub, min.pct = 0.1, logfc.threshold = 0.05, ident.1 = "mild", ident.2 = "post", test.use="wilcox"))
  de.uhA$cell= cell
  de.uhA$comp="Mild.vs.Post"
  de.uhA$gene= rownames(de.uhA)

  de.uhP <- ""
  try(de.uhP <- FindMarkers(sub, min.pct = 0.1, logfc.threshold = 0.05, ident.1 = "severe", ident.2 = "post", test.use="wilcox"))
  de.uhP$cell= cell
  de.uhP$comp="Severe.vs.Post"
  de.uhP$gene= rownames(de.uhP)

  de.uhQ <- ""
  try(de.uhQ <- FindMarkers(sub, min.pct = 0.1, logfc.threshold = 0.05, ident.1 = "severe", ident.2 = "mild", test.use="wilcox"))
  de.uhQ$cell= cell
  de.uhQ$comp="Severe.vs.Mild"
  de.uhQ$gene= rownames(de.uhQ)

  all.de <- rbind(all.de,de.uhA,de.uhP,de.uhQ)

#########

pbmc.pools.de <- all.de
pbmc.pools.de$method="wilcox"
saveRDS(pbmc.pools.de,file="DEuniq.wilcox.rds")

de.all <- pbmc.pools.de
de.all <- de.all[grep("^RPS|^RPL|^MT-", de.all$gene,invert = T),]
de.all$direction = "Down"
de.all$direction[de.all$avg_log2FC>0] <- "Up"
de.sig <- de.all[abs(de.all$avg_log2FC) > 1 & de.all$p_val_adj < 0.05, ]
de <- de.all[de.all$p_val_adj < 0.05, ]
write.table(de.sig, file="DEtable.uniq.log2fc1.txt",sep="\t",row.names=F,quote=F)
write.table(de, file="DEtable.uniq.txt",sep="\t",row.names=F,quote=F)



