library(ggplot2)

setwd("~/Desktop/covid19/MHH50/peak2geneLinks_Dm/")
test <- read.table("matched.DEG.txt",header=T,sep="\t",stringsAsFactors = F)

test.cMono <- test[test$cell=="cMono",]

plot <- test.cMono
plot$show <- "not"
plot[plot$avg_logFC > 0.7 & plot$DP_log2FC > 0.5, "show"] <- "yes"
plot[plot$avg_logFC > 0.3 & plot$DP_log2FC > 5, "show"] <- "yes"

ggplot(plot,aes(x=DP_log2FC,y=avg_logFC,color=comp))+geom_point()+theme_bw()+xlab("DiffPeak log2FC")+ylab("DEGs logFC")+geom_label_repel(aes(label = ifelse(show  == "yes", as.character(gene), '')), hjust = 1.25, vjust = 0, size = 3) 


test.ncMono <- test[test$cell=="ncMono",]
#ggplot(test.ncMono,aes(x=DP_log2FC,y=avg_logFC,color=comp))+geom_point()+theme_bw()+xlab("DiffPeak log2FC")+ylab("DEGs logFC")
plot <- test.ncMono
plot$show <- "not"
plot[plot$avg_logFC > 0.7 & plot$DP_log2FC > 0.5, "show"] <- "yes"
plot[plot$avg_logFC > 0.3 & plot$DP_log2FC > 5, "show"] <- "yes"

ggplot(plot,aes(x=DP_log2FC,y=avg_logFC,color=comp))+geom_point()+theme_bw()+xlab("DiffPeak log2FC")+ylab("DEGs logFC")+geom_label_repel(aes(label = ifelse(show  == "yes", as.character(gene), '')), hjust = 1.25, vjust = 0, size = 3) 


########## Figure 3C
test$dp="Up"
test$dp[test$DP_log2FC<0]="Down"
test$cell.dd <- paste(test$cell,test$direction,test$dp,sep=".")

plot2 <- as.data.frame(table(test$cell.dd,test$comp))

plot2$direct <- rep(c("Dw,Dw","Dw,Up","Up,Dw","Up,Up"),12)
ggplot(plot2, aes(x=Var1,y=Freq,fill=direct))+geom_bar(stat="identity")+facet_grid(.~Var2)+theme_bw()+theme(axis.text.x=element_text(angle = 70,hjust = 1))

t(table(test$cell.dd,test$comp))

stats.mild <- NULL
stats.seve <- NULL
#B
nums <- matrix(c(23, 17, 4, 14),2,2)
fisher.test(nums)$p
stats1$cell="B"
stats1$p=fisher.test(nums)$p
stats.mild <- as.data.frame(rbind(stats.mild, stats1))

nums <- matrix(c(1, 2, 13, 22),2,2)
fisher.test(nums)$p
stats2$cell="B"
stats2$p=fisher.test(nums)$p
stats.seve <- as.data.frame(rbind(stats.seve, stats2))

#CD4T
nums <- matrix(c(22, 48, 48, 82),2,2)
fisher.test(nums)$p
stats1$cell="CD4T"
stats1$p=fisher.test(nums)$p
stats.mild <- as.data.frame(rbind(stats.mild, stats1))

nums <- matrix(c(15, 20, 72, 94),2,2)
fisher.test(nums)$p
stats2$cell="CD4T"
stats2$p=fisher.test(nums)$p
stats.seve <- as.data.frame(rbind(stats.seve, stats2))

#CD8T
nums <- matrix(c(17, 9, 81, 103),2,2)
fisher.test(nums)$p
stats1$cell="CD8T"
stats1$p=fisher.test(nums)$p
stats.mild <- as.data.frame(rbind(stats.mild, stats1))

nums <- matrix(c(12, 13, 114, 122),2,2)
fisher.test(nums)$p
stats2$cell="CD8T"
stats2$p=fisher.test(nums)$p
stats.seve <- as.data.frame(rbind(stats.seve, stats2))

#cMono
nums <- matrix(c(43, 87, 248, 1130),2,2)
fisher.test(nums)$p
stats1$cell="cMono"
stats1$p=fisher.test(nums)$p
stats.mild <- as.data.frame(rbind(stats.mild, stats1))


nums <- matrix(c(79, 139, 300, 993),2,2)
fisher.test(nums)$p
stats2$cell="cMono"
stats2$p=fisher.test(nums)$p
stats.seve <- as.data.frame(rbind(stats.seve, stats2))



#ncMono
nums <- matrix(c(9, 25, 83, 444),2,2)
fisher.test(nums)$p
stats1$cell="ncMono"
stats1$p=fisher.test(nums)$p
stats.mild <- as.data.frame(rbind(stats.mild, stats1))

nums <- matrix(c(15, 3, 158, 104),2,2)
fisher.test(nums)$p
stats2$cell="ncMono"
stats2$p=fisher.test(nums)$p
stats.seve <- as.data.frame(rbind(stats.seve, stats2))


#NK
nums <- matrix(c(4, 12, 98, 378),2,2)
fisher.test(nums)$p
stats1$cell="NK"
stats1$p=fisher.test(nums)$p
stats.mild <- as.data.frame(rbind(stats.mild, stats1))

nums <- matrix(c(1, 4, 154, 209),2,2)
fisher.test(nums)$p
stats2$cell="NK"
stats2$p=fisher.test(nums)$p
stats.seve <- as.data.frame(rbind(stats.seve, stats2))



stats.mild$fdr <- p.adjust(stats.mild$p, method = "fdr")
stats.seve$fdr <- p.adjust(stats.seve$p, method = "fdr")

### results:
#> stats.mild#         cell            p          fdr#1      B   0.02171161 0.0651348427#2        CD4T    0.5343862 0.6412634482#3        CD8T    0.0575685 0.1151370057#4       cMono 0.0001005059 0.0006030352#5      ncMono    0.1466928 0.2200392375#6          NK    0.7529633 0.7529632925#> stats.seve#         cell            p          fdr#2      B            1 1.0000000000#2        CD4T            1 1.0000000000#3        CD8T            1 1.0000000000#4       cMono 6.800869e-05 0.0004080522#5      ncMono   0.07675458 0.2302637469#6          NK    0.4025796 0.8051591937




##########

spi1 <- read.table("Mono.SPI1.txt",header=T,sep="\t",stringsAsFactors = F)
gene1 <- unique(spi1$gene)
cebpd <- read.table("Mono.CEBPD.txt",header=T,sep="\t",stringsAsFactors = F)
gene2 <- unique(cebpd$gene)


library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)

gene1 <- bitr(gene1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene2 <- bitr(gene2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")




gene1$Label="SPI1"
gene2$Label="CEBPD"
de.list <- rbind(gene1,gene2)
de.list
kegg <- compareCluster(ENTREZID~Label, data=de.list, fun="enrichKEGG")

compareCluster(ENTREZID~Label, data= list.de, fun="enrichGO", OrgDb="org.Hs.eg.db", ont="BP")







tf <- read.table("/vol/projects/CIIM/Covid_50MHH/scATAC/analysis_new/gwas_meta_test/tf.txt", stringsAsFactors=F)[,1]

pbmc.uniq$cell.cond <- paste(pbmc.uniq$celltypeL0, pbmc.uniq$Severity2,sep=".")
pdf("test_tf_plot1.pdf", height=20, width=10)
DotPlot(pbmc.uniq,features = unique(tf),col.min= -1.5, col.max= 1.5,group.by="cell.cond")+
scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n =100, name = "RdBu")))+coord_flip()+xlab("")+ylab("")+theme(axis.text.x = element_text(angle = 45,hjust=1), axis.text.y=element_text(size = 8))
dev.off()


exp <- Get.Exp.Genes(subset(pbmc.uniq,Severity2 == "Severe"),pct=0.05)
intersect(tf,exp)


gene <- c("FYCO1", "CCR1", "CXCR6", "CCR1","CCR2","CCRL2", "DTX1", "OAS1", "OAS3", "OAS2")

pdf("test_gene_plot1.pdf", height=5, width=10)
DotPlot(pbmc.uniq,features = unique(gene),col.min= -1.5, col.max= 1.5,group.by="cell.cond")+
scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n =100, name = "RdBu")))+coord_flip()+xlab("")+ylab("")+theme(axis.text.x = element_text(angle = 45,hjust=1), axis.text.y=element_text(size = 8))
dev.off()

