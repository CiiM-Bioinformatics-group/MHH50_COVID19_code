
#### calculate AUCell scores for monocytes subclusters

library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(Matrix.utils)
library(reshape2)
library(pryr)
library(RColorBrewer)
require(AUCell)

options(future.globals.maxSize=20368709120)
mono <- readRDS("path to monocytes scRNAseq obj")  #mono.annoteV5.rds
DimPlot(mono)


####### AUCell scores from monocytes and macrophages markers of Covid-19 lungs ####
marker <- read.table("marker_leif.txt",sep="\t",stringsAsFactors = F)
colnames(marker) <- c("cell","gene")

Idents(mono) <- "new.id"
tmp<-as.data.frame(as.matrix(mono@assays$RNA@data))
filtered_matrix<-rowSums(tmp > 0)
filtered_matrix<-filtered_matrix[filtered_matrix>10]
tmp <- tmp[row.names(tmp) %in% names(filtered_matrix),]
cells_rankings <- AUCell_buildRankings(as.matrix(tmp))


markers_list <- list("FCN1.Mono" = marker$gene[marker$cell == "FCN1-Mono"], 'CD163.Mf' = marker$gene[marker$cell == "CD163/LGMN-Mf"], "mono.mac"=marker$gene[marker$cell == "Mono/Mf"])

cells_AUC1 <- AUCell_calcAUC(markers_list$FCN1.Mono, cells_rankings, aucMaxRank = ceiling(0.03 * nrow(cells_rankings)), normAUC = T)

tmp <- getAUC(cells_AUC1)
tmp <- as.data.frame(t(tmp))
  
  df <- data.frame(row.names = row.names(mono@meta.data), cluster = mono$new.id, stringsAsFactors = F)
  df <- as.data.frame(merge(df, tmp, by = 0))
  df$geneset<-df[,"geneSet"]   ## this can be cell type name when input a list
  
  df %>% group_by(cluster) %>% dplyr::summarize(Mean = mean(geneset, na.rm=TRUE)) -> df.tmp
  df <- merge(df, df.tmp, key = "cluster", all = T)
  df.th17 <- df
  #df$geneset<-df[,i]
  gg.th17 <- ggplot(df.th17, aes(y = geneset, x = cluster, group = cluster, fill = Mean)) + 
    geom_violin(scale = "area", trim = T, adjust = 3, draw_quantiles = c(0.25, 0.5, 0.75), color = "grey30")+
 	scale_fill_gradient2(low="#2C00F5", mid="#D084A2",high="#F5F700",midpoint =0.4)+
    #scale_fill_gradient2(low="#74add1",mid=,high="#d73027") +#geom_jitter(size=0.1,width=0.1)+
    ylab("AUC score") + xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(color="black", size = .3),
          axis.line.y = element_line(color="black", size = .3))+ylim(0,0.8) +ggtitle("FCN1-mono markers score")
          
  pdf("Leif.AUCell.fcn.pdf",width=4.5, height=3.3)            
 plot(gg.th17)
 dev.off()         
          
 cells_AUC2 <- AUCell_calcAUC(markers_list$CD163.Mf, cells_rankings, aucMaxRank = ceiling(0.03 * nrow(cells_rankings)), normAUC = T)       
 tmp <- getAUC(cells_AUC2)     
 tmp <- as.data.frame(t(tmp))
  
  df <- data.frame(row.names = row.names(mono@meta.data), cluster = mono$new.id, stringsAsFactors = F)
  df <- as.data.frame(merge(df, tmp, by = 0))
  df$geneset<-df[,"geneSet"]   ## this can be cell type name when input a list
  
  df %>% group_by(cluster) %>% dplyr::summarize(Mean = mean(geneset, na.rm=TRUE)) -> df.tmp
  df <- merge(df, df.tmp, key = "cluster", all = T)
  df.th17 <- df
  #df$geneset<-df[,i]
  gg.th17 <- ggplot(df.th17, aes(y = geneset, x = cluster, group = cluster, fill = Mean)) + 
    geom_violin(scale = "area", trim = T, adjust = 3, draw_quantiles = c(0.25, 0.5, 0.75), color = "grey30")+
 	scale_fill_gradient2(low="#2C00F5", mid="#D084A2",high="#F5F700",midpoint =0.4)+
    #scale_fill_gradient2(low="#74add1",mid=,high="#d73027") +#geom_jitter(size=0.1,width=0.1)+
    ylab("AUC score") + xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(color="black", size = .3),
          axis.line.y = element_line(color="black", size = .3))+ylim(0,0.8) +ggtitle("CD163-mac markers score")         
 
 pdf("Leif.AUCell.mac.pdf",width=4.5, height=3.3)            
 plot(gg.th17)
 dev.off()         
             
        
        
        cells_AUC1 <- AUCell_calcAUC(markers_list$mono.mac, cells_rankings, aucMaxRank = ceiling(0.03 * nrow(cells_rankings)), normAUC = T)


tmp <- getAUC(cells_AUC1)

tmp <- as.data.frame(t(tmp))
  
  df <- data.frame(row.names = row.names(mono@meta.data), cluster = mono$new.id, stringsAsFactors = F)
  df <- as.data.frame(merge(df, tmp, by = 0))
  df$geneset<-df[,"geneSet"]   ## this can be cell type name when input a list
  
  df %>% group_by(cluster) %>% dplyr::summarize(Mean = mean(geneset, na.rm=TRUE)) -> df.tmp
  df <- merge(df, df.tmp, key = "cluster", all = T)
  df.th17 <- df
  #df$geneset<-df[,i]
  gg.th17 <- ggplot(df.th17, aes(y = geneset, x = cluster, group = cluster, fill = Mean)) + 
    geom_violin(scale = "area", trim = T, adjust = 3, draw_quantiles = c(0.25, 0.5, 0.75), color = "grey30")+
 	scale_fill_gradient2(low="#2C00F5", mid="#D084A2",high="#F5F700",midpoint =0.4)+
    #scale_fill_gradient2(low="#74add1",mid=,high="#d73027") +#geom_jitter(size=0.1,width=0.1)+
    ylab("AUC score") + xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(color="black", size = .3),
          axis.line.y = element_line(color="black", size = .3))+ylim(0,0.8) +ggtitle("Mono/Mac markers score")
          
  pdf("Leif.AUCell.monomac.pdf",width=4.5, height=3.3)            
 plot(gg.th17)
 dev.off()     
          
          
          
          
          
          
 ######### Cell markers  from covid19  PBMC cohort1 (berlin) ###
 
 
 marker <- read.table("clustermarkers.cohort1.txt",sep="\t",stringsAsFactors = F,header=T)

marker.mono <- marker[marker$cluster %in% c(0,1,2,3,4),]

top30.mono <- marker.mono %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)

markers_list <- list("cMono" = top30.mono$gene[top30.mono$cluster == "0"], 'CD83hi.mono' = top30.mono$gene[top30.mono$cluster == "1"], 
                     'CD163hi.mono' = top30.mono$gene[top30.mono$cluster == "2"], 'S100Ahi.mono' = top30.mono$gene[top30.mono$cluster == "3"], 'ncMono'=top30.mono$gene[top30.mono$cluster == "4"] )
                     

dataframe<- data.frame("cMono" = top30.mono$gene[top30.mono$cluster == "0"], 'CD83hi.mono' = top30.mono$gene[top30.mono$cluster == "1"], 
                     'CD163hi.mono' = top30.mono$gene[top30.mono$cluster == "2"], 'S100Ahi.mono' = top30.mono$gene[top30.mono$cluster == "3"], 'ncMono'=top30.mono$gene[top30.mono$cluster == "4"] )

write.table(dataframe,"Cohort1.marker.txt",sep="\t",row.names=F,quote=F )

cells_AUC1 <- AUCell_calcAUC(markers_list, cells_rankings, aucMaxRank = ceiling(0.03 * nrow(cells_rankings)), normAUC = T)
#cells_AUC2 <- AUCell_calcAUC(markers_list$CD163.Mf, cells_rankings, aucMaxRank = ceiling(0.03 * nrow(cells_rankings)), normAUC = T)

tmp <- getAUC(cells_AUC1)
#tmp <- getAUC(cells_AUC2)
tmp <- as.data.frame(t(tmp))
 
  df <- data.frame(row.names = row.names(mono@meta.data), cluster = mono$new.id, stringsAsFactors = F)
  df <- as.data.frame(merge(df, tmp, by = 0))
  #df$geneset<-df[,"geneSet"]   ## this can be cell type name when input a list
  df$geneset<-df[,"S100Ahi.mono"]   ## this can be cell type name when input a list
  
  df %>% group_by(cluster) %>% dplyr::summarize(Mean = mean(geneset, na.rm=TRUE)) -> df.tmp
  df <- merge(df, df.tmp, key = "cluster", all = T)
  df.th17 <- df
  #df$geneset<-df[,i]
  gg.th17 <- ggplot(df.th17, aes(y = geneset, x = cluster, group = cluster, fill = Mean)) + 
    geom_violin(scale = "area", trim = T, adjust = 3, draw_quantiles = c(0.25, 0.5, 0.75), color = "grey30")+
 	scale_fill_gradient2(low="#2C00F5", mid="#D084A2",high="#F5F700",midpoint =0.4)+
    #scale_fill_gradient2(low="#74add1",mid=,high="#d73027") +#geom_jitter(size=0.1,width=0.1)+
    ylab("AUC score") + xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(color="black", size = .3),
          axis.line.y = element_line(color="black", size = .3))+ylim(0,0.8) +ggtitle("S100A-hi Mono. markers score")         
 pdf("cohort1.AUCell.3.pdf",width=4.5, height=3.3)            
 plot(gg.th17)
 dev.off()
 
   df <- data.frame(row.names = row.names(mono@meta.data), cluster = mono$new.id, stringsAsFactors = F)
  df <- as.data.frame(merge(df, tmp, by = 0))
 df$geneset<-df[,"cMono"]   ## this can be cell type name when input a list
   df %>% group_by(cluster) %>% dplyr::summarize(Mean = mean(geneset, na.rm=TRUE)) -> df.tmp
  df <- merge(df, df.tmp, key = "cluster", all = T)
  df.th17 <- df
  #df$geneset<-df[,i]
  gg.th17 <- ggplot(df.th17, aes(y = geneset, x = cluster, group = cluster, fill = Mean)) + 
    geom_violin(scale = "area", trim = T, adjust = 3, draw_quantiles = c(0.25, 0.5, 0.75), color = "grey30")+
 	scale_fill_gradient2(low="#2C00F5", mid="#D084A2",high="#F5F700",midpoint =0.4)+
    #scale_fill_gradient2(low="#74add1",mid=,high="#d73027") +#geom_jitter(size=0.1,width=0.1)+
    ylab("AUC score") + xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(color="black", size = .3),
          axis.line.y = element_line(color="black", size = .3))+ylim(0,0.8) +ggtitle("Classical Mono. markers score")             
 pdf("cohort1.AUCell.0.pdf",width=4.5, height=3.3)            
 plot(gg.th17)
 dev.off()
          
          
    df <- data.frame(row.names = row.names(mono@meta.data), cluster = mono$new.id, stringsAsFactors = F)
  df <- as.data.frame(merge(df, tmp, by = 0))        
  df$geneset<-df[,"CD83hi.mono"]   ## this can be cell type name when input a list
    df %>% group_by(cluster) %>% dplyr::summarize(Mean = mean(geneset, na.rm=TRUE)) -> df.tmp
  df <- merge(df, df.tmp, key = "cluster", all = T)
  df.th17 <- df
  #df$geneset<-df[,i]
  gg.th17 <- ggplot(df.th17, aes(y = geneset, x = cluster, group = cluster, fill = Mean)) + 
    geom_violin(scale = "area", trim = T, adjust = 3, draw_quantiles = c(0.25, 0.5, 0.75), color = "grey30")+
 	scale_fill_gradient2(low="#2C00F5", mid="#D084A2",high="#F5F700",midpoint =0.4)+
    #scale_fill_gradient2(low="#74add1",mid=,high="#d73027") +#geom_jitter(size=0.1,width=0.1)+
    ylab("AUC score") + xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(color="black", size = .3),
          axis.line.y = element_line(color="black", size = .3))+ylim(0,0.8) +ggtitle("CD83-hi Mono. markers score")             
 pdf("cohort1.AUCell.1.pdf",width=4.5, height=3.3)            
 plot(gg.th17)
 dev.off() 
  
  
    df <- data.frame(row.names = row.names(mono@meta.data), cluster = mono$new.id, stringsAsFactors = F)
  df <- as.data.frame(merge(df, tmp, by = 0))
  df$geneset<-df[,"CD163hi.mono"]   ## this can be cell type name when input a list
      df %>% group_by(cluster) %>% dplyr::summarize(Mean = mean(geneset, na.rm=TRUE)) -> df.tmp
  df <- merge(df, df.tmp, key = "cluster", all = T)
  df.th17 <- df
  #df$geneset<-df[,i]
  gg.th17 <- ggplot(df.th17, aes(y = geneset, x = cluster, group = cluster, fill = Mean)) + 
    geom_violin(scale = "area", trim = T, adjust = 3, draw_quantiles = c(0.25, 0.5, 0.75), color = "grey30")+
 	scale_fill_gradient2(low="#2C00F5", mid="#D084A2",high="#F5F700",midpoint =0.4)+
    #scale_fill_gradient2(low="#74add1",mid=,high="#d73027") +#geom_jitter(size=0.1,width=0.1)+
    ylab("AUC score") + xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(color="black", size = .3),
          axis.line.y = element_line(color="black", size = .3))+ylim(0,0.8) +ggtitle("CD163-hi markers score")             
  pdf("cohort1.AUCell.2.pdf",width=4.5, height=3.3)            
 plot(gg.th17)
 dev.off()
 
 
   df <- data.frame(row.names = row.names(mono@meta.data), cluster = mono$new.id, stringsAsFactors = F)
  df <- as.data.frame(merge(df, tmp, by = 0))
   df$geneset<-df[,"ncMono"]   ## this can be cell type name when input a list
      df %>% group_by(cluster) %>% dplyr::summarize(Mean = mean(geneset, na.rm=TRUE)) -> df.tmp
  df <- merge(df, df.tmp, key = "cluster", all = T)
  df.th17 <- df
  #df$geneset<-df[,i]
  gg.th17 <- ggplot(df.th17, aes(y = geneset, x = cluster, group = cluster, fill = Mean)) + 
    geom_violin(scale = "area", trim = T, adjust = 3, draw_quantiles = c(0.25, 0.5, 0.75), color = "grey30")+
 	scale_fill_gradient2(low="#2C00F5", mid="#D084A2",high="#F5F700",midpoint =0.4)+
    #scale_fill_gradient2(low="#74add1",mid=,high="#d73027") +#geom_jitter(size=0.1,width=0.1)+
    ylab("AUC score") + xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(color="black", size = .3),
          axis.line.y = element_line(color="black", size = .3))+ylim(0,0.8) +ggtitle("ncMono markers score")             
  pdf("cohort1.AUCell.4.pdf",width=4.5, height=3.3)            
 plot(gg.th17)
 dev.off()
 
 
 
 
 
 
 
 
 
 
        
 ######### Cell markers  from covid19  PBMC cohort2 (bonn) ###
 
 
 marker <- read.table("clustermarkers.cohort2.txt",sep="\t",stringsAsFactors = F,header=T)

marker.mono <- marker[marker$cluster %in% c(0,1,2,3,4),]

top30.mono <- marker.mono %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)

markers_list <- list("cMono" = top30.mono$gene[top30.mono$cluster == "0"], 'CD83hi.mono' = top30.mono$gene[top30.mono$cluster == "1"], 
                     'CD163hi.mono' = top30.mono$gene[top30.mono$cluster == "2"], 'S100Ahi.mono' = top30.mono$gene[top30.mono$cluster == "3"], 'ncMono'=top30.mono$gene[top30.mono$cluster == "4"] )
                       

cells_AUC1 <- AUCell_calcAUC(markers_list, cells_rankings, aucMaxRank = ceiling(0.03 * nrow(cells_rankings)), normAUC = T)
#cells_AUC2 <- AUCell_calcAUC(markers_list$CD163.Mf, cells_rankings, aucMaxRank = ceiling(0.03 * nrow(cells_rankings)), normAUC = T)

tmp <- getAUC(cells_AUC1)
#tmp <- getAUC(cells_AUC2)
tmp <- as.data.frame(t(tmp))
  
  df <- data.frame(row.names = row.names(mono@meta.data), cluster = mono$new.id, stringsAsFactors = F)
  df <- as.data.frame(merge(df, tmp, by = 0))


  #df$geneset<-df[,"geneSet"]   ## this can be cell type name when input a list
  

  df$geneset<-df[,"S100Ahi.mono"]   ## this can be cell type name when input a list
  
  df %>% group_by(cluster) %>% dplyr::summarize(Mean = mean(geneset, na.rm=TRUE)) -> df.tmp
  df <- merge(df, df.tmp, key = "cluster", all = T)
  df.th17 <- df
  #df$geneset<-df[,i]
  gg.th17 <- ggplot(df.th17, aes(y = geneset, x = cluster, group = cluster, fill = Mean)) + 
    geom_violin(scale = "area", trim = T, adjust = 3, draw_quantiles = c(0.25, 0.5, 0.75), color = "grey30")+
 	scale_fill_gradient2(low="#2C00F5", mid="#D084A2",high="#F5F700",midpoint =0.3)+
    #scale_fill_gradient2(low="#74add1",mid=,high="#d73027") +#geom_jitter(size=0.1,width=0.1)+
    ylab("AUC score") + xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(color="black", size = .3),
          axis.line.y = element_line(color="black", size = .3))+ylim(0,0.8) +ggtitle("S100A-hi Mono. markers score")         
 pdf("cohort2.AUCell.3.pdf",width=4.5, height=3.3)            
 plot(gg.th17)
 dev.off()
 
   df <- data.frame(row.names = row.names(mono@meta.data), cluster = mono$new.id, stringsAsFactors = F)
  df <- as.data.frame(merge(df, tmp, by = 0))
 df$geneset<-df[,"cMono"]   ## this can be cell type name when input a list
   df %>% group_by(cluster) %>% dplyr::summarize(Mean = mean(geneset, na.rm=TRUE)) -> df.tmp
  df <- merge(df, df.tmp, key = "cluster", all = T)
  df.th17 <- df
  #df$geneset<-df[,i]
  gg.th17 <- ggplot(df.th17, aes(y = geneset, x = cluster, group = cluster, fill = Mean)) + 
    geom_violin(scale = "area", trim = T, adjust = 3, draw_quantiles = c(0.25, 0.5, 0.75), color = "grey30")+
 	scale_fill_gradient2(low="#2C00F5", mid="#D084A2",high="#F5F700",midpoint =0.3)+
    #scale_fill_gradient2(low="#74add1",mid=,high="#d73027") +#geom_jitter(size=0.1,width=0.1)+
    ylab("AUC score") + xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(color="black", size = .3),
          axis.line.y = element_line(color="black", size = .3))+ylim(0,0.8) +ggtitle("Classical Mono. markers score")             
 pdf("cohort2.AUCell.0.pdf",width=4.5, height=3.3)            
 plot(gg.th17)
 dev.off()
          
          
    df <- data.frame(row.names = row.names(mono@meta.data), cluster = mono$new.id, stringsAsFactors = F)
  df <- as.data.frame(merge(df, tmp, by = 0))        
  df$geneset<-df[,"CD83hi.mono"]   ## this can be cell type name when input a list
    df %>% group_by(cluster) %>% dplyr::summarize(Mean = mean(geneset, na.rm=TRUE)) -> df.tmp
  df <- merge(df, df.tmp, key = "cluster", all = T)
  df.th17 <- df
  #df$geneset<-df[,i]
  gg.th17 <- ggplot(df.th17, aes(y = geneset, x = cluster, group = cluster, fill = Mean)) + 
    geom_violin(scale = "area", trim = T, adjust = 3, draw_quantiles = c(0.25, 0.5, 0.75), color = "grey30")+
 	scale_fill_gradient2(low="#2C00F5", mid="#D084A2",high="#F5F700",midpoint =0.3)+
    #scale_fill_gradient2(low="#74add1",mid=,high="#d73027") +#geom_jitter(size=0.1,width=0.1)+
    ylab("AUC score") + xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(color="black", size = .3),
          axis.line.y = element_line(color="black", size = .3))+ylim(0,0.8) +ggtitle("CD83-hi Mono. markers score")             
 pdf("cohort2.AUCell.1.pdf",width=4.5, height=3.3)            
 plot(gg.th17)
 dev.off() 
  
  
  
    df <- data.frame(row.names = row.names(mono@meta.data), cluster = mono$new.id, stringsAsFactors = F)
  df <- as.data.frame(merge(df, tmp, by = 0))
  df$geneset<-df[,"CD163hi.mono"]   ## this can be cell type name when input a list
      df %>% group_by(cluster) %>% dplyr::summarize(Mean = mean(geneset, na.rm=TRUE)) -> df.tmp
  df <- merge(df, df.tmp, key = "cluster", all = T)
  df.th17 <- df
  #df$geneset<-df[,i]
  gg.th17 <- ggplot(df.th17, aes(y = geneset, x = cluster, group = cluster, fill = Mean)) + 
    geom_violin(scale = "area", trim = T, adjust = 3, draw_quantiles = c(0.25, 0.5, 0.75), color = "grey30")+
 	scale_fill_gradient2(low="#2C00F5", mid="#D084A2",high="#F5F700",midpoint =0.3)+
    #scale_fill_gradient2(low="#74add1",mid=,high="#d73027") +#geom_jitter(size=0.1,width=0.1)+
    ylab("AUC score") + xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(color="black", size = .3),
          axis.line.y = element_line(color="black", size = .3))+ylim(0,0.8) +ggtitle("CD163-hi markers score")             
  pdf("cohort2.AUCell.2.pdf",width=4.5, height=3.3)            
 plot(gg.th17)
 dev.off()
 
 
   df <- data.frame(row.names = row.names(mono@meta.data), cluster = mono$new.id, stringsAsFactors = F)
  df <- as.data.frame(merge(df, tmp, by = 0))
   df$geneset<-df[,"ncMono"]   ## this can be cell type name when input a list
      df %>% group_by(cluster) %>% dplyr::summarize(Mean = mean(geneset, na.rm=TRUE)) -> df.tmp
  df <- merge(df, df.tmp, key = "cluster", all = T)
  df.th17 <- df
  #df$geneset<-df[,i]
  gg.th17 <- ggplot(df.th17, aes(y = geneset, x = cluster, group = cluster, fill = Mean)) + 
    geom_violin(scale = "area", trim = T, adjust = 3, draw_quantiles = c(0.25, 0.5, 0.75), color = "grey30")+
 	scale_fill_gradient2(low="#2C00F5", mid="#D084A2",high="#F5F700",midpoint =0.3)+
    #scale_fill_gradient2(low="#74add1",mid=,high="#d73027") +#geom_jitter(size=0.1,width=0.1)+
    ylab("AUC score") + xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(color="black", size = .3),
          axis.line.y = element_line(color="black", size = .3))+ylim(0,0.8) +ggtitle("ncMono markers score")             
  pdf("cohort2.AUCell.4.pdf",width=4.5, height=3.3)            
 plot(gg.th17)
 dev.off()