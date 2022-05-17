
library(ggplot2)

dp <- read.table("./DiffPeak/combine.cMono.txt",sep="\t",header = T)
#dp <- read.table("./DiffPeak/combine.ncMono.txt",sep="\t",header = T)



dp$sig = "not_sig"
dp$sig[dp$Pval1 < 0.05] = "sig_Mild"
dp$sig[dp$Pval2 < 0.05] = "sig_Severe"
dp$sig[dp$Pval1 < 0.05 & dp$Pval2 < 0.05] = "sig_Both"

dp.sig <- dp[dp$sig == "sig_Both",]
dp.sig <- dp[dp$sig != "not_sig",]


p <- ggplot(dp.sig, aes(x=log2FC_1, y=log2FC_2)) + geom_point(aes(color=factor(sig)), size=1, show.legend = T)
p2 <- p +  geom_hline(aes(yintercept=0),linetype="dashed") +
    geom_vline(aes(xintercept=0),linetype="dashed") + 
    theme_bw() + #xlim(c(-2,2))+ ylim(c(-2,2)) +
    geom_label_repel(aes(label = ifelse(sig  == "sig_Both", as.character(gene), '')),size = 3,alpha=0.7,label.size=NA, label.padding = 0.01,show.legend = F)
    p2 + geom_abline(intercept = 0, slope=1, color="darkgrey", linetype="dashed")+ theme(text=element_text(size=10)) +
    xlab("logFoldChange (Mild vs Post)") + ylab("logFoldChange (Severe vs Post)")
    
############    
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)

mild.up <- as.character(dp$gene[dp$Pval1 < 0.01 & dp$log2FC_1 > 1])
de1 <- bitr(unique(mild.up), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
nt1_kk <- enrichKEGG(gene = de1$ENTREZID,
                     organism     = "hsa",
                     pvalueCutoff = 0.1)
                     
both <- as.character(dp$gene[dp$sig == "sig_Both" & dp$log2FC_1 > 0 & dp$log2FC_2 > 0])
de1 <- bitr(unique(both), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
nt1_kk <- enrichKEGG(gene = de1$ENTREZID,
                     organism     = "hsa",
                     pvalueCutoff = 0.1)                     
                

m1 <- enrichGO(gene = de1$ENTREZID,
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",
               ont = "ALL",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.1,
               qvalueCutoff = 0.05,
               readable      = TRUE)

nt1_kk <- enrichKEGG(gene = de1$ENTREZID,
                     organism     = "hsa",
                     pvalueCutoff = 0.1)
                     
                     
                     
                     
dp <- read.table("./DiffPeak/match2DEG.cMono.txt",sep="\t",header = T)   


dp <- read.table("./DiffPeak/match2DEG.ncMono.txt",sep="\t",header = T)  
dp$sig = "not_sig"
dp$sig[dp$Pval1 < 0.05] = "sig_Mild"
dp$sig[dp$Pval2 < 0.05] = "sig_Severe"
dp$sig[dp$Pval1 < 0.05 & dp$Pval2 < 0.05] = "sig_Both"
dp$sig[dp$sig == "sig_Both" & abs(dp$log2FC_1) > 2] = "vsig_Both"
dp$sig[dp$sig == "sig_Both" & abs(dp$log2FC_2) > 2] = "vsig_Both"
#dp.sig <- dp[dp$sig == "sig_Both",]
dp.sig <- dp[dp$sig != "not_sig",]

dp.sig$RNA_direct <- factor(dp.sig$RNA_direct, levels = c("unknown","Up","Down"))
               
dp.sig$plot = "not"
dp.sig$plot[dp.sig$RNA_direct != "unknown" & dp.sig$sig == "sig_Both" ] = "yes"             
               
               
p <- ggplot(dp.sig, aes(x=log2FC_1, y=log2FC_2)) + geom_point(aes(color=RNA_direct), size=1, show.legend = T) +  
          geom_point(data = subset(dp.sig, RNA_direct != 'unknown'), aes(color=RNA_direct), size=1, show.legend = F) 
p2 <- p +  geom_hline(aes(yintercept=0),linetype="dashed") +  scale_color_manual(values=c("lightgrey","red","blue"))+
    geom_vline(aes(xintercept=0),linetype="dashed") + 
    theme_bw() + #xlim(c(-2,2))+ ylim(c(-2,2)) +
    geom_label_repel(aes(label = ifelse(plot  == "yes", as.character(gene), '')),size = 3,alpha=0.7,label.size=NA, label.padding = 0.01,show.legend = F)
    p2 + geom_abline(intercept = 0, slope=1, color="darkgrey", linetype="dashed") + theme(text=element_text(size=10)) +
    xlab("log2FC (Mild vs Post)") + ylab("log2FC (Severe vs Post)")
    
      