options(future.globals.maxSize = 10737418240)  # 10240*1024^2 = 10737418240  for 10GB

library(Seurat)
library(dplyr)
library(Matrix)
library(Matrix.utils)
library(ggplot2)
library(reshape2)
library(pryr)
library(RColorBrewer)
library(tidyr)
library(cowplot)


setwd("~/Desktop/covid19/MHH50")

covid <- readRDS("pbmc.pools.harmony.rds")
###plots

FeaturePlot(pbmc,features=c("CD3D","CCR7","IL7R","CD8A","GNLY",
"NKG7","GZMB","CD14","LYZ","FCGR3A",
"CD163","CD83","CD68","CD74","HLA-DRA","PPBP"),
order = F)

mks <- c("IL7R","CD3E","CD8A","GNLY","NKG7","GZMB","CD14","LYZ","FCGR3A","CST3","CD163",
         "CD83","CD68","IFI6","IFITM3","HLA-DRA","CD74","CD79A","CD27","SDC1","HBA1","PPBP","SOX4")
genes <- intersect(mks,row.names(covid@assays$RNA@data))

violin.plot.data <- get.violin.data1(covid, genes)
ggplot(violin.plot.data, aes(factor(ident),value)) +
  geom_violin(scale="width",adjust=1,trim=TRUE,aes(fill=factor(ident)),show.legend = F) + 
  ylab("") + xlab("") +
  coord_flip() +facet_wrap(~ gene,scales = "free_x", ncol = length(levels(violin.plot.data$gene))) +
  theme(strip.text.x = element_text(size=12, angle=-90),
        axis.text.y = element_text(size=15),
        strip.background = element_blank(),
        panel.spacing.x = unit(c(-0.2), "lines"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0.5,0.4,0.5,1.3), "cm")) +
  scale_x_discrete(limits = rev(levels(violin.plot.data$ident)), position = "top")
  
  
  tt.markers <- readRDS("mhh.cluster.rds")
  
#tt.markers <- FindAllMarkers(object = pbmc.pools, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)

#tt.markers <- tt.markers[grep("^RPS|^RPL",tt.markers$gene,invert = T),]
top10 <- tt.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

#top10$cluster <- factor(top10$cluster, levels=c(0,3,5,16,6,13,15,7,14,2,1,8,4,12,17,9,10,18,11))
top10 <- top10[order(top10$cluster),]
#pdf("marker_plots.newid2.pdf",width=12,height=28)
DotPlot(covid,features = rev(unique(top10$gene)),group.by="seurat_clusters",cols="RdBu")+coord_flip()+xlab("")+ylab("")
#dev.off()
  
  
  
  
  
Idents(covid)<-"seurat_clusters"
covid <- RenameIdents(covid,
                      `0` = "0",
                      `5` = "1",
                      `12`= "2",
                      `13`= "3",
                      
                      `10`= "4",
                      `8` = "5",
                      `9` = "6",
                      
                      `17`= "7" , 
                                   
                      `2` = "8",
                      `15`= "9",
                      `1` = "10",
                      `3` = "11",
                      `4` = "12",
                      `11`= "13",
                      
                      `6` = "14" ,     
                      `16`= "15",
                      
                      `14`= "16",
                      `7` = "17",
                      `18`= "18"
)
covid$new.id <- Idents(covid)

Idents(covid)<-"RNA_snn_res.0.4"
covid <- RenameIdents(covid,
                      `0` = "0: cMono",
                      `5` = "1: cMono",
                      `12`= "2: cMono",
                      `13`= "3: CD163+ cMono",
                      
                      `10`= "4: CD163+ cMono",
                      `8` = "5: ncMono",
                      `9` = "6: mDC",
                      
                      `17`= "7: pDC" , 
                                   
                      `2` = "8: CD4 T",
                      `15`= "9: CD4 T",
                      `1` = "10: CD8 T",
                      `3` = "11: CD4 T",
                      `4` = "12: NK",
                      `11`= "13: NK",
                      
                      `6` = "14: B" ,     
                      `16`= "15: Plasmablast",
                      
                      `14`= "16: Erythrocytes",
                      `7` = "17: Megakaryocyte",
                      `18`= "18: undefined"
)
covid$celltype.idL0 <- Idents(covid)

Idents(covid)<-"seurat_clusters"
covid <- RenameIdents(covid,
                      `0` = "cMono",
                      `5` = "cMono",
                      `12`= "cMono",
                      `13`= "CD163.cMono",
                      
                      `10`= "CD163.cMono",
                      `8` = "ncMono",
                      `9` = "mDC",
                      
                      `17`= "pDC" , 
                                
                      `2` = "CD4.T",
                      `15`= "CD4.T",
                      `1` = "CD8.T",
                      `3` = "CD4.T",
                      `4` = "NK",
                      `11`= "NK",
                      
                      `6` = "B" ,     
                      `16`= "Plasmablast",
                      
                      `14`= "Erythrocytes",
                      `7` = "Megakaryocyte",
                      `18`= "undefined"
)
covid$celltypeL0 <- Idents(covid)
cell.essential <- c("cMono","CD163.cMono","ncMono","mDC","pDC","CD4.T","CD8.T","NK","B","Plasmablast")
covid.sub <- subset(covid, celltypeL0 %in% cell.essential)

Idents(pbmc)<-"celltypeL0"

pbmc <- RenameIdents(pbmc,
`cMono`="merged.cMonos",
`CD163.cMono`="merged.cMonos")
pbmc$celltypeL1<- Idents(pbmc)
pbmc$celltypeL1 <- factor(pbmc$celltypeL1, levels = c("merged.cMonos", "ncMono","mDC" ,"CD4.T", "CD8.T", "NK", "pDC","B","Plasmablast" ))


pbmc <- RenameIdents(pbmc,
`cMono`="cMonos",
`CD163.cMono`="cMonos",
`Plasmablast`= "B cells",
`CD4.T`="T cells",
`CD8.T`="T cells"
)
pbmc$celltypeL1<- Idents(pbmc)
pbmc$celltypeL1 <- factor(pbmc$celltypeL1, levels = c("cMonos", "ncMono","T cells", "NK", "mDC", "pDC","B cells"))

Idents(pbmc)<-"Severity"
pbmc <- RenameIdents(pbmc,
`severe`="active",
`mild`="active",
`post`="convalescent")
pbmc$Disease <- Idents(pbmc)


#saveRDS(covid, "pbmc.harmony.annotev3.rds")
#saveRDS(covid.sub, "pbmc.filter.annotev3.rds")

umap.colors.id<-c(
  "#66c2a4",  #0 cmono
  "#41ae76",  #1
  "#006d2c",  #2
  "#8c96c6",  #3   cd163 
  "#88419d",  #4
  "#e7298a",  #5  ncmono
  "#3288bd",  #6  mdc
  "#2166ac",  #7  pdc
  
  "#d6604d",  #8  Ts
  "#f4a582",  #9
  "#fddbc7",  #10
  "#b2182b",  #11
  
  "#92c5de",  #12 NK
  "#67001f",  #13 cycT
  
  "#bf812d",  #14
  "#543005",  #15
  
  "#fee391",  #16
  "#737373",  #17
  "#bdbdbd" 
)


umap.colors.ct <-c(
"cMono"= "#66c2a4",  #1
"CD163.cMono"=  "#8c96c6",  #4
"ncMono"=  "#e7298a",
"mDC"=  "#3288bd",  #5
"pDC"=  "#2166ac",  #6
"CD4.T"=  "#FDDBC7",  #7
"NK"=  "#92c5de",  #9
"CD8.T"=  "#67001f",  #8

"B"=  "#bf812d",  #13
"Plasmablast"=  "#543005",  #14

"Erythrocytes"=  "#fee391",  #16
"Megakaryocyte" =     "#737373",  #15
"undefined"=  "#bdbdbd"  #17
)

color.severity <- c("Severe"="#9e0142","Mild"="#fee08b","Post"="#4393c3")


covid.sub$Severity2 <- factor(covid.sub$Severity2, levels = c("Severe","Mild","Post"))


pdf("umap_celltypeL0.pdf",width=6,height=5)
DimPlot(covid.sub, group.by="celltypeL0", cols= umap.colors.ct)
dev.off()

pdf("umap_celltypeL0_split.pdf",width=15,height=5)
DimPlot(covid.sub, group.by="celltypeL0", split.by="Severity",cols= umap.colors.ct)
dev.off()

pdf("umap_celltypeL0_split.pdf",width=10,height=5)
DimPlot(covid, group.by="celltypeL0", split.by="Disease",cols= umap.colors.ct)
dev.off()


pdf("umap_celltypeidL0.pdf",width=6,height=5)
DimPlot(covid.sub, group.by="celltype.idL0", cols= umap.colors.id)
dev.off()

pdf("umap_celltypeidL0_split.pdf",width=15,height=5)
DimPlot(covid.sub, group.by="celltype.idL0", split.by="Severity2",cols= umap.colors.id)
dev.off()

unique <- read.table("./unique_RNA.txt",stringsAsFactors=F)[,1]
covid.uniq <- subset(covid.sub, SampleID %in% unique)
saveRDS(covid.uniq, "pbmc.uniq.annotev3.rds")



Idents(covid) <- "new.id"
mks <- c("CD14","LYZ","FCGR3A","CST3","CD163",
         "CD83","CD68","IFI6","IFITM3","HLA-DRA","CD74","IL7R","CD3E","CD8A","GNLY","NKG7","GZMB","CD79A","CD27","SDC1","HBA1","PPBP","SOX4")
genes <- intersect(mks,row.names(covid@assays$RNA@data))

violin.plot.data <- get.violin.data1(covid, rev(genes))
ggplot(violin.plot.data, aes(factor(ident),value)) +
  geom_violin(scale="width",adjust=1,trim=TRUE,aes(fill=factor(ident)),show.legend = F) + 
  ylab("") + xlab("") +
  coord_flip() +facet_wrap(~ gene,scales = "free_x", ncol = length(levels(violin.plot.data$gene))) +
  theme(strip.text.x = element_text(size=12, angle=-90),
        axis.text.y = element_text(size=15),
        strip.background = element_blank(),
        panel.spacing.x = unit(c(-0.2), "lines"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0.5,0.4,0.5,1.3), "cm")) +
  scale_x_discrete(limits = rev(levels(violin.plot.data$ident)), position = "top")+
  scale_fill_manual(values= umap.colors.id)
  
  
tt.markers <- readRDS("mhh.cluster.rds")
tt.markers <- FindAllMarkers(object = covid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
tt.markers <- tt.markers[grep("^RPS|^RPL",tt.markers$gene,invert = T),]
top10 <- tt.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

top10$cluster <- factor(top10$cluster, levels=c(0,5,12,13,10,8,9,17,2,15,1,3,4,11,6,16,14,7,18))
top10 <- top10[order(top10$cluster),]

DotPlot(covid,features = rev(unique(top10$gene)),group.by="seurat_clusters",cols="RdBu")+coord_flip()+xlab("")+ylab("")
#dev.off()


meta <- read.table("MHH_smaple_ID.txt",sep="\t",header=T,stringsAsFactors = F)


row.names(meta) <- paste(meta$pool.id, meta$soup.id,sep="_")

#loadingID SampleID gender patient genoID pool.id soup.id Age
aa <- meta$loadingID
names(aa) <- row.names(meta)
Idents(covid) = "sample_assign"
covid <- RenameIdents(covid,aa)
covid$loadingID <- Idents(covid)

aa <- meta$gender
names(aa) <- row.names(meta)
Idents(covid) = "sample_assign"
covid <- RenameIdents(covid,aa)
covid$gender <- Idents(covid)


aa <- meta$patient
names(aa) <- row.names(meta)
Idents(covid) = "sample_assign"
covid <- RenameIdents(covid,aa)
covid$patient <- Idents(covid)

aa <- meta$Age
names(aa) <- row.names(meta)
Idents(covid) = "sample_assign"
covid <- RenameIdents(covid,aa)
covid$Age <- Idents(covid)


aa <- meta$Severity
names(aa) <- row.names(meta)
Idents(covid) = "sample_assign"
covid <- RenameIdents(covid,aa)
covid$Severity <- Idents(covid)

aa <- meta$Condition
names(aa) <- row.names(meta)
Idents(covid) = "sample_assign"
covid <- RenameIdents(covid,aa)
covid$Condition <- Idents(covid)

aa <- meta$SampleID
names(aa) <- row.names(meta)
Idents(covid) = "sample_assign"
covid <- RenameIdents(covid,aa)
covid$SampleID <- Idents(covid)


covid$Severity2 <- as.character(covid$Severity)
covid$Severity2[covid$sample_assign == "P1b_1"]="Mild"
covid$Severity2[covid$sample_assign == "P3b_1"]="Mild"
covid$Severity2[covid$sample_assign == "P1c_3"]="Mild"
covid$Severity2[covid$Severity2 == "light_or_post"]="Post"



library(DirichletReg)
library(ggsignif)

cells <- c(
"cMono",  #1
"CD163.cMono",  #4
"ncMono",
"mDC",  #5
"pDC",  #6
"CD4.T",  #7
"CD8.T",  #8
"NK",  #9
"B",  #13
"Plasmablast"  #14
)

covid <- subset(covid.uniq , celltypeL0 %in% cells)
covid$donor.cond <- paste(covid$SampleID, covid$Severity2, sep=".")

#Idents(covid) <- "celltype.idL0"

Idents(covid) <- "celltypeL0"

dirichlet_regression = function(counts, covariates){
# Calculate regression
counts = as.data.frame(counts)
counts$counts = DR_data(counts)
data = cbind(counts, covariates)
fit = DirichReg(counts ~ condition, data)
# Get p-values
u = summary(fit)
pvals = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 4]
v = names(pvals)
pvals = matrix(pvals, ncol=length(u$varnames))
rownames(pvals) = gsub('condition', '', v[1:nrow(pvals)])
colnames(pvals) = u$varnames
fit$pvals = pvals
fit
}

######
covid.sub <- subset(covid,Severity2 %in% c("mild","post"))
covid.sub$donor.cond <- as.factor(as.character(covid.sub$donor.cond))
counts = do.call(cbind, tapply(covid.sub@meta.data$donor.cond, covid.sub@active.ident, table))
covariates = data.frame(condition=gsub('[0-9a-c]+.', '', rownames(counts)))

fit <- dirichlet_regression(counts,covariates)
dat.test <- data.frame("Comp1"="mild","Comp2"="post","P"=t(fit$pvals))
names(dat.test)<-c("Comp1","Comp2","P")
dat.test$fdr <- p.adjust(dat.test$P,method="fdr")
dat.test$sig <- "ns"
dat.test$sig[dat.test$fdr < 0.05]="*"
dat.test$sig[dat.test$fdr < 0.01]="**"
dat.test$sig[dat.test$fdr < 0.001]="***"
dat.test$cell_type = row.names(dat.test)
dat.test.1 <- dat.test

######
covid.sub <- subset(covid,Severity %in% c("severe","post"))
covid.sub$donor.cond <- as.factor(as.character(covid.sub$donor.cond))
counts = do.call(cbind, tapply(covid.sub@meta.data$donor.cond, covid.sub@active.ident, table))
covariates = data.frame(condition=gsub('[0-9a-c]+.', '', rownames(counts)))

fit <- dirichlet_regression(counts,covariates)
dat.test <- data.frame("Comp1"="severe","Comp2"="post","P"=t(fit$pvals))
names(dat.test)<-c("Comp1","Comp2","P")
dat.test$fdr <- p.adjust(dat.test$P,method="fdr")
dat.test$sig <- "ns"
dat.test$sig[dat.test$fdr < 0.05]="*"
dat.test$sig[dat.test$fdr < 0.01]="**"
dat.test$sig[dat.test$fdr < 0.001]="***"
dat.test$cell_type = row.names(dat.test)
dat.test.2 <- dat.test

######
covid.sub <- subset(covid,Severity %in% c("severe","mild"))
covid.sub$donor.cond <- as.factor(as.character(covid.sub$donor.cond))
counts = do.call(cbind, tapply(covid.sub@meta.data$donor.cond, covid.sub@active.ident, table))
covariates = data.frame(condition=gsub('[0-9a-c]+.', '', rownames(counts)))

fit <- dirichlet_regression(counts,covariates)
dat.test <- data.frame("Comp1"="severe","Comp2"="mild","P"=t(fit$pvals))
names(dat.test)<-c("Comp1","Comp2","P")
dat.test$fdr <- p.adjust(dat.test$P,method="fdr")
dat.test$sig <- "ns"
dat.test$sig[dat.test$fdr < 0.05]="*"
dat.test$sig[dat.test$fdr < 0.01]="**"
dat.test$sig[dat.test$fdr < 0.001]="***"
dat.test$cell_type = row.names(dat.test)
dat.test.3 <- dat.test

######
dat.test <- rbind(dat.test.1, dat.test.2, dat.test.3)

#Idents(covid) <- "celltype.idL0"
Idents(covid) <- "celltypeL0"
tmp<-data.frame(cell_type= covid@active.ident,sample= covid@meta.data$SampleID,donor= covid@meta.data$Severity2)
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

y2$other <- 1-y2$f

y2$Condition <- factor(y2$Condition,levels=c("severe","mild","post"))

#dat.sig <- dat.test
dat.sig <- dat.test[dat.test$sig != "ns",]
maxs <- NULL
for(i in 1:nrow(dat.sig)){
	c <- dat.sig[i,"cell_type"]
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



pdf("cpp_merge_v4.pdf",width=15,height=4.5)
p <- ggplot(data=y2, aes(x=Condition,y=f))+geom_boxplot(aes(fill=Condition),outlier.shape=NA)
p2 <- p + facet_grid(.~cell_type)+geom_jitter(width = 0.2)+ theme_classic()+scale_fill_manual(values=c("#9e0142","#fee08b","#4393c3"))+xlab("")+ylab("percent")
p2 + geom_signif(data=dat.sig, aes(xmin=Comp1,xmax=Comp2, annotations=sig, y_position=maxs), manual = TRUE, vjust=0.8, tip_length=0.005)+theme(axis.text.x=element_text(angle = 70,hjust = 1))
#filename <- paste("BoxPlot_cohort2_cluster",i,".pdf",sep="")
dev.off()




#########################

Idents(covid)<-"Severity"
covid <- RenameIdents(covid,
`severe`="Active",
`mild`="Active",
`post`="Convalescent")
covid$Disease <- Idents(covid)

covid$donor.cond <- paste(covid$SampleID, covid$Disease, sep=".")

Idents(covid) <- "celltypeL0"

dirichlet_regression = function(counts, covariates){
counts = as.data.frame(counts)
counts$counts = DR_data(counts)
data = cbind(counts, covariates)
fit = DirichReg(counts ~ condition, data)
# Get p-values
u = summary(fit)
pvals = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 4]
v = names(pvals)
pvals = matrix(pvals, ncol=length(u$varnames))
rownames(pvals) = gsub('condition', '', v[1:nrow(pvals)])
colnames(pvals) = u$varnames
fit$pvals = pvals
fit
}

covid$donor.cond <- as.factor(as.character(covid$donor.cond))
counts = do.call(cbind, tapply(covid@meta.data$donor.cond, covid@active.ident, table))
covariates = data.frame(condition=gsub('[0-9a-c]+.', '', rownames(counts)))

fit <- dirichlet_regression(counts,covariates)
dat.test <- data.frame("Comp1"="Active","Comp2"="Convalescent","P"=t(fit$pvals))
names(dat.test)<-c("Comp1","Comp2","P")
dat.test$fdr <- p.adjust(dat.test$P,method="fdr")
dat.test$sig <- "ns"
dat.test$sig[dat.test$fdr < 0.05]="*"
dat.test$sig[dat.test$fdr < 0.01]="**"
dat.test$sig[dat.test$fdr < 0.001]="***"
dat.test$cell_type = row.names(dat.test)
dat.test.1 <- dat.test


Idents(covid) <- "celltypeL0"
tmp<-data.frame(cell_type= covid@active.ident,sample= covid@meta.data$SampleID,donor= covid@meta.data$Disease)
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
y2$other <- 1-y2$f
y2$Condition <- factor(y2$Condition,levels=c("Active","Convalescent"))


dat.sig <- dat.test[dat.test$sig != "ns",]
maxs <- NULL
for(i in 1:nrow(dat.sig)){
	c <- dat.sig[i,"cell_type"]
	if(dat.sig[i, "Comp1"]=="Active" & dat.sig[i,"Comp2"]=="Convalescent"){
	        maxs <- c(maxs,max(y2$f[y2$cell_type == c])+0.03)
	}
}
dat.sig$maxs <- maxs


pdf("cpp_2cond_v5.pdf",width=15,height=4.5)
p <- ggplot(data=y2, aes(x=Condition,y=f))+geom_boxplot(aes(fill=Condition),outlier.shape=NA)
p2 <- p + facet_grid(.~cell_type)+geom_jitter(width = 0.2)+ theme_classic()+scale_fill_manual(values=c("#E35C2E","#4393c3"))+xlab("")+ylab("percent")
p2 + geom_signif(data=dat.sig, aes(xmin=Comp1,xmax=Comp2, annotations=sig, y_position=maxs), manual = TRUE, vjust=0.8, tip_length=0.005)+theme(axis.text.x=element_text(angle = 70,hjust = 1))

dev.off()








dat.test <- NULL
for(i in levels(as.factor(y2$cell_type))){
#for(i in c("PLA1","PLA2","PLA3")){  
  print(i)
  y2.i <- y2[y2$cell_type==i,]
  y2.i$Smp <- DR_data(y2.i[,5:6])
  res <- DirichReg(Smp~Condition,y2.i,model="alternative",base=2)
  x <- summary(res)
  sc <- x$coef.mat["ConditionConvalescent",3:4]

  dat.test <- rbind(dat.test, c(i,sc,"Active","Convalescent"))
}

dat.test<- as.data.frame(dat.test)
names(dat.test) <- c("cell_type","Z","Pr","Comp1","Comp2")
dat.test$Pr <- as.numeric(as.character(dat.test$Pr))
dat.test$Z <- as.numeric(as.character(dat.test$Z))
dat.test$sig <- "ns"
dat.test$sig[dat.test$Pr < 0.05]="*"
dat.test$sig[dat.test$Pr < 0.01]="**"
dat.test$sig[dat.test$Pr < 0.001]="***"