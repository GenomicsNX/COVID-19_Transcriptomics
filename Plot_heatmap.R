
# Author: Komudi Singh
# Contact: komudi.singh@nih.gov

rm(list=ls())
setwd("set wrokign directory")

#read DE result file
DataActi<-read.csv("A549_series16_Mock_Sars_highMOI_ALLgenes_DEresults.csv") 
DataActi<-DataActi[!duplicated(DataActi$external_gene_name),]
count_data<-DataActi[,10:15]
rownames(count_data)<-DataActi$external_gene_name
dge = DGEList(counts = count_data)
dge = calcNormFactors(dge, method="TMM")
tmmCounts<-cpm(dge)
tmmlog<-log(tmmCounts+1,2)
tmmlog2<-as.data.frame(tmmlog)
tmmlog2$external_gene_name<-rownames(tmmlog2)
geneOfInterest<-read.csv("pathwayGenes_series16.csv", header = F)

SARS<-tmmlog2[tmmlog2$external_gene_name%in%geneOfInterest$V1,]
names(SARS)[7]<-"external_gene_name"
SARSinfo<-merge(DataActi[,c(1:9)],SARS, by="external_gene_name")


SARSinfo2<-SARSinfo[order(SARSinfo$logFC, decreasing = T),]
SARSinfo2<-SARSinfo2[,c(10:15)]

######for heatmap
library(gplots)
Data<-as.matrix(SARSinfo2)

pdf("HEATMAP_file.pdf", width=11, hei=12)

par(oma=c(4,0,0,0))
par(oma=c(10,4,4,2))
heatmap.2(Data, margins = c(5,20), main="heatmap", trace='none', scale = 'row', 
          col = colorpanel(300, low='green', mid = 'black', high = 'red'), 
          labRow='', dendrogram='none', Colv=FALSE, Rowv=FALSE,cexRow=0.8)


cData <- t(apply(Data, 1, function(x) x - median(x)))
scale.range <- c(-3,3)
scale.brakes <- seq(scale.range[1], scale.range[2], by = 0.01)
n.brakes <- length(scale.brakes) - 1

dev.off()
