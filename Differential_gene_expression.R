
# Author: Komudi Singh
# Contact: komudi.singh@nih.gov

rm(list=ls())

setwd("~/Documents/covid_cellLine_RNAseq/GSE147507_redo")
#formating Merged featureCounts Data
CountData<-read.delim(file = "GSE147507_RawReadCounts_Human.tsv", header =T)
head(CountData)
names(CountData)[1]<-"GENE"
CountData<-CountData[!duplicated(CountData$GENE),]
rownames(CountData)<-CountData$GENE
phenoData<-read.delim("SampleInfo_fromGSE147507_formergeV2.txt", header = F)

pheno_5<-subset(phenoData,phenoData$V2=="Series5") #or "Series6" or "Series16
pheno_6<-phenoData[grep("Series6",phenoData$V2),]
####subsetting coundData to phenoFile
Count_comp_16<-CountData[,c(71:76)] ###select appropriate columns
names(Count_comp_16)
isexpr = rowSums(Count_comp_16>1) >= 2
table(isexpr)
Count_comp_16 = Count_comp_16[isexpr,]

####experiment design to label control and infected
ExpDesign<-c(rep("cont",3),rep("Infection",3))
grouptype<-factor(ExpDesign)
designs<-model.matrix(~grouptype)

library(limma)
library(edgeR)

dge=DGEList(countData)
dge=calcNormFactors(dge)
v2 <- voom(dge, designs)
countdataobo2<-v2$E
genecount2<-as.data.frame(countdataobo2)
genecount2$external_gene_name<-rownames(genecount2)

fit2 <- lmFit(v2, designs)
fit2 <- eBayes(fit2)
topTable(fit2, coef=ncol(designs))
DGEResults_voom22<- topTable(fit2, coef=ncol(designs), n=61000)
DGEResults_voom22$X<-rownames(DGEResults_voom22)
Count_comp_16$X<-rownames(Count_comp_16)
library(dplyr)
MergeResults<-merge(DGEResults_voom22,Count_comp_16, by = "X") #####DGE result
