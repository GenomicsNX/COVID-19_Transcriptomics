
# Author: Komudi Singh
# Contact: komudi.singh@nih.gov

rm(list=ls())

setwd("set working directory")
 
#read DE result file
DataActi<-read.csv("A549_series16_Mock_Sars_highMOI_ALLgenes_DEresults.csv")
library(ggplot2)
library(ggrepel)
table(geneOfInterest$V1%in%DataActi$external_gene_name)

###choose relevant columns
data1<-DataActi[,c(1,4,7,8)] 
head(data1)

####set pvalue and logfc thresholds
data1$pval_threshold = as.factor(data1$adj.P.Val > 0.05)
data1$Increased<-as.factor(data1$adj.P.Val<0.05 & data1$logFC>0)
data1$Decreased<-as.factor(data1$adj.P.Val<0.05 & data1$logFC<0)

###lable up, down, and non significant genes
data1$cumul<-"none"
data1$cumul[data1$pval_threshold=='TRUE']<-"nonsignificant"
data1$cumul[data1$Increased=='TRUE']<-"UP"
data1$cumul[data1$Decreased=='TRUE']<-"DOWN"

ans<-order(data1$cumul, decreasing = F)
data1<-data1[c(ans),]
head(data1)

table(data1$cumul)
cols <- c("blue","lightgrey","red")
p=ggplot(data=data1, 
         aes(x=logFC, y =-log10(adj.P.Val))) +
  scale_color_manual(values=cols) +
  geom_point(alpha=0.8, size=1.75, aes(col=c(cumul))) +
  xlim(c(-3.5, 5.5)) +
  xlab("log2 fold change") + ylab("-log10 adjusted p value") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black")+
  theme_classic() 

pdf('Volcano_plot.pdf', width=11, height=11)
plot(p) 
dev.off()
