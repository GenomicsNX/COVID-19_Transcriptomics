
# Author: Komudi Singh
# Contact: komudi.singh@nih.gov

rm(list=ls())
library("biomaRt")

####read DGE result
geneExpr<-read.csv("A549_series16_Mock_Sars_highMOI_ALLgenes_DEresults.csv") 

####subset significant genes
sigGenes<-subset(geneExpr,geneExpr$adj.P.Val<0.05) 

####run pathway enrichment
library(clusterProfiler) 

library(org.Hs.eg.db)

ego2 <- enrichGO(gene         = sigGenes$entrezgene_id,
                 OrgDb         = org.Hs.eg.db,
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)
head(ego2)

dim(ego2)

####save files
write.csv(ego2, "ClusterProfiler_A549_Mock_Sars_series16_DEgenes.csv") 

