
# Author: Komudi Singh
# Contact: komudi.singh@nih.gov

rm(list=ls())
setwd("set working directory")
ClusterFile<-read.csv("ClusterProfiler_A549_Mock_Sars_series16_DEgenes.csv")
Cfile2<-as.data.frame(ClusterFile[,c(4,9)])
Cfile<-Cfile2[order(Cfile2$qvalue,decreasing = F),]
Cfile<-Cfile[1:25,]
Cfile<-Cfile[order(Cfile$qvalue,decreasing = T),]
Cfile$log10<--log10(Cfile$qvalue)
col=Cfile$log10

A=min(col)
B=max(col)

#keeps the same order as in the cluster file
Pathways<-Cfile[,1] 

p<-ggplot(Cfile,aes(x=Pathways, y=log10, fill=log10)) + 
  scale_x_discrete(limits=Pathways)+
  scale_fill_gradient(low="red",high="blue", name="-log10(q value)",limits=c(A,B))+
  geom_bar(width=0.5, stat = "identity")+
  ylab("-log10 q value") +
  theme( axis.text.y  = element_text(face="bold", vjust=0.5, size=40))+
  coord_flip()

pdf("ClusterProfiler_A549_series16_DEgenesv2_DOWNREGULATED.pdf", width=36, hei=16)
par(oma=c(4,0,0,0))
plot(p)
dev.off()
