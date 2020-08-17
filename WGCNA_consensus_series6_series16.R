
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
names(CountData)
phenoData<-read.delim("SampleInfo_fromGSE147507_formergeV3.txt", header = F)

PhenoSplit<-colsplit(phenoData$V2,"_",c("series","treatment"))
phenoData2<-cbind(phenoData,PhenoSplit)

##select low viral load and high viral load datasets
pheno_25<-subset(phenoData,phenoData2$series=="Series6" | phenoData2$series=="Series16") 

####format pheno file
pheno_25$V2<-gsub("-","\\.",pheno_25$V2)

###select count matrix
Count_comp_25<-CountData[,names(CountData)%in%pheno_25$V2]
Count_comp_25<-Count_comp_25[,c(1:3,7:9,4:6,10:12)]
names(Count_comp_25)

#####filter low expression genes, cut off count > 5 in all samples
Count_comp_25$obo5<-rowSums(Count_comp_25[,c(1:12)]>5)
Count_comp_25_final<-subset(Count_comp_25,Count_comp_25$obo5>10)

###run batch corrections using combat
library(sva)

my_DGEList <- DGEList(counts=Count_comp_25_final[,-c(13)])
ExpDesign<-c(rep("cont",6),rep("Infection",6))
batch<-c(rep("series6",3),rep("series16",3),rep("series6",3),rep("series16",3))
my_batch<-as.factor(batch)
my_mod = model.matrix(~as.factor(ExpDesign))
All_count_norm <- calcNormFactors(my_DGEList)
my_data = cpm(All_count_norm, log=TRUE, prior.count=1)
combat <- ComBat(dat=my_data, batch=my_batch, mod=my_mod)
head(combat)


library(WGCNA)
library(dplyr)
library(dynamicTreeCut)
library(fastcluster)
library(stringi)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()


####RUN WGCNA
#########################################################################################################

##Data formatting for WGCNA

#########################################################################################################
#write.table(combat, "Corrected_cpms.txt", row.names = T, quote = F)

All.data=combat
#All.data<-All.data[rownames(All.data)%in%rownames(Count_comp_25_final),]
names(All.data)
dim(All.data)
grp1Data<-data.frame(All.data[,c(1:6)])
grp2Data<-data.frame(All.data[,c(7:12)])


annotID<-as.data.frame(as.character(rownames(Count_comp_25_final)))
names(annotID)<-"external_gene_name"

library(biomaRt)####get gene names for all ensembl ids
mart <- useMart( "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "useast.ensembl.org" )


Gene_list <- getBM(filters= "external_gene_name", 
                   attributes= c("ensembl_gene_id",
                                 "external_gene_name",
                                 "entrezgene_id"),values=annotID$external_gene_name,mart= mart)


finalannot<-merge(annotID,Gene_list,by="external_gene_name")
finalannot<-finalannot[,-c(3)]
write.csv(finalannot,"gene_annot.csv", row.names = F)



grp1.Data=as.data.frame(t(grp1Data))
head(grp1.Data)[1:5,1:5]
grp2.Data=as.data.frame(t(grp2Data))

# We work with two sets:
nSets = 2;
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("grp1", "grp2")
shortLabels = c("grp1", "grp2")
# Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = vector(mode = "list", length = nSets)

multiExpr[[1]] = list(data = grp1.Data);
names(multiExpr[[1]]$data) = rownames(grp1Data);
rownames(multiExpr[[1]]$data) = names(grp1Data);
multiExpr[[2]] = list(data = grp2.Data);
names(multiExpr[[2]]$data) = rownames(grp2Data);
rownames(multiExpr[[2]]$data) = names(grp2Data);
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)

###check if the list have gene names:
DataExpr1=multiExpr[[2]]
Expr1DataFrame=as.data.frame(DataExpr1["data"])
names(Expr1DataFrame)[1:5]

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK


#####save data input
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;
save(multiExpr, nGenes, nSamples, setLabels, shortLabels, exprSize, 
     file = "Consensus-dataInput_series6_series16.RData");


#############################################################################################

####NETWORK CONSTRUCTION use the data input file from previos step

#############################################################################################
rm(list = ls())
options(stringsAsFactors = FALSE);
allowWGCNAThreads() 
setwd("set working directory")
lnames = load(file = "Consensus-dataInput_series6_series16.RData");
lnames

#=====================================================================================
#
#  NETWORK CONSTRUCTION Code 1
#
#=====================================================================================
nSets=length(multiExpr)
# Choose a set of soft-thresholding powers
powers = c(seq(4,12,by=1), seq(12,30, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]]);
collectGarbage();
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
pdf(file = "scaleFreeAnalysis_series6_series16_A549.pdf", wi = 12, he = 6)

par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off();
#=====================================================================================
#
#  NETWORK CONSTRUCTION Code 2
#
#=====================================================================================
net = blockwiseConsensusModules(maxBlockSize = 15000,
                                multiExpr, power = 11, minModuleSize = 30, deepSplit = 2,
                                pamRespectsDendro = FALSE, 
                                mergeCutHeight = 0.25, numericLabels = TRUE,
                                minKMEtoStay = 0, 
                                saveTOMs = FALSE, saveTOMFileBase = "none",verbose = 5)
#=====================================================================================
#
#  NETWORK CONSTRUCTION Code 3
#
#=====================================================================================
consMEs = net$multiMEs;
moduleLabels = net$colors;
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
consTree = net$dendrograms[[1]]; 

##ksedits--------------------------------------------
#good to have the module numbers and the number of genes present in each module information store.
#note: columns number 0 refers to the genes that were unassigned to any module
moduleList=as.data.frame(table(moduleLabels))
moduleColorList=as.data.frame(table(moduleColors))
moduleList=as.data.frame(table(moduleLabels))
colorList=as.data.frame(table(moduleColors))
ModuleListColor=cbind(moduleList,colorList)
write.csv(ModuleListColor, file="ModuleNumbers_moduleSize_series6_series16_p11_A549.csv")
#=====================================================================================
#
#  NETWORK CONSTRUCTION Code 4
#
#=====================================================================================
sizeGrWindow(12,10);
pdf(file = "ConsensusDendrogram-auto_series6_series16_p11_A549.pdf", wi = 12, he = 10)
plotDendroAndColors(consTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")

dev.off()
#=====================================================================================
#
#  NETWORK CONSTRUCTION Code 4 Save!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#=====================================================================================
consMEs = net$multiMEs;
moduleLabels = net$colors;
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
consTree = net$dendrograms[[1]]; 


save(consMEs, moduleLabels, moduleColors, consTree, file = "Consensus-NetworkConstruction-auto_series6_series16_p11_A549.RData")


#############################################################################################

####EXPORT NETWORK

#############################################################################################
rm(list=ls())
setwd("set working directory")
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
#######################################1_4 comparison###################################################
lnames = load(file = "Consensus-dataInput_series6_series16.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Consensus-NetworkConstruction-auto_series6_series16_p11_A549.RData");
lnames


#=====================================================================================
#  export NETWORK  Code 3
#=====================================================================================

# We work with two sets:
nSets = 2;
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("grp1", "grp2")
shortLabels = c("grp1", "grp2")
softPower = 11;
# Initialize an appropriate array to hold the adjacencies
adjacencies = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate adjacencies in each individual data set
for (set in 1:nSets)
  adjacencies[set, , ] = abs(cor(multiExpr[[set]]$data, use = "p"))^softPower;
# Recalculate topological overlap if needed
TOM=array(0,dim=c(nSets,nGenes,nGenes));
#calculate TOM for individual data set
for (set in 1:nSets)
  TOM[set,,] = TOMsimilarityFromExpr(adjacencies[set,,], power = softPower);
consensusTOM = pmin(TOM[1, , ], TOM[2, , ]);
save(consensusTOM, nGenes, setLabels, shortLabels, nSets, 
     file = "consensusTOM_series6_series16_A549.RData");
# Read in the annotation file
annotMain = read.csv("Mock_A549_log2.csv"); 
annotMain<-as.data.frame(annotMain$external_gene_name)
#table(annot$gene_id==All.data$gene_id)
names(annotMain)<-"external_gene_name"

annot<-read.csv("gene_annot.csv")

annot<-annot[,c(1,2)]
head(annot)


moduleList<-read.csv("ModuleNumbers_moduleSize_series6_series16_p11_A549.csv") ##total 67 modules
moduleList1<-subset(moduleList,moduleList$Freq.1>50) ##47 modules>50
intModules<-moduleList1$moduleColors
# Select module probes
#first get one of the datasest from multiExpr###ks edits-------------------
DataExpr1=multiExpr[[2]]
Expr1DataFrame=as.data.frame(DataExpr1["data"])
names(Expr1DataFrame)[1:5]
Preprobes=names(Expr1DataFrame)
probes = substring(Preprobes,6)
head(probes)[1:5]
dir.create("series6_series16_consensus_network")
modules='blue'
for (modules in intModules)
{
  inModule = is.finite(match(moduleColors, modules));
  modProbes = probes[inModule];
  #ks edits---------------------------------------------------
  annotMod=match(modProbes,annot$external_gene_name)
  annotMod2=annot[annotMod,]
  modGenes=annotMod2[,2]
  #modGenes = annot$gene_name[match(modProbes, annot$gene_id)];
  #ks edits ends---------------------------------------------
  # Select the corresponding Topological Overlap
  modTOM = consensusTOM[inModule, inModule];
  #modTOMtop=modTOM[top, top];
  #dimnames(modTOM) = list(modProbes, modProbes)
  # Export the network into edge and node list files Cytoscape can read
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("series6_series16_consensus_network/Consensus_Cytoscape-edges_", paste(modules, collapse="-"), "_series6_series16_A549.txt", sep=""),
                                 nodeFile = paste("series6_series16_consensus_network/Consensus_Cytoscape-nodes_", paste(modules, collapse="-"), "_series6_series16_A549.txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0,
                                 nodeNames = modProbes,
                                 altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule]);
  
}

