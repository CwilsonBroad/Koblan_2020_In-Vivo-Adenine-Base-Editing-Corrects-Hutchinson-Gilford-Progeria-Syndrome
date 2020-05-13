library(limma)
library(Glimma)
library(edgeR)
library(tidyr)
library(statmod)
library(ggplot2)
library(readr)
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(utility)
library(philentropy)
library(tidyverse)
library(pheatmap)

setwd("/Volumes/Extreme/Luke_Progeria/Human/Differential_Experession/Human /Final_Set/")

#Open RSEM matrix:
rsem_gene_data_full<-read.table("count_matrix.txt",header=TRUE, row.names = 1)

#Trim Matrix:
#rsem_gene_data <- dplyr::select(rsem_gene_data_full, matches("yr."))
rsem_gene_data <- rsem_gene_data_full

##Filter low counts:

# First create the filter:
filter=rowSums(cpm(rsem_gene_data)>=1)>=6

# Second filter the data:
cpm_filtered_matrix=rsem_gene_data[filter,]

#Create a genelist:

d <- cbind(rownames(cpm_filtered_matrix), data.frame(cpm_filtered_matrix, row.names = NULL))

#Limma requires a Digital Expression Matrix (DGE):
DGE<-DGEList(cpm_filtered_matrix)

#Get Gene name list:
DGE$genes<-as.data.frame(d$`rownames(cpm_filtered_matrix)`)

#Take the samplenames from the DGE object.

samplenames <- colnames(DGE)

#Create the design table inside the DGE object
DGE$samples$group <- as.factor(c("Disease", "Disease", "Disease", 
                                  "Disease", "Disease", "Disease", 
                                  "Treated", "Treated", "Treated", 
                                  "Treated", "Treated", "Treated", 
                                  "WT", "WT","WT","WT","WT","WT"))

DGE$samples$individual <- as.factor(c("Patient1", "Patient1", "Patient1",
                                      "Patient2", "Patient2", "Patient2",
                                      "Patient1", "Patient1", "Patient1",
                                      "Patient1", "Patient1", "Patient1",
                                      "Patient3", "Patient3", "Patient3",
                                      "Patient3", "Patient3", "Patient3"))
                                      

group <- as.factor(c("Disease", "Disease", "Disease", 
                     "Disease", "Disease", "Disease", 
                     "Treated", "Treated", "Treated", 
                     "Treated", "Treated", "Treated", 
                     "WT", "WT","WT","WT","WT","WT"))

individual <- as.factor(c("Patient1", "Patient1", "Patient1",
                          "Patient2", "Patient2", "Patient2",
                          "Patient1", "Patient1", "Patient1",
                          "Patient1", "Patient1", "Patient1",
                          "Patient3", "Patient3", "Patient3",
                          "Patient3", "Patient3", "Patient3"))

#Convert raw counts to cpm and logcpm
cpm <- cpm(DGE)
lcpm <- cpm(DGE, log=TRUE)

#get the mean and median library size.
L <- mean(DGE$samples$lib.size) * 1e-6
M <- median(DGE$samples$lib.size) * 1e-6
c(L, M)

#Plot the coverage:
lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(DGE)
col <- brewer.pal(nsamples, "Paired")
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

##Normalizing gene expression distributions:

#Calculate the normalization factors using TMM:
DGE <- calcNormFactors(DGE, method = "TMM")
DGE$samples$norm.factors

DGE2 <- DGE
DGE2$samples$norm.factors <- 1


#Create plot comparing samples:
#Unormalized samples
par(mfrow=c(1,2))
lcpm <- cpm(DGE2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Unnormalised data",ylab="Log-cpm")
DGE2 <- calcNormFactors(DGE2)  
DGE2$samples$norm.factors

#Normalized samples
lcpm <- cpm(DGE2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Normalised data",ylab="Log-cpm")

#Plot MDS to compare samples:
lcpm <- cpm(DGE, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.individual <- individual
levels(col.individual) <-  brewer.pal(nlevels(col.individual), "Set2")
col.individual <- as.character(col.individual)
title(main="A. Sample groups")
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")
plotMDS(lcpm, labels=individual, col=col.individual, dim=c(3,4))
title(main="B. individual")

#You can create an active MDS plot with Glimma:
glMDSPlot(lcpm, labels=paste(group, individual, sep="_"), 
          groups=DGE$samples[,c(1,4)], launch=TRUE)

##Differential Expression Analysis

#Create design matrix
design <- model.matrix(~0+group+individual)
colnames(design) <- gsub("group", "", colnames(design))
design

#Create contrasts:
contr.matrix <- makeContrasts(
  WTvsDisease = WT - Disease,
  TreatedtvsDisease = Treated - Disease,
  WTvsTreated = WT - Treated,
  levels = colnames(design))
contr.matrix

#Removing heteroscedascity from count data:
par(mfrow=c(1,2))
v <- voom(DGE, design, plot=TRUE)
v

#fit linear model
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

#Examine the # of differentially expressed genes in the empirical bayes fit:
summary(decideTests(efit))

#Use the treat method can be used to calculate p-values from empirical Bayes
#moderated t-statistics with a minimum log-FC requirement. 

tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

#Genes that are DE in multiple comaprisons can be extracted using the results from
#decide tests where 0's represent that are not 1's upregulated -1 downregulated. 

de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)

head(tfit$de.common, n=20)
vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))
vennDiagram(dt[,1:3], circle.col=c("turquoise", "salmon"))

#Get the top genes in each condition:
WT.vs.Disease <- topTreat(tfit, coef=1, n=Inf)
Treated.vs.Disease <- topTreat(tfit, coef=2, n=Inf)
WT.vs.Treated <- topTreat(tfit, coef=3, n=Inf)
head(Treated.vs.Disease)

#Get the top tables
WTvsDisease_all_genes<-topTable(efit, adjust="BH",coef="WTvsDisease", p.value=1, number=Inf ,resort.by="P")
TreatedvsDisease_all_genes<-topTable(efit, adjust="BH",coef="TreatedtvsDisease", p.value=1, number=Inf ,resort.by="P")
WTvsTreated_all_genes<-topTable(efit, adjust="BH",coef="TreatedtvsDisease", p.value=1, number=Inf ,resort.by="P")

write.csv(WT.vsDisease_all_genes, file = "MA_Plot/WTvsDiseaseAllGene.csv")
write.csv(TreatedvsDisease_all_genes, file = "MA_Plot/TreatedvdDiseaseAllGene.csv")
write.csv(WTvsTreated_all_genes, file = "MA_Plot/WTvsTreatedAllGene.csv")

#apply threshold:
all_genes$threshold <- as.factor(all_genes$adj.P.Val < 0.0000000000000000000000000000000000005 | abs(all_genes$logFC) > 2)


#Mean difference plot:
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[2], 
       xlim=c(-8,13))

glMDPlot(tfit, coef=2, status=dt, main=colnames(tfit)[2],
         side.main="genes", counts=lcpm, groups=group, launch=TRUE)

#Heatmap
Treated.vs.Disease.topgenes <- Treated.vs.Disease$d..rownames.cpm_filtered_matrix..[1:300]
i <- which(v$genes[,1] %in% Treated.vs.Disease.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,], scale="row",
          labRow=v$genes$SYMBOL[i], labCol=group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")

WT.vs.Disease.topgenes <- WT.vs.Disease$d..rownames.cpm_filtered_matrix..[1:300]
i <- which(v$genes[,1] %in% WT.vs.Disease.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,], scale="row",
          labRow=v$genes$SYMBOL[i], labCol=group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column",
          keysize=1)
               

#Write CSV files for normalized counts:
write.csv(DGE2$counts, file = "WTvsDisease_counts.csv")

#Camera and Barcode Plot:
load("human_c5_v5p2.rdata")

#Jennsen-Shannon Divergence plot:

##Full plot:
Jcounts <- t(DGE2$counts) 
JSDMatrix <- JSD(Jcounts)

heatmap(JSDMatrix)

##Top 500 plot:

geneList <- c(as.character(WT.vs.Disease.topgenes))

#Select for that genelist:
Normalizedcounts <- read_csv("CountPlots/WTvsDisease_counts.csv")

selected_Counts <-  Normalizedcounts %>%
  filter(Normalizedcounts$X1 %in% geneList)

Jcounts <- t(selected_Counts)
Jcounts <- as.data.frame(Jcounts)
Jcounts <- dplyr::slice(Jcounts, 2:18)
Jcounts <- data.matrix(Jcounts)
JSDMatrix <- JSD(Jcounts)


pheatmap(JSDMatrix, cutree_rows = 17)

#Camera 
idx <- ids2indices(Hs.c5,id=v$genes$ENTREZID)
cam.Control.vs.Disease <- camera(v,idx,design,contrast=contr.matrix[,1])
head(cam.Control.vs.Disease,5)  
