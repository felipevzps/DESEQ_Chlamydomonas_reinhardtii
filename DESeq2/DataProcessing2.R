if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("tximport")
BiocManager::install("DESeq2")
BiocManager::install("impute")
BiocManager::install("preprocessCore")
BiocManager::install("GO.db")

#install.packages("jsonlite")
#install.packages("WGCNA")
#install.packages("pheatmap")


library(DESeq2)
library(tximport)
library(jsonlite)
library(pheatmap)
#library(WGCNA)

rm(list=ls())
setwd("/home/felipevzps/Documentos/DESEQ_Chlamydomonas_reinhardtii/DESeq2")
#targets<-read.csv("target.csv",header=TRUE)

##testanto gambiarra
targets<-read.csv("target_gambiarra.csv",header=TRUE)

rownames(targets)<-targets$SampleName
targets

files <- paste("/home/felipevzps/Documentos/DESEQ_Chlamydomonas_reinhardtii/DESeq2/", targets$SampleName, "/quant.sf",sep="")
tx2gene<-read.delim("transcript2gene.txt",header=FALSE)
names(files) <- targets$SampleName
all(file.exists(files))

#txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut=FALSE)
txi.salmon <- tximport(files, type = "salmon", txOut= TRUE, countsFromAbundance = 'lengthScaledTPM')
head(txi.salmon$abundance)

Arabtpm<-txi.salmon$abundance

############# FUNCIONANDO
##dds <- DESeqDataSetFromTximport(txi.salmon, colData = targets, design = ~Condition)
##dds$Condition <-relevel(dds$Condition, ref='SALT')
#############

#Testando novos tratamentos -> TEMPO
targets$Time <- factor(targets$Time, levels = c("24","48"))
dds <- DESeqDataSetFromTximport(txi.salmon, colData = targets, design = ~Time)
dds$Condition <-relevel(dds$Time, ref='24')

head(dds)
dds <- DESeq(dds)

vsd <- vst(dds, blind=FALSE)
vsd

#Check distribution of samples in a PCA, showing the top 500 varying genes
plotPCA(vsd, intgroup=c("Time"),ntop=500)
plotPCA(vsd, intgroup=c("Condition"),ntop=500)
plotPCA(vsd, intgroup=c("Condition","Time"),ntop=500)

#Show clustering of samples, should be in agreement with PCA
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Condition, vsd$Time, sep="-")
colnames(sampleDistMatrix) <- paste(vsd$Condition, vsd$Time, sep="-")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)
#Getting differentially expressed genes, comparisons against control (MOCK)
res_SALT_vs_CONTROL <- results(dds, lfcThreshold=1.5, altHypothesis="greaterAbs", contrast = c('Condition','CONTROL','SALT'), alpha=0.05)
summary(res_SALT_vs_CONTROL)

sig_SALT_vs_CONTROL<-res_SALT_vs_CONTROL[which(res_SALT_vs_CONTROL$padj<0.05),]
dim(sig_SALT_vs_CONTROL)
head(sig_SALT_vs_CONTROL,200)

####
res_24_vs_48 <- results(dds, lfcThreshold=1.5, altHypothesis="greaterAbs", contrast = c('Time','24','48'), alpha=0.05)
summary(res_24_vs_48)

sig_24_vs_48<-res_24_vs_48[which(res_24_vs_48$padj<0.05),]
dim(sig_24_vs_48)
head(sig_24_vs_48,200)
####



#sig_ABA_vs_MOCK<-res_ABA_vs_MOCK[which(res_ABA_vs_MOCK$padj<0.05),]
#dim(sig_ABA_vs_MOCK)
#head(sig_ABA_vs_MOCK)

#sig_DEHYD_vs_MOCK<-res_DEHYD_vs_MOCK[which(res_DEHYD_vs_MOCK$padj<0.05),]
#dim(sig_DEHYD_vs_MOCK)
#head(sig_DEHYD_vs_MOCK)

#Exploring DGE results
drawLines <- function() abline(h=c(-1,1),col="dodgerblue",lwd=2)
#plotMA(res_ABA_vs_MOCK); drawLines()
plotMA(res_SALT_vs_CONTROL); drawLines()
plotMA(res_24_vs_48); drawLines()

#plotMA(res_DEHYD_vs_MOCK); drawLines()

TPM_SALT_vs_CONTROL<-Arabtpm[which(rownames(Arabtpm) %in% rownames(sig_SALT_vs_CONTROL)),
                          as.character(targets[which(targets$Condition %in% c('SALT','CONTROL')),'SampleName'])]

##################
TPM_24_vs_48<-Arabtpm[which(rownames(Arabtpm) %in% rownames(sig_24_vs_48)),
                             as.character(targets[which(targets$Time %in% c('24','48')),'SampleName'])]

TPM_24_vs_48
##PLOT 
annot_col= data.frame(Time=factor(targets$Time),
                      Condition=factor(targets$Condition))
annot_col
rownames(annot_col)<-targets$SampleName

pheatmap(TPM_24_vs_48, scale='row', annotation_col = annot_col)
#################

#TPM_ABA_vs_MOCK<-Arabtpm[which(rownames(Arabtpm) %in% rownames(sig_ABA_vs_MOCK)),
#                          as.character(targets[which(targets$Condition %in% c('ABA','MOCK')),'SampleName'])]

#TPM_DEHYD_vs_MOCK<-Arabtpm[which(rownames(Arabtpm) %in% rownames(sig_DEHYD_vs_MOCK)),
#                          as.character(targets[which(targets$Condition %in% c('DEHYD','MOCK')),'SampleName'])]

annot_col= data.frame(Time=factor(targets$Time),
                      Condition=factor(targets$Condition))
annot_col
rownames(annot_col)<-targets$SampleName
pheatmap(TPM_SALT_vs_CONTROL, scale='row', annotation_col = annot_col)
#pheatmap(TPM_ABA_vs_MOCK, scale='row', annotation_col = annot_col)
#pheatmap(TPM_DEHYD_vs_MOCK, scale='row', annotation_col = annot_col)

