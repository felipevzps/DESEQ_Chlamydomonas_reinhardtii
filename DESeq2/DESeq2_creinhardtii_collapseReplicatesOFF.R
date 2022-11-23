if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("tximport")
BiocManager::install("DESeq2")
BiocManager::install("impute")
BiocManager::install("preprocessCore")
BiocManager::install("GO.db")

install.packages("jsonlite")
install.packages("WGCNA")
install.packages("pheatmap")

library(DESeq2)
library(tximport)
library(jsonlite)
library(pheatmap)

# Reset environment variables
rm(list=ls())

setwd("/home/felipevzps/Documentos/DESEQ_Chlamydomonas_reinhardtii/DESeq2")

targets<-read.csv("target.csv",header=TRUE)
rownames(targets)<-targets$SampleName
targets

files <- paste("/home/felipevzps/Documentos/DESEQ_Chlamydomonas_reinhardtii/DESeq2/", targets$SampleName, "/quant.sf",sep="")
tx2gene<-read.delim("transcript2gene.txt",header=FALSE)
names(files) <- targets$SampleName
all(file.exists(files))

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut=FALSE)
head(txi.salmon$abundance)

Chla_TPM<-txi.salmon$abundance

####
max(txi.salmon$abundance)
min(txi.salmon$abundance)
mean(txi.salmon$abundance)
####

dds <- DESeqDataSetFromTximport(txi.salmon, colData = targets, design = ~Condition)
dds$Condition <-relevel(dds$Condition, ref='CONTROL')

dds <- DESeq(dds)
dds

vsd <- vst(dds, blind=FALSE)
vsd
#############imagem pra apresentacao -- ignorar codigos

dds <- DESeqDataSetFromTximport(txi.salmon, colData = targets, design = ~Condition)
dds$Condition <-relevel(dds$Condition, ref='CONTROL')

dds <- DESeq(dds)

vsd <- vst(dds, blind=FALSE)

plotPCA(vsd, intgroup=c("Condition","Time"),ntop=500)


#############imagem pra apresentacao -- ignorar codigos

#Check distribution of samples in a PCA, showing the top 500 varying genes
plotPCA(vsd, intgroup=c("Time"),ntop=500)
plotPCA(vsd, intgroup=c("Condition"),ntop=500)
plotPCA(vsd, intgroup=c("Condition","Time"),ntop=500)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Condition, vsd$Time, sep="-")
colnames(sampleDistMatrix) <- paste(vsd$Condition, vsd$Time, sep="-")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)

res_SALT_vs_CONTROL <- results(dds, lfcThreshold=1, altHypothesis="greaterAbs", contrast = c('Condition','CONTROL','SALT'), alpha=0.05)

summary(res_SALT_vs_CONTROL)

sig_SALT_vs_CONTROL<-res_SALT_vs_CONTROL[which(res_SALT_vs_CONTROL$padj<0.05),]
summary(sig_SALT_vs_CONTROL)
sig_SALT_vs_CONTROL
dim(sig_SALT_vs_CONTROL)
head(sig_SALT_vs_CONTROL,20)
#write.csv(sig_SALT_vs_CONTROL, "/home/felipevzps/Documentos/DESEQ_Chlamydomonas_reinhardtii/DESeq2/sig_SALT_vs_CONTROL.csv")

# Get the significant up-regulated genes
up = subset(sig_SALT_vs_CONTROL, log2FoldChange > 0)
up

# Get the significant down-regulated genes
down = subset(sig_SALT_vs_CONTROL, log2FoldChange < 0)
down

# Order genes by log2FC
up_ordered = up[order(up$log2FoldChange, decreasing = TRUE),]
down_ordered = down[order(down$log2FoldChange, decreasing = TRUE),]
up_ordered
#write.csv(up_ordered, "/home/felipevzps/Documentos/DESEQ_Chlamydomonas_reinhardtii/DESeq2/sig_up_ordered_SALT_vs_CONTROL.csv")
down_ordered
#write.csv(down_ordered, "/home/felipevzps/Documentos/DESEQ_Chlamydomonas_reinhardtii/DESeq2/sig_down_ordered_SALT_vs_CONTROL.csv")

########3
up_ordered<-read.csv("/home/felipevzps/Documentos/DESEQ_Chlamydomonas_reinhardtii/DESeq2/Suplementar/Tabela1/sig_up_ordered_converted_SALT_vs_CONTROL.csv", header = TRUE)
down_ordered<-read.csv("/home/felipevzps/Documentos/DESEQ_Chlamydomonas_reinhardtii/DESeq2/Suplementar/Tabela1/sig_down_ordered_converted_SALT_vs_CONTROL.csv", header = TRUE)

rownames(up_ordered)<-up_ordered$GeneleName
rownames(down_ordered)<-down_ordered$GeneName

up_ordered$GeneName
down_ordered$GeneName

gp_up_ordered = gost(row.names(up_ordered), organism = "creinhardtii",  significant = FALSE ,ordered_query = FALSE)
gp_down_ordered = gost(row.names(down_ordered), organism = "creinhardtii",  significant = FALSE ,ordered_query = FALSE)
head(gp_up_ordered$result, 20)
head(gp_down_ordered$result, 20)

#Visualization of GO and KEGG Enrichment Analysis
gostplot(gp_up_ordered, interactive = TRUE)
gostplot(gp_down_ordered, interactive = TRUE)

########3


#Exploring DGE results
drawLines <- function() abline(h=c(-1,1),col="red",lwd=2)
plotMA(res_SALT_vs_CONTROL); drawLines()

TPM_SALT_vs_CONTROL<-Chla_TPM[which(rownames(Chla_TPM) %in% rownames(sig_SALT_vs_CONTROL)),
                              as.character(targets[which(targets$Condition %in% c('SALT','CONTROL')),'SampleName'])]

annot_col= data.frame(Time=factor(targets$Time),
                      Condition=factor(targets$Condition))
annot_col
rownames(annot_col)<-targets$SampleName
pheatmap(TPM_SALT_vs_CONTROL, scale='row', annotation_col = annot_col)

