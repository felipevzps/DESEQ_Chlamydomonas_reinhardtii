library(DESeq2)
library(airway)
library(gprofiler2)
library(tximport)
library(pheatmap)

rm(list=ls())

setwd("/home/felipevzps/Documentos/DESEQ_Chlamydomonas_reinhardtii/DESeq2")
targets<-read.csv("target.csv",header=TRUE)
targets
rownames(targets)<-targets$SampleName

files <- paste("/home/felipevzps/Documentos/DESEQ_Chlamydomonas_reinhardtii/DESeq2/", targets$SampleName, "/quant.sf",sep="")
tx2gene<-read.delim("transcript2gene_atualizado.txt",header=FALSE)
names(files) <- targets$SampleName
all(file.exists(files))

#txi.salmon <- tximport(files, type = "salmon", txOut= TRUE, countsFromAbundance = 'lengthScaledTPM')
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut=FALSE)
head(txi.salmon$abundance)

####### construct the DESeqDataSet object
ddsMat = DESeqDataSetFromMatrix(countData = round(txi.salmon$abundance),
                                colData = targets,
                                design = ~ Condition)