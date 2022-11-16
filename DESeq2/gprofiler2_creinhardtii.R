library(DESeq2)
library(airway)
library(gprofiler2)
library(tximport)
library(pheatmap)

setwd("/home/felipevzps/Documentos/DESEQ_Chlamydomonas_reinhardtii/DESeq2")
targets<-read.csv("target.csv",header=TRUE)
rownames(targets)<-targets$SampleName

files <- paste("/home/felipevzps/Documentos/DESEQ_Chlamydomonas_reinhardtii/DESeq2/", targets$SampleName, "/quant.sf",sep="")
tx2gene<-read.delim("transcript2gene.txt",header=FALSE)
names(files) <- targets$SampleName
all(file.exists(files))

#txi.salmon <- tximport(files, type = "salmon", txOut= TRUE, countsFromAbundance = 'lengthScaledTPM')
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut=FALSE)
head(txi.salmon$abundance)

####### construct the DESeqDataSet object
ddsMat = DESeqDataSetFromMatrix(countData = round(txi.salmon$abundance),
                                colData = targets,
                                design = ~ Condition)
# run DESeq2 pipeline
dds = DESeq(ddsMat)
# get the results
results = results(dds, contrast = c("Condition", "CONTROL", "SALT"),
                  alpha = 0.05, lfcThreshold = 1)
# keep only the significant genes
results_sig = subset(results, padj < 0.05)
# get the significant up-regulated genes
up = subset(results_sig, log2FoldChange > 0)
# get the significant down-regulated genes
down = subset(results_sig, log2FoldChange < 0)
# enrichment analysis
gp_up = gost(row.names(up), organism = "creinhardtii", significant = FALSE)
gp_down = gost(row.names(down), organism = "creinhardtii") 



head(gp_up$result)