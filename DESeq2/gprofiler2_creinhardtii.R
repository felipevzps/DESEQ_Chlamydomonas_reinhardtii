library(DESeq2)
library(airway)
library(gprofiler2)
library(tximport)
library(pheatmap)

# Reset environment variables
rm(list=ls())

setwd("/home/felipevzps/Documentos/DESEQ_Chlamydomonas_reinhardtii/DESeq2")

# Opening experiment design
targets<-read.csv("target.csv",header=TRUE)
rownames(targets)<-targets$SampleName

# Opening reads quantification
files <- paste("/home/felipevzps/Documentos/DESEQ_Chlamydomonas_reinhardtii/DESeq2/", targets$SampleName, "/quant.sf",sep="")
tx2gene<-read.delim("transcript2gene_atualizado.txt",header=FALSE)
names(files) <- targets$SampleName
all(file.exists(files))

# Import reads quantification
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut=FALSE)
head(txi.salmon$abundance)

# Building the DESeqDataSet object
ddsMat = DESeqDataSetFromMatrix(countData = round(txi.salmon$abundance),
                                colData = targets,
                                design = ~ Condition)
# Run DESeq2 pipeline
dds = DESeq(ddsMat)

# Get the results
results = results(dds, contrast = c("Condition", "CONTROL", "SALT"),
                  alpha = 0.05, lfcThreshold = 1)
results

sig_SALT_vs_CONTROL
# Keep only the significant genes
results_sig = subset(results, padj < 0.05) #0.05
results_sig

# Get the significant up-regulated genes
up = subset(results_sig, log2FoldChange > 0)
up

# Get the significant down-regulated genes
down = subset(results_sig, log2FoldChange < 0)
down

# GO enrichment analysis
#gp_up = gost(row.names(up), organism = "creinhardtii", significant = FALSE)
#gp_down = gost(row.names(down), organism = "creinhardtii", significant = FALSE) 
#head(gp_up$result)
#head(gp_down$result)

# Order genes by log2FC
up_ordered = up[order(up$log2FoldChange, decreasing = TRUE),]
down_ordered = down[order(down$log2FoldChange, decreasing = TRUE),]
up_ordered
down_ordered

# GO ordered enrichment analysis
gp_up_ordered = gost(row.names(up_ordered), organism = "creinhardtii",  significant = TRUE ,ordered_query = TRUE)
gp_down_ordered = gost(row.names(down_ordered), organism = "creinhardtii",  significant = TRUE ,ordered_query = TRUE)
head(gp_up_ordered$result, 20)
head(gp_down_ordered$result, 20)

#Visualization of GO and KEGG Enrichment Analysis
gostplot(gp_up_ordered, interactive = TRUE)
gostplot(gp_down_ordered, interactive = TRUE)

p1 = gostplot(gp_up_ordered, interactive = FALSE)                                                 
p2 = gostplot(gp_down_ordered, interactive = FALSE)

# Beautiful plot with table
publish_gostplot(p1, highlight_terms = c("GO:0071704", "KEGG:00970"))
publish_gostplot(p2, highlight_terms = c("GO:0033573", "GO:1905862"))

