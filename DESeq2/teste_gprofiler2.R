## From https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7859841/

# installing Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2","airway"))

library(DESeq2)
library(airway)
library(gprofiler2)

## Functional enrichment of differentially expressed genes

# load the airway data
data(airway)

# construct the DESeqDataSet object
ddsMat = DESeqDataSetFromMatrix(countData = assay(airway),
                                colData = colData(airway),
                                design = ~ cell + dex)
assay(airway)
colData(airway)
# run DESeq2 pipeline
dds = DESeq(ddsMat)
# get the results
results = results(dds, contrast = c("dex", "trt", "untrt"),
                  alpha = 0.05, lfcThreshold = 1)
# keep only the significant genes
results_sig = subset(results, padj < 0.05)
# get the significant up-regulated genes
up = subset(results_sig, log2FoldChange > 0)
# get the significant down-regulated genes
down = subset(results_sig, log2FoldChange < 0)
# enrichment analysis
gp_up = gost(row.names(up), organism = "hsapiens")
gp_down = gost(row.names(down), organism = "hsapiens") 

head(gp_up$result)


# order genes by log2FC
up_ordered = up[order(up$log2FoldChange, decreasing = TRUE),]
# ordered enrichment analysis
gp_up_ordered = gost(row.names(up_ordered), organism = "hsapiens",
                     ordered_query = TRUE)
head(gp_up_ordered$result, 8)

## Visualisation of functional enrichment results
gostplot(gp_up, interactive = TRUE)

p1 = gostplot(gp_up, interactive = FALSE)
publish_gostplot(p1, highlight_terms = c("GO:0050896", "KEGG:04978",
                                         "REAC:R-HSA-5661231", "WP:WP3286"))

