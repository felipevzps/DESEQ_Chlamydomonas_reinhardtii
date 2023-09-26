# topGO Chlamydomonas reinhardtii 24h

library(topGO)

setwd("/home/felipe/Documents/DESEQ_Chlamydomonas_reinhardtii/GO_enrichment/NOISeq_24h")

dados <- read.table("DEGs_NOISeqTech_SALT24H.mod.mod.csv", header = FALSE, sep = ",")
dados

go_annotations <- read.table("GO_annotations.txt", header = FALSE, sep = "\t")
go_annotations

go_annotations_separados <- strsplit(go_annotations[, 2], ",")
go_annotations_separados

novo_go_annotations <- data.frame(gene = rep(go_annotations[, 1], sapply(go_annotations_separados, length)),
                                  GO = unlist(go_annotations_separados))
novo_go_annotations

genes_ranking <- dados[, c(1, ncol(dados))]
genes_ranking
genes_ranking[, 2] <- as.numeric(genes_ranking[, 2])
genes_ranking[, 2]

allGenes <- setNames(genes_ranking[, 2], genes_ranking[, 1])
allGenes

limiar <- 100

topGOdata <- new("topGOdata", ontology="BP", allGenes=allGenes, 
                 annot = annFUN.gene2GO, gene2GO = novo_go_annotations,
                 geneSelectionFun = function(x) x >= limiar)

