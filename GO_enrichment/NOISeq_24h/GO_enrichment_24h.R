# topGO Chlamydomonas reinhardtii 24h

rm(list=ls())

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("topGO")

library(topGO)

# Notebook
#setwd("/home/felipe/Documents/DESEQ_Chlamydomonas_reinhardtii/GO_enrichment/NOISeq_24h")

# CENA
setwd("/home/felipevzps/Documentos/DESEQ_Chlamydomonas_reinhardtii/GO_enrichment/NOISeq_24h")

dados <- read.table("DEGs_NOISeqTech_SALT24H.mod.mod.csv", header = FALSE, sep = ",")
dados

go_annotations <- readMappings("GO_annotations.txt", sep = "\t", IDsep=",")
go_annotations

genes_ranking <- dados[, c(1, ncol(dados))]
genes_ranking

genes_ranking[, 2] <- as.numeric(genes_ranking[, 2])
genes_ranking[, 2]

allGenes <- setNames(genes_ranking[, 2], genes_ranking[, 1])
allGenes

# Função de seleção de genes simples (selecione todos os genes como significativos)
geneSelectionFun <- function(allGenes) {
  return(rep(TRUE, length(allGenes)))
}

# Crie o objeto topGOdata com o vetor numérico
topGOdata <- new("topGOdata", ontology = "BP", allGenes = allGenes,
                 annot = annFUN.gene2GO, gene2GO = go_annotations,
                 geneSel = geneSelectionFun)

fishers_result <- runTest(topGOdata, algorithm = "classic", statistic = "fisher")
#			 -- Classic Algorithm -- 
#the algorithm is scoring 1496 nontrivial nodes
#parameters: 
#  test statistic: fisher

allRes <- GenTable(topGOdata, classic = fishers_result, ranksOf = "classic", topNodes=1496)
allRes

topGO_all_table <- allRes[order(allRes$classic),]
topGO_all_table

write.table(topGO_all_table, file = paste0("topGO_Chlamydomonas_reinhardtii_24h", ".BP", ".csv"), sep =",", quote= F, col.names = T)