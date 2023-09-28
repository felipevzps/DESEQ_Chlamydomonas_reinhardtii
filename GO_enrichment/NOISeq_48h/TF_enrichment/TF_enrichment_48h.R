# topGO Chlamydomonas reinhardtii 48h (enrichment of TFs)

rm(list=ls())

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("topGO")

library(topGO)

# Notebook
setwd("/home/felipe/Documents/DESEQ_Chlamydomonas_reinhardtii/GO_enrichment/NOISeq_48h/TF_enrichment")

# CENA
#setwd("/home/felipevzps/Documentos/DESEQ_Chlamydomonas_reinhardtii/GO_enrichment/NOISeq_24h")

dados <- read.table("DEGs_NOISeq_CRE_v5.6.hmmsearch_PlnTFDVPfam32.families.48H_scores.txt", header = FALSE, sep = ",")
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

# Create a map of geneIDs to GO terms
ann.genes <- genesInTerm(topGOdata)
str(ann.genes)

fishers_result <- runTest(topGOdata, algorithm = "classic", statistic = "fisher")
#			 -- Classic Algorithm -- 
#the algorithm is scoring 58 nontrivial nodes
#parameters: 
#  test statistic: fisher

allRes <- GenTable(topGOdata, classic = fishers_result, ranksOf = "classic", topNodes=58)
allRes

topGO_all_table <- allRes[order(allRes$classic),]
topGO_all_table

write.table(topGO_all_table, file = paste0("topGO_TFs_Chlamydomonas_reinhardtii_48h", ".BP", ".csv"), sep =",", quote= F, col.names = T)

# Create a map of geneIDs/GO terms on the GO terms from the Fisher analysis
fisher.go <- names(score(fishers_result))
fisher.ann.genes <- genesInTerm(topGOdata, whichGO = fisher.go)
fisher.ann.genes

# Use lapply para formatar os elementos da lista como linhas de texto
formatted_lines <- lapply(names(fisher.ann.genes), function(go_id) {
  genes <- fisher.ann.genes[[go_id]]
  return(paste(go_id, paste(genes, collapse = ", "), sep = "\t"))
})
formatted_lines

# Escreva todas as linhas formatadas em um arquivo de texto
write.table(formatted_lines, file = paste0("topGO_Chlamydomonas_reinhardtii_24h_ann.genes.BP.csv"), quote= F, col.names = F, row.names = F, sep="\n")
