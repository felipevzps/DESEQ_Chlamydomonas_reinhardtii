# topGO Chlamydomonas reinhardtii 48h

rm(list=ls())

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("topGO")

library(topGO)

# Notebook
setwd("/home/felipe/Documents/DESEQ_Chlamydomonas_reinhardtii/GO_enrichment/NOISeq_48h")

# CENA
#setwd("/home/felipevzps/Documentos/DESEQ_Chlamydomonas_reinhardtii/GO_enrichment/NOISeq_48h")

dados <- read.table("DEGs_NOISeqTech_SALT48H.mod.mod.csv", header = FALSE, sep = ",")
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
#the algorithm is scoring 836 nontrivial nodes
#parameters: 
#  test statistic: fisher

allRes <- GenTable(topGOdata, classic = fishers_result, ranksOf = "classic", topNodes=836)
allRes

topGO_all_table <- allRes[order(allRes$classic),]
topGO_all_table

write.table(topGO_all_table, file = paste0("topGO_Chlamydomonas_reinhardtii_48h", ".BP", ".csv"), sep =",", quote= F, col.names = T)

library(ggplot2)

ntop <- 30
ggdata <- topGO_all_table[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
ggdata

# Todos valores da coluna classic são 1
# Adicionar 0.001 para o log10(1) ser diferente de 0
ggdata$classic <- as.numeric(ggdata$classic) + 0.001
ggdata

ggplot(ggdata,
       aes(x = Term, y = -log10(classic), size = -log10(classic), fill = -log10(classic))) +
  
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +
  
  xlab('') + ylab('Enrichment score') +
  labs(
    title = 'GO Analysis',
    #subtitle = 'Top 50 terms ordered by Kolmogorov-Smirnov p-value',
    subtitle = 'Top 30 terms ordered by Fisher Exact p-value',
    caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001') +
  
  geom_hline(yintercept = c(-log10(0.01), -log10(0.001), -log10(0.0001)),
             linetype = c("dotted", "longdash", "solid"),
             colour = c("black", "black", "black"),
             size = c(0.5, 1.5, 3)) +
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")) +
  
  coord_flip()

ggplot2::ggsave("GOTerms_P30_Cortex_Fisher.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
