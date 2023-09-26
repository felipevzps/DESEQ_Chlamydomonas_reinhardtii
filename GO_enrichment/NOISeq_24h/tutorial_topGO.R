#tutorial topGO (enrichment of GO terms with differentially expressed genes)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("topGO")
BiocManager::install("ALL")
BiocManager::install("hgu95av2.db")
BiocManager::install("Rgraphviz")

install.packages("Rgraphviz")

library(topGO)

# Creating R object of class topGOdata
# The user needs to provide:
# - Gene universe and a criteria for selecting interesting genes (e.g. differentially expressed genes)
# - GO annotations

library(ALL)
data(ALL)
data("geneList")
head(geneList)
length(geneList)
columns(geneList)

affyLib <- paste(annotation(ALL), "db", sep = ".")
library(package = affyLib, character.only = TRUE)

sum(topDiffGenes(geneList))

sampleGOdata <- new("topGOdata", description = "Simple session", ontology = "BP", allGenes = geneList, geneSel = topDiffGenes, nodeSize = 10, annot = annFUN.db, affyLib = affyLib)
sampleGOdata

resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
resultFisher

resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(sampleGOdata, classicFisher = resultFisher, classicKS = resultKS, elimKS = resultKS.elim, orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)
allRes

showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = "all")

# How custom annotations can be used for building a topGOdata object
geneID2GO <- readMappings(file = system.file("examples/geneid2go.map", package = "topGO"))
str(head(geneID2GO))
#gene_ID  GO_ID1, GO_ID2, GO_ID3, GO_ID4, ...

# Predefined list of interesting genes
# Gene universe
geneNames <- names(geneID2GO)
head(geneNames)
