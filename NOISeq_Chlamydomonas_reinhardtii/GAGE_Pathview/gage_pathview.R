setwd('/home/felipe/Documents/DESEQ_Chlamydomonas_reinhardtii/NOISeq_Chlamydomonas_reinhardtii/GAGE_Pathview')


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("gage")
BiocManager::install("gageData")

library(gage)
#library(help=gage)

#Example: Basic Analysis - https://bioconductor.org/packages/2.13/bioc/vignettes/gage/inst/doc/gage.pdf

#Matriz com count
data(gse16873)
cn=colnames(gse16873)
hn=grep('HN', cn, ignore.case = T)
adh=grep('ADH', cn, ignore.case = T)
dcis=grep('DCIS', cn, ignore.case = T)

#Posicao das samples na Matriz com count
print(dcis)

data(kegg.gs)
print(kegg.gs)
data(go.gs)
lapply(kegg.gs[1:3], head)

#UHUU! Tem anotacao da chlamy
#kegg.gsets("cre")

kg.hsa=kegg.gsets("hsa")
kegg.gs=kg.hsa$kg.sets[kg.hsa$sigmet.idx]
#save(kegg.gs, file="kegg.sigmet.gsets.RData")

gse16873.kegg.p <- gage(gse16873, gsets = kegg.gs,
                        ref = hn, samp=dcis)

go.gs #only yhe first 1000 entries
gse16873.go.p <- gage(gse16873, gsets = go.gs,
                      ref = hn, samp=dcis)

str(gse16873.kegg.p, strict.width = 'wrap')

head(gse16873.kegg.p$greater[, 1:5],4)
head(gse16873.kegg.p$less[, 1:5],4)
head(gse16873.kegg.p$stats[, 1:5],4)

#capture pathways perturbed towards both directions - only for KEGGs
gse16873.kegg.2d.p <- gage(gse16873, gsets = kegg.gs,
                        ref = hn, samp=dcis, same.dir = F)

head(gse16873.kegg.2d.p$greater[,1:5],4)
head(gse16873.kegg.2d.p$stats[,1:5],4)

###
# precisamos:
#   - gse16873: Matriz com counts de todas samples
#   - kegg.gs: Anotação dos genes (gene ID, desc, K)

head(gse16873)
head(kegg.gs)
