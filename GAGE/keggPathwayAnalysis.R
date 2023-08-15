# Generally Applicable Gene-set/Pathway Analysis Pipeline
# Example: Basic Analysis - https://bioconductor.org/packages/2.13/bioc/vignettes/gage/inst/doc/gage.pdf

rm(list=ls())

setwd('/home/felipe/Documents/DESEQ_Chlamydomonas_reinhardtii/GAGE')

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("gage")
library(gage)

matrix_tpm=read.table('cre_quantmerge_mod_ID.txt', header = T, sep = "\t")
# Remover linhas duplicadas com base na coluna "Name"
matrix_tpm <- matrix_tpm[!duplicated(matrix_tpm$Name), ]
matrix_tpm

# Definir a coluna "GENE_NAMES" como os índices das linhas
rownames(matrix_tpm) <- matrix_tpm$Name
matrix_tpm$Name <- NULL  # Remover a coluna "GENE_NAMES" dos dados
matrix_tpm

cn=colnames(matrix_tpm)
hn=grep('CONTROLE', cn, ignore.case = T)
dcis=grep('SAL', cn, ignore.case = T)

# Posicao das samples na Matriz com count
print(cn)
print(hn) # Controle
print(dcis) # Sal

# Load gene set data for Chlamydomonas reinhardtii
kg.cre=kegg.gsets("cre")
kegg.gs=kg.cre$kg.sets[kg.cre$sigmet.idx]
#save(kegg.gs, file="kegg.sigmet.gsets.RData")
print(kegg.gs)

# If the samples are not one-on-one paired in the experiment design, we take samples as unpaired
gseCre.kegg.p <- gage(matrix_tpm, gsets = kegg.gs, ref = hn, samp=dcis, compare = "unpaired")

str(gseCre.kegg.p, strict.width = 'wrap')

# Significant gene sets based on p-value cutoffs
head(gseCre.kegg.p$greater[, 1:5],4)
head(gseCre.kegg.p$less[, 1:5],4)
head(gseCre.kegg.p$stats[, 1:5],400)

# Capture pathways perturbed towards both directions - only for KEGGs
gseCre.kegg.p.kegg.2d.p <- gage(matrix_tpm, gsets = kegg.gs,ref = hn, samp=dcis, same.dir = F, compare = "unpaired")

head(gseCre.kegg.p.kegg.2d.p$greater[,1:5],4)
head(gseCre.kegg.p.kegg.2d.p$stats[,1:5],4)

# Heatmaps -> arquivo gseCre.kegg.teste.gs.heatmap.pdf
gseCre.kegg.sig <- sigGeneSet(gseCre.kegg.p, outname = "gseCre.kegg.sig") 

#tutorial pagina 11 -- continuar estudando.. nao sei direito o que isso faz
gseCre.kegg.esg.up <- esset.grp(gseCre.kegg.p$greater, matrix_tpm, gsets = kegg.gs, ref=hn, samp=dcis,
                                test4up = T, output = F, make.plot = T, compare = "unpaired")
gseCre.kegg.esg.up