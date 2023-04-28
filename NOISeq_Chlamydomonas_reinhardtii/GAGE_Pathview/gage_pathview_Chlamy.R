setwd('/home/felipevzps/Documentos/DESEQ_Chlamydomonas_reinhardtii/NOISeq_Chlamydomonas_reinhardtii/GAGE_Pathview')


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("gage")
#BiocManager::install("gageData")

library(gage)
#library(help=gage)

rm(list=ls())

#Example: Basic Analysis - https://bioconductor.org/packages/2.13/bioc/vignettes/gage/inst/doc/gage.pdf

#Matriz com count
matrix_tpm=read.table('gseChlamy_counts.salmon_quantmerge_modID.txt', header=TRUE)
cn=colnames(matrix_tpm)
hn=grep('CONTROLE', cn, ignore.case = T)
#adh=grep('ADH', cn, ignore.case = T)
dcis=grep('SAL', cn, ignore.case = T)

#Posicao das samples na Matriz com count
print(cn)
print(hn)
print(dcis)

#data(kegg.gs)
#head(kegg.gs)
#data(go.gs)

kegg.gs=read.csv('gse_onlyKvalues_formatado_todasLinhasPreenchidas.Creinhardtii_281_v5.6.annotation_info.txt', sep="\t", header = FALSE)
kegg.gs
lapply(kegg.gs[1:3], head)
## TEMOS ISSO
#$V1
#[1] "CHLRE_01g000300v5" "CHLRE_01g000600v5" "CHLRE_01g000650v5" "CHLRE_01g002200v5" "CHLRE_01g002300v5" "CHLRE_01g002350v5"

#$V2
#[1] "alpha/beta-Hydrolases superfamily protein"                           "WD-40 repeat family protein / notchless protein, putative"          
#[3] "Copper amine oxidase family protein"                                 "RNA polymerase Rpb6"                                                
#[5] "Cyclophilin-like peptidyl-prolyl cis-trans isomerase family protein" "NAD(P)-bindig Rossmann-fold superfamily protein"                    

#$V3
#[1] "13535" "14855" "276"   "3014"  "1802"  "4532" 

## PRECISA FICAR ASSIM
#$`hsa00010 Glycolysis / Gluconeogenesis`
#[1] "10327" "124"   "125"   "126"   "127"   "128"  

#$`hsa00020 Citrate cycle (TCA cycle)`
#[1] "1431" "1737" "1738" "1743" "2271" "3417"

#$`hsa00030 Pentose phosphate pathway`
#[1] "2203"   "221823" "226"    "229"    "22934"  "230" 



#UHUU! Tem anotacao da chlamy
kegg.gsets("cre")

kg.hsa=kegg.gsets("cre")

kegg.gs=kg.hsa$kg.sets[kg.hsa$sigmet.idx]

#save(kegg.gs, file="kegg.sigmet.gsets.RData")

print(kegg.gs)
gseCre.kegg.p <- gage(matrix_tpm, gsets = kegg.gs,ref = hn, samp=dcis)


# go.gs #only yhe first 1000 entries
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
