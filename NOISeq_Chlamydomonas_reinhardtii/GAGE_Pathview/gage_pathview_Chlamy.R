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
<<<<<<< HEAD
##matrix_tpm=read.table('gseChlamy_counts.salmon_quantmerge_modID_head_and_uniq_matrizquadrada.txt', header = F, sep = "\t")

#working matriz quadrada
matrix_tpm=read.table('gseChlamy_counts.salmon_quantmerge_modID_matrizQuadrada_uniqrows.txt', header = T, row.names = 1,sep = "\t")

#working matriz nao quadrada -> mais pra frente precisa usar o compare = unpaired na linha do gseCre.kegg.p
matrix_tpm=read.table('gseChlamy_counts.salmon_quantmerge_modID_matrizNAOQuadrada_uniqrows.txt', header = T, row.names = 1,sep = "\t")
matrix_tpm

### transformando matriz de entrada: removendo genes duplicados e transformando coluna GENE_NAMES nos indices das linhas
matrix_tpm=read.table('gseChlamy_counts.salmon_quantmerge_modID.txt', header = T, sep = "\t")

# Remover linhas duplicadas com base na coluna "GENE_NAMES"
matrix_tpm <- matrix_tpm[!duplicated(matrix_tpm$GENE_NAMES), ]
matrix_tpm

# Definir a coluna "GENE_NAMES" como os índices das linhas
rownames(matrix_tpm) <- matrix_tpm$GENE_NAMES
matrix_tpm$GENE_NAMES <- NULL  # Remover a coluna "GENE_NAMES" dos dados
matrix_tpm


###

#Verificar se todos valores na matriz são numéricos
all_numeric <- all(sapply(matrix_tpm, is.numeric))
if (all_numeric) {
  print("Todos os valores na matriz matrix_tpm são numéricos.")
} else {
  print("A matriz matrix_tpm contém valores não numéricos.")
}

### Existem valores nao numericos: Investigando isso

# Função para verificar se um valor é numérico
is_numeric <- function(x) {
  is.numeric(x) && !is.na(x)
}

# Encontre as posições dos valores não numéricos na matriz
non_numeric_positions <- which(!sapply(matrix_tpm, is_numeric), arr.ind = TRUE)


# Verifique se há valores não numéricos
if (length(non_numeric_positions) > 0) {
  print("Posições dos valores não numéricos na matriz:")
  print(non_numeric_positions)
 } else {
  print("Não foram encontrados valores não numéricos na matriz.")
}

###

=======
matrix_tpm=read.table('gseChlamy_counts.salmon_quantmerge_modID.txt', header=TRUE)
>>>>>>> d2af2509868d2a699c03c66c55ac94bddc9ea98f
cn=colnames(matrix_tpm)
hn=grep('CONTROLE', cn, ignore.case = T)
#cname=grep('GENE_NAMES', cn, ignore.case = T)
dcis=grep('SAL', cn, ignore.case = T)

#Posicao das samples na Matriz com count
print(cn)
print(hn)
#print(cname)
print(dcis)

#data(kegg.gs)
#head(kegg.gs)
#data(go.gs)

<<<<<<< HEAD
=======
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



>>>>>>> d2af2509868d2a699c03c66c55ac94bddc9ea98f
#UHUU! Tem anotacao da chlamy
kegg.gsets("cre")

kg.hsa=kegg.gsets("cre")

kegg.gs=kg.hsa$kg.sets[kg.hsa$sigmet.idx]

#save(kegg.gs, file="kegg.sigmet.gsets.RData")

print(kegg.gs)
<<<<<<< HEAD

#apenas para matriz quadrada, onde experimentos tem 2 controles e 2 condições experimentais (2vs2)
##gseCre.kegg.p <- gage(matrix_tpm, gsets = kegg.gs, ref = hn, samp=dcis)

##if the samples are not one-on-one paired in the experiment design, we take samples as unpaired
gseCre.kegg.p <- gage(matrix_tpm, gsets = kegg.gs, ref = hn, samp=dcis, compare = "unpaired")
=======
gseCre.kegg.p <- gage(matrix_tpm, gsets = kegg.gs,ref = hn, samp=dcis)


# go.gs #only yhe first 1000 entries
gse16873.go.p <- gage(gse16873, gsets = go.gs,
                      ref = hn, samp=dcis)
>>>>>>> d2af2509868d2a699c03c66c55ac94bddc9ea98f

# Trabalhando com termos GO
##go.gs #only yhe first 1000 entries
##gse16873.go.p <- gage(gse16873, gsets = go.gs,
##                      ref = hn, samp=dcis)

str(gseCre.kegg.p, strict.width = 'wrap')

head(gseCre.kegg.p$greater[, 1:5],400)
head(gseCre.kegg.p$less[, 1:5],400)
head(gseCre.kegg.p$stats[, 1:5],400)

#capture pathways perturbed towards both directions - only for KEGGs
gseCre.kegg.p.kegg.2d.p <- gage(matrix_tpm, gsets = kegg.gs,
                        ref = hn, samp=dcis, same.dir = F, compare = "unpaired")

head(gseCre.kegg.p.kegg.2d.p$greater[,1:5],4)
head(gseCre.kegg.p.kegg.2d.p$stats[,1:5],4)

#Heatmaps -> arquivo gseCre.kegg.teste.gs.heatmap.pdf
gseCre.kegg.sig <- sigGeneSet(gseCre.kegg.p, outname = "gseCre.kegg.teste") 

#tutorial pagina 11 -- continuar estudando.. nao sei direito o que isso faz
gseCre.kegg.esg.up <- esset.grp(gseCre.kegg.p$greater, matrix_tpm, gsets = kegg.gs, ref=hn, samp=dcis, 
                                test4up = T, output = F, make.plot = T, compare = "unpaired")
gseCre.kegg.esg.up
