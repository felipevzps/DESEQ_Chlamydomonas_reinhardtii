if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("NOISeq")

library(tximport)
library(NOISeq)

#tutorial: https://rstudio-pubs-static.s3.amazonaws.com/525119_64c1fe6e1a514b89a1ef26d23bf4aae3.html
getwd()
setwd("/home/felipe/Documents/DESEQ_Chlamydomonas_reinhardtii/NOISeq_Chlamydomonas_reinhardtii/Chlamy_24_48")
#PC CENA
#setwd("/home/felipevzps/Documentos/DESEQ_Chlamydomonas_reinhardtii/NOISeq_Chlamydomonas_reinhardtii")

experimentalMatrix<-read.csv("target_without_LOM1.24.1_LOM1.24.2_LOM2.24.1_LOM2.24.3.csv",header=TRUE)
rownames(experimentalMatrix)<-experimentalMatrix$SampleName
experimentalMatrix

#PC CENA
#files <- paste("/home/felipevzps/Documentos/DESEQ_Chlamydomonas_reinhardtii/DESeq2/", experimentalMatrix$SampleName, "/quant.sf",sep="")
files <- paste("/home/felipe/Documents/DESEQ_Chlamydomonas_reinhardtii/DESeq2/", experimentalMatrix$SampleName, "/quant.sf",sep="")
#PC CENA
#tx2gene<-read.delim("transcript2gene.txt",header=FALSE)
tx2gene<-read.delim("/home/felipe/Documents/DESEQ_Chlamydomonas_reinhardtii/DESeq2/transcript2gene.txt",header=FALSE)

names(files) <- experimentalMatrix$SampleName
all(file.exists(files))

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut=FALSE)
head(txi.salmon$abundance)

Chla_TPM<-txi.salmon$abundance
Chla_TPM #mycounts

rownames(experimentalMatrix)
colnames(Chla_TPM)

#Use ?readData to see all options and descriptions of the parameters that you can include.
?readData

#Convert your generic data frame into a NOISeq object
myDataMyFactors <- readData(data=Chla_TPM, factors=experimentalMatrix)
myDataMyFactors

#Data exploration and quality
#JUST CAN BE DONE IF WE HAVE INCLUDED OPTIONAL ARGUMENTS TO NOISEQ OBJECT

myPCA = dat(myDataMyFactors, type = "PCA")
explo.plot(myPCA, factor = "Condition")

#Looking at the NOISeq object we created with our own factors, let's also look at `Run` and `Lane`
#Function to plot multiple figures (1,3) for 3 figures
par(mfrow = c(1, 3))

explo.plot(myPCA, factor = "Condition")
explo.plot(myPCA, factor = "Time")
explo.plot(myPCA, factor = "Sample")

#Create and save a data exploration QC report

#to see in console
#QCreport(myDataMyFactors, samples = NULL, factor = "Condition", norm = FALSE)

#to export as a pdf with a custom name
#QCreport(myDataMyFactors, file="QCreportForRNASeqWithNOISeq.pdf", samples = NULL, factor = "Condition", norm=FALSE)

#Normalization - Some way to address and normalize the following:
# - Within sample: For technical biases that affect detected counts
# - Between sample: For technical biases like lib preparation
# - Between batches: For systemic sources of variation and technical noise due to factors that cluster data independent of your experimental design.

#Normalize between samples for sequencing depth bias

#RPKM = reads per kilobase of transcript, per million mapped reads.
myRPKM = rpkm(assayData(myDataMyFactors)$exprs, k = 0, lc = 1)
head(myRPKM[,1:4])

#compare this to the original raw count data we started with
head(Chla_TPM[,1:4])

#Remove zero count and low count features (CPM or Wilcoxon test or Proportion test)

#Example with method 1 (CPM):

# depth = sequencing depth of samples (column totals before normalizing), used in low count feature filter methods 2 & 3
# cv.cutoff = cutoff for coefficient of variation used in method 1
# cpm = cutoff for low expression in counts per million (used in methods 1 & 3)
# p.adj = method for multiple testing correction

#mycounts
myfilt = filtered.data(Chla_TPM, factor = experimentalMatrix$Condition, norm=FALSE, depth=NULL, method=1, cv.cutoff=100, cpm=1, p.adj="fdr")
myfilt

#Differential Expression

#NOISeq with technical replicates
mynoiseqTech = noiseq(myDataMyFactors, k = 0.5, norm = "rpkm", factor = "Condition", pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "technical")

head(mynoiseqTech@results[[1]])
#Note: the output myresults@results[[1]]$prob gives the estimated probability of differential expression for each feature. These probabilities are not equivalent to p-values. The higher the probability, the more likely that the difference in expression is due to the change in the experimental condition and not to chance.

mynoiseqTech.deg = degenes(mynoiseqTech, q = 0.8, M = NULL)
#"1200 differentially expressed features"

#The probability of differential expression is not equivalent to 1 − pvalue. 
#Developers recommend using q values around 0.8 for `NOISeq` with technical replicates. If no technical replicates are available and `NOISeq-sim` is used, a more stringent threshold such as q = 0.9 is preferable.

#Plot expression with DEGs highlighted in red
par(mfrow = c(1, 1))
DE.plot(mynoiseqTech, q = 0.8, graphic = "expr", log.scale = TRUE)
DE.plot(mynoiseqTech, q = 0.8, graphic = "MD", log.scale = TRUE, title(main = "teste", xlab = "log-fold change", ylab = "absolute value of the difference in expression between conditions"))
?DE.plot


#export results to archive or for downstream analyses
#recall that for this particular output, the object mynoiseq.deg was filtered for q>= 0.8. For the complete gene list, use the `degenes` function above and set q = 0.
write.csv(mynoiseqTech.deg, file="DEGs_NOISeqTech.csv")
