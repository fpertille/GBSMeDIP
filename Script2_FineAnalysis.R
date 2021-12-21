###############
# Authors: Fabio Pertille;
# Description: Fine analysis to functional annotation
# Date: 10th of oct, 2018
# Input Files: .rda
# Outpu Files: .txt
# Code: GBSMEDIP002

module load HPC/R-3.5.1
R

#libraries
library("GenomicRanges")
library("BSgenome.Hsapiens.UCSC.hg38")
library(seqinr)
library("rGADEM")
library("Biostrings")
library(stats)
library(ade4)
library(ape)
library(adegenet)
library(phangorn)
library(hierfstat)
library(vcfR)
library(pegas)
library(poppr)
library(ggplot2)


#set a workspace where your .rda files were created
setwd("/.")
#load them
load(file="mr.edgeR_.rda")
load(file="mr.edgeR_ROI.rda")

#create a genomic range object to work using mr.edegeR object, for example
library("GenomicRanges")
GR_Human<-with(mr.edgeR, GRanges(chr, IRanges(start, stop), strand, edgeR.logFC, edgeR.p.value, edgeR.adj.p.value))
#extract the sequence of the ranges
Seq_GR = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, GR_Human)
#seq the number of CpGs per DMR
seq_CG_Human <- vcountPattern("CG", Seq_GR)

#add CpG collumn to the Grange objects
CG<- seq_CG_Human
Human <- cbind(Human, CG)

#Genomic Ranges Object including CG counts
GR_Human <-with(Human, GRanges(chr, IRanges(start, stop), strand, edgeR.logFC, edgeR.p.value, edgeR.adj.p.value, CG))

#subset by p-value: Ex: All low p-value 0.05
GR_Human05 <- subset(GR_Human, edgeR.p.value<.05)

#For annotation
supportedUCSCtables(genome = "hg38")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

library('TxDb.Hsapiens.UCSC.hg38.knownGene')
hs_txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#hs_txdb <- makeTxDbFromUCSC(genome="hg38", tablename="ensGene")
hs_txdb <- makeTxDbFromUCSC(genome="hg38", tablename="xenoRefGene")

View(hs_txdb)

peak_file <- as.data.frame(GR_Human)
peak_file <- subset(peak_file, select=c(seqnames, start,end))
save(peak_file, file = "peak_file.rda")
names(peak_file) <- c("CHR", "BP", "BP2")

library ("org.Hs.eg.db")
peak_anno_function <- function(input_data){
  b2 <- input_data[,c("CHR" , "BP" , "BP2")]
  write.table(b2 , "annotate_file.txt" , sep = "\t" , col.names = TRUE , row.names = FALSE)
  peakAnno <- annotatePeak("somewhereInYourHD/annotate_file.txt",
                           TxDb=hs_txdb, annoDb="org.Hs.eg.db")
  return(peakAnno)
}

Human <- as.data.frame (GR_Human)
peak <- peak_anno_function(peak_file)
plotAnnoPie(peak)
info <- peak@anno@elementMetadata
all <- cbind(Human, info , all=T)

write.table(Human, "Annotate_Human.txt" , sep="/t" , row.names = F , quote = F)


######the fine analysis depends a lot on the researcher's interest. Please send an email to the correspondence author if you have any questions.#########