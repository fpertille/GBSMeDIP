###############
# Authors: Fabio Pertille;
# Description: MeDIPS
# Date: 15th of jun, 2018
# Files input: .bam files from each one of the demultiplex individuals.
# Files output: .rda file with DMR statistics per treatment.
# Code: GBSMEDIP001

#necessary libraries

library("MEDIPS")
library("BSgenome.Hsapiens.UCSC.hg38")


#set wd where the .bam files are located per experiment
setwd("/.")

genome="BSgenome.Hsapiens.UCSC.hg38"


#set the chromosomes are going to be analysed
chr_all=paste("chr", c(1:22, "X", "Y"), sep="")#eliminating chr32, avoide the erro in XStringSet in evaluating the argument 'x' in selecting a method for function 'XStringSet': Error in ans[] <- x : replace


#MeDIP general parameters
#To avoid artefacts caused by PCR over amplification MEDIPS determines a maximal allowed number of stacked reads per genomic position by a poisson distribution of stacked reads genome wide and by a given p-value:
uniq=1e-3
#The smaller the p-value, the more reads at the same genomic position are potentially allowed. Alternatively, all reads mapping to exactly the same genomic position can be maintained (uniq = 0) or replaced by only one representative
#(uniq = 1).All reads will be extended to a length of 300nt according to the given strand information:
#extend=300
extend=100
#As an alternative to the extend parameter, the shift parameter can be used. Here, the reads are not extended but
#shifted by the specified number of nucleotides with respect to the given strand infomation. One of the two parameters
#extend or shift has to be 0.
shift=0
#The genome will be divided into adjacent windows of length 100nt and all further calculations (short read coverage,
#differential coverage between conditions etc.) will be applied to these windows.
ws=100


###########PREPARING DATA#############FOR MEDIPS###########START############

#create a list with the names of the individuals that will be analyzed based on the .bam files.
filenames <- list.files(pattern="*.bam$")
names <-sub('*.bam', '',filenames)

#run MEDIPS CpG enrichment, saving under sample_ID
g <- list()
for(i in seq_len(length(filenames))){
  richnames <-sub("","Enrich_",names)
  g[[i]] <-   MEDIPS.CpGenrich(file= filenames[i], BSgenome=genome, chr.select=chr_all, extend = extend, shift = shift, uniq = uniq)
  
}

#Create a list with all sets
f <- list()
for(i in seq_len(length(filenames))){
  setnames <-sub("","Set_",names)
  f[[i]] <- MEDIPS.createSet(file= filenames[i], BSgenome=genome, extend=extend, shift=shift, uniq=uniq, window_size=ws, chr.select=chr_all, paired=T)
}

#subset the list according to some pattern of the file name that differentiates the individuals from the treated group (_T) versus the control group (_C).
x <- grep("_T", unlist(setnames))
subsetList <- function(myList, elementNames) {
  lapply(elementNames, FUN=function(x) myList[[x]])
}
###new name list
MeDIP_T <-subsetList(f, x)

x <- grep("_C", unlist(setnames))
subsetList <- function(myList, elementNames) {
  lapply(elementNames, FUN=function(x) myList[[x]])
}
###new name list
MeDIP_C <- subsetList(f, x)


#For CpG density dependent normalization of MeDIP-seq data, we need to generate a coupling set. The coupling set
#must be created based on the same reference genome, the same set of chromosomes, and with the same window size
#used for the MEDIPS SETs. For this, we specify the first MEDIPS SET in the hESCs object as reference, but any
#of the other MEDIPS SETs would be fine as well, because all of them consist of the same set of chromosomes (here
#chr22 only, hg19) and have been generated with the same window size.
#CS = MEDIPS.couplingVector(pattern = "CG", refObj = MeDIP_allRJFGBS)
CS = MEDIPS.couplingVector(pattern = "CG", refObj = MeDIP_C[[1]]) #a sample that have a very good coverage compare with the other ones

###########PREPARING DATA#############FOR ORDINARY MEDIPS###########END############



#############################################################################################
#############################################################################################
######Alternatively PREPARE DATA to build libraries based on MACS2 peak calling coordinates##
#############################################################################################
#############################################################################################
#set worspace where the output of MACS2 is located
setwd(/.)
ROI<- read.delim("xxxx.broadPeak", sep="", header=F, stringsAsFactors = FALSE)
newname<-list("chr", "start", "end", "name")
ROI<- data.frame(ROI$V1, ROI$V2, ROI$V3, ROI$V4, stringsAsFactors = FALSE)
names(ROI)<-newname

filenames_T <- list.files (pattern="*_T_sorted.bam$")
filenames_C <- list.files (pattern="*_C_sorted.bam$")
###############################################################
###############################################################
r_T <- list()
for(i in seq_len(length(filenames_T))){
  setnames <-sub("","ROISet_",names)
  r_T[[i]] <- MEDIPS.createROIset(file= filenames_T[i], BSgenome=genome, extend=extend, shift=shift, uniq=uniq, ROI=ROI, chr.select=chr_all, paired=T)
}

r_C <- list()
for(i in seq_len(length(filenames_C))){
  setnames <-sub("","ROISet_",names)
  r_C[[i]] <- MEDIPS.createROIset(file= filenames_C[i], BSgenome=genome, extend=extend, shift=shift, uniq=uniq, ROI=ROI, chr.select=chr_all, paired=T)
}

CST = MEDIPS.couplingVector(pattern = "CG", refObj = r_T[[1]])


#####################RUNNING MEDIP CALL ON THE PREPARED DATA#######################################################

#It is possible to calculate genome wide coverage and methylation profiles for only one MEDIPS SET or for only one
#group of MEDIPS SETs using the function MEDIPS.meth. However, in this case study we also want to calculate
#differential coverage (i.e. differential methylation) between two conditions. Whenever two groups of MEDIPS SETs
#are provided to the MEDIPS.meth function differential coverage will be calculated.

######for MEDIPS based on adjacent windows (ADWs)
mr.edgeR = MEDIPS.meth(MSet1 = MeDIP_T, MSet2 = MeDIP_C, CSet = CS, p.adj = "fdr", diff.method = "edgeR", MeDIP = F, CNV = F, minRowSum = 10)

######for MEDIPS based on reagions of interest (ROIs) identified using MACS2
mr.edgeROI = MEDIPS.meth(MSet1 = r_T, MSet2 = r_C, CSet = CST, p.adj = "fdr", diff.method = "edgeR", MeDIP = F, CNV = F, minRowSum = 10)

#omitting NAS
mr.edgeR <- na.omit(mr.edgeR, cols=c("edgeR.p.value"))
mr.edgeROI <- na.omit(mr.edgeROI, cols=c("edgeR.p.value"))

#saving R objects for further analysis
save(mr.edgeR, file="mr.edgeR_.rda")
save(mr.edgeROI, file="mr.edgeR_ROI.rda")






