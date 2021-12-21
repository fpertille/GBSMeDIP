# GBSMeDIP
GBS-MeDIP: A novel method for parallel identification of genetic and epigenetic variation in the same reduced fraction of genomes across individuals 

We developed the GBS-MeDIP, a method that combines two previously described techniques, Genotype-by-Sequencing (GBS) and Methylated DNA Immunoprecipitation (MeDIP). Our  method allows for parallel and cost-efficient interrogation of genetic and methylomic variants in large sample numbers of reduced genomes. Our protocol takes advantage of the barcoding of DNA samples performed in the GBS and the subsequent creation of DNA pools that are then used as an input for the MeDIP. The GBS-MeDIP allows to identify genetic and methylomic biomarkers when many individual samples are available and there is lack of resources for whole genome interrogation.


In this repository we are providing 4 files with instructions for data analysis and a file with the methodological description of this study. 

*Bioinformatic Processing and Analysis.docx is the methodological description. 

*BasicCommandsInBash_forRAWDATAanalysis.txt: is a .txt file where we provided basic bash commands to run the main programs for rawdata analysis. This programs use the .fastq files (with index not determined by Illumina's CASAVA program), it also provide commands for demultiplexing of the reads, read screening and trimming, read alignment against the reference genome, the merge commands for the .bam files according to the contrast they belong, and the call peaking calling between the merged files.

*Script1_MEDIPS.R: is a script that uses .bam files and peak coordinates obtained through bash commands to generate an .rda file with statistics for differentially methylated regions between a control group and a treated group of individuals.

*Script2_FineAnalysis.R: is a script that uses the .rda object obtained in the previous script to explore some features of the DMRs. The script provides commands to establish the cut-offline for defining the DMR, the sequence of nucleotides of each one of the regions analyzed, CpGs counting by analyzed region, and annotating the analyzed regions. All this information can be extract in a .txt file.
