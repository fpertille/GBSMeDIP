# GBSMeDIP
GBS-MeDIP: A novel method for parallel identification of genetic and epigenetic variation in the same reduced fraction of genomes across individuals
We developed the GBS-MeDIP, a method that combines two previously described techniques, Genotype-by-Sequencing (GBS) and Methylated DNA Immunoprecipitation (MeDIP). Our method allows for parallel and cost-efficient interrogation of genetic and methylomic variants in large sample numbers of reduced genomes. Our protocol takes advantage of the barcoding of DNA samples performed in the GBS and the subsequent creation of DNA pools that are then used as an input for the MeDIP. The GBS-MeDIP allows to identify genetic and methylomic biomarkers when many individual samples are available and there is lack of resources for whole genome interrogation.

In this repository we are providing one file with the methodological description of this study, and 4 files with instructions for data analysis.

*Bioinformatic Processing and Analysis.docx is a file that contains the methodological description of the data used in the referred paper.

*BasicCommandsInBash_forRAWDATAanalysis.txt: is a .txt file where we provide basic bash commands and the main programs for raw data analysis. These programs use the .fastq files (with index not determined by Illumina's CASAVA program). Moreover, it also provide commands for demultiplexing of the reads, read screening and trimming, read alignment against the reference genome, the merge commands for the .bam files according to the contrast they belong and, the peaking calling between the merged files.

*Script1_MEDIPS.R: is a script with bash commands, which uses .bam files and peak coordinates as input and generate an .rda file as output. This output contains statistics for differentially methylated regions between a control group and a treated group of individuals which underwent the GBS-MeDIP approach .

*Script2_FineAnalysis.R: is a script that uses the .rda object obtained in the previous script (Script 1) to explore some features of the analyzed regions. The script provides commands to subset significant DMRs, the sequence of nucleotides of each one of the regions analyzed, CpG counting by each one of the analyzed regions, and annotating of each one of the analyzed regions against the available reference genome. Finally, we provide commands to extract these merged information’s in a .txt file.

PS: We are  providing here a practical driving guide and not a suggested analysis. There are several programs and pipelines published in recognized journals that should be considered when choosing the approach according to the objectives of each researcher.
