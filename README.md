# GBSMeDIP

GBS-MeDIP: A novel method for parallel identification of genetic and epigenetic variation in the same reduced fraction of genomes across individuals

Short summary

We developed the GBS-MeDIP, a method that combines two previously described techniques, Genotype-by-Sequencing (GBS) and Methylated DNA Immunoprecipitation (MeDIP). Our method allows for parallel and cost-efficient interrogation of genetic and methylomic variants in large sample numbers of reduced genomes. Our protocol takes advantage of the barcoding of DNA samples performed in the GBS and the subsequent creation of DNA pools that are then used as an input for the MeDIP. The GBS-MeDIP allows to identify genetic and methylomic biomarkers when many individual samples are available and there is lack of resources for whole genome interrogation.

Short description

The genomic library undergoes PstI enzymatic selection which represents the GBS reduced-representation method firstly described by Elshire, R (2011)(1) in maize. This method, named GBS-MeDIP, is based on the combination of Genotype-by-Sequencing (GBS, initially described in maize(1), and employed in chickens for the first time by our group(2)) and methylated DNA immunoprecipitation (MeDIP)(3). GBS-MeDIP addresses DNA methylation variations in reduced genomes of large number of animals in a cost-effective manner. This combination of techniques was done because current methods assessing DNA methylation in reduced genomes are biased towards CpG sites located in CpG-rich regions of the genome(4–7). GBS, instead, reduces the genome by digesting the genome independently of CpG density. This enzymatic fragmentation reduces the genome to approximately 2% of its original size(8). Based on the GBS, the fragments (200-500 bp) are then barcoded to later identify individual samples in a pool. Sequencing libraries are then created with DNA pooled from several individuals. The methylated fraction is then captured from the pooled and barcoded DNA by MeDIP(9). In fact, this representation is unique, since we found regions never associated before with performance traits(10), and feather pecking behavior (non-published data) in chickens. Moreover, our method proved to be applicable, in this sense, in a wide range of studies and species models(8,10–18). As the template for this article is protocol-focused, we have decided to provide methodological information as supplementary material. Therefore, we have uploaded a file called "Overview and development of the GBS-MeDIP protocol:" in the GitHub repository. Please see the "Data and code availability" section.


In this repository we are providing one file with the Overview and development of the GBS-MeDIP protocol, one file with the methodological description of this study, two files with the "adapter templates", and 4 files with instructions for data analysis.

*Bioinformatic Processing and Analysis.docx is a file that contains the methodological description of the data used in the referred paper.

*BasicCommandsInBash_forRAWDATAanalysis.txt: is a .txt file where we provide basic bash commands and the main programs for raw data analysis. These programs use the .fastq files (with index not determined by Illumina's CASAVA program). Moreover, it also provide commands for demultiplexing of the reads, read screening and trimming, read alignment against the reference genome, the merge commands for the .bam files according to the contrast they belong and, the peaking calling between the merged files.

*Script1_MEDIPS.R: is a script with bash commands, which uses .bam files and peak coordinates as input and generate an .rda file as output. This output contains statistics for differentially methylated regions between a control group and a treated group of individuals which underwent the GBS-MeDIP approach .

*Script2_FineAnalysis.R: is a script that uses the .rda object obtained in the previous script (Script 1) to explore some features of the analyzed regions. The script provides commands to subset significant DMRs, the sequence of nucleotides of each one of the regions analyzed, CpG counting by each one of the analyzed regions, and annotating of each one of the analyzed regions against the available reference genome. Finally, we provide commands to extract these merged information’s in a .txt file.

PS: We are  providing here a practical driving guide and not a suggested analysis. There are several programs and pipelines published in recognized journals that should be considered when choosing the approach according to the objectives of each researcher.


Lead contact

Fábio Pértille
fabiopertille@gmail.com
Carlos Guerrero Bosagna
carlos.guerrero.bosagna@ebc.uu.se 

Data and code availability

The dataset generated with the samples employed in this article is available from the European Nucleotide Archive (ENA) repository (EMBL-EBI), under accession number PRJEB35669 (www.ebi.ac.uk/ena/data/view/PRJEB35669). All the programs and R packages used in this protocol are open access and all the scripts are available online in the manual of the programs or packages used  . For more details, refer  to our previous publications employing the GBS-MeDIP (Pértille, Brantsæter, et al., 2017; Pertille et al., 2019; Guerrero-Bosagna et al., 2020; Pértille et al., 2020, 2021).

References

01. 	Elshire RJ, Glaubitz JC, Sun Q, Poland JA, Kawamoto K, Buckler ES, et al. A Robust, Simple Genotyping-by-Sequencing (GBS) Approach for High Diversity Species. Orban L, editor. PLoS One. 2011 May 4;6(5):e19379. 
02. 	Pértille F, Guerrero-Bosagna C, Silva VH da, Boschiero C, Nunes J de R da S, Ledur MC, et al. High-throughput and Cost-effective Chicken Genotyping Using Next-Generation Sequencing. Sci Rep. 2016 May 25;6(January):26929. 
03. 	Bock C, Tomazou EM, Brinkman AB, Müller F, Simmer F, Gu H, et al. Quantitative comparison of genome-wide DNA methylation mapping technologies. Nat Biotechnol. 2010 Oct 19;28(10):1106–14. 
04. 	Habermann F a., Cremer M, Walter J, Kreth G, von Hase J, Bauer K, et al. Arrangements of macro- and microchromosomes in chicken cells. Chromosome Res. 2001;9(7):569–84. 
05. 	McQueen HA, Fantes J, Cross SH, Clark VH, Archibald AL, Bird AP. CpG islands of chicken are concentrated on microchromosomes. Nat Genet. 1996;12(3):321–4. 
06. 	Smith J, Burt DW. Parameters of the chicken genome (Gallus gallus). Anim Genet. 1998;29(4):290–4. 
07. 	Smith J, Bruley CK, Paton IR, Dunn I, Jones CT, Windsor D, et al. Differences in gene density on chicken macrochromosomes and microchromosomes. Anim Genet. 2000;31(2):96–103. 
08. 	Pértille F, Brantsæter M, Nordgreen J, Coutinho LL, Janczak AM, Jensen P, et al. DNA methylation profiles in red blood cells of adult hens correlate with their rearing conditions. J Exp Biol. 2017 Oct 1;220(19):3579–87. 
09. 	Guerrero-Bosagna C, Jensen P. Optimized method for methylated DNA immuno-precipitation. MethodsX. 2015;2:432–9. 
10. 	Pértille F, Moreira GCM, Zanella R, Nunes J de R da S, Boschiero C, Rovadoscki GA, et al. Genome-wide association study for performance traits in chickens using genotype by sequencing approach. Sci Rep. 2017 Feb 9;7(September 2016):41748. 
11. 	Pértille F, Guerrero-Bosagna C, Silva VH da, Boschiero C, Nunes J de R da S, Ledur MC, et al. High-throughput and Cost-effective Chicken Genotyping Using Next-Generation Sequencing. Sci Rep. 2016 May 25;Major Ri(January):26929. 
12. 	Nunes J de R da S, Liu S, Pértille F, Perazza CA, Villela PMS, de Almeida-Val VMF, et al. Large-scale SNP discovery and construction of a high-density genetic map of Colossoma macropomum through genotyping-by-sequencing. Sci Rep. 2017 Apr 7;7(March):46112. 
13. 	Pértille F, Da Silva VH, Johansson AM, Lindström T, Wright D, Coutinho LL, et al. Mutation dynamics of CpG dinucleotides during a recent event of vertebrate diversification. Epigenetics. 2019 May 9;14(1):1–23. 
14. 	Guerrero-Bosagna C, Pértille F, Gomez Y, Rezaei S, Gebhardt-Henrich SG, Vögeli S, et al. DNA methylation variation in the brain of laying hens in relation to differential behavioral patterns. Comp Biochem Physiol - Part D Genomics Proteomics. 2020;35(May):100700. 
15. 	Pértille F, Alvarez-Rodriguez M, Silva AN da, Barranco I, Roca J, Guerrero-Bosagna C, et al. Sperm Methylome Profiling Can Discern Fertility Levels in the Porcine Biomedical Model. Int J Mol Sci. 2021 Mar 6;22(5):2679. 
16. 	Sundman A, Pértille F, Lehmann Coutinho L, Jazin E, Guerrero-Bosagna C, Jensen P. DNA methylation in canine brains is related to domestication and dog-breed formation. Bartos L, editor. PLoS One. 2020 Oct 29;15(10):e0240787. 
17. 	Nunes JRS, Pértille F, Andrade SCS, Perazza CA, Villela PMS, Almeida-Val VMF, et al. Genome-wide association study reveals genes associated with the absence of intermuscular bones in tambaqui (Colossoma macropomum). Anim Genet. 2020; 
18. 	Pértille F, Ibelli AMG, Sharif M El, Poleti MD, Fröhlich AS, Rezaei S, et al. Putative Epigenetic Biomarkers of Stress in Red Blood Cells of Chickens Reared Across Different Biomes. Front Genet. 2020 Nov 2;11(Llc):1–44. 
