# RNA-Seq-Analysis
Pipeline for analysing RNA sequencing data from fastq to adjusted p-values and box plots
The bash script informs how to trim reads for higher quality, using Trimmomatic. The quality of the reads is checked using fastQC and multiQC. Salmon is used for quantification of transcript expression and to correct for GC bias.
The R script uses EdgeR for differential expression using Fisher's exact test, principal component analysis and False discovery rates.
