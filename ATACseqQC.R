#!/usr/bin/env Rscript

###########################
## DARRAGH NIMMO
## JANUARY 2021
## ATACseqQC plots for ATAC-seq quality assessment
##########################


###################################################################################################################################################
#####LOAD PACKAGES
###################################################################################################################################################
library(ATACseqQC)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#library(BSgenome.Hsapiens.UCSC.hg38)

######################################################################################################################################################
#####SET UP VARIABLES AND PATHS
######################################################################################################################################################
## Getting the BAM file name and sample ID
args = commandArgs(trailingOnly=TRUE)
bamfile = args[1]
bamfile.sample.ID = gsub(".bam", "", basename(bamfile))


#######################################################################################################################################################
#####LIBRARY COMPLEXITY PLOT
#######################################################################################################################################################
histFile <- readsDupFreq(bamFile=bamfile, index =bamfile)
pdf(paste0(bamfile,".library.complexity.pdf"), width=8, height=5)
complexity <- estimateLibComplexity(histFile=histFile, times=1, interpolate.sample.sizes=seq(0.1, 1, by=0.01), extrapolate.sample.sizes=seq(5, 20, by=5))
dev.off()
write.table(complexity, paste0(bamfile, ".complexity.txt"), sep ="\t", quote=F, row.names=FALSE)


#######################################################################################################################################################
#####FRAGMENT SIZE DISTRIBUTION PLOT
#######################################################################################################################################################
pdf(paste0(bamfile.sample.ID, ".fragment.size.distribution.pdf"),width =10, height=8)
fragSize <- fragSizeDist(bamFiles=bamfile, bamFiles.labels=bamfile.sample.ID)
dev.off()
