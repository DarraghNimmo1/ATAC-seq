#!/usr/bin/env Rscript

###########################
## DARRAGH NIMMO
## JANUARY 2021
## ATACseqQC plots for ATAC-seq quality assessment
## To run on command line type 'Rscript ATACseqQC_plots.R my.bam chroms.txt'
##########################


###################################################################################################################################################
#####LOAD PACKAGES
###################################################################################################################################################
library(motifStack)
library(ATACseqQC)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(MotifDb)
