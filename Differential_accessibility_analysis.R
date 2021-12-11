```{r}

#Differential accessibility analysis - HPEC vs CAPAN2

library(GenomicRanges)
library(csaw)
library(edgeR)
library(dplyr)
library(data.table)
library(ggplot2)

pe.bams <- c("~/ATAC-Seq/Round2/differential/Pancreas/ATAC_J/results/Bam/ATAC_J_markedDuplicates.bam", "~/ATAC-Seq/Round2/differential/Pancreas/ATAC_K/results/Bam/ATAC_K_markedDuplicates.bam", "~/ATAC-Seq/Round2/ATAC_L/results/Bam/ATAC_L_markedDuplicates.bam", "~/ATAC-Seq/Round2/differential/Pancreas/ATAC_M/results/Bam/ATAC_M_markedDuplicates.bam")

##############################
# read mm10 blacklist
blacklist <- read.table("~/ENCODE_resources/ENCFF356LFX_hg38.bed", sep="\t")
colnames(blacklist) <- c("chrom", "start", "end")
blacklist <- GRanges(blacklist)

# define read parameters
standard.chr <- paste0("chr", c(1:19, "X", "Y")) # only use standard chromosomes
param <- readParam(max.frag=1000, pe="both", discard=blacklist, restrict=standard.chr)


##############################
# count BAM reads in, e.g. 300 bp windows
counts <- windowCounts(pe.bams, width=300, param=param) # set width as desired from the fragment length distribution analyses

# filter uninteresting features (windows) by local enrichment
# local background estimator: 2kb neighborhood
neighbor <- suppressWarnings(resize(rowRanges(counts), width=2000, fix="center")) # change width parameter as desired
wider <- regionCounts(pe.bams, regions=neighbor, param=param) # count reads in neighborhoods
# filter.stat <- filterWindows(counts, wider, type="local") # the filterWindows() function is deprecated and has been replaced by filterWindowsLocal(). This is an archived step.
filter.stat <- filterWindowsLocal(counts, wider)
counts.local.filt <- counts[filter.stat$filter > log2(3),] # threshold of 3-fold increase in enrichment over 2kb neighborhood abundance; change as desired

# method 4: csaw de novo peaks by local enrichment, csaw loess-normalization
counts.local.loess <- counts.local.filt
counts.local.loess <- normOffsets(counts.local.loess, method="loess", se.out=TRUE)

# DIFFERENTIAL ACCESSIBILITY ANALYSIS
working.windows <- counts.local.loess

# setup design matrix
# see edgeR manual for more information
y <- asDGEList(working.windows)
colnames(y$counts) <- c("HPEC_1", "HPEC_2", "CAPAN2_1", "CAPAN2_2")
rownames(y$samples) <- c("HPEC_1", "HPEC_2", "CAPAN2_1", "CAPAN2_2")
y$samples$group <- c("HPEC", "HPEC", "CAPAN2", "CAPAN2")
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- c("HPEC", "CAPAN2")
# design

# stabilize dispersion estimates with empirical bayes
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)

# testing for differentially-accessible windows
results <- glmQLFTest(fit, contrast=makeContrasts(CAPAN2-HPEC, levels=design))
# head(results$table)
rowData(working.windows) <- cbind(rowData(working.windows), results$table) # combine GRanges rowdata with differential statistics
# working.windows@rowRanges

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case; max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(working.windows), tol=500L, max.width=5000L)
# summary(width(merged.peaks$region))
# should merge some peaks; change as desired

# use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)

final.merged.peaks <- merged.peaks$region
final.merged.peaks@elementMetadata <- cbind(final.merged.peaks@elementMetadata, tab.best[,-1])

# sort by FDR
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ]

# filter by FDR threshold
FDR.thresh <- 0.05 # set as desired
final.merged.peaks.sig <- final.merged.peaks[final.merged.peaks@elementMetadata$FDR < FDR.thresh, ]
final.merged.peaks.sig # significant differentially-accessible windows

write.table(final.merged.peaks, "~/CAPAN2_vs_HPEC_csaw_DA-windows_all.txt", sep="\t", quote=F, col.names=T, row.names=F)
write.table(final.merged.peaks.sig, "~/CAPAN2_vs_HPEC_csaw_DA-windows_significant.txt", sep="\t", quote=F, col.names=T, row.names=F)

final.merged.peaks = as.data.frame(final.merged.peaks)
  
New = filter(final.merged.peaks, logFC > 0)
Lost = filter(final.merged.peaks, logFC < 0)
HPEC_naive_overlap = as.data.frame(read.table("~/ATAC-Seq/Round2/differential/Pancreas/naive/HPEC_naive_overlap.filt.broadPeak"))[,1:3]
colnames(HPEC_naive_overlap) = c("seqnames", "start" ,"end")
HPEC_naive_overlap_gr = makeGRangesFromDataFrame(HPEC_naive_overlap, keep.extra.columns=F)

#######
##Define what is significant for New regions
#######

New_gr = makeGRangesFromDataFrame(New, keep.extra.columns=T)
overlaps = countOverlaps(New_gr, HPEC_naive_overlap_gr, minoverlap = 50L)
New = as.data.frame(New_gr)
New$overlap = overlaps
New$sig = ifelse(( (New$FDR < 0.05 ) & (New$logFC > 3)),"Significant","N.S")
New$overlapping = ifelse(( (New$overlap > 0)) ,"Yes","No")

#a = filter(New, sig == "Significant")
#dim(a) # 2,519

#######
##Define what is significant for Lost regions
#######

Lost_gr = makeGRangesFromDataFrame(Lost, keep.extra.columns=T)
overlaps = countOverlaps(Lost_gr, HPEC_naive_overlap_gr, minoverlap = 50L)
Lost = as.data.frame(Lost_gr)
Lost$overlap = overlaps
Lost$sig = ifelse(( (Lost$FDR < 0.05 ) & (Lost$logFC < -3)),"Significant","N.S")
Lost$overlapping = ifelse(( (Lost$overlap > 0)) ,"Yes","No")

a = filter(Lost, overlapping == "No")
dim(a) # 6,805

DA_HPEC_CAPAN = rbind(New, Lost)

nrow(New %>% filter(sig == "Significant" & overlapping == "Yes")) # 5931
nrow(New %>% filter(sig == "Significant" & overlapping == "No")) # 874

nrow(Lost %>% filter(sig == "Significant" & overlapping == "Yes")) # 136
nrow(Lost %>% filter(sig == "Significant" & overlapping == "No")) # 2383 


#########
##Volcano plot of results
#########

pdf("~/tmp.pdf", 7, 5)
ggplot(data=data.frame(DA_HPEC_CAPAN)) +
        geom_point(aes(x=logFC, y=-log10(FDR), col = factor(sig, levels=c("N.S", "Significant")))) +
        ggtitle("Differential accessibility between HPEC and CAPAN-2") +
        xlab("log2 fold change") + 
        ylab("-log10 false discovery rate") +
        #scale_y_continuous(limits = c(0,50)) +
        theme(legend.position = "none",
              plot.title = element_text(size = rel(1.5), hjust = 0.5),
              axis.title = element_text(size = rel(1.25))) + scale_color_manual(values = c("black", "red"))
dev.off()




Differential_accessibility_regions <- c(rep("Gained" , 2) , rep("Lost" , 2) )
Overlapping_HPEC_peak <- rep(c("Yes" , "No") , 2)
Number <- c(5931,874,136,2383)
data <- data.frame(specie,Overlapping_HPEC_peak,value)
 
# Stacked
ggplot(data, aes(fill=Overlapping_HPEC_peak, y=Number, x=Differential_accessibility_regions)) + 
    geom_bar(position="stack", stat="identity")





```
