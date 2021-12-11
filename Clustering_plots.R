```{r}
library(ggplot2)

count_matrix = as.data.frame(read.table("~/ATAC-Seq/Round2/clustering/bigwigs/raw_counts.bed", header = FALSE, sep="\t"))

head(count_matrix)

colnames(count_matrix) = c("chr", "start", "end", "MDAMB468_replicate1", "MDAMB468_replicate2", "MCF7_replicate1", "MCF7_replicate2", "MCF10A_replicate1", "MCF10A_replicate2", "HPEC_replicate1", "HPEC_replicate2", "CAPAN_2_replicate1", "CAPAN_2_replicate2", "BxPC_3_replicate1", "BxPC_3_replicate2","HMEC_replicate1", "HMEC_replicate2" )

rownames(count_matrix) = paste0(count_matrix$chr, '_', count_matrix$start, '_', count_matrix$end)

count_matrix = count_matrix[count_matrix$chr %in% c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY" ), ]

count_matrix = count_matrix[,4:ncol(count_matrix)]

sample = c("MDAMB468_replicate1", "MDAMB468_replicate2", "MCF7_replicate1", "MCF7_replicate2", "MCF10A_replicate1", "MCF10A_replicate2", "HPEC_replicate1", "HPEC_replicate2", "CAPAN_2_replicate1", "CAPAN_2_replicate2", "BxPC_3_replicate1", "BxPC_3_replicate2","HMEC_replicate1", "HMEC_replicate2")

celltype = c("MDAMB468","MDAMB468", "MCF7", "MCF7", "MCF10A", "MCF10A", "HPEC", "HPEC", "CAPAN2", "CAPAN2", "BxPC3", "BxPC3", "HMEC", "HMEC")
donor = c("NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL")
si = data.frame(sample,celltype,donor )

pca = prcomp(t(count_matrix))
print(summary(pca))

pcaData = as.data.frame(pca$x)
head(pcaData)

pcaData$sample=rownames(pcaData)

pcaData=merge(pcaData, si)

percentVar = round(100 * (pca$sdev^2 / sum( pca$sdev^2 ) ))

p=ggplot(data=pcaData, aes(x = PC1, y = PC2, color=celltype)) + geom_point(size=3)
p=p+xlab(paste0("PC1: ", percentVar[1], "% variance"))
p=p+ylab(paste0("PC2: ", percentVar[2], "% variance"))

pdf("~/pca.pdf", 7, 5)
print(p)
dev.off()

```

```{r}
install.packages("Rtsne")
library(Rtsne)

count_matrix

count_matrix_transposed = t(count_matrix)

dim(count_matrix_transposed)



tsne_realData <- Rtsne(count_matrix_transposed,perplexity=4, check_duplicates = FALSE)
par(mfrow=c(1,2))
plot(tsne_realData$Y, col = "blue", pch = 19, cex = 1)

plot(tsne_realData$Y, col = "blue", bg= as.factor(si$celltype), pch = 21, cex = 1)

tsne_plot <- data.frame(x = tsne_realData$Y[,1], y = tsne_realData$Y[,2], col = si$celltype)

pdf("~/tsne.pdf", 7, 5)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col),size = 4)+ xlab("t-SNE Dimension 1") + ylab("t-SNE Dimension 2")
dev.off()

```
