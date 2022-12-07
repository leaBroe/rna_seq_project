
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library("DESeq2")

count_matrix <- read.csv("/Users/leabroennimann/Desktop/Master_Bioinformatik/1._Semester/RNA_Seq/feature_counts_2nd_output/gene_matrix_count.csv", header = T, row.names = 1)
colnames(count_matrix) <- c("HER21", "Normal2", "NonTNBC3", "TNBC1", "NonTNBC2", "Normal3", "TNBC2", "HER22", "Normal1", "HER23", "TNBC3", "NonTNBC1")


sample_info <- data.frame(c("HER21", "Normal2", "NonTNBC3", "TNBC1", "NonTNBC2", "Normal3", "TNBC2", "HER22", "Normal1", "HER23", "TNBC3", "NonTNBC1"), 
                          c("HER2", "Normal", "NonTNBC", "TNBC", "NonTNBC", "Normal", "TNBC", "HER2", "Normal", "HER2", "TNBC", "NonTNBC" ))
colnames(sample_info) <- c("sample", "type")

countDataMatrix <- as.matrix(count_matrix[ , ])

col.order <- c("HER21", "HER22", "HER23", "NonTNBC1", "NonTNBC2", "NonTNBC3", "TNBC1", "TNBC2", "TNBC3", "Normal1", "Normal2", "Normal3")
countDataMatrix <- countDataMatrix[ , col.order]

dds <- DESeqDataSetFromMatrix(countData = countDataMatrix, colData = sample_info, design = ~ type)
dds2 <- DESeq(dds)

vsd <- vst(dds2, blind=TRUE)
rlog <- rlog(dds2, blind=TRUE)

library("pheatmap")

rlog <- rlog(dds)
select <- order(rowMeans(counts(dds2,normalized=FALSE)), decreasing=TRUE)[1:20]

df <- as.data.frame(colData(dds2)[,"type"])

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, main="Heatmap of the count matrix")
?pheatmap

plotPCA(vsd, intgroup="type")
plotPCA(rlog, intgroup="type")

resultsNames(dds2)

(results_HER2 <- results(dds2, contrast=c("type", "HER2", "Normal")))

(results_TNBC <- results(dds2, contrast=c("type", "TNBC", "Normal")))

(results_NonTNBC <- results(dds2, contrast=c("type", "NonTNBC", "Normal")))


(p_values_HER2 <- results_HER2$pvalue)

p_values_HER2[which(p_values_HER2 < 0.05)]







