#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# RNA Sequencing Project
# Steps 5 and 6 (exploratory and differential expression analysis)
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# load required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library("DESeq2")
library("pheatmap")
library("ggplot2")
library("RColorBrewer")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano", force=TRUE)
library("EnhancedVolcano")

# Read count matrix and change column names
count_matrix <- read.csv("/Users/leabroennimann/Desktop/Master_Bioinformatik/1._Semester/RNA_Seq/feature_counts_2nd_output/gene_matrix_count.csv", header = T, row.names = 1)
colnames(count_matrix) <- c("HER21", "Normal2", "NonTNBC3", "TNBC1", "NonTNBC2", "Normal3", "TNBC2", "HER22", "Normal1", "HER23", "TNBC3", "NonTNBC1")

# Create colData for DESeq2 object creation
sample_info <- data.frame(c("HER21", "Normal2", "NonTNBC3", "TNBC1", "NonTNBC2", "Normal3", "TNBC2", "HER22", "Normal1", "HER23", "TNBC3", "NonTNBC1"), 
                          c("HER2", "Normal", "NonTNBC", "TNBC", "NonTNBC", "Normal", "TNBC", "HER2", "Normal", "HER2", "TNBC", "NonTNBC" ))
colnames(sample_info) <- c("sample", "type")

sample_info <- sample_info[match(colnames(countDataMatrix), as.character(sample_info$sample)), ]

# Turn count matrix from dataframe to matrix
countDataMatrix <- as.matrix(count_matrix[ , ])

# change the order of the samples in matrix
col.order <- c("HER21", "HER22", "HER23", "NonTNBC1", "NonTNBC2", "NonTNBC3", "TNBC1", "TNBC2", "TNBC3", "Normal1", "Normal2", "Normal3")
countDataMatrix <- countDataMatrix[ , col.order]


# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = countDataMatrix, colData = sample_info, design = ~ type)
dds2 <- DESeq(dds)

# Transformations
vsd <- vst(dds2, blind=TRUE)
rlog <- rlog(dds2, blind=TRUE)

# Create heatmap
select <- order(rowMeans(counts(dds2,normalized=FALSE)), decreasing=TRUE)[1:20]

df <- as.data.frame(colData(dds2)[,"type"])

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, main="Heatmap of the count matrix")
?pheatmap

# Heatmap of the sample-to-sample distances

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Oranges")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Create PCA plots with different transformations
pcaData <- plotPCA(vsd, intgroup="type", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  labs(title="PCA plot of vst transformed data", colour="Group") +
  theme_light()


plotPCA(rlog, intgroup="type")





# Differential expression analysis
resultsNames(dds2)

# COMPARISON HER VS NORMAL
(results_HER2 <- results(dds2, contrast=c("type", "HER2", "Normal")))
(results_HER2_df <- as.data.frame(results_HER2))

# COMPARISON TNBC VS NORMAL
(results_TNBC <- results(dds2, contrast=c("type", "TNBC", "Normal")))
results_TNBC_df <- as.data.frame(results_TNBC)

# COMPARISON NonTNBC VS NORMAL
(results_NonTNBC <- results(dds2, contrast=c("type", "NonTNBC", "Normal")))
results_NonTNBC_df <- as.data.frame(results_NonTNBC)

# COMPARISON TNBC VS NonTNBC
(results_TNBC_Non_TNBC <- results(dds2, contrast=c("type", "TNBC", "NonTNBC")))
results_TNBC_Non_TNBC_df <- as.data.frame(results_TNBC_Non_TNBC)

# COMPARISON TNBC VS HER2
(results_TNBC_HER2 <- results(dds2, contrast=c("type", "TNBC", "HER2")))
results_TNBC_HER2_df <- as.data.frame(results_TNBC_HER2)


# How many genes are differentially expressed (DE) in the pairwise comparison you selected (e.g. padj < 0.05)
sum(results_TNBC_Non_TNBC$padj < 0.05, na.rm=TRUE)

# How many of the DE genes are down-regulated?
sum(results_TNBC_Non_TNBC$log2FoldChange < 0 & results_TNBC_Non_TNBC$padj < 0.05, na.rm=TRUE)

# How many of the DE genes are up-regulated?
sum(results_TNBC_Non_TNBC$log2FoldChange > 0 & results_TNBC_Non_TNBC$padj < 0.05, na.rm=TRUE)


(p_values_HER2 <- results_HER2$pvalue)

p_values_HER2[which(p_values_HER2 < 0.05)]

sum(results_HER2$padj < 0.05, na.rm=TRUE)

sum(results_TNBC$padj < 0.05, na.rm=TRUE)

sum(results_NonTNBC$padj < 0.05, na.rm=TRUE)

EnhancedVolcano(results_HER2_df,
                lab = rownames(results_HER2_df),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "HER2 vs WT",
                subtitle = element_blank(),
                colAlpha = 0.5,
                legendPosition = 'bottom',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                max.overlaps = 10,
                pointSize = 1,
                labSize = 4)


EnhancedVolcano(results_TNBC_df,
                lab = rownames(results_TNBC_df),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "TNBC vs WT",
                subtitle = element_blank(),
                colAlpha = 0.5,
                legendPosition = 'bottom',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                max.overlaps = 10,
                pointSize = 1,
                labSize = 4)

EnhancedVolcano(results_NonTNBC_df,
                lab = rownames(results_NonTNBC_df),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "NonTNBC vs WT",
                subtitle = element_blank(),
                colAlpha = 0.5,
                legendPosition = 'bottom',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                max.overlaps = 10,
                pointSize = 1,
                labSize = 4)


EnhancedVolcano(results_TNBC_Non_TNBC_df,
                lab = rownames(results_TNBC_Non_TNBC_df),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "TNBC vs NonTNBC",
                subtitle = element_blank(),
                colAlpha = 0.5,
                legendPosition = 'bottom',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                pointSize = 1,
                labSize = 4)


# Based on the original publication, select 2-3 genes that are of particular interest and investigate their expression level. 
# You could use, for example, the normalised counts (see DESeq2::counts) where the effect of between-sample differences in 
# sequencing depth has been removed.


counts <- counts(dds2, normalized = TRUE)

# SPARC
counts["ENSG00000113140", ]

# RACK1
counts["ENSG00000204628", ]

# B2M
counts["ENSG00000166710", ]






