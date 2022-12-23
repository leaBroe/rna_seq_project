if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)

# Get rid of NAs 
results_TNBC_Non_TNBC_df_omit_NA <- na.omit(results_TNBC_Non_TNBC_df)

results_TNBC_df_omit_NAs <- na.omit(results_TNBC_df)

results_TNBC_HER2_df_omit_NAs <- na.omit(results_TNBC_HER2_df)

results_NonTNBC_HER2_df_omit_NAs <- na.omit(results_NonTNBC_HER2_df)

results_HER2_df_omit_NAs <- na.omit(results_HER2_df)

results_NonTNBC_df_omit_NAs <- na.omit(results_NonTNBC_df)

# Select all Gene IDS for TNBC vs NonTNBC DE genes
gene_IDs_TNBC_NonTNBC <- row.names(results_TNBC_Non_TNBC_df_omit_NA[results_TNBC_Non_TNBC_df_omit_NA$padj<0.05 & abs(results_TNBC_Non_TNBC_df_omit_NA$log2FoldChange) > 1,])

# Select all Gene IDS for TNBC vs normal DE genes
gene_IDs_TNBC <- row.names(results_TNBC_df_omit_NAs[results_TNBC_df_omit_NAs$padj<0.05 &  abs(results_TNBC_df_omit_NAs$log2FoldChange) > 1,])

# Select all Gene IDs for TNBC vs HER2 DE genes
gene_IDs_TNBC_HER2 <- row.names(results_TNBC_HER2_df_omit_NAs[results_TNBC_HER2_df_omit_NAs$padj<0.05 &  abs(results_TNBC_HER2_df_omit_NAs$log2FoldChange) > 1,])

# Select all Gene IDs for NonTNBC vs HER2 DE genes
gene_IDs_NonTNBC_HER2 <- row.names(results_NonTNBC_HER2_df_omit_NAs[results_NonTNBC_HER2_df_omit_NAs$padj<0.05 &  abs(results_NonTNBC_HER2_df_omit_NAs$log2FoldChange) > 1,])

# Select all Gene IDs for HER2 vs normal DE genes
gene_IDs_HER2_Normal <- row.names(results_HER2_df_omit_NAs[results_HER2_df_omit_NAs$padj<0.05 &  abs(results_HER2_df_omit_NAs$log2FoldChange) > 1,])

# Select all Gene IDs for NonTNBC vs normal DE genes
gene_IDs_NonTNBC_Normal <- row.names(results_NonTNBC_df_omit_NAs[results_NonTNBC_df_omit_NAs$padj<0.05 &  abs(results_NonTNBC_df_omit_NAs$log2FoldChange) > 1,])


# Check with sum(results_TNBC_Non_TNBC$padj < 0.05, na.rm=TRUE) = 1683
length(gene_IDs_TNBC_NonTNBC)
sum(results_TNBC_Non_TNBC$padj < 0.05, na.rm=TRUE)

length(gene_IDs_TNBC)
sum(results_TNBC_df$padj < 0.05, na.rm=TRUE)

length(gene_IDs_TNBC_HER2)
sum(results_TNBC_HER2_df$padj < 0.05, na.rm=TRUE)

length(gene_IDs_NonTNBC_HER2)
sum(results_NonTNBC_HER2_df$padj < 0.05, na.rm=TRUE)

# Get GeneIDS of all genes

all_gene_IDs <- row.names(results_TNBC_Non_TNBC_df)

# Check if Gene IDs are the same for all samples
row.names(results_TNBC_Non_TNBC_df) == row.names(results_TNBC_df)

length(all_gene_IDs)

# ---------------------------------------------------------------------------------------------------------------------------------------
#
#  GO enrichment analysis
#
#  --------------------------------------------------------------------------------------------------------------------------------------

# EnrichGO TNBC vs NonTNBC
enrichGO_TNBC_NonTNBC <-
  enrichGO(
    gene = gene_IDs_TNBC_NonTNBC,
    universe = all_gene_IDs,
    OrgDb = 'org.Hs.eg.db',
    ont = "ALL",
    keyType = "ENSEMBL",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = T
  )

# EnrichGO TNBC vs Normal
enrichGO_TNBC_Normal <-
  enrichGO(
    gene = gene_IDs_TNBC,
    universe = all_gene_IDs,
    OrgDb = 'org.Hs.eg.db',
    ont = "ALL",
    keyType = "ENSEMBL",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = T
  )
# EnrichGO TNBC vs HER2
enrichGO_TNBC_HER2 <-
  enrichGO(
    gene = gene_IDs_TNBC_HER2,
    universe = all_gene_IDs,
    OrgDb = 'org.Hs.eg.db',
    ont = "ALL",
    keyType = "ENSEMBL",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = T
  )

# EnrichGO NonTNBC vs HER2
enrichGO_NonTNBC_HER2 <-
  enrichGO(
    gene = gene_IDs_NonTNBC_HER2,
    universe = all_gene_IDs,
    OrgDb = 'org.Hs.eg.db',
    ont = "ALL",
    keyType = "ENSEMBL",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = T
  )


# EnrichGO HER2 vs Normal
enrichGO_HER2_Normal <-
  enrichGO(
    gene = gene_IDs_HER2_Normal,
    universe = all_gene_IDs,
    OrgDb = 'org.Hs.eg.db',
    ont = "ALL",
    keyType = "ENSEMBL",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = T
  )



# EnrichGO NonTNBC vs Normal
enrichGO_NonTNBC_Normal <-
  enrichGO(
    gene = gene_IDs_NonTNBC_Normal,
    universe = all_gene_IDs,
    OrgDb = 'org.Hs.eg.db',
    ont = "ALL",
    keyType = "ENSEMBL",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = T
  )


# ---------------------------------------------------------------------------------------------------------------------------------------
#
#  Visualization of GO enrichment analysis
#
#  --------------------------------------------------------------------------------------------------------------------------------------
library(DOSE)
library(ggplot2)

barplot(enrichGO_TNBC_NonTNBC, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") + ggtitle("TNBC vs. NonTNBC Bar plot of enriched terms")
barplot(enrichGO_TNBC_Normal, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")+ ggtitle("TNBC vs. Normal Bar plot of enriched terms")
barplot(enrichGO_TNBC_HER2, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") + ggtitle("TNBC vs. HER2 Bar plot of enriched terms")
barplot(enrichGO_NonTNBC_HER2, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") + ggtitle("NonTNBC vs. HER2 Bar plot of enriched terms")
barplot(enrichGO_HER2_Normal, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") + ggtitle("HER2 vs. WT Bar plot of enriched terms")
barplot(enrichGO_NonTNBC_Normal, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") + ggtitle("NonTNBC vs. WT Bar plot of enriched terms")



dotplot(enrichGO_TNBC_NonTNBC, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") + ggtitle("TNBC vs. NonTNBC Enriched GO Terms")
dotplot(enrichGO_TNBC_Normal, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") + ggtitle("TNBC vs. Normal Enriched GO Terms")
dotplot(enrichGO_TNBC_HER2, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") + ggtitle("TNBC vs. HER2 Enriched GO Terms")
dotplot(enrichGO_NonTNBC_HER2, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") + ggtitle("NonTNBC vs. HER2 Enriched GO Terms")
dotplot(enrichGO_HER2_Normal, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") + ggtitle("HER2 vs. WT Enriched GO Terms")
dotplot(enrichGO_NonTNBC_Normal, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") + ggtitle("NonTNBC vs. WT Bar plot of enriched terms")


