if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)

# Get rid of NAs 
results_TNBC_Non_TNBC_df_omit_NA <- na.omit(results_TNBC_Non_TNBC_df)

results_TNBC_df_omit_NAs <- na.omit(results_TNBC_df)

results_TNBC_HER2_df_omit_NAs <- na.omit(results_TNBC_HER2_df)

# Select all Gene IDS for TNBC vs NonTNBC DE genes
gene_IDs_TNBC_NonTNBC <- row.names(results_TNBC_Non_TNBC_df_omit_NA[results_TNBC_Non_TNBC_df_omit_NA$padj<0.05,])

# Select all Gene IDS for TNBC vs normal DE genes
gene_IDs_TNBC <- row.names(results_TNBC_df_omit_NAs[results_TNBC_df_omit_NAs$padj<0.05,])

# Select all Gene IDs for TNBC vs HER2 DE genes
gene_IDs_TNBC_HER2 <- row.names(results_TNBC_HER2_df_omit_NAs[results_TNBC_HER2_df_omit_NAs$padj<0.05,])

# Check with sum(results_TNBC_Non_TNBC$padj < 0.05, na.rm=TRUE) = 1683
length(gene_IDs_TNBC_NonTNBC)
sum(results_TNBC_Non_TNBC$padj < 0.05, na.rm=TRUE)

length(gene_IDs_TNBC)
sum(results_TNBC_df$padj < 0.05, na.rm=TRUE)

length(gene_IDs_TNBC_HER2)
sum(results_TNBC_HER2_df$padj < 0.05, na.rm=TRUE)

all_gene_IDs <- row.names(results_TNBC_Non_TNBC_df)
row.names(results_TNBC_Non_TNBC_df) == row.names(results_TNBC_df)

length(all_gene_IDs)

# ---------------------------------------------------------------------------------------------------------------------------------------
#
#  GO enrichment analysis
#
#  --------------------------------------------------------------------------------------------------------------------------------------

# EnrichGO TNBC vs NonTNBC
enrichGO_TNBC_NonTNBC <- enrichGO(gene=gene_IDs_TNBC_NonTNBC, universe=all_gene_IDs, OrgDb='org.Hs.eg.db', ont = "ALL", keyType = "ENSEMBL")

# EnrichGO TNBC vs Normal
enrichGO_TNBC_Normal <- enrichGO(gene=gene_IDs_TNBC, universe=all_gene_IDs, OrgDb='org.Hs.eg.db', ont = "ALL", keyType = "ENSEMBL")

# EnrichGO TNBC vs HER2
enrichGO_TNBC_HER2 <- enrichGO(gene=gene_IDs_TNBC_HER2, universe=all_gene_IDs, OrgDb='org.Hs.eg.db', ont = "ALL", keyType = "ENSEMBL")
enrichGO_TNBC_HER2_BP <- enrichGO(gene=gene_IDs_TNBC_HER2, universe=all_gene_IDs, OrgDb='org.Hs.eg.db', ont = "BP", keyType = "ENSEMBL")


# ---------------------------------------------------------------------------------------------------------------------------------------
#
#  Visualization of GO enrichment analysis
#
#  --------------------------------------------------------------------------------------------------------------------------------------
library(DOSE)


barplot(enrichGO_TNBC_NonTNBC, showCategory = 11)
barplot(enrichGO_TNBC_Normal, showCategory = 11)
barplot(enrichGO_TNBC_HER2, showCategory = 20) + ggtitle("TNBC vs. HER2 Bar plot of enriched terms")
barplot(enrichGO_TNBC_HER2_BP, showCategory = 11)

dotplot(enrichGO_TNBC_NonTNBC) + ggtitle("TNBC vs. NonTNBC Enriched GO")
dotplot(enrichGO_TNBC_Normal) + ggtitle("TNBC vs. Normal Enriched GO")
dotplot(enrichGO_TNBC_HER2, showCategory = 20) + ggtitle("TNBC vs. HER2 Enriched GO")
dotplot(enrichGO_TNBC_HER2_BP)




