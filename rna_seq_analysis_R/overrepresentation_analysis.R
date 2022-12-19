if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)

results_TNBC_Non_TNBC_df_omit_NA <- na.omit(results_TNBC_Non_TNBC_df)

# Select all Gene IDS for TNBC vcNonTNBC DE genes
gene_IDs_TNBC_NonTNBC <- row.names(results_TNBC_Non_TNBC_df_omit_NA[results_TNBC_Non_TNBC_df_omit_NA$padj<0.05,])

# Check with sum(results_TNBC_Non_TNBC$padj < 0.05, na.rm=TRUE) = 1683
length(gene_IDs_TNBC_NonTNBC)

all_gene_IDs <- row.names(results_TNBC_Non_TNBC_df)

length(all_gene_IDs)

# ---------------------------------------------------------------------------------------------------------------------------------------
#
#  GO enrichment analysis
#
#  --------------------------------------------------------------------------------------------------------------------------------------

enrichGO_TNBC_NonTNBC <- enrichGO(gene=gene_IDs_TNBC_NonTNBC, universe=all_gene_IDs, OrgDb='org.Hs.eg.db', ont = "ALL", keyType = "ENSEMBL")





# ---------------------------------------------------------------------------------------------------------------------------------------
#
#  Visualization of GO enrichment analysis
#
#  --------------------------------------------------------------------------------------------------------------------------------------
library(DOSE)


barplot(enrichGO_TNBC_NonTNBC, showCategory = 11)

dotplot(enrichGO_TNBC_NonTNBC)





