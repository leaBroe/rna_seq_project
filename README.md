# RNA Sequencing Project HS2022 

This is the code used for the RNA sequencing project where the goal was to detect differentially expressed genes in breast cancer subtypes based on data of the paper "Transcriptomic landscape of breast cancers through mRNA sequencing" by Eswaran et al. (2012, [https://doi.org/10.1038/srep00264](https://doi.org/10.1038/srep00264)). 

The script for step 1-3 of the analysis can be found in the bash_scripts folder, the subsequent analysis with R in the rna_seq_project/rna_seq_analysis_R/ folder.

## Steps of the analysis:

1. ***Quality Control***: For the quality control, the rna_seq_project/bash_scripts/quality_checks.sh as well as the rna_seq_project/bash_scripts/trimmomatic.sh script was used 

2. ***Mapping reads to the reference genome***: The reference genome used for this step ("Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz") can be found at [https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/](https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/), the annotation file "Homo_sapiens.GRCh38.108.gtf.gz" at [https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/](https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/).  
2.1 **Production of Index files using Hisat2**: rna_seq_project/bash_scripts/produce_index_files_hisat2.sh  
2.2 **Mapping of reads to reference genome for each sample separately**: rna_seq_project/bash_scripts/map_reads_to_genome_hisat2.sh
2.3 **Convert the resulting sam files to bam format using Samtools**: rna_seq_project/bash_scripts/convert_sam_to_bam_samtools.sh
2.4 **Sort the bam files by genomic coordinates using Samtools**: rna_seq_project/bash_scripts/sort_bam_files_samtools.sh
2.5 **Index the coordinate sorted bam files using Samtools**: rna_seq_project/bash_scripts/index_samtools.sh

3. ***Counting the number of genes using FeatureCounts***: rna_seq_project/bash_scripts/feature_counts.sh
Since feature counts was run on each sample separately, the count matrix had to be created using the bioinfokit python module (rna_seq_project/bash_scripts/generate_count_matrix.py). 

4. ***Exploratory data and differential expression analysis using DESeq2***: rna_seq_project/rna_seq_analysis_R/diff_expression_analysis_final.R

5. ***Overrepresentation analysis using ClusterProfiler***: rna_seq_project/rna_seq_analysis_R/overrepresentation_analysis_final.R
