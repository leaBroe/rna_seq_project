# Generates count matrix from featureCounts output if run bam files separately
# run this Python code from a folder where all files are present (featureCounts_output.txt files)!

from bioinfokit.analys import HtsAna

# Creates count matrix
HtsAna.merge_featureCount()






