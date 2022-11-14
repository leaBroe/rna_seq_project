#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=20:00:00
#SBATCH --job-name=index_bam_files
#SBATCH --mail-user=lea.broennimann@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/lbroennimann/feature_counts_%j.o
#SBATCH --error=/data/users/lbroennimann/feature_counts_%j.e

cd /data/users/lbroennimann/rna_seq_project/map_reads

module load UHTS/Analysis/subread/2.0.1;

for i in HER21 HER22 HER23 NonTNBC1 NonTNBC2 NonTNBC3 Normal1 Normal2 Normal3 TNBC1 TNBC2 TNBC3;
do featureCounts -p -T 4 -a Homo_sapiens.GRCh38.108.gtf -G Homo_sapiens.GRCh38.dna.primary_assembly.fa -o ${i}_featureCounts_output.txt sorted_${i}.bam;
done


