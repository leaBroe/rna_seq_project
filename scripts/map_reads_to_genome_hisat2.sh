#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=17:00:00
#SBATCH --job-name=mapping_reads
#SBATCH --mail-user=lea.broennimann@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/lbroennimann/mapping_files_%j.o
#SBATCH --error=/data/users/lbroennimann/error_mapping_reads_%j.e

cd /data/users/lbroennimann/rna_seq_project/map_reads

module add UHTS/Aligner/hisat/2.2.1

for i in HER21 HER22 HER23 NonTNBC1 NonTNBC2 NonTNBC3 Normal1 Normal2 Normal3 TNBC1 TNBC2 TNBC3;
do hisat2 -x index_file -1 /data/courses/rnaseq_course/breastcancer_de/reads/${i}_R1.fastq.gz -2 /data/courses/rnaseq_course/breastcancer_de/reads/${i}_R2.fastq.gz -S mapped_reads_${i}.sam -p 4;
done




