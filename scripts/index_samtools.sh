#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=03:00:00
#SBATCH --job-name=index_bam_files
#SBATCH --mail-user=lea.broennimann@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/lbroennimann/index_bam_files_%j.o
#SBATCH --error=/data/users/lbroennimann/index_bam_files_%j.e

cd /data/users/lbroennimann/rna_seq_project/map_reads

module add UHTS/Analysis/samtools/1.10

for i in HER21 HER22 HER23 NonTNBC1 NonTNBC2 NonTNBC3 Normal1 Normal2 Normal3 TNBC1 TNBC2 TNBC3;
do samtools index sorted_${i}.bam;
done
