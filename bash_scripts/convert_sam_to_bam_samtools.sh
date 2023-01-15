#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=42000M
#SBATCH --time=05:00:00
#SBATCH --job-name=convert_sam_to_bam
#SBATCH --mail-user=lea.broennimann@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/lbroennimann/convert_to_bam_%j.o
#SBATCH --error=/data/users/lbroennimann/convert_to_bam_%j.e

cd /data/users/lbroennimann/rna_seq_project/map_reads

module add UHTS/Analysis/samtools/1.10

for i in HER21 HER22 HER23 NonTNBC1 NonTNBC2 NonTNBC3 Normal1 Normal2 Normal3 TNBC1 TNBC2 TNBC3;
do samtools view -hbS mapped_reads_${i}.sam > mapped_reads_${i}.bam;
done


