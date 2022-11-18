#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=20:00:00
#SBATCH --job-name=feature_counts
#SBATCH --mail-user=lea.broennimann@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/lbroennimann/feature_counts2_%j.o
#SBATCH --error=/data/users/lbroennimann/feature_counts2_%j.e
#SBATCH --array=0-11%12

cd /data/users/lbroennimann/rna_seq_project/map_reads

module load UHTS/Analysis/subread/2.0.1;

samples=('HER21' 'HER22' 'HER23' 'NonTNBC1' 'NonTNBC2' 'NonTNBC3' 'Normal1' 'Normal2' 'Normal3' 'TNBC1' 'TNBC2' 'TNBC3')

featureCounts -p -C -s 0 -Q 25 -T 4 --tmpDir $SCRATCH -a Homo_sapiens.GRCh38.108.gtf -G Homo_sapiens.GRCh38.dna.primary_assembly.fa -o ${samples[$SLURM_ARRAY_TASK_ID]}_featureCounts_output.txt /data/users/lbroennimann/rna_seq_project/map_reads/sorted_files/sorted_${samples[$SLURM_ARRAY_TASK_ID]}.bam



