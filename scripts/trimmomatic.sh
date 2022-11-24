#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=05:00:00
#SBATCH --job-name=trimmomatic
#SBATCH --mail-user=lea.broennimann@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/lbroennimann/trimmomatic_%j.o
#SBATCH --error=/data/users/lbroennimann/trimmomatic_%j.e

# perform trimmomatic on HER22 for quality improvement

module load UHTS/Analysis/trimmomatic/0.36;
module load UHTS/Quality_control/fastqc/0.11.9;

cd /data/courses/rnaseq_course/breastcancer_de/reads

trimmomatic PE -threads 1 HER22_R1.fastq.gz HER22_R2.fastq.gz /data/users/lbroennimann/rna_seq_project/quality_checks/trimmomatic/HER22_R1_trim.fastq.gz /data/users/lbroennimann/rna_seq_project/quality_checks/trimmomatic/HER22_R1_unpaired.fastq.gz /data/users/lbroennimann/rna_seq_project/quality_checks/trimmomatic/HER22_R2_trim.fastq.gz /data/users/lbroennimann/rna_seq_project/quality_checks/trimmomatic/HER22_R2_unpaired.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40

cd /data/users/lbroennimann/rna_seq_project/quality_checks/trimmomatic

fastqc -t 1 HER22_*_trim.fastq.gz