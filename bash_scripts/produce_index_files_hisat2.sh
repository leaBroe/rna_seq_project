#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=05:00:00
#SBATCH --job-name=indexing_hisat2
#SBATCH --mail-user=lea.broennimann@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/lbroennimann/index_files_%j.o
#SBATCH --error=/data/users/lbroennimann/error_index_files_%j.e

cd /data/users/lbroennimann/rna_seq_project/map_reads

module add UHTS/Aligner/hisat/2.2.1

hisat2-build -f Homo_sapiens.GRCh38.dna.primary_assembly.fa index_file
