#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=02:30:00
#SBATCH --job-name=quality_checks
#SBATCH --mail-user=lea.broennimann@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/courses/rnaseq_course/breastcancer_de/lea_broe/quality_checks/quality_checks_%j.o
#SBATCH --error=/data/courses/rnaseq_course/breastcancer_de/lea_broe/quality_checks/error_quality_checks_%j.e

## Change directories to where the fastq files are located
cd /data/courses/rnaseq_course/breastcancer_de/reads

## Load modules required for script commands
module load UHTS/Quality_control/fastqc/0.11.9

## Run FASTQC
fastqc -o /data/courses/rnaseq_course/breastcancer_de/lea_broe/quality_checks/ -t 6 *.gz
