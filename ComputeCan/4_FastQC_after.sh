#!/bin/bash
#SBATCH --account=def-gagnoned
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --output=FastQC_after.out

### This script is designed to check the quality of reads after trimming
### Author: Mia Pei
### Last modified: July 10, 2024

# Load required module
module load StdEnv/2023
module load fastqc/0.12.1

# Quality check of reads using FastQC
mkdir ./FastQC_output_after
fastqc -o ./FastQC_output_after ./TrimmedSequence_20240611/*clean.fastq.gz

