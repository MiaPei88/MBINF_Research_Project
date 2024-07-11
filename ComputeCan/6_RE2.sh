#!/bin/bash
#SBATCH --account=def-gagnoned
#SBATCH --mail-user=qpei@uoguelph.ca
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=16
#SBATCH --time=96:00:00
#SBATCH --mem=128G
#SBATCH --output=RE2_%a.out
#SBATCH --array=0-14


### This script is designed to analyze the repeat profile of all samples using RepeatExplore2
### Author: Mia Pei
### Last modified: July 10, 2024

cd /home/miapei/projects/def-gagnoned/miapei/BINF6999/MergedSequence_20240619
FILES=(*.fasta)
filename=${FILES[$SLURM_ARRAY_TASK_ID]}

# Trim the suffix of the file for directory name
dir_name="${filename%_S*_merged.fasta}"

# Create output directory for each job
output_dir="/home/miapei/projects/def-gagnoned/miapei/BINF6999/RE2_output_20240619/${dir_name}"
mkdir -p ${output_dir}

# Load required module
module load apptainer

# Unpack the sandbox images to the temperary directory for slurm this job
tar -xvf /home/miapei/projects/def-gagnoned/miapei/repex_tarean.tar.gz -C $SLURM_TMPDIR

# Back to the main directory for the project
cd /home/miapei/projects/def-gagnoned/miapei/BINF6999

# Run RepeatExplorer2
apptainer exec -e --bind ${PWD}:/data/ $SLURM_TMPDIR/repex_tarean \
seqclust  -p -v /data/RE2_output_20240619/${dir_name} \
-c 16 -tax VIRIDIPLANTAE3.0 \
-C \
/data/MergedSequence_20240619/${filename}
