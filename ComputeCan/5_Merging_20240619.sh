#!/bin/bash
#SBATCH --account=def-gagnoned
#SBATCH --cpus-per-task=2
#SBATCH --time=00:15:00
#SBATCH --mem=16G
#SBATCH --output=Merging_20240619.out
#SBATCH --mail-user=qpei@uoguelph.ca
#SBATCH --mail-type=ALL

# Modified on June 19, 2024
# This script is to merge all the cleaned fastq files into a paired file

# Load required modules
module load StdEnv/2023 nixpkgs/16.09  gcc/7.3.0
module load seqtk/1.3

cd /home/miapei/projects/def-gagnoned/miapei/BINF6999/TrimmedSequence_20240611

# Loop through all R1 files
for r1_file in *_R1_clean.fastq.gz
do

  # Extract the nmae of the file
  r1_name="$(basename $r1_file)"
  # Infer the corresponding R2 file name by replacing R1 with R2
  r2_file="${r1_file/_R1_/_R2_}"
  # Print the process
  echo "Merging $r1_file and $r2_file"
  # Create a new name for the merged file
  merged_name="${r1_name%_L001_R1_clean.fastq.gz}_merged"
  # Merge r1 and r2 files and convert to a fasta file
  seqtk mergepe $r1_file $r2_file > ../MergedSequence_20240619/${merged_name}.fastq
  seqtk seq -A ../MergedSequence_20240619/${merged_name}.fastq > ../MergedSequence_20240619/${merged_name}.fasta

done


