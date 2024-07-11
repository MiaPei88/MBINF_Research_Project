#!/bin/bash
#SBATCH --account=def-gagnoned
#SBATCH --cpus-per-task=2
#SBATCH --time=0:30:00
#SBATCH --mem=16G
#SBATCH --output=Trimming.out

### This script is designed to trim the low quality reads and adaptors
### Author: Mia Pei
### Last modified: July 10, 2024

# Load modules
module load StdEnv/2023
module load trimmomatic/0.39
# Downlaod TruSeq3-SE.fa from Trimmomatic github and upload it to directory TrimmedSequence in computecanda
# Trim the adaptor sequence
# Change the path to the directory to save the trimmed sequences
cd /home/miapei/projects/def-gagnoned/miapei/BINF6999_ComputeCan/TrimmedSequence_20240611

# Loop through all R1 files
for r1_file in ../RawData/MiaPei_BINF/*_R1_001.fastq.gz; do

	# Save the basename of the forward reads file to a vector
	r1_name="$(basename $r1_file)"
	r1_base="${r1_name%_001.fastq.gz}"
	# Infer the corresponding R2 file name by replacing R1 with R2
	r2_file="${r1_file/_R1_/_R2_}"
	# Save the basename of the backward reads file to a vector
	r2_name="$(basename $r2_file)"
	r2_base="${r2_name%_001.fastq.gz}"
	echo "Processing pair: $r1_name and $r2_name ..."
	# Running Trimmomatic
	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -phred33 \
	$r1_file $r2_file \
	${r1_base}_clean.fastq.gz ${r1_base}_unpaired.fastq.gz \
	${r2_base}_clean.fastq.gz ${r2_base}_unpaired.fastq.gz \
	ILLUMINACLIP:./TruSeq3-SE.fa:2:40:15 SLIDINGWINDOW:4:15 CROP:245 HEADCROP:10 MINLEN:90;

done
