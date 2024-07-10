#!/bin/bash
#SBATCH --account=def-gagnoned
#SBATCH --time=0:10:00
#SBATCH --mem=8G
#SBATCH --output=seqkit_stats.out

# Load module
module load StdEnv/2023
module load seqkit/2.5.1

cd /home/miapei/projects/def-gagnoned/miapei/BINF6999_ComputeCan/TrimmedSequence_20240611
# Check statistics of fastq files:
seqkit stats ./*fastq.gz
