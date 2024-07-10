# This script is for renaming RepeatExplorer2 output files for manually checking annotion
# and extract the sampleID and number of analyzed reads into a txt file
# Mia Pei 06/19/2024

# Rename CLUSTER_TABLE.csv to the following format: SampleID_CLUSTER_TABLE.csv, etc.
cd /home/miapei/projects/def-gagnoned/miapei/BINF6999/RE2_output_20240619

for d in XAS*; do

  # Add prefix of the SampleID to the file
  cp ${d}/sorted_CLUSTER_TABLE.csv ${d}/${d}_CLUSTER_TABLE.csv

  # Extract the number of analyzed reads and print it along with the directory name
  num_reads=$(grep "Number_of_analyzed_reads" ${d}/${d}_CLUSTER_TABLE.csv | awk '{print $2}')
  echo -e "$d\t$num_reads" >> N_reads_analyzed_MiaPei.txt

  # Move all the renamed files to a new directory
  mv ${d}/${d}_CLUSTER_TABLE.csv ./RE2_mannual_annotation
done

### After all these steps, add species names manually into the second column of the txt file
### then add genome size information manually into the last column fo the txt file
