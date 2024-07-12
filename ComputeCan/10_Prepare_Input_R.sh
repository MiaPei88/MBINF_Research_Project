### This script is designed to prepare proper input format for R
### Author: Mia Pei
### Last Modified: July 11, 2024

# Remove the first block of text from each annotated file
cd /home/miapei/projects/def-gagnoned/miapei/BINF6999_ComputeCan/RE2_output_20240619/RE2_manual_annotation
for f in annotated_*.csv;
do (tail -n +7 $f > ${f%_CLUSTER_TABLE.csv}_MiaPei.csv);
done

# Then remove all the notes that have been added in the first two columns 
# (these notes are to show which clusters' annotations are changed)
