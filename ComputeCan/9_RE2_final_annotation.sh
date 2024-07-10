### This script is to copy the content of the 5th column "Automatic_annotation" to 7th column "Final_annotation" in all CLUSTER_TABLE.csv files
### In order to save time copy paste right automated annotation manually.
### After this step, we can manually change the specific row on column "Final_annotation" when the automatic annotation is wrong.
### Wrote by Mia Pei
### 2024/06/24

cd /home/miapei/projects/def-gagnoned/miapei/BINF6999/RE2_output_20240619/RE2_manual_annotation

for f in *.csv; do

  echo "Processing ${f}" # Show the process
  annotated_file="annotated_${f}" # Make a new file to save the output
  head -n 7 "${f}" > "${annotated_file}" # Copy paste the first 7 lines to the new file

  # Copy paste the 5th column to the 7th column and append to the new file
  tail -n +8 "$f" | awk -F'\t' 'BEGIN {OFS="\t"} {$7=$5; print}' >> "$annotated_file"

done
