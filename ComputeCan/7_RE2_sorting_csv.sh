### This script is to sort all the CLUSTER_TABLE.csv based on the second column
### Author: Mia Pei
### Last modified: July 10, 2024

cd /home/miapei/projects/def-gagnoned/miapei/BINF6999_ComputeCan/RE2_output_20240619

# Sort CLUSTER_TABLE.csv based on column B (the number of the supercluster)
for d in XAS*; do
  if [ -f "$d/CLUSTER_TABLE.csv" ]; then
    echo "Processing $d/CLUSTER_TABLE.csv"
    sorted_file=$d/sorted_CLUSTER_TABLE.csv
    head -n 7 "$d/CLUSTER_TABLE.csv" > "$sorted_file"
    tail -n +8 "$d/CLUSTER_TABLE.csv" | sort -t$'\t' -k2,2n >> "$sorted_file"
  else
    echo "File $d/CLUSTER_TABLE.csv does not exist."
  fi
done

