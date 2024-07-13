import pandas as pd
import os
import re
import csv

# Load the "N_reads_analyzed_MiaPei.csv" file
species_file_path = '/Users/miapei/Documents/GitHub/BINF_Research_Project/RE2_output/N_reads_analyzed_MiaPei.csv'
species_df = pd.read_csv(species_file_path, sep='\t')

# Print column names to verify
print("Columns in species_df:", species_df.columns)

# Specify the directory containing the annotated CLUSTER_TABLE.csv files
annotated_files_path = '/Users/miapei/Documents/GitHub/BINF_Research_Project/RE2_output/RE2_manual_annotation'

example_file_path = os.path.join(annotated_files_path, "annotated_XAS130_MiaPei.csv")
example_df = pd.read_csv(example_file_path, sep='\t')
print(example_df.dtypes)

# Iterate over all files in the specified directory matching the pattern
for filename in os.listdir(annotated_files_path):
    if filename.endswith('_MiaPei.csv'):

        # Extract sample ID from the filename
        sample_id = re.match(r'annotated_(.*)_MiaPei.csv', filename).group(1)

        # Find the corresponding species name, number of analyzed reads, and genome size
        species_name = species_df.loc[species_df['Sample_ID'] == sample_id, 'Species_Name'].iloc[0]

        # Load the current CSV file (the original file is seperated by tab)
        file_path = os.path.join(annotated_files_path, filename)
        current_df = pd.read_csv(file_path, sep='\t')

        # Add the species column
        current_df['Species_Name'] = species_name

        # Save the modified DataFrame back to CSV
        # quoting=csv.QUOTE_NONNUMERIC is to only quote nonnumeric values, to keep the format consistent
        current_df.to_csv(file_path, index=False, sep='\t', quoting=csv.QUOTE_NONNUMERIC)

print("All files have been updated with columns of species")

# Make a list to save all data frame from each file
dfs = []

# Iterate over all the files in the specified directory matching the pattern
for filename in os.listdir(annotated_files_path):
    if filename.endswith('_MiaPei.csv'):
        file_path = os.path.join(annotated_files_path, filename)
        # Read in the csv file (the original file is seperated by tab, use the first line as the header)
        df = pd.read_csv(file_path, sep='\t', header=0)
        dfs.append(df)

# Combine multiple DataFrame objects into one, ignore_index=True is used to avoid duplicate index values
combined_df = pd.concat(dfs, ignore_index=True)

# Save the dataframe into a csv file
combined_df.to_csv('/Users/miapei/Documents/GitHub/BINF_Research_Project/RE2_output/All_Clusters_RE2_MiaPei.csv', index=False)

print("All files have been combined as one")