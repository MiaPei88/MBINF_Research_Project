import pandas as pd
import os
import re

# Load the "N_reads_analyzed_MiaPei.csv" file
species_file_path = '~/Documents/GitHub/BINF_Research_Project/RE2_output/N_reads_analyzed_MiaPei.csv'
species_df = pd.read_csv(species_file_path)

# Specify the directory containing the annotated CLUSTER_TABLE.csv files
annotated_files_path = '~/Documents/GitHub/BINF_Research_Project/RE2_output/RE2_manual_annotation'

# Iterate over all files in the specified directory matching the pattern
for filename in os.listdir(annotated_files_path):
    if filename.endswith('_MiaPei.csv'):

        # Extract sample ID from the filename
        sample_id = re.match(r'annotated_(.*)_MiaPei.csv', filename).group(1)

        # Find the corresponding species name, number of analyzed reads, and genome size
        species_name = species_df.loc[species_df['Sample_Id'] == sample_id, 'Species_Name'].iloc[0]
        n_reads = species_df.loc[species_df['Sample_Id'] == sample_id, 'Number_of_analyzed_reads'].iloc[0]
        genome_size = species_df.loc[species_df['Sample_Id'] == sample_id, 'Genome_Size(mean2C)'].iloc[0]

        # Load the current CSV file
        file_path = os.path.join(annotated_files_path, filename)
        current_df = pd.read_csv(file_path)

        # Add the species column
        current_df['Species_Name'] = species_name
        current_df['N_Reads'] = n_reads
        current_df['Genome_Size'] = genome_size

        # Save the modified DataFrame back to CSV
        current_df.to_csv(file_path, index=False)

print("All files have been updated with columns of species, n_reads and genome size.")