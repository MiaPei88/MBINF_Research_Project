### Analysing repeat profiles from RepeatExplorer2 for genome size correlation
### Author: Mia Pei
### Last day modified: July 11, 2024

##### Repeat quantification #####
# Load required packages
library(dplyr)
library(tidyverse)

# Set the current working directory
setwd("~/Documents/GitHub/BINF_Research_Project/RE2_output")

# Read the combined cluster table
All_Eggplant_Repeat <- read.csv("All_Clusters_RE2_MiaPei.csv")

# Read the table with number of analyzed reads and genome size
N_reads_analyzed_GS <- read.csv("N_reads_analyzed_MiaPei.csv", sep = "\t")

# Check the number of species
unique(All_Eggplant_Repeat$Species_Name)

# Calculate the total 'Size_adjusted' for each species for rows containing 'organelle'
number_of_reads_by_species <- All_Eggplant_Repeat %>%
  filter(grepl("organelle", Final_annotation)) %>%
  group_by(Species_Name) %>%
  summarise(Total_organelle_size_adjusted = sum(Size_adjusted, na.rm = TRUE))

# Append the column of total_organelle size
N_reads_analyzed_GS <- N_reads_analyzed_GS %>%
  left_join(number_of_reads_by_species, by = "Species_Name")

# Add a new column for the number of reads representative of nuclear sequences as:
# the number of analyzed reads - the number of reads annotated as organelle
# (Because in the output from RE2 of this project do not have contamination, 
# we don't need to minus the number of reads of contamination here.)
N_reads_analyzed_GS$Number_of_reads_nuclear <- N_reads_analyzed_GS[,3] - N_reads_analyzed_GS[,5]

# Delete rows that are annotated as organelle
All_Eggplant_Repeat <- All_Eggplant_Repeat[-which(str_detect(All_Eggplant_Repeat$Final_annotation,"organelle")),]


All_Eggplant_Repeat_ag <- aggregate(Size_adjusted~Species_Name+Final_annotation, data=All_Eggplant_Repeat, FUN=sum) 

# Calculate the genomic proportion of individual clusters by dividing the size 
# (number of reads; column 4 in All clusters df) of the clusters by Number_of_reads_nuclear
All_Eggplant_Repeat_ag <- left_join(All_Eggplant_Repeat_ag, N_reads_analyzed_GS %>%
                                   select(Species_Name, Genome_Size.mean2C., Number_of_reads_nuclear), 
                                 by = "Species_Name") %>%
  mutate(Genomic_proportion = Size_adjusted/Number_of_reads_nuclear)

unique(All_Eggplant_Repeat$Final_annotation)
