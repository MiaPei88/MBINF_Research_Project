### Analysing repeat profiles from RepeatExplorer2 for genome size correlation
### Author: Mia Pei
### Last day modified: July 11, 2024

# Load required packages
if (!require(vegan)) install.packages("vegan")
if (!require(reshape2)) install.packages("reshape2")
library(dplyr)
library(tidyverse)
library(vegan)
library(ggplot2)
library(reshape2)
library(stringr)

############# Repeat Quantification #############
#################################################

# Set the current working directory
setwd("~/Documents/GitHub/BINF_Research_Project/RE2_output")

# Read the combined cluster table
All_Eggplant_Repeat <- read.csv("All_Clusters_RE2_MiaPei.csv")

# Read the table with number of analyzed reads and genome size
N_reads_analyzed_GS <- read.csv("N_reads_analyzed_MiaPei.csv", sep = "\t")

# Check the number of species
unique(All_Eggplant_Repeat$Species_Name)

# Calculate the total 'Size_adjusted' for each species for rows containing 'organelle'
number_of_organelle_reads_by_species <- All_Eggplant_Repeat %>%
  filter(grepl("organelle", Final_annotation)) %>%
  group_by(Species_Name) %>%
  summarise(Total_organelle_size_adjusted = sum(Size_adjusted, na.rm = TRUE))

# Append the column of total_organelle size
N_reads_analyzed_GS <- N_reads_analyzed_GS %>%
  left_join(number_of_organelle_reads_by_species, by = "Species_Name")

# Add a new column for the number of reads representative of nuclear sequences as:
# the number of analyzed reads - the number of reads annotated as organelle
# (Because in the output from RE2 of this project do not have contamination, 
# we don't need to minus the number of reads of contamination here.)
N_reads_analyzed_GS$Number_of_reads_nuclear <- N_reads_analyzed_GS[,3] - N_reads_analyzed_GS[,5]

# Delete rows that are annotated as organelle
All_Eggplant_Repeat <- All_Eggplant_Repeat[-which(str_detect(All_Eggplant_Repeat$Final_annotation,"organelle")),]

# Check the repeat classification
unique(All_Eggplant_Repeat$Final_annotation)

# Aggregate counts of each repeat classification in the annotation for each species 
All_Eggplant_Repeat_ag <- aggregate(Size_adjusted~Species_Name+Final_annotation, data=All_Eggplant_Repeat, FUN=sum) 

# Calculate the genomic proportion of individual clusters by dividing the size 
# (number of reads, column called Size_adjusted) of the clusters by Number_of_reads_nuclear
All_Eggplant_Repeat_ag <- left_join(All_Eggplant_Repeat_ag, N_reads_analyzed_GS %>%
                                   select(Species_Name, Genome_Size.mean2C., Number_of_reads_nuclear), 
                                 by = "Species_Name") %>%
  mutate(Genomic_proportion = (Size_adjusted/Number_of_reads_nuclear) * 100) %>%
  mutate(Final_annotation = case_when(
    Final_annotation == "All"  ~ "Other",
    Final_annotation == "All/repeat" | Final_annotation == "All/repeat/mobile_element" | Final_annotation == "All/repeat/mobile_element/Class_I" ~ "All/repeat_unclassified",
    Final_annotation == "All/repeat/mobile_element/Class_I/LTR" ~ "All/repeat/mobile_element/Class_I/LTR/nonspecific",
    Final_annotation == "All/repeat/mobile_element/Class_I/LTR/Ty3_gypsy" ~ "All/repeat/mobile_element/Class_I/LTR/Ty3_gypsy/nonspecific",
    Final_annotation == "All/repeat/rDNA/45S_rDNA" ~ "All/repeat/rDNA/45S_rDNA/nonspecific",
    TRUE ~ Final_annotation
  ))


# Check the repeat classification after further annotation
unique(All_Eggplant_Repeat_ag$Final_annotation)


############### Data Exploration ################
#################################################

# Make a new table for input to vegan, with column names as repeat annotation, row names as species names
# values as reads
All_Eggplant_Repeat_ag_reads <- dcast(All_Eggplant_Repeat_ag, Species_Name ~ Final_annotation, value.var = "Size_adjusted", fun.aggregate = sum)
row.names(All_Eggplant_Repeat_ag_reads) <- All_Eggplant_Repeat_ag_reads$Species_Name
All_Eggplant_Repeat_ag_reads <- All_Eggplant_Repeat_ag_reads[,-1]

# Make a new table for input to vegan, with column names as repeat annotation, row names as species names
# values as genomic proportion
All_Eggplant_Repeat_ag_pct <- dcast(All_Eggplant_Repeat_ag, Species_Name ~ Final_annotation, value.var = "Genomic_proportion", fun.aggregate = sum)
row.names(All_Eggplant_Repeat_ag_pct) <- All_Eggplant_Repeat_ag_pct$Species_Name
All_Eggplant_Repeat_ag_pct <- All_Eggplant_Repeat_ag_pct[,-1]


# Replace all the NA values as 0
All_Eggplant_Repeat_ag_reads[is.na(All_Eggplant_Repeat_ag_reads)] <- 0
All_Eggplant_Repeat_ag_pct[is.na(All_Eggplant_Repeat_ag_pct)] <- 0


# 
All_Eggplant_Repeat_ag_pct %>%
  as.data.frame() %>%
  rownames_to_column("Species_Name") %>%
  pivot_longer(-c(Species_Name), names_to = "Repeat_type", values_to = "Genomic_proportion") %>%
  ggplot(aes(x=Repeat_type, y=Species_Name, fill=Genomic_proportion)) + 
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  scale_fill_viridis_c()


# Let's try plotting a stacked bar chart for each species and the abundances of different repeats
# Repeat % of all repeat types
ggplot(All_Eggplant_Repeat_ag, aes(fill=Final_annotation, y=Size_adjusted, x=Species_Name)) + 
  theme(axis.text.x = element_text(angle = 60, size= 10, hjust = 1)) + 
  theme(legend.text = element_text(size=6)) + 
  xlab("Species") +
  ylab("Percent of all repeats") + 
  geom_bar(position="fill", stat="identity")


# Absolute repeat number
ggplot(All_Eggplant_Repeat_ag, aes(fill=Final_annotation, y=Size_adjusted, x=Species_Name)) + 
  theme(axis.text.x = element_text(angle = 90, size= 10)) + 
  theme(legend.text = element_text(size=6)) + 
  ylab("Number of repeat reads") + 
  geom_bar(stat="identity")


# Repeat % of whole genome
ggplot(All_Eggplant_Repeat_ag, aes(fill=Final_annotation, y=Genomic_proportion, x=Species_Name)) + 
  theme(axis.text.x = element_text(angle = 90, size= 10)) + 
  theme(legend.text = element_text(size=6)) + 
  ylab("Percentage occupied within the genome") + 
  ylim(0,50) + 
  geom_bar(position="stack", stat="identity")


##### Total genome proportion and diversity #####
#################################################
### Calculate the total percentage of the genome occupied by repeats in each species
Genome_occupancy <- aggregate(Genomic_proportion~Species_Name, data=All_Eggplant_Repeat_ag, FUN=sum) 

### Calculate the Shannon's diversity index
Shannon_Eggplant_Repeat <- diversity(All_Eggplant_Repeat_ag_reads, index = "shannon")
# Compare Shannon indices
Shannon_df <- Shannon_Eggplant_Repeat %>%
  as.data.frame() %>%
  rownames_to_column("Species_Name")%>%
  rename(Shannon_Index = names(.)[2]) # Change the 2nd column name to what we want

### Richness (Mehinick's Index) - Closer to 1 = higher richness
n <- apply(All_Eggplant_Repeat_ag_reads > 0, 1, sum) # number of repeat types for each species
N <- apply(All_Eggplant_Repeat_ag_reads, 1, sum) # total number of reads for each species
# Calculate the Mehinick's Index
MI_df <- n/sqrt(N)
MI_df <- MI_df %>%
  as.data.frame() %>%
  rownames_to_column("Species_Name")%>%
  rename(Mehinicks_index = names(.)[2]) # Change the 2nd column name to what we want

#### Subset the genome size information
GS_df <- subset(N_reads_analyzed_GS[,c(2,4)])
colnames(GS_df)[2] <- "Genome_Size_gbp"
# convert pg to Gbp
GS_df$Genome_Size_gbp <- GS_df$Genome_Size_gbp*0.978

### Combine the data frames of GS, richness, diversity together with the aggregated repeat profile percentage
Cmbd_All_Eggplant_Repeat_ag_pct <- All_Eggplant_Repeat_ag_pct %>%
  rownames_to_column("Species_Name")%>%
  left_join(Genome_occupancy, by = "Species_Name") %>%
  left_join(Shannon_df, by = "Species_Name") %>%
  left_join(MI_df, by = "Species_Name") %>%
  left_join(GS_df, by = "Species_Name")
row.names(Cmbd_All_Eggplant_Repeat_ag_pct) <- Cmbd_All_Eggplant_Repeat_ag_pct$Species_Name 
Cmbd_All_Eggplant_Repeat_ag_pct <- subset(Cmbd_All_Eggplant_Repeat_ag_pct, select = -Species_Name)

##~~~~~~~~~~Exploratory plotting~~~~~~~~~~## 

# Richness and genome size
plot(Cmbd_All_Eggplant_Repeat_ag_pct$Mehinicks_index ~ Cmbd_All_Eggplant_Repeat_ag_pct$Genome_Size_gbp, ylab = "Repeat type richness (Mehinicks Index)", xlab = "Genome size (Gbp)")

ggplot(Cmbd_All_Eggplant_Repeat_ag_pct, aes(x=as.numeric(Genome_Size_gbp), y= Mehinicks_index)) + 
  geom_point() +
  geom_smooth(method="loess", span = 1) +
  labs(y="Repeat type richness (Mehinicks Index)", x = "Genome size (Gbp)") + 
  theme(text=element_text(size = 14),axis.text = element_text(size = 12))


# Genomic proportion of repeats and genome size
plot(Cmbd_All_Eggplant_Repeat_ag_pct$Genomic_proportion ~ Cmbd_All_Eggplant_Repeat_ag_pct$Genome_Size_gbp, ylab = "Genomic proportion of repeats (%)", xlab="Genome size (Gbp)")

ggplot(Diploid_Eggplant_Repeat_ag_pct, aes(x=as.numeric(Genome_Size_gbp), y= Genomic_proportion)) + 
  geom_point() +
  geom_smooth(method = "loess", span = 1.5) +
  labs(y = "Genomic proportion of repeats (%)", x = "Genome size (Gbp)") + 
  theme(text=element_text(size = 14),axis.text = element_text(size = 12))

# Diversity and genome size
plot(Diploid_Eggplant_Repeat_ag_pct$Shannon_Index~Diploid_Eggplant_Repeat_ag_pct$Genome_Size_gbp, ylab = "Repeat diversity (Shannon Index)", xlab="Genome size (Gbp)")

ggplot(Diploid_Eggplant_Repeat_ag_pct, aes(x=as.numeric(Genome_Size_gbp), y= Shannon_Index)) + 
  geom_point() +
  geom_smooth(method="loess", span = 1.5) +
  labs(y="Repeat type diversity (Shannon Index)", x = "Genome size (Gbp)") + 
  theme(text=element_text(size = 14),axis.text = element_text(size = 12))
