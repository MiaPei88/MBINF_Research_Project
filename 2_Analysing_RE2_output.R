### Analysing repeat profiles from RepeatExplorer2 for genome size correlation
### Author: Mia Pei
### Last day modified: July 11, 2024

# Load required packages
#if (!require(vegan)) install.packages("vegan")
#if (!require(reshape2)) install.packages("reshape2")
#if (!require(car)) install.packages("car")
#if (!require(caper)) install.packages("caper")
#if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("ggtree")
#if (!require(phangorn)) install.packages("phangorn")
#if (!require(paletteer)) install.packages("paletteer")

library(dplyr) #version 1.1.4
library(tidyverse) #version 2.0.0
library(vegan) #version 2.6.6.1
library(ggplot2) #version 3.5.1
library(reshape2) #version 1.4.4
library(stringr) #version 1.5.1
library(car) #version 3.1.2
library(caper) #version 1.0.3
library(ape) #version 5.8
library(muscle) #version 3.46.0
library(BiocManager) #version 1.30.23
library(ggtree) #version 3.12.0
library(phytools) #version 2.3.0
library(phangorn) #version 2.11.1
library(RColorBrewer)
library(paletteer)
library(grid)

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

sub_N_reads_analyzed_GS <- N_reads_analyzed_GS %>%
  select(Species_Name, `Genome_Size.mean2C.`, Number_of_reads_nuclear)

# Calculate the genomic proportion of individual clusters by dividing the size 
# (number of reads, column called Size_adjusted) of the clusters by Number_of_reads_nuclear
All_Eggplant_Repeat_ag <- left_join(All_Eggplant_Repeat_ag, N_reads_analyzed_GS %>%
                                   select(Species_Name, Genome_Size.mean2C.,
                                          Number_of_reads_nuclear), 
                                 by = "Species_Name") %>%
  mutate(Genomic_proportion = (Size_adjusted/Number_of_reads_nuclear) * 100) %>%
  mutate(Final_annotation = case_when(
    Final_annotation == "All"  ~ "Other",
    TRUE ~ Final_annotation
  ))


# Check the repeat classification after further annotation
unique(All_Eggplant_Repeat_ag$Final_annotation)


############### Data Exploration ################
#################################################

# Make a new table for input to vegan, with column names as repeat annotation, row names as species names, values as reads
All_Eggplant_Repeat_ag_reads <- dcast(All_Eggplant_Repeat_ag, Species_Name ~ Final_annotation, value.var = "Size_adjusted", fun.aggregate = sum)
row.names(All_Eggplant_Repeat_ag_reads) <- All_Eggplant_Repeat_ag_reads$Species_Name
All_Eggplant_Repeat_ag_reads <- All_Eggplant_Repeat_ag_reads[,-1]

# Make a new table for input to vegan, with column names as repeat annotation, row names as species names, values as genomic proportion
All_Eggplant_Repeat_ag_pct <- dcast(All_Eggplant_Repeat_ag, Species_Name ~ Final_annotation, value.var = "Genomic_proportion", fun.aggregate = sum)
row.names(All_Eggplant_Repeat_ag_pct) <- All_Eggplant_Repeat_ag_pct$Species_Name
All_Eggplant_Repeat_ag_pct <- All_Eggplant_Repeat_ag_pct[,-1]

# Replace all the NA values as 0
All_Eggplant_Repeat_ag_reads[is.na(All_Eggplant_Repeat_ag_reads)] <- 0
All_Eggplant_Repeat_ag_pct[is.na(All_Eggplant_Repeat_ag_pct)] <- 0

# Calculate the total percentage of the genome occupied by repeats in each species
Genome_occupancy <- aggregate(Genomic_proportion~Species_Name, data=All_Eggplant_Repeat_ag, FUN=sum)
Genome_occupancy <- Genome_occupancy %>%
  rename(Genomic_proportion = "Total_Genomic_Proportion")

# Subset the genome size information
GS_df <- subset(N_reads_analyzed_GS[,c(2,4)])
colnames(GS_df)[2] <- "Genome_Size_gbp"
# Convert pg to Gbp
GS_df$Genome_Size_gbp <- GS_df$Genome_Size_gbp*0.978

# Add the total genomic proportion and genome size columns
sum_All_Eggplant_Repeat_ag <- All_Eggplant_Repeat_ag %>%
  left_join(Genome_occupancy, by = "Species_Name") %>%
  left_join(GS_df, by = "Species_Name") %>%
  mutate(Species_Name = reorder(Species_Name, Total_Genomic_Proportion))

##### Heatmap of genomic proportion for different repeat type
All_Eggplant_Repeat_ag_pct %>%
  as.data.frame() %>%
  rownames_to_column("Species_Name") %>%
  pivot_longer(-c(Species_Name), names_to = "Repeat_type", values_to = "Genomic_proportion") %>%
  ggplot(aes(x = Species_Name, y = Repeat_type, fill = Genomic_proportion)) + 
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  xlab("Species Name") +
  ylab("Repeat Type") +
  scale_fill_viridis_c("Genomic Proportion (%)")

##### A stacked bar chart for each species and the abundances of different repeats
### Repeat % of all repeat types
ggplot(All_Eggplant_Repeat_ag, aes(fill=Final_annotation, y=Species_Name, x=Size_adjusted)) + 
  theme(axis.text.x = element_text(size= 10, hjust = 1)) + 
  theme(legend.text = element_text(size=6)) + 
  xlab("Percentage of all repeats") + 
  ylab("Species Name") +
  geom_bar(position="fill", stat="identity")


### Absolute repeat number
ggplot(All_Eggplant_Repeat_ag, aes(fill=Final_annotation, y=Species_Name, x=Size_adjusted)) + 
  theme(axis.text.x = element_text(size= 10)) + 
  theme(legend.text = element_text(size=6)) + 
  xlab("Number of repeat reads") + 
  ylab("Species Name") +
  geom_bar(stat="identity")


### Repeat % of whole genome with Genome Size data as a line (proportion goes from low to high)
## Plot without custom colors
RptProp_GS_bySpecies <- ggplot(sum_All_Eggplant_Repeat_ag, aes(x = Species_Name)) +
  geom_bar(aes(y = Genomic_proportion, fill = Final_annotation), stat = "identity", position="stack") +
  geom_line(aes(y = Genome_Size_gbp * 10), group=1, color="blue4") +  # Adjust the scaling factor (* 10) as needed
  scale_y_continuous(sec.axis = sec_axis(~ (. / 10), name = "Average Genome Size (Gbp/2C)")) +
  labs(x = "Species Name", y = "Percentage of each repeat occupied within the genome") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), 
        legend.text = element_text(size = 7),
        legend.location = "plot") +
  guides(fill = guide_legend(title = "Repeat Annotation", ncol = 1))

print(RptProp_GS_bySpecies)

## Plot with custom colors
# Custom the colors of the repeat types
custom_palette <- c(
  "#A5A79EFF", "#A16928FF", "#BD925AFF", "#D6BD8DFF", "#D0D3A2FF", 
  "#B5C8B8FF", "#79A7ACFF", "#2887A1FF", "#A7DBD8FF",
  "#8FA3ABFF","#76A1CDFF", "#4E79A5FF", "#1B3A6BFF",
  paletteer::paletteer_d("khroma::sunset"),
  paletteer::paletteer_d("PrettyCols::RedBlues"),
  "#486078FF"
)

RptProp_GS_bySpecies_custm <- ggplot(sum_All_Eggplant_Repeat_ag, aes(x = Species_Name)) +
  geom_bar(aes(y = Genomic_proportion, fill = Final_annotation), stat = "identity", position="stack") +
  geom_line(aes(y = Genome_Size_gbp * 10), group=1, color="blue4") +  # Adjust the scaling factor (* 10) as needed
  scale_y_continuous(sec.axis = sec_axis(~ (. / 10), name = "Average Genome Size (Gbp/2C)")) +
  labs(title = "Genomic Repeat Type Proportions and Genome Size Across Different Species", 
       x = "Species Name", y = "Percentage of each repeat occupied within the genome") +
  scale_fill_manual(values = custom_palette) +
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20, 
                                    margin = margin(r = 15, unit = "pt")), 
        axis.text.x = element_text(angle = 30, hjust = 1, size = 15), 
        axis.text.y = element_text(hjust = 1, size = 15), 
        legend.text = element_text(size = 13),
        legend.location = "plot",
        plot.margin = margin(t = 15, r = 15, b = 10, l = 60, unit = "pt"),
        plot.title = element_text(size = 25, face = "bold")) +
  guides(fill = guide_legend(title = "Repeat Annotation", ncol = 1,
                             theme(legend.title = element_text(size = 18))))

print(RptProp_GS_bySpecies_custm)

# Save the figure as a pdf file
ggsave("../R_Plots/RptProp_GS_bySpecies.pdf", RptProp_GS_bySpecies_custm, width = 20, height = 10)

##### Total Genome Proportion and Diversity #####
#################################################

### Calculate the Shannon's diversity index
Shannon_Eggplant_Repeat <- vegan::diversity(All_Eggplant_Repeat_ag_reads, index = "shannon")
# Compare Shannon indices
Shannon_df <- Shannon_Eggplant_Repeat %>%
  as.data.frame() %>%
  rownames_to_column("Species_Name")

Shannon_df <- Shannon_df %>%
  rename(`.`="Shannon_Index") # Change the 2nd column name to what we want

### Richness (Mehinick's Index) - Closer to 1 = higher richness
n <- apply(All_Eggplant_Repeat_ag_reads > 0, 1, sum) # number of repeat types for each species
N <- apply(All_Eggplant_Repeat_ag_reads, 1, sum) # total number of reads for each species
# Calculate the Mehinick's Index
MI_df <- n/sqrt(N)
MI_df <- MI_df %>%
  as.data.frame() %>%
  rownames_to_column("Species_Name")

MI_df <- MI_df %>%
  rename(`.`= "Mehinicks_index") # Change the 2nd column name to what we want

### Combine the data frames of GS, richness, diversity together with the aggregated repeat profile percentage
Cmbd_All_Eggplant_Repeat_ag_pct <- All_Eggplant_Repeat_ag_pct %>%
  rownames_to_column("Species_Name")%>%
  left_join(Genome_occupancy, by = "Species_Name") %>%
  left_join(Shannon_df, by = "Species_Name") %>%
  left_join(MI_df, by = "Species_Name") %>%
  left_join(GS_df, by = "Species_Name")
row.names(Cmbd_All_Eggplant_Repeat_ag_pct) <- Cmbd_All_Eggplant_Repeat_ag_pct$Species_Name 
Cmbd_All_Eggplant_Repeat_ag_pct <- subset(Cmbd_All_Eggplant_Repeat_ag_pct, select = -Species_Name)

# Remove the polyploid Solanum_campylacanthum
Diploid_Eggplant_Repeat_ag_pct <- Cmbd_All_Eggplant_Repeat_ag_pct[-3,]

##~~~~~~~~~~Exploratory plotting~~~~~~~~~~## 

# Richness and genome size
plot(Cmbd_All_Eggplant_Repeat_ag_pct$Mehinicks_index ~ Cmbd_All_Eggplant_Repeat_ag_pct$Genome_Size_gbp, ylab = "Repeat type richness (Mehinicks Index)", xlab = "Genome size (Gbp/2C)")

# All samples
All_Richness_GS <- ggplot(Cmbd_All_Eggplant_Repeat_ag_pct, aes(x = as.numeric(Genome_Size_gbp), 
                                            y = Mehinicks_index)) + 
  geom_point() +
  geom_smooth(method="loess", span = 1, formula = y ~ x) +
  labs(title = "Genome Size vs. Repeat Type Richness", 
       y = "Repeat type richness (Mehinicks Index)", x = "Genome size (Gbp/2C)") + 
  theme(text = element_text(size = 14), axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 25, face = "bold"))

print(All_Richness_GS)

ggsave("../R_Plots/All_Richness_GS.pdf", All_Richness_GS, width = 10, height = 8)        

# Only diploid samples
Diploid_Richness_GS <- ggplot(Diploid_Eggplant_Repeat_ag_pct, 
                              aes(x = as.numeric(Genome_Size_gbp), 
                                  y = Mehinicks_index)) + 
  geom_point() +
  geom_smooth(method = "loess", span = 1, formula = y ~ x) +
  labs(title = "Genome Size vs. Repeat Type Richness", 
       y = "Repeat type richness (Mehinicks Index)", x = "Genome size (Gbp/2C)") + 
  theme(text = element_text(size = 14), axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 25, face = "bold"))

print(Diploid_Richness_GS)

ggsave("../R_Plots/Diploid_Richness_GS.pdf", Diploid_Richness_GS, width = 10, height = 8)

### Genomic proportion of repeats and genome size
plot(Cmbd_All_Eggplant_Repeat_ag_pct$Total_Genomic_Proportion ~ Cmbd_All_Eggplant_Repeat_ag_pct$Genome_Size_gbp, ylab = "Genomic proportion of repeats (%)", xlab="Genome size (Gbp/2C)")

# All samples
All_GP_GS <- ggplot(Cmbd_All_Eggplant_Repeat_ag_pct, aes(x = as.numeric(Genome_Size_gbp), 
                                                         y = Total_Genomic_Proportion)) + 
  geom_point() +
  geom_smooth(method = "loess", span = 1, formula = y ~ x) +
  labs(title = "Genome Size vs. Genomic Proportion of Repeats", 
       y = "Genomic proportion of repeats (%)", x = "Genome size (Gbp/2C)") + 
  theme(text = element_text(size = 14), axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 25, face = "bold"))

print(All_GP_GS)

ggsave("../R_Plots/All_GP_GS.pdf", All_GP_GS, width = 10, height = 8)        

# Only diploid samples
Diploid_GP_GS <- ggplot(Diploid_Eggplant_Repeat_ag_pct, 
                        aes(x = as.numeric(Genome_Size_gbp), 
                            y = Total_Genomic_Proportion)) + 
  geom_point() +
  geom_smooth(method = "loess", span = 1, formula = y ~ x) +
  labs(title = "Genome Size vs. Genomic Proportion of Repeats", 
       y = "Genomic proportion of repeats (%)", x = "Genome size (Gbp/2C)") + 
  theme(text = element_text(size = 14), axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 25, face = "bold"))

print(Diploid_GP_GS)

ggsave("../R_Plots/Diploid_GP_GS.pdf", Diploid_GP_GS, width = 10, height = 8)

### Diversity and genome size
plot(Cmbd_All_Eggplant_Repeat_ag_pct$Shannon_Index ~ Cmbd_All_Eggplant_Repeat_ag_pct$Genome_Size_gbp, ylab = "Repeat diversity (Shannon Index)", xlab = "Genome size (Gbp/2C)")

# All samples
All_Diversity_GS <- ggplot(Cmbd_All_Eggplant_Repeat_ag_pct, 
                           aes(x = as.numeric(Genome_Size_gbp),
                               y = Shannon_Index)) + 
  geom_point() +
  geom_smooth(method = "loess", span = 1, formula = y ~ x) +
  labs(title = "Genome Size vs. Repeat Type Diversity", 
       y = "Repeat type diversity (Shannon Index)", x = "Genome size (Gbp/2C)") + 
  theme(text = element_text(size = 14), axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 25, face = "bold"))

print(All_Diversity_GS)

ggsave("../R_Plots/All_Diversity_GS.pdf", All_Diversity_GS, width = 10, height = 8)

# Only diploid samples
Diploid_Diversity_GS <- ggplot(Diploid_Eggplant_Repeat_ag_pct, 
                               aes(x = as.numeric(Genome_Size_gbp), 
                                   y = Shannon_Index)) + 
  geom_point() +
  geom_smooth(method = "loess", span = 1, formula = y ~ x) +
  labs(title = "Genome Size vs. Repeat Type Diversity", 
       y = "Repeat type diversity (Shannon Index)", x = "Genome size (Gbp/2C)") + 
  theme(text = element_text(size = 14), axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 25, face = "bold"))

print(Diploid_Diversity_GS)

ggsave("../R_Plots/Diploid_Diversity_GS.pdf", Diploid_Diversity_GS, width = 10, height = 8)




############### Phylogenetic Tree ###############
#################################################


######### Repeat Types and Genome Size ##########
#################################################
