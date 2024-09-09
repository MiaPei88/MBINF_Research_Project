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
library(readxl)
library(rgbif)
library(countrycode)
library(CoordinateCleaner)
library(rnaturalearth)
library(maps)
library(ggmap)
library(raster)
library(sf)
library(sp)

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

## Calculate the genome coverage of the sequencing data (Change the unit of Number of Reads Analyzed and Genome Size to Mbp)
N_reads_analyzed_GS$Genome_Coverage <- (N_reads_analyzed_GS[,3] * 235 / 1000000) / (N_reads_analyzed_GS[,4] * 97.8)
# When the genome coverage is 0.1-0.5, the repeat analysis is more reliable

# Delete rows that are annotated as organelle
All_Eggplant_Repeat <- All_Eggplant_Repeat[-which(str_detect(All_Eggplant_Repeat$Final_annotation,"organelle")),]

# Check the repeat classification
unique(All_Eggplant_Repeat$Final_annotation)

# Aggregate counts of each repeat classification in the annotation for each species 
All_Eggplant_Repeat_ag <- aggregate(Size_adjusted~Species_Name+Final_annotation, data=All_Eggplant_Repeat, FUN=sum) 

sub_N_reads_analyzed_GS <- N_reads_analyzed_GS %>%
  dplyr::select(Species_Name, `Genome_Size.mean2C.`,
                Number_of_reads_nuclear)

# Calculate the genomic proportion of individual clusters by dividing the size 
# (number of reads, column called Size_adjusted) of the clusters by Number_of_reads_nuclear
All_Eggplant_Repeat_ag <- left_join(All_Eggplant_Repeat_ag, N_reads_analyzed_GS %>%
                                      dplyr::select(Species_Name, Genome_Size.mean2C.,
                                                    Number_of_reads_nuclear), 
                                    by = "Species_Name") %>%
  mutate(Genomic_proportion = (Size_adjusted/Number_of_reads_nuclear)) %>%
  mutate(Final_annotation = case_when(
    Final_annotation == "All"  ~ "Other",
    TRUE ~ Final_annotation))


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

### Repeat % of all repeat types
ggplot(All_Eggplant_Repeat_ag, aes(fill=Final_annotation, 
                                   y=Species_Name, 
                                   x=Size_adjusted)) + 
  scale_fill_manual(values = custom_palette) +
  theme(axis.text.x = element_text(size= 10, hjust = 1)) + 
  theme(legend.text = element_text(size=6)) + 
  xlab("Percentage of all repeats") + 
  ylab("Species Name") +
  geom_bar(position="fill", stat="identity")


### Absolute repeat number
ggplot(All_Eggplant_Repeat_ag, aes(fill = Final_annotation, 
                                   y = Species_Name, 
                                   x = Size_adjusted)) + 
  scale_fill_manual(values = custom_palette) +
  theme(axis.text.x = element_text(size= 10)) + 
  theme(legend.text = element_text(size=6)) + 
  xlab("Number of repeat reads") + 
  ylab("Species Name") +
  geom_bar(stat="identity")


### Repeat % of whole genome with Genome Size data as a line (proportion goes from low to high)

RptProp_GS_bySpecies_custm <- ggplot(sum_All_Eggplant_Repeat_ag, aes(x = reorder(Species_Name, Genome_Size_gbp))) +
  geom_bar(aes(y = Genomic_proportion, fill = Final_annotation), stat = "identity", 
           position = "stack") +
    labs(x = "Species name", y = "Total genomic proportion of repeats") +
  scale_fill_manual(values = custom_palette) +
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20, 
                                    margin = margin(r = 15, unit = "pt")), 
        axis.text.x = element_text(angle = 30, hjust = 1, size = 15), 
        axis.text.y = element_text(hjust = 1, size = 15), 
        legend.text = element_text(size = 13),
        legend.location = "plot",
        legend.margin = margin(6, 6, 6, 10),
        plot.margin = margin(t = 15, r = 15, b = 10, l = 60, unit = "pt")) +
  guides(fill = guide_legend(title = "Repeat Annotation", ncol = 1,
                             theme(legend.title = element_text(size = 18))))

print(RptProp_GS_bySpecies_custm)

# Save the figure as a pdf file
ggsave("../R_Plots/RptProp_GS_bySpecies.pdf", RptProp_GS_bySpecies_custm, width = 20, height = 10)

RptProp_GS_bySpecies_horizontal <- ggplot(sum_All_Eggplant_Repeat_ag, aes(y = Species_Name)) +
  geom_bar(aes(x = Genomic_proportion, fill = Final_annotation), stat = "identity", 
           position = "stack") +
  labs(title = "Figure. 1", 
       x = "Genomic proportion of total repeat", y = "Species Name") +
  scale_fill_manual(values = custom_palette) +
  theme(axis.title.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20, 
                                    margin = margin(r = 15, unit = "pt")), 
        axis.text.y = element_text(hjust = 1, size = 15), 
        axis.text.x = element_text(hjust = 1, size = 15), 
        legend.text = element_text(size = 13),
        legend.location = "plot",
        plot.title = element_text(size = 25, face = "bold")) +
  guides(fill = guide_legend(title = "Repeat Annotation", ncol = 1,
                             theme(legend.title = element_text(size = 18))))

print(RptProp_GS_bySpecies_horizontal)
ggsave("../R_Plots/RptProp_GS_bySpecies_horizontal.pdf", RptProp_GS_bySpecies_horizontal, width = 20, height = 10)

##### Total Genome Proportion and Diversity #####
#################################################
colnames(All_Eggplant_Repeat_ag_reads)

# Exclude repeat types that are not specific enough 
Specific_Eggplant_Repeat_ag_reads <- All_Eggplant_Repeat_ag_reads[, -c(1, 2, 3, 6, 15, 16, 29)]
colnames(Specific_Eggplant_Repeat_ag_reads)

### Calculate the Shannon's diversity index
Shannon_Eggplant_Repeat <- vegan::diversity(Specific_Eggplant_Repeat_ag_reads, index = "shannon")
# Compare Shannon indices
Shannon_df <- Shannon_Eggplant_Repeat %>%
  as.data.frame() %>%
  rownames_to_column("Species_Name")

Shannon_df <- Shannon_df %>%
  rename(`.`="Shannon_Index") # Change the 2nd column name to what we want

### Richness (Mehinick's Index) - Closer to 1 = higher richness
n <- apply(Specific_Eggplant_Repeat_ag_reads > 0, 1, sum) # number of repeat types for each species
N <- apply(Specific_Eggplant_Repeat_ag_reads, 1, sum) # total number of reads for each species
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
All_Richness_GS <- ggplot(Cmbd_All_Eggplant_Repeat_ag_pct, 
                          aes(x = as.numeric(Genome_Size_gbp), 
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
       y = "Genomic proportion of repeats", x = "Genome size (Gbp/2C)") + 
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
       y = "Genomic proportion of repeats", x = "Genome size (Gbp/2C)") + 
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

# Linear Regression Modeling between Genome size and Total Repeat Proportion
lm_GS_TGP <- lm(Genome_Size_gbp ~ Total_Genomic_Proportion, 
               data = Diploid_Eggplant_Repeat_ag_pct)
summary(lm_GS_TGP)

# Linear Regression Modeling between Genome size and Repeat Richness
lm_GS_Richness <- lm(Genome_Size_gbp ~ Mehinicks_index, 
                     data = Diploid_Eggplant_Repeat_ag_pct)
summary(lm_GS_Richness)

lm_GS_Richness_log <- lm(log(Genome_Size_gbp) ~ Mehinicks_index, 
                data = Diploid_Eggplant_Repeat_ag_pct)
summary(lm_GS_Richness_log)

# Linear Regression Modeling between Genome size and Repeat Diversity
lm_GS_Diversity <- lm(Genome_Size_gbp ~ Shannon_Index, 
                     data = Diploid_Eggplant_Repeat_ag_pct)
summary(lm_GS_Diversity)


## Ploting the linear correlation between Genome size and Total Repeat Proportion
# Extract coefficients
intercept <- coef(lm_GS_TGP)[1]
slope <- coef(lm_GS_TGP)[2]

# Create equation text
equation <- paste("y = ", round(slope, 4), "x", " + ", round(intercept, 4), sep = "")

# Create the plot with the linear model line
lm_GS_TGP_Plot <- ggplot(Diploid_Eggplant_Repeat_ag_pct, 
                         aes(x = Total_Genomic_Proportion, 
                             y = Genome_Size_gbp)) +
  labs(y = "Genome size (Gbp/2C)", x = "Total genomic proportion of repeats") + 
  geom_point(cex = 2) +
  geom_smooth(method = lm, formula = y ~ x) +
  theme(text = element_text(size = 14), axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20, margin = margin(r = 10, unit = "pt"))) +
  # Add the equation to the plot
  annotate("text", x = 0.44, y = 2.71, label = equation, hjust = 0, size = 7, color = "blue")

print(lm_GS_TGP_Plot)

ggsave("../R_Plots/lm_GS_TGP_Plot.pdf", lm_GS_TGP_Plot, width = 12, height = 8)


############### Phylogenetic Tree ###############
#################################################
# Read in phylo tree from Xavier
Xavier_Tree <- read.nexus("/Users/miapei/Desktop/MBinf/Research_Project/Data/Xavier_phylo/FINALConsensus_MAFFT.tre")
plot(Xavier_Tree)

# Names of tips to be dropped
tips_to_drop <- setdiff(Xavier_Tree$tip.label, N_reads_analyzed_GS$Sample_ID) 
# Drop the tips
Mia_tree <- drop.tip(Xavier_Tree, tips_to_drop)
plot(Mia_tree)

# Make a name map to correspond the species name to each sample id
name_map <- setNames(N_reads_analyzed_GS$Species_Name, N_reads_analyzed_GS$Sample_ID)

# Change the tip labels using the name map
Mia_tree$tip.label <- name_map[Mia_tree$tip.label]
plot(Mia_tree)

################# Occurrence Data ###############
#################################################
### Read in the Occurrence data from Edeline's previous study
Edeline_Occ <- read.csv("/Users/miapei/Desktop/MBinf/Research_Project/Data/Occurence.csv")

### Add the original country of each species we have for this project
NextGen <- read_excel("/Users/miapei/Desktop/MBinf/Research_Project/Data/Melongena_NextGen.xlsx")
names(NextGen)[1] <- "Sample_ID" # Change column name
NextGen <- NextGen %>%
  dplyr::select(Sample_ID, country) # Select columns we need

# Join the two dataframe to add the country column
Species_Country <- N_reads_analyzed_GS %>%
  inner_join(NextGen, by = "Sample_ID") %>%
  dplyr::select(Species_Name, country) %>%
  dplyr::rename(genus.sp = Species_Name) %>%
  dplyr::rename(COUNTRY = country)

Species_Country[11,1] <- "Solanum_incanum"
Species_Country[13,1] <- "Solanum_incanum"

####### Species from one country in this project
# Filter only the species and its original country we have for this project
Mia_Occ_onectry <- inner_join(Edeline_Occ, Species_Country, by = c("genus.sp", "COUNTRY"))

# Check the number of occurrence data points of each species
table(Mia_Occ_onectry$genus.sp)
unique(Mia_Occ_onectry$genus.sp)
unique(Species_Country$genus.sp)

# Check which species do not have occurence data
setdiff(Species_Country$genus.sp, Mia_Occ_onectry$genus.sp) 
## The three species without occurence data are cultivated crops and hard to observe their natural distribution

# Check the number of occurrence data points of Solanum incanum from the two different countries
nrow(Mia_Occ[Mia_Occ_onectry$genus.sp == "Solanum_incanum" & Mia_Occ_onectry$COUNTRY == "Burkina Faso",])
nrow(Mia_Occ[Mia_Occ_onectry$genus.sp == "Solanum_incanum" & Mia_Occ_onectry$COUNTRY == "Kenya",])

#### Add occurrence data points of Solanum incanum in Burkina Faso from GBIF
S_incanum_BFaso <- read.csv("/Users/miapei/Desktop/MBinf/Research_Project/Data/Solanum_incanum_BFaso.csv", sep = "\t")

# Remove data points that are observed by human or the uncertainty in meters is too high
S_incanum_BFaso_cl <- S_incanum_BFaso %>%
  filter(basisOfRecord != "HUMAN_OBSERVATION") %>%
  filter(coordinateUncertaintyInMeters / 1000 <= 100 | 
           is.na (coordinateUncertaintyInMeters)) 

### Wrangle the columns we want and combine them with the data points from Edeline's occurrence data
S_incanum_BFaso_sub <- S_incanum_BFaso_cl %>%
  dplyr::rename(LATDEC = decimalLatitude, 
                LONGDEC = decimalLongitude, 
                COUNTRY = countryCode) %>%
  mutate(COUNTRY = "Burkina Faso") %>%
  mutate(SP1 = word(species, 2)) %>%
  mutate(genus.sp = str_c(word(species, 1), word(species, 2), sep = "_")) %>%
  dplyr::select(LATDEC, LONGDEC, COUNTRY, SP1, genus.sp)

Mia_BFaso <- Mia_Occ_onectry %>% 
  dplyr::filter(genus.sp == "Solanum_incanum" & COUNTRY == "Burkina Faso")

S_incanum_BFaso_cmbd <- bind_rows(Mia_BFaso, S_incanum_BFaso_sub)


# Change the countryCode into the format CoordinateCleaner accepts
S_incanum_BFaso_cmbd$COUNTRY <-  countrycode(S_incanum_BFaso_cmbd$COUNTRY, 
                                            origin =  'country.name',
                                            destination = 'iso3c')

# Flag the problematic coordinates using CoordinateCleaner
S_incanum_BFaso_flags <- clean_coordinates(x = S_incanum_BFaso_cmbd, 
                                           lon = "LONGDEC", 
                                           lat = "LATDEC",
                                           countries = "COUNTRY",
                                           species = "SP1",
                                           tests = c("capitals", "centroids",
                                                     "equal", "zeros", "countries",
                                                     "gbif", "institutions", 
                                                     "outliers", "seas", "duplicates"),
                                           inst_rad = 1000)

# Check and plot the flags
summary(S_incanum_BFaso_flags)
plot(S_incanum_BFaso_flags, lon = "LONGDEC", lat = "LATDEC")

# Clean the occurrence data to exclude the coordinates with flags
S_incanum_BFaso_cl <- S_incanum_BFaso_cmbd[S_incanum_BFaso_flags$.summary,]

### Remove coordinates that are spatially autocorrelated
## A function to find coordinates that are spatially autocorrelated - based on script from Edeline Gagnon: https://github.com/edgagnon/Geophytes_Solanum-/blob/main/00_Data/02_Occurrence_Data/003_Spatial_filtering_scripts_SOLANUM.R
filterByProximity <- function(xy, dist, mapUnits = FALSE) {
  require(sp)
  if (!mapUnits) {
    d <- spDists(xy, longlat = TRUE)
  } else {
    d <- spDists(xy, longlat = FALSE)
  }
  diag(d) <- NA
  close <- (d <= dist)
  diag(close) <- NA
  close_indices <- which(close, arr.ind = TRUE)
  
  if (length(close_indices) > 0) {
    # Get unique indices to remove (any row participating in a close pair)
    to_remove <- unique(as.vector(close_indices))
    to_keep <- setdiff(1:nrow(xy), to_remove)
  } else {
    to_keep <- 1:nrow(xy)
  }
  
  return(to_keep)  # Return indices of the rows to keep
}

# Set the distance for spatial filtering (e.g., 10 km)
distance_threshold <- 10  # distance in kilometers if mapUnits = FALSE

# Convert DataFrame coordinates to a matrix
S_incanum_BFaso_coordmatrix <- cbind(S_incanum_BFaso_cl$LONGDEC,
                                     S_incanum_BFaso_cl$LATDEC)

# Get indices of the data points to keep
S_incanum_BFaso_indices_keep <- filterByProximity(S_incanum_BFaso_coordmatrix, distance_threshold, mapUnits = FALSE)
S_incanum_BFaso_indices_keep

# Now extract those rows from the original DataFrame
S_incanum_BFaso_filtered <- S_incanum_BFaso_cl[S_incanum_BFaso_indices_keep, ]


### Mapping the coordinates of S.incanum in Burkina Faso from GBIF
# Get map data of Burkina Faso
BFaso_map <- map_data("world", region = "Burkina Faso")

# Plot the distribution of S.incanum in the map of Burkina Faso
ggplot() +
  geom_polygon(data = BFaso_map, aes(x = long, y = lat, group = group), fill = "lightblue", color = "black") +
  geom_point(data = S_incanum_BFaso_filtered, aes(x = LONGDEC, y = LATDEC), color = "red", size = 2) +
  labs(title = "Distribution of Solanum incanum in Burkina Faso (GBIF)") +
  theme_minimal()


#### Add occurrence data points of Solanum incanum in Kenya from GBIF(already exclude data points without coordinates)
S_incanum_Kenya <- read.csv("/Users/miapei/Desktop/MBinf/Research_Project/Data/Solanum_incanum_Kenya.csv", sep = "\t")

# Remove data points that are observed by human or the uncertainty in meters is too high
S_incanum_Kenya <- S_incanum_Kenya %>%
  filter(basisOfRecord != "HUMAN_OBSERVATION") %>%
  filter(coordinateUncertaintyInMeters / 1000 <= 100 | 
           is.na (coordinateUncertaintyInMeters))

### Wrangle the columns we want to combine with the major occurrence df
S_incanum_Kenya_sub <- S_incanum_Kenya %>%
  dplyr::rename(LATDEC = decimalLatitude, 
                LONGDEC = decimalLongitude, 
                COUNTRY = countryCode) %>%
  mutate(COUNTRY = "Kenya") %>%
  mutate(SP1 = word(species, 2)) %>%
  mutate(genus.sp = str_c(word(species, 1), word(species, 2), sep = "_")) %>%
  dplyr::select(LATDEC, LONGDEC, COUNTRY, SP1, genus.sp)

Mia_Kenya <- Mia_Occ_onectry %>% 
  dplyr::filter(genus.sp == "Solanum_incanum" & COUNTRY == "Kenya")

S_incanum_Kenya_cmbd <- bind_rows(Mia_Kenya, S_incanum_Kenya_sub)

# Change the countryCode into the format CoordinateCleaner accepts
S_incanum_Kenya_cmbd$COUNTRY <-  countrycode(S_incanum_Kenya_cmbd$COUNTRY, 
                                            origin =  'country.name',
                                            destination = 'iso3c')

# Flag the problematic coordinates using CoordinateCleaner
S_incanum_Kenya_flags <- clean_coordinates(x = S_incanum_Kenya_cmbd, 
                                           lon = "LONGDEC", 
                                           lat = "LATDEC",
                                           countries = "COUNTRY",
                                           species = "SP1",
                                           tests = c("capitals", "centroids",
                                                     "equal", "zeros", "countries",
                                                     "gbif", "institutions", 
                                                     "outliers", "seas", "duplicates"),
                                           inst_rad = 1000)

# Check and plot the flags
summary(S_incanum_Kenya_flags)
plot(S_incanum_Kenya_flags, lon = "LONGDEC", lat = "LATDEC")

# Clean the occurrence data to exclude the coordinates with flags
S_incanum_Kenya_cl <- S_incanum_Kenya_cmbd[S_incanum_Kenya_flags$.summary,]

### Remove coordinates that are spatially autocorrelated
# Convert DataFrame coordinates to a matrix
S_incanum_Kenya_coordmatrix <- cbind(S_incanum_Kenya_cl$LONGDEC,
                                     S_incanum_Kenya_cl$LATDEC)

# Get indices of the data points to keep
S_incanum_Kenya_indices_keep <- filterByProximity(S_incanum_Kenya_coordmatrix, distance_threshold, mapUnits = FALSE)

# Now extract those rows from the original DataFrame
S_incanum_Kenya_filtered <- S_incanum_Kenya_cl[S_incanum_Kenya_indices_keep, ]
nrow(S_incanum_Kenya_filtered)


### Mapping the coordinates of S.incanum in Kenya from GBIF
# Get map data of Kenya
Kenya_map <- map_data("world", region = "Kenya")

# Plot the distribution of S.incanum in the map of Kenya
ggplot() +
  geom_polygon(data = Kenya_map, aes(x = long, y = lat, group = group), fill = "lightblue", color = "black") +
  geom_point(data = S_incanum_Kenya_filtered, aes(x = LONGDEC, y = LATDEC), color = "red", size = 2) +
  labs(title = "Distribution of Solanum incanum in Kenya (GBIF)") +
  theme_minimal()


#### Check the distribution of the 3 cultivated crops
Origin_Occ <- read_excel("/Users/miapei/Desktop/MBinf/Research_Project/Data/Solanum_all_COLLEXTRACT_21-04-2021_at_08-39-16_v2.xlsx")

## Filtering
# Filter only the 3 cultivated species we need
Cultiv_Species_Occ <- Origin_Occ %>%
  filter(SP1 %in% c("melongena", "macrocarpon", "aethiopicum"))

rm(Origin_Occ) #Clear out some space

Cultiv_Species_Occ_cl_onectry <- Cultiv_Species_Occ %>%
  filter(SP1 == "melongena" & COUNTRY == "Thailand"|
           SP1 == "macrocarpon" & COUNTRY == "Tanzania"|
           SP1 == "aethiopicum" & COUNTRY == "Tanzania") %>%
  filter(CULTIVATED == "No") %>% # keep only the non-cultivated specimens
  filter(!grepl("Cultivated|cultivated|cultivation", HABITATTXT)) %>%
  filter(!grepl("cultivated|Planted", HABITATTXT))# keep only the non-cultivated specimens

Cultiv_Species_Occ_cl_allctry <- Cultiv_Species_Occ %>%
  filter(SP1 == "melongena"|
           SP1 == "macrocarpon"|
           SP1 == "aethiopicum" ) %>%
  filter(CULTIVATED == "No") %>% # keep only the non-cultivated specimens
  filter(!grepl("Cultivated|cultivated|cultivation", HABITATTXT)) %>%
  filter(!grepl("cultivated|Planted", HABITATTXT))# keep only the non-cultivated specimens

# Change the countryCode into the format CoordinateCleaner accepts
Cultiv_Species_Occ_cl_onectry$COUNTRY <-  countrycode(Cultiv_Species_Occ_cl_onectry$COUNTRY,
                                                    origin =  'country.name',
                                                    destination = 'iso3c')

Cultiv_Species_Occ_cl_allctry$COUNTRY <-  countrycode(Cultiv_Species_Occ_cl_allctry$COUNTRY,
                                                      origin =  'country.name',
                                                      destination = 'iso3c')

# Flag the problematic coordinates using CoordinateCleaner
Cultiv_Species_flags_onectry <- clean_coordinates(x = Cultiv_Species_Occ_cl_onectry, 
                                           lon = "LONGDEC", 
                                           lat = "LATDEC",
                                           countries = "COUNTRY",
                                           species = "SPECIES",
                                           tests = c("capitals", "centroids",
                                                     "equal", "zeros", "countries",
                                                     "gbif", "institutions", 
                                                     "outliers", "seas", "duplicates"),
                                           inst_rad = 1000)

Cultiv_Species_flags_allctry <- clean_coordinates(x = Cultiv_Species_Occ_cl_allctry, 
                                                  lon = "LONGDEC", 
                                                  lat = "LATDEC",
                                                  countries = "COUNTRY",
                                                  species = "SPECIES",
                                                  tests = c("capitals", "centroids",
                                                            "equal", "zeros", "countries",
                                                            "gbif", "institutions", 
                                                            "outliers", "seas", "duplicates"),
                                                  inst_rad = 1000)

# Check and plot the flags
summary(Cultiv_Species_flags_onectry)
plot(Cultiv_Species_flags_onectry, lon = "LONGDEC", lat = "LATDEC")

summary(Cultiv_Species_flags_allctry)
plot(Cultiv_Species_flags_allctry, lon = "LONGDEC", lat = "LATDEC")


# Clean the occurrence data to exclude the coordinates with flags and select columns we need
Cultiv_Species_cl_onectry <- Cultiv_Species_Occ_cl_onectry[Cultiv_Species_flags_onectry$.summary,] %>%
  mutate(genus.sp = str_c(word(SPECIES, 1), word(SPECIES, 2), sep = "_")) %>%
  dplyr::select(intersect(colnames(Cultiv_Species_cl_onectry),colnames(Mia_Occ_onectry)))

Cultiv_Species_cl_allctry <- Cultiv_Species_Occ_cl_allctry[Cultiv_Species_flags_allctry$.summary,] %>%
  mutate(genus.sp = str_c(word(SPECIES, 1), word(SPECIES, 2), sep = "_")) %>%
  dplyr::select(intersect(colnames(Cultiv_Species_cl_allctry),colnames(Mia_Occ_allctry)))


# Check the number of each species left after filtering
table(Cultiv_Species_cl_onectry$SP1)
table(Cultiv_Species_cl_allctry$SP1)

### Remove coordinates that are spatially autocorrelated
# Convert DataFrame coordinates to a matrix
Cultiv_Species_onectry_coordmatrix <- cbind(Cultiv_Species_cl_onectry$LONGDEC,
                                     Cultiv_Species_cl_onectry$LATDEC)
Cultiv_Species_allctry_coordmatrix <- cbind(Cultiv_Species_cl_allctry$LONGDEC,
                                            Cultiv_Species_cl_allctry$LATDEC)

# Get indices of the data points to keep
Cultiv_Species_allctry_indices_keep <- filterByProximity(Cultiv_Species_allctry_coordmatrix, distance_threshold, mapUnits = FALSE)
Cultiv_Species_allctry_indices_keep <- filterByProximity(Cultiv_Species_allctry_coordmatrix, distance_threshold, mapUnits = FALSE)

# Now extract those rows from the original DataFrame
Cultiv_Species_onectry_filtered <- Cultiv_Species_cl_onectry[Cultiv_Species_onectry_indices_keep, ]
nrow(Cultiv_Species_onectry_filtered)
Cultiv_Species_allctry_filtered <- Cultiv_Species_cl_allctry[Cultiv_Species_allctry_indices_keep, ]
nrow(Cultiv_Species_allctry_filtered)

# Check if the columns are the same between Cultivated species data frame and the data frame for the rest of the species
ncol(Cultiv_Species_onectry_filtered)
ncol(Mia_Occ_onectry)
setdiff(colnames(Mia_Occ_onectry), colnames(Cultiv_Species_onectry_filtered))

ncol(Cultiv_Species_allctry_filtered)
ncol(Mia_Occ_allctry)
setdiff(colnames(Mia_Occ_allctry), colnames(Cultiv_Species_allctry_filtered))


#### Combine the occurrence data of S.incanum from GBIF to the major occurrence df (one country)
# Exclude S.cerasferum (only one occurrence) and S.incanum (already moved to the two data frames of this species from two different countries)
Mia_Occ_onectry <- Mia_Occ_onectry %>%
  filter(genus.sp %in% Species_Country$genus.sp) %>%
  filter(genus.sp != "Solanum_cerasiferum" & genus.sp != "Solanum_incanum")

Mia_Occ_onectry_cmbd <- rbind(Mia_Occ_onectry, S_incanum_BFaso_filtered, S_incanum_Kenya_filtered, Cultiv_Species_onectry_filtered)

# Check the number of occurrence data points of each species after combination
table(Mia_Occ_onectry_cmbd$SP1)

####### Species from all country in this project
# Filter the species we have for this project
Mia_Occ_allctry <- Edeline_Occ %>%
  filter(genus.sp %in% Species_Country$genus.sp) %>%
  filter(genus.sp != "Solanum_cerasiferum" & genus.sp != "Solanum_incanum")

# Check the number of occurrence data points of each species
table(Mia_Occ_allctry$genus.sp)
unique(Mia_Occ_allctry$genus.sp)

# Check which species do not have occurence data
setdiff(Species_Country$genus.sp, Mia_Occ_allctry$genus.sp) 

# Combine the data points of S.incanum to the occurrence data frame
Mia_Occ_allctry_cmbd <- rbind(Mia_Occ_allctry, S_incanum_BFaso_filtered, S_incanum_Kenya_filtered, Cultiv_Species_allctry_filtered)

# Check the number of occurrence data points of each species after combination
table(Mia_Occ_allctry_cmbd$genus.sp)


################## Climate Data #################
#################################################

###### List the climate data files 

data_path <- "/Users/miapei/Desktop/MBinf/Research_Project/Data/Env_Data/"
files <- list.files(path = data_path, pattern = '.tif', full.names=TRUE)
list(files)

###### Chelsa Bio 1-19

extract_env <- function(filenumber, occurrence_df) {
  predictors <- stack(files[[filenumber]])

  # Extract raw values which you'll need for running SDMs (12 equals long, 11 equal lat)
  bio_values <- extract(predictors, occurrence_df[,12:11])
  bio_values1 <- as.data.frame(bio_values)
  
  return(bio_values1)
}

# Find the column numbers of lat and long
dim(Mia_Occ_allctry_cmbd)
names(Mia_Occ_allctry_cmbd) #11-12 are Lat and Long

dim(Mia_Occ_onectry_cmbd)
names(Mia_Occ_onectry_cmbd) #11-12 are Lat and Long


Mia_Env_allctry <- Mia_Occ_allctry_cmbd %>%
  dplyr::select(COUNTRY, LONGDEC, LATDEC, genus.sp, SP1)

Mia_Env_onectry <- Mia_Occ_onectry_cmbd %>%
  dplyr::select(COUNTRY, LONGDEC, LATDEC, genus.sp, SP1)

# Iterate through all tif files exclude the annual MI one for occurrence data from all countries
print("Iterate through df of all country")
for (i in 2:20) {
  # Extract the variable name from the filename, e.g. bio1
  variable_name <- str_extract(files[[i]], "bio\\d+")
  
  # Print the variable name
  print(variable_name)
  
  # Use the function extract_env
  biox <- extract_env(i, Mia_Occ_allctry_cmbd)
  
  # Print the column name of the extracted environmental data
  print(colnames(biox))
  
  # Adding extracted values to the data frame
  Mia_Env_allctry[[variable_name]] <- biox[,1]
  
  # Check and print NA values count
  print(table(is.na(Mia_Env_allctry[[variable_name]]))) 
  
  print("------------------------")
}

# Iterate through all tif files exclude the annual MI one for occurrence data from one countries
print("Iterate through df of one country")
for (i in 2:20) {
  # Extract the variable name from the filename, e.g. bio1
  variable_name <- str_extract(files[[i]], "bio\\d+")
  
  # Print the variable name
  print(variable_name)
  
  # Use the function extract_env
  biox <- extract_env(i, Mia_Occ_onectry_cmbd)
  
  # Print the column name of the extracted environmental data
  print(colnames(biox))
  
  # Adding extracted values to the data frame
  Mia_Env_onectry[[variable_name]] <- biox[,1]
  
  # Check and print NA values count
  print(table(is.na(Mia_Env_onectry[[variable_name]]))) 
  
  print("------------------------")
}


###### Annual Moisture Index (MI)

predictors <- stack(files[[1]]) #MI
predictors

## All countries
# extract raw values which you'll need for running SDMs (12 equals long, 11 equal lat)
AMI_values <- raster::extract(predictors, Mia_Occ_allctry_cmbd[,12:11])
AMI_values1 <- as.data.frame(AMI_values)
colnames(AMI_values1)

Mia_Env_allctry$mi<-AMI_values1$X04_AnnualMI_world
table(is.na(Mia_Env_allctry$mi))#2 occurrence records don't have GDD values

## One country
# extract raw values which you'll need for running SDMs (12 equals long, 11 equal lat)
AMI_values_1 <- raster::extract(predictors, Mia_Occ_onectry_cmbd[,12:11])
AMI_values2 <- as.data.frame(AMI_values_1)
colnames(AMI_values2)

Mia_Env_onectry$mi<-AMI_values2$X04_AnnualMI_world
table(is.na(Mia_Env_onectry$mi))#2 occurrence records don't have GDD values


