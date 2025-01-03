library(dplyr)
library(phyloseq)
library(vegan)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(Matrix)

# Use corrected data after
bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/vaginal_bacteria_cleaned.rds")
# samples.data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_samples.csv")
##########################################################################################

#### Exploratory Data Analysis

# Relative abundances
vaginal_relative_abundances <- transform_sample_counts(bacterial.data, function(x) x/sum(x))
# vaginal_relative_abundances <- transform_sample_counts(bacteria_subset_vaginal, function(x) x/sum(x))
relative_abundance_otu <- as.data.frame(otu_table(vaginal_relative_abundances))

relative_abundance_otu_t <- t(relative_abundance_otu) %>% as.data.frame()
relative_abundance_otu_t$SampleID <- rownames(relative_abundance_otu_t)

otu_with_participant <- relative_abundance_otu_t %>%
  left_join(bacteria_metadata_df, by="SampleID")

# bray_curtis_dist <- vegdist(t(relative_abundance_otu), method = "bray")
specnumber.ra <- specnumber(relative_abundance_otu_t)

###### Five Community State Type (CSTs) Clustering

# species to CST mapping
species_to_cst <- data.frame(
  Species = c("Lactobacillus crispatus", "Lactobacillus gasseri", "Lactobacillus iners", "Lactobacillus jensenii"),
  CST = c("I", "II", "III", "V") # IV is anaerobic/diverse cluster
)

# Get most abundant OTU per sample
max_taxa <- apply(relative_abundance_otu, 2, function(sample) {
  taxa_idx <- which.max(sample)
  taxa_names(vaginal_relative_abundances)[taxa_idx]
})

# Map most abundant OTU to the sample data
bacteria_metadata_df$max_taxa <- max_taxa
# sample_data(bacteria_subset_vaginal)$max_taxa <- max_taxa
# samples.data <- sample_data(bacteria_subset_vaginal) %>%
#   group_by(biome_id) %>%
#   summarise(max_taxa=names(which.max(table(max_taxa)))) %>%
#   ungroup()
# samples.data$OTU <- bacteria_taxa_table[samples.data$max_taxa, "Species"]
# samples.data <- as.data.frame(samples.data)
# samples.data$OTU <- as.character(samples.data$OTU)
# Join most abundant OTU for each participant to the metadata
# bacteria_metadata_df <- bacteria_metadata_df %>% 
#   left_join(samples.data, by="biome_id")
# rownames(bacteria_metadata_df) <- bacteria_metadata_df$SampleID
# sample_data(bacteria_subset_Species) <- bacteria_metadata_df

bacteria_taxa_table <- tax_table(bacteria_subset_vaginal)
# Map max_taxa to name
bacteria_metadata_df$OTU <- bacteria_taxa_table[bacteria_metadata_df$max_taxa, "Species"]
# Set sample data in phyloseq obj
sample_data(bacteria_subset_vaginal) <- bacteria_metadata_df

# Add CST col
bacteria_taxa_table$CST <- ifelse(
  bacteria_taxa_table[,"Species"] %in% species_to_cst$Species,
  species_to_cst$CST[match(bacteria_taxa_table[,"Species"], species_to_cst$Species)],
  # Assign to CST IV if not in I, II, III, V
  "IV"  
)
table(bacteria_taxa_table$CST)
# Set taxa data in phyloseq obj

########################################################

# Aggregate data by participants by mean relative abundance for a given OTU
participant_otu <- tapply(sample_names(bacteria_subset_vaginal), 
                          sample_data(bacteria_subset_vaginal)$biome_id, 
                          function(samples) rowMeans(otu_table(vaginal_relative_abundances)[, samples, drop = FALSE]))
participant_otu <- do.call(cbind, participant_otu)
rownames(participant_otu) <- taxa_names(vaginal_relative_abundances)

########################################################
## Corr with diet?

# max rel OTU per sample
# table it - for cleaning purposes


