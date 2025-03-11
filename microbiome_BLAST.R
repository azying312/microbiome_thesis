library(vegan)
library(pheatmap)
library(tidyverse)

source("~/Microbiome Thesis/functions.R")

# RELABELED DATA
bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/relabeled_data/vaginal_bacteria_cleanedv3.rds")
fecal.bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/relabeled_data/fecal_bacteria_cleanedv3.rds")

# bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/vaginal_bacteria_cleanedv3.rds")
# BLAST of max taxa
BLAST_taxa <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/max taxa BLAST.csv")

###############################################################################################

#### Exploratory Data Analysis

# Relative abundances
vaginal_relative_abundances <- transform_sample_counts(bacterial.data, function(x) x/sum(x))
relative_abundance_otu <- as.data.frame(otu_table(vaginal_relative_abundances))
relative_abundance_otu_t <- t(relative_abundance_otu) %>% as.data.frame()

# Add sample ID
relative_abundance_otu$SampleID <- rownames(relative_abundance_otu)
bacteria_metadata_df <- sample_data(bacterial.data)
bacteria_metadata_df <- bacteria_metadata_df[,-2]
otu_with_participant <-  relative_abundance_otu %>%
  left_join(bacteria_metadata_df, by="SampleID")
rownames(otu_with_participant) <- otu_with_participant$SampleID
relative_abundance_otu <- relative_abundance_otu %>% 
  select(!SampleID)

###### Five Community State Type (CSTs) Clustering

# species to CST mapping
species_to_cst <- data.frame(
  Species = c("crispatus", "gasseri", "iners", "jensenii"),
  CST = c("I", "II", "III", "V") # IV is anaerobic/diverse cluster
)

# Get most abundant OTU per sample
# BEFORE RELABELED
# max_taxa <- apply(relative_abundance_otu_t, 2, function(sample) {
# RELABELED
max_taxa <- apply(relative_abundance_otu_t, 1, function(sample) {
  taxa_idx <- which.max(sample)
  taxa_names(vaginal_relative_abundances)[taxa_idx]
})

# Map most abundant OTU to the sample data
bacteria_metadata_df$max_taxa <- max_taxa

###################################################################################################

bacteria_taxa_table <- tax_table(bacterial.data)
bacteria_taxa_df <- as.data.frame(bacteria_taxa_table)

################## FASTA
# fasta.file <- "/Volumes/T7/microbiome_data/sequenced_data/otu_seq.fasta"
fasta.file <- "/Volumes/T7/microbiome_data/sequenced_data/relabeled_data/otu_seq.fasta"
# Get unique max taxa to BLAST
blast_taxa <- unique(bacteria_metadata_df$max_taxa)
blast_sequences <- bacteria_taxa_df[blast_taxa, "sequence", drop=FALSE]
sink(fasta.file)
for (otu_id in rownames(blast_sequences)) {
  cat(paste0(">", otu_id, "\n", blast_sequences[otu_id, "sequence"], "\n"))
}
sink()
############################################################

# Set taxa to taxa from BLAST
tax.22 <- bacteria_taxa_df %>% 
  mutate(OTU_name=rownames(bacteria_taxa_df))
tax.22 <- tax.22 %>% 
  left_join(BLAST_taxa, by=c("OTU_name", "sequence"))
tax.22 <- tax.22 %>% 
  mutate(BLAST_species = ifelse(is.na(BLAST_species), Species_exact, BLAST_species))
rownames(tax.22) <- taxa_names(bacterial.data)
tax.22 <- tax.22 %>% 
  select(-OTU_name)
tax_table(bacterial.data) <- as.matrix(tax.22)

vaginal_max_taxa <- phyloseq(otu_table(bacterial.data), sample_data(bacterial.data), tax_table(bacterial.data))

# Save new obj
# saveRDS(vaginal_max_taxa, file = "/Volumes/T7/microbiome_data/sequenced_data/vaginal_cleaned_max_taxa.rds")

# RELABELED VERSION
saveRDS(vaginal_max_taxa, file = "/Volumes/T7/microbiome_data/sequenced_data/relabeled_data/vaginal_cleaned_max_taxa.rds")

###############################################################################################

## Set fecal/gut data

#### Exploratory Data Analysis

# Relative abundances
fecal_relative_abundances <- transform_sample_counts(fecal.bacterial.data, function(x) x/sum(x))
relative_abundance_otu <- as.data.frame(otu_table(fecal_relative_abundances))
relative_abundance_otu_t <- t(relative_abundance_otu) %>% as.data.frame()

# Add sample ID
relative_abundance_otu$SampleID <- rownames(relative_abundance_otu)
bacteria_metadata_df <- sample_data(fecal.bacterial.data)
bacteria_metadata_df <- bacteria_metadata_df[,-2]
otu_with_participant <-  relative_abundance_otu %>%
  left_join(bacteria_metadata_df, by="SampleID")
rownames(otu_with_participant) <- otu_with_participant$SampleID
relative_abundance_otu <- relative_abundance_otu %>% 
  select(!SampleID)

max_taxa <- apply(relative_abundance_otu_t, 1, function(sample) {
  taxa_idx <- which.max(sample)
  taxa_names(fecal.bacterial.data)[taxa_idx]
})

# Map most abundant OTU to the sample data
bacteria_metadata_df$max_taxa <- max_taxa

bacteria_taxa_table <- tax_table(fecal.bacterial.data)
bacteria_taxa_df <- as.data.frame(bacteria_taxa_table)

# Set taxa to taxa from BLAST
tax.22 <- bacteria_taxa_df %>% 
  mutate(OTU_name=rownames(bacteria_taxa_df))
tax.22 <- tax.22 %>% 
  left_join(BLAST_taxa, by=c("OTU_name", "sequence"))
tax.22 <- tax.22 %>% 
  mutate(BLAST_species = ifelse(is.na(BLAST_species), Species_exact, BLAST_species))
rownames(tax.22) <- taxa_names(fecal.bacterial.data)
tax.22 <- tax.22 %>% 
  select(-OTU_name)
tax_table(fecal.bacterial.data) <- as.matrix(tax.22)

fecal_max_taxa <- phyloseq(otu_table(fecal.bacterial.data), sample_data(fecal.bacterial.data), tax_table(fecal.bacterial.data))

# Save new obj
saveRDS(fecal_max_taxa, file = "/Volumes/T7/microbiome_data/sequenced_data/relabeled_data/fecal_cleaned_max_taxa.rds")


