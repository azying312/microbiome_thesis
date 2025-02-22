library(phyloseq)
library(decontam)
library(tidyverse)
library(Matrix)


bacteria_physeq <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/02_11/bacteria_intermediary2.rds")
samples.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_samples.csv")


bacteria_physeq_otu <- otu_table(bacteria_physeq)
bacteria_physeq_tax <- tax_table(bacteria_physeq)
bacteria_physeq_meta <- sample_data(bacteria_physeq)

## Metadata
bacteria_metadata_df <- as.data.frame(as.matrix(sample_data(bacteria_physeq))) %>% 
  select(SampleID, is_blank, qr)

################################################################################
# Get fecal swabs
fecal_metadata <- bacteria_metadata_df %>% 
  left_join(samples.data, by="qr") %>% 
  filter(!is.na(biome_id)) %>%
  filter(sampleType=="fecal")


# Subset samples - fecal
bacteria_subset_fecal <- subset_samples(
  bacteria_physeq,
  SampleID %in% fecal_metadata$SampleID
)

bacteria_physeq_otu <- otu_table(bacteria_subset_fecal)
bacteria_physeq_tax <- tax_table(bacteria_subset_fecal)
bacteria_physeq_meta <- sample_data(bacteria_subset_fecal)

# Identify contaminants based on prevalence
contam_prev <- isContaminant(bacteria_subset_fecal, method = "prevalence", neg = "is_blank")
contaminants <- contam_prev$contaminant
table(contam_prev$contaminant)


# Filter the OTUs in the phyloseq object
bacteria_physeq_no_contam <- prune_taxa(!contaminants, bacteria_subset_fecal)
bacteria_physeq_no_contam <- phyloseq(otu_table(bacteria_physeq_no_contam), bacteria_physeq_meta, tax_table(bacteria_physeq_no_contam))

# Save new obj
saveRDS(bacteria_physeq_no_contam, file = "/Volumes/T7/microbiome_data/sequenced_data/fecal_bacteria_decontam.rds")


## Data Processing

bacteria_physeq_no_contam <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/fecal_bacteria_decontam.rds")


# Filter ASVs w/o Phylum asst
bacteria_physeq_subset <- subset_taxa(bacteria_physeq_no_contam, !is.na(Phylum) & Phylum != "")
# Filter ASVs w/o Genus asst
bacteria_physeq_subset <- subset_taxa(bacteria_physeq_subset, !is.na(Genus) & Genus != "")

# bacteria_metadata_df <- as.data.frame(as.matrix(sample_data(bacteria_subset_vaginal)))
bacteria_physeq_clean <- phyloseq(otu_table(bacteria_physeq_subset), bacteria_physeq_meta, tax_table(bacteria_physeq_subset))

# Save new obj
saveRDS(bacteria_physeq_clean, file = "/Volumes/T7/microbiome_data/sequenced_data/fecal_bacteria_cleaned.rds")

