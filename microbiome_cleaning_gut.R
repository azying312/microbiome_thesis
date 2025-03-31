library(phyloseq)
library(decontam)
library(tidyverse)
# library(Matrix)

# after vaginal_gut_sample_check.R
bacteria_physeq <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/02_11/bacteria_cleanedv2.rds")
samples.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/relabeled_data/cleaned_samplesv2.csv")

# FIRST RUN (NO RELABELING)
# bacteria_physeq <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/02_11/bacteria_intermediary2.rds")
# samples.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_samples.csv")

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
# saveRDS(bacteria_physeq_no_contam, file = "/Volumes/T7/microbiome_data/sequenced_data/fecal_bacteria_decontam.rds")

## Data Processing
# bacteria_physeq_no_contam <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/fecal_bacteria_decontam.rds")

# Filter ASVs w/o Phylum asst
bacteria_physeq_subset <- subset_taxa(bacteria_physeq_no_contam, !is.na(Phylum) & Phylum != "")
# Filter ASVs w/o Genus asst
bacteria_physeq_subset <- subset_taxa(bacteria_physeq_subset, !is.na(Genus) & Genus != "")

# bacteria_metadata_df <- as.data.frame(as.matrix(sample_data(bacteria_subset_vaginal)))
bacteria_physeq_clean <- phyloseq(otu_table(bacteria_physeq_subset), bacteria_physeq_meta, tax_table(bacteria_physeq_subset))

# Save new obj
# saveRDS(bacteria_physeq_clean, file = "/Volumes/T7/microbiome_data/sequenced_data/fecal_bacteria_cleanedv2.rds")

################################################################################

# SKIP FILTERING
bacterial.data.filtered <- bacteria_physeq_clean
# 
# # Filtering
# bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/fecal_bacteria_cleanedv2.rds")
# 
# otu.22 <- otu_table(bacterial.data)
# tax.22 <- tax_table(bacterial.data)
# meta.22 <- sample_data(bacterial.data)
# 
# # Initial Reads
# total_reads_initial <- sum(otu.22) # 138697737 reads
# 
# # number of reads per taxa
# taxa_total_reads <- taxa_sums(bacterial.data)
# taxa_sample_presence <- colSums((otu.22) > 0) # 0 - taxon absent in sample; 1 - taxon in sample -- then sum
# 
# ## Filter data
# # Identify non-rare taxa to keep (present in at least 3 samples)
# keep_taxa <- (taxa_sample_presence >= 2) # 2 samples
# # keep_taxa <- (taxa_sample_presence >= 3)
# bacterial.data.filtered <- prune_taxa(keep_taxa, bacterial.data)
# 
# # Save new obj
# saveRDS(bacterial.data.filtered, file = "/Volumes/T7/microbiome_data/sequenced_data/fecal_bacteria_filtered.rds")
# 
# total_reads_after <- sum(otu_table(bacterial.data.filtered))
# total_reads_after # 138222649 (475088 reads removed)
# 
# # Percent removed
# total_reads_initial-total_reads_after # 475088
# 100*(total_reads_initial-total_reads_after)/total_reads_initial # 0.343% removed
# 
# # Number of taxa removed
# num_taxa_before <- ntaxa(bacterial.data)
# num_taxa_after <- ntaxa(bacterial.data.filtered)
# taxa_removed <- num_taxa_before - num_taxa_after
# taxa_removed # 76222
# 100*(taxa_removed)/num_taxa_before # 82.403%

############################################################

### Get Species variable cleaned in new data
bacteria_taxa_table <- tax_table(bacterial.data.filtered)
bacteria_taxa_df <- as.data.frame(bacteria_taxa_table)
dim(bacteria_taxa_df) # 16277 10

length(bacteria_taxa_df$Species) # 16277
# number of missing Species/Species_exact
sum(bacteria_taxa_df$Species=="") # 6940
sum(bacteria_taxa_df$Species_exact=="") # 15262

# fill species_exact with species when blank
bacteria_taxa_df <- bacteria_taxa_df %>% 
  mutate(Species_exact=ifelse(Species_exact=="", Species, Species_exact))
sum(bacteria_taxa_df$Species_exact=="") # 6696 missing
sum(bacteria_taxa_df$Species_exact!="") # 9581 not missing

# Subset taxa; Filter ASVs w/o Species_exact assignment
bacterial.data_subset <- subset_taxa(
  bacterial.data.filtered,
  Species_exact != ""
)

bacteria_taxa_table <- tax_table(bacterial.data_subset)
bacteria_taxa_df <- as.data.frame(bacteria_taxa_table)
dim(bacteria_taxa_df) # 1015 OTU

# After Species filtering
total_reads_filtered <- sum(otu_table(bacterial.data.filtered))
total_reads_Species <- sum(otu_table(bacterial.data_subset))

# Percent removed
total_reads_filtered-total_reads_Species # 50122517 removed
100*(total_reads_filtered-total_reads_Species)/total_reads_filtered # 36.262% removed

bacteria_metadata <- sample_data(bacterial.data_subset)
bacteria_metadata_df <- as.data.frame(bacteria_metadata)
dim(bacteria_metadata_df)
head(bacteria_metadata_df)
length(unique(bacteria_metadata_df$biome_id))

# Save new obj
# BEFORE RE-LABEL DATA
# saveRDS(bacterial.data_subset, file = "/Volumes/T7/microbiome_data/sequenced_data/fecal_bacteria_filteredv2.rds")

# NEW LABELED DATA
saveRDS(bacterial.data_subset, file = "/Volumes/T7/microbiome_data/sequenced_data/relabeled_data/fecal_bacteria_cleanedv3.rds")
