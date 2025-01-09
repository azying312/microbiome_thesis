library(phyloseq)
library(decontam)
library(tidyverse)
library(Matrix)

# bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/Walther-Antonio_Project_022_16S.rds")
bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/Walther-Antonio_Project_022_16S_SILVA138.rds")
uminn_data <- read.csv("/Volumes/T7/microbiome_data/Swabs with blood - Sheet1.csv", header=TRUE)
samples.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_samples.csv")

bacterial_otu_table <- otu_table(bacterial.data)
bacterial_tax_table <- tax_table(bacterial.data)

## Metadata

uminn_data <- uminn_data %>% 
  select(Sample.ID, Special.Notes) %>% 
  mutate(qr=sub("_.*", "", Sample.ID))

# Data Cleaning
# uminn_data_keep <- uminn_data %>% 
#   filter(str_detect(Special.Notes, "error") | str_detect(Special.Notes, "No swab in tube")) %>% 
#   select(Sample.ID, Special.Notes, qr)

metadata_bacteria <- data.frame(
  SampleID = colnames(bacterial_otu_table),
  is_blank = grepl("BLANK", colnames(bacterial_otu_table)),
  stringsAsFactors = FALSE
) %>%
  mutate(qr = sub("\\_.*", "", SampleID))

metadata_bacteria <- metadata_bacteria %>% 
  left_join(samples.data, by="qr")
# metadata_bacteria <- metadata_bacteria %>% 
#   mutate(sampleType=ifelse(is_blank, "BLANK", sampleType))
rownames(metadata_bacteria) <- metadata_bacteria$SampleID

# Set sample data
sample_data_obj <- sample_data(metadata_bacteria)

## Saving as new obj
bacteria_physeq <- phyloseq(bacterial_otu_table, sample_data_obj, bacterial_tax_table)

### Data Cleaning
bacteria_physeq_otu <- otu_table(bacteria_physeq)

## Filter out error data - filter samples
uminn_data_keep <- uminn_data %>% 
  filter(str_detect(Special.Notes, "error") | str_detect(Special.Notes, "No swab in tube")) %>% 
  select(Sample.ID, Special.Notes, qr)
otu_table_bacteria<- as.data.frame(t(bacteria_physeq_otu))
dim(otu_table_bacteria) # 3659 122403

sample_ids <- rownames(otu_table_bacteria)
sample_ids <- sub("_.*", "", sample_ids)
length(sample_ids) # 3659

# Prune samples that had errors
samples_to_keep <- !sample_ids %in% uminn_data_keep$qr
t_bacterial_otu_table <- t(bacterial_otu_table)
physeq_no_error <- prune_samples(samples_to_keep, t_bacterial_otu_table)

# Prep for decontam
physeq_no_error_otu <- otu_table(physeq_no_error)
# physeq_no_error_otu_t <- t(physeq_no_error_otu)
physeq_no_error_otu_df <- as.data.frame(physeq_no_error_otu)

# physeq_no_error_tax <- tax_table(physeq_no_error)

# Sample data
metadata_bacteria <- data.frame(
  SampleID = rownames(physeq_no_error_otu_df),
  is_blank = grepl("BLANK", rownames(physeq_no_error_otu_df)),
  stringsAsFactors = FALSE
) %>%
  mutate(qr = sub("\\_.*", "", SampleID))
metadata_bacteria <- metadata_bacteria %>% 
  left_join(samples.data, by="qr")
rownames(metadata_bacteria) <- metadata_bacteria$SampleID

# Convert back to phyloseq obj
otu_table_obj <- otu_table(physeq_no_error_otu_df, taxa_are_rows = FALSE)
sample_data_obj <- sample_data(metadata_bacteria)

otu_names <- taxa_names(otu_table_obj)
tax_names <- rownames(bacterial_tax_table)

## Saving as new obj
bacteria_physeq <- phyloseq(otu_table_obj, sample_data_obj, bacterial_tax_table)

# Save new obj
saveRDS(bacteria_physeq, file = "/Volumes/T7/microbiome_data/sequenced_data/vaginal_bacteria_intermediary.rds")

#
#
#
#

### Filter vaginal data
bacteria_physeq <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/vaginal_bacteria_intermediary.rds")
samples.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_samples.csv")

bacteria_physeq_otu <- otu_table(bacteria_physeq)
bacteria_physeq_tax <- tax_table(bacteria_physeq)
bacteria_physeq_meta <- sample_data(bacteria_physeq)

## Metadata
bacteria_metadata_df <- as.data.frame(as.matrix(sample_data(bacteria_physeq))) %>% 
  select(SampleID, is_blank, qr)

# Get vaginal swabs
vaginal_metadata <- bacteria_metadata_df %>% 
  left_join(samples.data, by="qr") %>% 
  filter(!is.na(biome_id)) %>%
  filter(sampleType=="vaginal")

# Subset samples
bacteria_subset_vaginal <- subset_samples(
  bacteria_physeq,
  SampleID %in% vaginal_metadata$SampleID
)

bacteria_physeq_otu <- otu_table(bacteria_subset_vaginal)
bacteria_physeq_tax <- tax_table(bacteria_subset_vaginal)
bacteria_physeq_meta <- sample_data(bacteria_subset_vaginal)

# Identify contaminants based on prevalence
contam_prev <- isContaminant(bacteria_subset_vaginal, method = "prevalence", neg = "is_blank")
contaminants <- contam_prev$contaminant
table(contam_prev$contaminant)

# bacteria_subset_vaginal_otu <- otu_table(bacteria_subset_vaginal)
# bacteria_subset_vaginal_tax <- tax_table(bacteria_subset_vaginal)
# bacteria_subset_vaginal_meta <- sample_data(bacteria_subset_vaginal)
# ntaxa(bacteria_subset_vaginal)

# Filter the OTUs in the phyloseq object
bacteria_physeq_no_contam <- prune_taxa(!contaminants, bacteria_subset_vaginal)
bacteria_physeq_no_contam <- phyloseq(otu_table(bacteria_physeq_no_contam), bacteria_physeq_meta, tax_table(bacteria_physeq_no_contam))

# Save new obj
saveRDS(bacteria_physeq_no_contam, file = "/Volumes/T7/microbiome_data/sequenced_data/vaginal_bacteria_decontam.rds")

## Data Processing

# Subset taxa; Filter ASVs w/o Species assignment
bacteria_physeq_no_contam <- subset_taxa(
  bacteria_physeq_no_contam,
  !is.na(tax_table(bacteria_physeq_no_contam)[, "Species"]) & 
    tax_table(bacteria_physeq_no_contam)[, "Species"] != ""
)

# Filter ASVs w/o Phylum asst
bacteria_physeq_subset <- subset_taxa(bacteria_physeq_no_contam, !is.na(Phylum) & Phylum != "")
# Filter ASVs w/o Genus asst
bacteria_physeq_subset <- subset_taxa(bacteria_physeq_subset, !is.na(Genus) & Genus != "")

# bacteria_metadata_df <- as.data.frame(as.matrix(sample_data(bacteria_subset_vaginal)))
bacteria_physeq_clean <- phyloseq(otu_table(bacteria_physeq_subset), bacteria_physeq_meta, tax_table(bacteria_physeq_subset))

# Save new obj
saveRDS(bacteria_physeq_clean, file = "/Volumes/T7/microbiome_data/sequenced_data/vaginal_bacteria_cleaned.rds")

# # ## Blanks
# # bacteria_metadata_df <- as.data.frame(as.matrix(sample_data(bacteria_physeq))) %>% 
# #   select(SampleID, is_blank, qr)
# # 
# # blanks_df <- bacteria_metadata_df %>% 
# #   filter(is_blank==TRUE)
# 
# ### Vaginal Data
# 
# # bacteria_metadata_df <- as.data.frame(as.matrix(sample_data(bacteria_physeq))) %>% 
# #   select(SampleID, is_blank, qr)
# 
# # Get vaginal swabs
# vaginal_metadata <- bacteria_metadata_df %>% 
#   left_join(samples.data, by="qr") %>% 
#   filter(!is.na(biome_id)) %>%
#   filter(sampleType=="vaginal")
# 
# # Subset samples
# bacteria_subset_vaginal <- subset_samples(
#   bacteria_physeq,
#   SampleID %in% vaginal_metadata$SampleID
# )
# 
# ### Data Cleaning
# 
# uminn_data_keep <- uminn_data %>% 
#   filter(str_detect(Special.Notes, "error") | str_detect(Special.Notes, "No swab in tube")) %>% 
#   select(Sample.ID, Special.Notes, qr)
# otu_table_vaginal <- as.data.frame(t(otu_table(bacteria_subset_vaginal)))
# dim(otu_table_vaginal)
# 
# sample_ids <- rownames(otu_table_vaginal)
# sample_ids <- sub("_.*", "", sample_ids)
# length(sample_ids)
# 
# samples_to_keep <- !sample_ids %in% uminn_data_keep$qr
# physeq_no_error <- prune_samples(samples_to_keep, t(otu_table(bacteria_subset_vaginal)))
# 
# # Subset sample data
# noerr_phy_otu <- otu_table(physeq_no_error)
# physeq_sample_names <- rownames(noerr_phy_otu)
# vaginal_metadata_subset <- vaginal_metadata %>% 
#   filter(SampleID %in% physeq_sample_names)
# 
# # Set sample data
# rownames(vaginal_metadata_subset) <- vaginal_metadata_subset$SampleID
# sample_data_obj <- sample_data(vaginal_metadata_subset)
# 
# # Taxa table
# # noerr_phy_tax <- tax_table(bacteria_subset_vaginal)
# 
# ##################
# 
# ## Saving as new obj
# physeq_no_error <- phyloseq(noerr_phy_otu, sample_data_obj, tax_table(bacteria_subset_vaginal))
# 
# ## Saving as new obj
# # physeq_no_error <- phyloseq(physeq_no_error, sample_data(vaginal_metadata_subset))
# 
# # Filter OTUs/ASVs - decontam pkg
# 
# # no_error_otu <- otu_table(physeq_no_error)
# # no_error_otu <- t(no_error_otu)
# # otu_table_bacteria <- as.data.frame(no_error_otu)
# 
# ##########################################################################################
# 
# meta_data <- sample_data(physeq_no_error)
# 
# # Identify contaminants based on prevalence
# contam_prev <- isContaminant(physeq_no_error, method = "prevalence", neg = "is_blank")
# contaminants <- contam_prev$contaminant
# table(contam_prev$contaminant)
# 
# # Filter the OTUs in the phyloseq object
# bacteria_physeq_no_contam <- prune_taxa(!contaminants, bacteria_subset_vaginal)
# bacteria_physeq_no_contam <- phyloseq(otu_table(bacteria_physeq_no_contam), sample_data_obj, tax_table(bacteria_physeq_no_contam))
# 
# ## Data Processing
# 
# # Subset taxa; Filter ASVs w/o Species assignment
# bacteria_subset_vaginal <- subset_taxa(
#   bacteria_subset_vaginal,
#   !is.na(tax_table(bacteria_subset_vaginal)[, "Species"]) & 
#     tax_table(bacteria_subset_vaginal)[, "Species"] != ""
# )
# 
# bacteria_metadata_df <- as.data.frame(as.matrix(sample_data(bacteria_subset_vaginal)))
# 
# ## Saving as new obj
# # bacteria_physeq <- phyloseq(bacteria_subset_vaginal, sample_data_obj)
# 
# # Filter ASVs w/o Phylum asst
# bacteria_physeq_subset <- subset_taxa(bacteria_physeq_no_contam, !is.na(Phylum) & Phylum != "")
# # Filter ASVs w/o Genus asst
# bacteria_physeq_subset <- subset_taxa(bacteria_physeq_subset, !is.na(Genus) & Genus != "")
# 
# bacteria_physeq_clean <- phyloseq(otu_table(bacteria_physeq_subset), sample_data_obj, tax_table(bacteria_physeq_subset))
# 
# # Save new obj
# saveRDS(bacteria_physeq_clean, file = "/Volumes/T7/microbiome_data/sequenced_data/vaginal_bacteria_cleaned.rds")
# 
# ### Data Cleaning
# 
# ## Filter out error data - filter samples
# # uminn_data_keep <- uminn_data %>% 
# #   filter(str_detect(Special.Notes, "error") | str_detect(Special.Notes, "No swab in tube")) %>% 
# #   select(Sample.ID, Special.Notes, qr)
# # otu_table_bacteria<- as.data.frame(t(bacterial_otu_table))
# # dim(otu_table_bacteria)
# # 
# # sample_ids <- rownames(otu_table_bacteria)
# # sample_ids <- sub("_.*", "", sample_ids)
# # length(sample_ids)
# # 
# # samples_to_keep <- !sample_ids %in% uminn_data_keep$qr
# # physeq_no_error <- prune_samples(samples_to_keep, t(bacterial_otu_table))
# 
# ## Filter OTUs/ASVs - decontam pkg
# # no_error_otu <- otu_table(physeq_no_error)
# # no_error_otu <- t(no_error_otu)
# # otu_table_bacteria<- as.data.frame(no_error_otu)
# 
# # metadata_bacteria <- data.frame(
# #   SampleID = colnames(otu_table_bacteria),
# #   is_blank = grepl("BLANK", colnames(otu_table_bacteria)),
# #   stringsAsFactors = FALSE
# # ) %>%
# #   mutate(qr = sub("\\_.*", "", SampleID))
# # metadata_bacteria <- metadata_bacteria %>% 
# #   left_join(samples.data, by="qr")
# # rownames(metadata_bacteria) <- metadata_bacteria$SampleID
# 
# # # Convert back to phyloseq obj
# # otu_table_obj <- otu_table(otu_table_bacteria, taxa_are_rows = TRUE)
# # sample_data_obj <- sample_data(metadata_bacteria)
# # 
# # ## Saving as new obj
# # bacteria_physeq <- phyloseq(otu_table_obj, sample_data_obj, bacterial_tax_table)
# 
# # ### Filter vaginal data
# # 
# # ## Metadata
# # bacteria_metadata_df <- as.data.frame(as.matrix(sample_data(bacteria_physeq))) %>% 
# #   select(SampleID, is_blank, qr)
# # 
# # # Get vaginal swabs
# # vaginal_metadata <- bacteria_metadata_df %>% 
# #   left_join(samples.data, by="qr") %>% 
# #   filter(!is.na(biome_id)) %>%
# #   filter(sampleType=="vaginal")
# # 
# # # Subset samples
# # bacteria_subset_vaginal <- subset_samples(
# #   bacteria_physeq,
# #   SampleID %in% vaginal_metadata$SampleID
# # )
# # 
# # # Subset taxa; Filter ASVs w/o Species assignment
# # bacteria_subset_vaginal <- subset_taxa(
# #   bacteria_subset_vaginal,
# #   !is.na(tax_table(bacteria_subset_vaginal)[, "Species"]) & 
# #     tax_table(bacteria_subset_vaginal)[, "Species"] != ""
# # )
# # 
# # bacteria_metadata_df <- as.data.frame(as.matrix(sample_data(bacteria_subset_vaginal)))
# # 
# # ## Saving as new obj
# # # bacteria_physeq <- phyloseq(otu_table_obj, sample_data_obj, bacterial_tax_table)
# 
# # # Identify contaminants based on prevalence
# # contam_prev <- isContaminant(bacteria_subset_vaginal, method = "prevalence", neg = "is_blank")
# # contaminants <- contam_prev$contaminant
# # table(contam_prev$contaminant)
# # 
# # # Filter the OTUs in the phyloseq object
# # bacteria_physeq_no_contam <- prune_taxa(!contaminants, bacteria_subset_vaginal)
# # bacteria_physeq_no_contam <- phyloseq(otu_table(bacteria_physeq_no_contam), sample_data_obj, tax_table(bacteria_physeq_no_contam))
# # 
# # ## Data Processing
# # 
# # # Filter ASVs w/o Phylum asst
# # bacteria_physeq_subset <- subset_taxa(bacteria_physeq_no_contam, !is.na(Phylum) & Phylum != "")
# # # Filter ASVs w/o Genus asst
# # bacteria_physeq_subset <- subset_taxa(bacteria_physeq_subset, !is.na(Genus) & Genus != "")
# # 
# # bacteria_physeq_clean <- phyloseq(otu_table(bacteria_physeq_subset), sample_data_obj, tax_table(bacteria_physeq_subset))
# # 
# # # Save new obj
# # saveRDS(bacteria_physeq_clean, file = "/Volumes/T7/microbiome_data/sequenced_data/vaginal_bacteria_cleaned.rds")
