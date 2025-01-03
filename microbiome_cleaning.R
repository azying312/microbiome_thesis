library(phyloseq)
library(decontam)
library(tidyverse)
library(Matrix)

bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/Walther-Antonio_Project_022_16S.rds")
bacterial_otu_table <-otu_table(bacterial.data)
bacterial_tax_table<-tax_table(bacterial.data)
# uminn_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_uminn_data.csv", header=TRUE)
uminn_data <- read.csv("/Volumes/T7/microbiome_data/Swabs with blood - Sheet1.csv", header=TRUE)
samples.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_samples.csv")

uminn_data <- uminn_data %>% 
  select(Sample.ID, Special.Notes) %>% 
  mutate(qr=sub("_.*", "", Sample.ID))

### Data Cleaning

## Filter out error data - filter samples
uminn_data_keep <- uminn_data %>% 
  filter(str_detect(Special.Notes, "error") | str_detect(Special.Notes, "No swab in tube")) %>% 
  select(Sample.ID, Special.Notes, qr)
otu_table_bacteria<- as.data.frame(t(bacterial_otu_table))
dim(otu_table_bacteria)

sample_ids <- rownames(otu_table_bacteria)
sample_ids <- sub("_.*", "", sample_ids)
length(sample_ids)

samples_to_keep <- !sample_ids %in% uminn_data_keep$qr
physeq_no_error <- prune_samples(samples_to_keep, t(bacterial_otu_table))

## Filter OTUs/ASVs - decontam pkg
otu_table_bacteria<- as.data.frame(t(otu_table(physeq_no_error)))

metadata_bacteria <- data.frame(
  SampleID = colnames(otu_table_bacteria),
  is_blank = grepl("BLANK", colnames(otu_table_bacteria)),
  stringsAsFactors = FALSE
) %>%
  mutate(qr = sub("\\_.*", "", SampleID))
metadata_bacteria <- metadata_bacteria %>% 
  left_join(samples.data, by="qr")
rownames(metadata_bacteria) <- metadata_bacteria$SampleID

# Convert back to phyloseq obj
otu_table_obj <- otu_table(otu_table_bacteria, taxa_are_rows = TRUE)
sample_data_obj <- sample_data(metadata_bacteria)

## Saving as new obj
bacteria_physeq <- phyloseq(otu_table_obj, sample_data_obj, bacterial_tax_table)

### Filter vaginal data

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

# Subset taxa; Filter ASVs w/o Species assignment
bacteria_subset_vaginal <- subset_taxa(
  bacteria_subset_vaginal,
  !is.na(tax_table(bacteria_subset_vaginal)[, "Species"]) & 
    tax_table(bacteria_subset_vaginal)[, "Species"] != ""
)

bacteria_metadata_df <- as.data.frame(as.matrix(sample_data(bacteria_subset_vaginal)))

## Saving as new obj
# bacteria_physeq <- phyloseq(otu_table_obj, sample_data_obj, bacterial_tax_table)

# Identify contaminants based on prevalence
contam_prev <- isContaminant(bacteria_subset_vaginal, method = "prevalence", neg = "is_blank")
contaminants <- contam_prev$contaminant
table(contam_prev$contaminant)

# Filter the OTUs in the phyloseq object
bacteria_physeq_no_contam <- prune_taxa(!contaminants, bacteria_subset_vaginal)
bacteria_physeq_no_contam <- phyloseq(otu_table(bacteria_physeq_no_contam), sample_data_obj, tax_table(bacteria_physeq_no_contam))

## Data Processing

# Filter ASVs w/o Phylum asst
bacteria_physeq_subset <- subset_taxa(bacteria_physeq_no_contam, !is.na(Phylum) & Phylum != "")
# Filter ASVs w/o Genus asst
bacteria_physeq_subset <- subset_taxa(bacteria_physeq_subset, !is.na(Genus) & Genus != "")

bacteria_physeq_clean <- phyloseq(otu_table(bacteria_physeq_subset), sample_data_obj, tax_table(bacteria_physeq_subset))

# Save new obj
saveRDS(bacteria_physeq_clean, file = "/Volumes/T7/microbiome_data/sequenced_data/vaginal_bacteria_cleaned.rds")
