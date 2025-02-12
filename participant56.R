library(phyloseq)
library(tidyverse)
library(vegan)

# load("/Volumes/T7/microbiome_data/R_environments/microbiome_cleaningv1.RData")

bacteria_physeq <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/02_11/bacteria_intermediary2.rds")
samples_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_samples.csv", header=TRUE)

samples.data.56 <- samples_data %>% 
  filter(biome_id==56) %>% 
  filter(sampleType != "saliva")

otu56 <- otu_table(bacteria_physeq)
sample_data56 <- sample_data(bacteria_physeq)
tax_table56 <- tax_table(bacteria_physeq)

# Remove BLANK samples
sample_data56 <- sample_data56[!grepl("^BLANK", rownames(sample_data56)), ]
otu56 <- otu56[, !grepl("^BLANK", colnames(otu56))]

new_sample_names <- sub("\\_.*", "", rownames(sample_data56))
non_unique_samples <- names(table(new_sample_names))[table(new_sample_names) > 1]
# Remove these samples from both the sample metadata and OTU table
sample_data56 <- sample_data56[!new_sample_names %in% non_unique_samples, ]
otu56 <- otu56[, !new_sample_names %in% non_unique_samples]

# Rename sample IDs to QR
rownames(otu56) <- sub("\\_.*", "", rownames(otu56))
rownames(sample_data56) <- sub("\\_.*", "", rownames(sample_data56))

# Create a new phyloseq object with updated names
physeq_updated <- phyloseq(otu_table(otu56), sample_data(sample_data56), tax_table56)

## Prune
samples_to_keep <- samples.data.56$qr
physeq56 <- prune_samples(samples_to_keep, physeq_updated)

otu56.df <- as.data.frame(otu_table(physeq_updated))
bacteria_metadata_df <- sample_data(physeq_updated)

# Save R environment
save.image("/Volumes/T7/microbiome_data/R_environments/participant56.RData")

load("/Volumes/T7/microbiome_data/R_environments/participant56.RData")

## Alpha Div - Shannon Index
shannon.24 <- vegan::diversity(otu56.df, "shannon")

# bacteria_metadata_df <- as(bacteria_metadata_df, "data.frame")
shannon.df.24 <- data.frame("SampleID"=names(shannon.24), "shannon"=shannon.24)

# Add participant IDs from sample data | Merge the calculated Shannon diversity values with metadata
bacteria_metadata_df <- as(bacteria_metadata_df, "data.frame")
shannon.df.24 <- data.frame("qr"=names(shannon.24), "shannon"=shannon.24)
shannon.qr.merged.24 <- merge(shannon.df.24, bacteria_metadata_df, by="qr")
shannon.cst.qr.merged.24 <- shannon.qr.merged.24 %>% 
  filter(biome_id==56) %>% 
  # Filter the data to be within study days: 10-14 to 12-14
  filter(logDate > "2022-10-13" & logDate < "2022-12-15") # 2038 to 1971

write.csv(shannon.cst.qr.merged.24,
          file = "/Volumes/T7/microbiome_data/cleaned_data/person56.csv",
          row.names = FALSE)
