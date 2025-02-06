library(dplyr)
library(phyloseq)
library(vegan)
library(pheatmap)
library(tidyverse)
library(Matrix)

library(ggplot2)

library(cluster)
library("igraph")
library("markovchain")
library("RColorBrewer")
library("gridExtra")

# Use corrected data after
bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/vaginal_bacteria_cleaned.rds")

otu.22 <- otu_table(bacterial.data)
tax.22 <- tax_table(bacterial.data)
meta.22 <- sample_data(bacterial.data)

head(otu.22[1:10,1:10])
dim(otu.22)

# Initial Reads
total_reads_initial <- sum(otu.22) # 173324711

# number of reads per taxa
taxa_total_reads <- taxa_sums(bacterial.data)
taxa_sample_presence <- colSums((otu.22) > 0) # 0 - taxon absent in sample; 1 - taxon in sample -- then sum

## Filter data
# Identify non-rare taxa to keep (Callahan paper: present in at least 3 samples and ≥100 total reads)
keep_taxa <- (taxa_total_reads >= 100) & (taxa_sample_presence >= 3)
bacterial.data.filtered <- prune_taxa(keep_taxa, bacterial.data)

# Save new obj
saveRDS(bacterial.data.filtered, file = "/Volumes/T7/microbiome_data/sequenced_data/vaginal_bacteria_cleanedv2.rds")

total_reads_after <- sum(otu_table(bacterial.data.filtered))
total_reads_after # 172889277 (435434 reads removed)

# Percent removed
total_reads_initial-total_reads_after
100*(total_reads_initial-total_reads_after)/total_reads_initial

# Number of taxa removed
num_taxa_before <- ntaxa(bacterial.data)
num_taxa_after <- ntaxa(bacterial.data.filtered)
taxa_removed <- num_taxa_before - num_taxa_after
taxa_removed
100*(taxa_removed)/num_taxa_before

# Figure: Visualize removed taxa
reads_before <- taxa_sums(bacterial.data)
reads_after <- taxa_sums(bacterial.data.filtered)

taxa.df <- data.frame(
  reads = c(reads_before, reads_after),
  status = rep(c("Before Filtering", "After Filtering"), 
               c(length(reads_before), length(reads_after)))
)

ggplot(taxa.df, aes(x = reads, fill = status)) +
  geom_histogram(data = subset(taxa.df, status == "Before Filtering"), bins = 50, alpha = 0.3) + 
  geom_histogram(data = subset(taxa.df, status == "After Filtering"), bins = 50, alpha = 0.7) +  
  # geom_histogram(bins = 50, alpha = 0.4, position = "identity") +
  scale_x_log10() +
  theme_minimal() +
  scale_fill_manual(values = c("Before Filtering" = "cyan", "After Filtering" = "green")) +
  labs(title = "Distribution of Read Counts Before and After Filtering",
       x = "Number of Reads per Taxon (log scale)",
       y = "Number of Taxa",
       fill = "Status")

############################################################

### Get Species variable cleaned in new data
bacteria_taxa_table <- tax_table(bacterial.data.filtered)
bacteria_taxa_df <- as.data.frame(bacteria_taxa_table)
dim(bacteria_taxa_df) # 3830 10

length(bacteria_taxa_df$Species) # 3830
# number of missing Species/Species_exact
sum(bacteria_taxa_df$Species=="") # 1871
sum(bacteria_taxa_df$Species_exact=="") # 3162

bacteria_taxa_df <- bacteria_taxa_df %>% 
  mutate(Species_exact=ifelse(Species_exact=="", Species, Species_exact))
sum(bacteria_taxa_df$Species_exact=="") # 1716 missing
sum(bacteria_taxa_df$Species_exact!="") # 2114 not missing

# Subset taxa; Filter ASVs w/o Species_exact assignment
bacterial.data_subset <- subset_taxa(
  bacterial.data.filtered,
  Species_exact != ""
)

bacteria_taxa_table <- tax_table(bacterial.data_subset)
bacteria_taxa_df <- as.data.frame(bacteria_taxa_table)
dim(bacteria_taxa_df) # 668 OTU

# After Species filtering
total_reads_filtered <- sum(otu_table(bacterial.data.filtered))
total_reads_Species <- sum(otu_table(bacterial.data_subset))

# Percent removed
total_reads_filtered-total_reads_Species
100*(total_reads_filtered-total_reads_Species)/total_reads_filtered

# Save new obj
saveRDS(bacterial.data_subset, file = "/Volumes/T7/microbiome_data/sequenced_data/vaginal_bacteria_cleanedv3.rds")