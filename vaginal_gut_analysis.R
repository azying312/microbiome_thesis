########################
#
# Vaginal and Gut Microbiome Analysis with Lifestyle Factors
# Last updated: 03/23/2025
#
#########################

library(tidyverse)
library(vegan)
library(viridis)
library(phyloseq)

source("~/Microbiome Thesis/functions.R")

gut.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/relabeled_data/fecal_cleaned_max_taxa.rds")
vag.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/relabeled_data/vaginal_cleaned_max_taxa.rds")

###############################################################################################

gut.metadata <- sample_data(gut.data)
vag.metadata <- sample_data(vag.data)

# Extract OTU tables and metadata
gut_otu <- otu_table(gut.data)
vaginal_otu <- otu_table(vag.data)

# Get one sample per day for each participant
gut_otu_with_date <- as.data.frame(t(gut_otu)) %>% 
  mutate(logDate = gut.metadata$logDate, 
         biome_id = gut.metadata$biome_id,
         rownames = gut.metadata$SampleID)

gut_otu_distinct <- gut_otu_with_date %>% 
  group_by(biome_id, logDate) %>% 
  slice(1) %>% 
  ungroup() %>% 
  column_to_rownames("rownames")

vaginal_otu_with_date <- as.data.frame(t(vaginal_otu)) %>% 
  mutate(logDate = vag.metadata$logDate, 
         biome_id = vag.metadata$biome_id,
         rownames = vag.metadata$SampleID)

vag_otu_distinct <- vaginal_otu_with_date %>% 
  group_by(biome_id, logDate) %>% 
  slice(1) %>% 
  ungroup() %>% 
  column_to_rownames("rownames")

# Filter otu data
overlap_ids <- 
  intersect(paste(gut_otu_distinct$biome_id, gut_otu_distinct$logDate, sep="_"), 
            paste(vag_otu_distinct$biome_id, vag_otu_distinct$logDate, sep="_"))

gut_otu_distinct$key <- paste(gut_otu_distinct$biome_id, gut_otu_distinct$logDate, sep = "_")
vag_otu_distinct$key <- paste(vag_otu_distinct$biome_id, vag_otu_distinct$logDate, sep = "_")
gut_otu_matched <- gut_otu_distinct[gut_otu_distinct$key %in% overlap_ids, ]
vag_otu_matched <- vag_otu_distinct[vag_otu_distinct$key %in% overlap_ids, ]

dim(gut_otu_matched)
dim(vag_otu_matched)

# Filter metadata
gut_metadata_filtered <- gut.metadata[
  paste(gut.metadata$biome_id, gut.metadata$logDate, sep = "_") %in% overlap_ids, 
]

vag_metadata_filtered <- vag.metadata[
  paste(vag.metadata$biome_id, vag.metadata$logDate, sep = "_") %in% overlap_ids, 
]

# Compute cor matrix
gut_matrix <- as.matrix(gut_otu_matched[, !(colnames(gut_otu_matched) %in% c("logDate", "biome_id", "key"))])
vag_matrix <- as.matrix(vag_otu_matched[, !(colnames(vag_otu_matched) %in% c("logDate", "biome_id", "key"))])

gut_matrix_rownames <- rownames(gut_matrix)
vag_matrix_rownames <- rownames(vag_matrix)

rownames(gut_matrix) <- paste(gut_otu_matched$biome_id, gut_otu_matched$logDate, sep = "_")
rownames(vag_matrix) <- paste(vag_otu_matched$biome_id, vag_otu_matched$logDate, sep = "_")

cor_matrix <- cor(gut_matrix, vag_matrix, method = "spearman", use = "pairwise.complete.obs")

gut_matrix[1, 15]
cor(gut_matrix[1:15, 1:15], vag_matrix[1:15, 1:15])

##

# Heatmap of cor matrix
library(pheatmap)

cor_matrix <- as.matrix(cor_matrix)
pheatmap(cor_matrix, 
         clustering_distance_rows = "euclidean",  
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         display_numbers = TRUE,
         number_format = "%.2f",  
         main = "Correlation Heatmap of Gut vs Vaginal Microbiome OTUs")

