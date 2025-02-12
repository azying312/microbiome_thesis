# install.packages("BiocManager")
# library(BiocManager)
# BiocManager::install("phyloseq")
# BiocManager::install("decontam")

library(dplyr)
library(phyloseq)
library(decontam)
library(vegan)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(Matrix)

# id_mapping <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Original Study Mapping - Sheet3.csv", header = TRUE)
fungal.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/Walther-Antonio_Project_022_ITS2.rds")
uminn_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_uminn_data.csv", header=TRUE)
samples.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_samples.csv")

# uminn_data 
uminn_data <- uminn_data %>% 
  select(Sample.ID, Special.Notes, qr)

# fungal.data

## Accessors
ntaxa(fungal.data)
nsamples(fungal.data)
sample_names(fungal.data)[1:5]
rank_names(fungal.data)
sample_variables(fungal.data, errorIfNULL=FALSE) # no colnames
otu_table(fungal.data)[1:5, 1:5]
tax_table(fungal.data)[1:5, 1:4]
phy_tree(fungal.data, errorIfNULL=FALSE) # no phylo tree
taxa_names(fungal.data)[1:10]

## Taxa table
taxa_fungal.data <- as.data.frame(tax_table(fungal.data))
dim(taxa_fungal.data)
summary(as.numeric(taxa_fungal.data$confidence))
table(as.numeric(taxa_fungal.data$confidence))

## OTU table
# otu_table: Operational taxonomic unit (OTU) or abundance table.
otu_table_fungal.data <- as.data.frame(otu_table(fungal.data))

length(rownames(otu_table(fungal.data)))
head(colnames(otu_table(fungal.data)))

### Data Preprocessing

otu_table_fungal.data <- as.data.frame(t(otu_table_fungal.data))

metadata_fungal <- otu_table_fungal.data %>% 
  mutate(SampleID=rownames(otu_table_fungal.data),
         is_blank=as.logical(ifelse(str_detect(SampleID, "BLANK"), "TRUE", "FALSE")),
         SampleID= sub("\\..*", "", SampleID)) %>%
  select(SampleID, is_blank)

# Convert back to phyloseq obj
otu_table_obj <- otu_table(otu_table_fungal.data, taxa_are_rows = FALSE)
sample_data_obj <- sample_data(metadata_fungal)
fungal_physeq <- phyloseq(otu_table_obj, sample_data_obj, tax_table(fungal.data))
# pck Decontam
# Identify contaminants based on prevalence
contam_prev <- isContaminant(fungal_physeq, method = "prevalence", neg = "is_blank")
contaminants <- contam_prev$contaminant

# Filter the OTUs in the phyloseq object
fungal_physeq_no_contam <- prune_taxa(!contaminants, fungal_physeq)
fungal_physeq_subset <- subset_taxa(fungal_physeq_no_contam, Kingdom == "Fungi" & !is.na(Phylum) & Phylum != "")

# Filter ASVs - for now didn't
fungal_otu_table <- otu_table(fungal_physeq_subset)
asv_presence <- colSums(fungal_otu_table > 0)
fungal_taxa_table <- tax_table(fungal_physeq_subset)

# get df
# fungal_otu_table_df <- as.data.frame(otu_table(fungal_physeq_subset))
# fungal_taxa_table_df <- as.data.frame(tax_table(fungal_physeq_subset))
fungal_metadata_df <- as.data.frame(as.matrix(sample_data(fungal_physeq_subset)))

## Map uminn_data
# samples.data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Report 8-Sample.csv")
samples.data$SampleID <- sub("_.*", "", samples.data$qr)
# samples.data <- samples.data %>% 
#   rename(SampleID=qr)
# samples.data <- samples.data %>% 
#   rename(biome_id=uid)
# full_sample_data <- uminn_data %>% 
#   left_join(samples.data, by="qr") %>% 
#   rename(SampleID=Sample.ID) %>% 
#   select(biome_id, SampleID, qr, logDate, timestamp, sampleType, Special.Notes)

# Filter samples
# fungal_metadata_df_IDs <- fungal_metadata_df %>% 
#   full_join(full_sample_data, by="SampleID") %>% 
#   filter(!is.na(biome_id))
fungal_metadata_df_IDs <- fungal_metadata_df %>% 
  left_join(samples.data, by="SampleID") %>% 
  filter(!is.na(biome_id))

# Update OTU and taxa table
rownames(fungal_otu_table) <- sub("\\..*", "", rownames(fungal_otu_table))
# rownames(fungal_taxa_table) <- sub("\\..*", "", rownames(fungal_taxa_table))
fungal_metadata_df_IDs_rownames <- fungal_metadata_df_IDs$SampleID
rownames(fungal_metadata_df_IDs) <- fungal_metadata_df_IDs$SampleID
otu_table_filtered <- fungal_otu_table[fungal_metadata_df_IDs_rownames, , drop = FALSE]
fungal_taxa_table_filtered <- fungal_taxa_table[colnames(otu_table_filtered), , drop = FALSE]

# Filtered phyloseq obj
fungal_physeq_filtered <- phyloseq(otu_table(otu_table_filtered), sample_data(fungal_metadata_df_IDs), tax_table(fungal_taxa_table_filtered))

### Get OTU table and Taxa table
otu_table_fungal.data <- (otu_table(fungal_physeq_filtered))
taxa_table_fungal.data <- ((tax_table(fungal_physeq_filtered)))
sample_data_phyloseq <- as.data.frame(as.matrix(sample_data(fungal_physeq_filtered)))

### Basic Analysis

ntaxa(fungal_physeq_filtered)
nsamples(fungal_physeq_filtered)

# OTU Confidence
# summary(as.numeric(taxa_table_fungal.data$confidence))

# filter for just fecal
sample_data_phyloseq_fec <- sample_data_phyloseq %>% 
  filter(sampleType=="fecal")
fungal_otu_fec <- otu_table_fungal.data[rownames(sample_data_phyloseq_fec), , drop = FALSE]
fungal_taxa_fec <- taxa_table_fungal.data[colnames(fungal_otu_fec), , drop = FALSE]

# filter for just vaginal
sample_data_phyloseq_vag <- sample_data_phyloseq %>% 
  filter(sampleType=="vaginal")
fungal_otu_vag <- otu_table_fungal.data[rownames(sample_data_phyloseq_vag), , drop = FALSE]
fungal_taxa_vag <- taxa_table_fungal.data[colnames(fungal_otu_vag), , drop = FALSE]

# vaginal phyloseq obj
fungal_physeq_vag <- phyloseq(otu_table(fungal_otu_vag, taxa_are_rows=FALSE), 
                              sample_data(sample_data_phyloseq_vag), tax_table(fungal_taxa_vag))

tax_table_vag <- as.data.frame(tax_table(fungal_physeq_vag))
# dominant Species
dominant_spec <- apply(otu_table(fungal_physeq_vag), 1, function(x) {
  spec <- tax_table_vag$Species[which.max(x)]
  ifelse(is.na(spec), "Unknown", spec)
})
sample_data(fungal_physeq_vag)$DominantSpecies <- dominant_spec

# Alpha Diversity
alpha_div <- estimate_richness(fungal_physeq_vag, measures = c("Shannon"))
plot_richness(fungal_physeq_vag, x = "sampleType", measures = c("Shannon")) + theme_minimal()

plot_richness(fungal_physeq_vag, x = "sampleType", color = "DominantSpecies", measures = c("Shannon")) +
  theme_minimal() +
  ggtitle("Alpha Diversity Colored by Dominant Species")

# Beta Diversity
beta_div <- ordinate(fungal_physeq_vag, method = "NMDS", distance = "bray")

plot_ordination(fungal_physeq_vag, beta_div, color = "DominantSpecies") +
  theme_minimal() +
  ggtitle("Beta Diversity Colored by Dominant Species")


# Normalize using relative abundance
fungal_physeq_filtered_rel <- transform_sample_counts(fungal_physeq_vag, function(x) x / sum(x))

######################################################

otu_table_data <- as.data.frame(otu_table(fungal_physeq_filtered))
taxa_data <- as.data.frame(tax_table(fungal_physeq_filtered))

# Calculate total abundance per species
total_abundance <- colSums(otu_table_data, na.rm = TRUE)
species_abundance <- data.frame(Species = taxa_data$Species, Total_Abundance = total_abundance)
# works but there might be some separation I need to do here
species_abundance <- species_abundance %>% 
  group_by(Species) %>% 
  summarise(Total_Abundance=sum(Total_Abundance, na.rm=TRUE)) %>% 
  arrange(desc(Total_Abundance)) %>% 
  ungroup() %>% 
  filter(Total_Abundance != 0) %>% 
  filter(Species != "")  %>% 
  mutate(Relative_Abundance=
           round((species_abundance$Total_Abundance / sum(species_abundance$Total_Abundance)) * 100, 2)) %>% 
  mutate(Rank=row_number()) 

summary(species_abundance$Total_Abundance)

## Boxplot
q1 <- quantile(species_abundance$Total_Abundance, 0.25, na.rm = TRUE)
q3 <- quantile(species_abundance$Total_Abundance, 0.75, na.rm = TRUE)
iqr <- q3 - q1
lower_bound <- 0
upper_bound <- q3 + 1.5 * iqr

outliers <- species_abundance$Total_Abundance[species_abundance$Total_Abundance < lower_bound |
                                                species_abundance$Total_Abundance > upper_bound]
ggplot(species_abundance, aes(y = Total_Abundance)) +
  geom_boxplot(fill = "orange", color = "black") +
  # Set axis limits
  coord_cartesian(ylim = c(lower_bound, upper_bound)) +  
  labs(title = " ",
       y = "Total Abundance") +
  theme_minimal()


# Relative Abundance (percentage)
# species_abundance <- species_abundance %>% 
#   mutate(Relative_Abundance=
#            round((species_abundance$Total_Abundance / sum(species_abundance$Total_Abundance)) * 100, 3))

summary(species_abundance$Relative_Abundance)

## Boxplot
ggplot(species_abundance, aes(y = Relative_Abundance)) +
  geom_boxplot(fill = "yellow", color = "black") +
  # Set axis limits
  coord_cartesian(ylim = c(0, 80)) +  
  labs(title = " ",
       y = "Relative Abundance (%)") +
  theme_minimal()

## Rank Abundance
# ggplot(species_abundance, aes(x = Rank, y = Total_Abundance)) +
#   geom_line() +
#   geom_point() +
#   geom_text(
#     data = species_abundance[c(1:3),], 
#     aes(label = Species), 
#     nudge_y = 0.01*max(species_abundance$Total_Abundance),
#     nudge_x = 8,
#     size = 2, hjust = 0
#   ) +
#   theme_minimal() +
#   labs(title = "Rank Abundance Curve", x = "Rank", y = "Total Abundance")

ggplot(species_abundance[c(1:20),], aes(x = Rank, y = Total_Abundance)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  geom_text(
    data = species_abundance[c(1:3),],
    aes(label = Species),
    nudge_y = 0.02*max(species_abundance$Total_Abundance),
    nudge_x = 0.05,
    size = 5, hjust = 0
  ) +
  labs(title = " ", x = "Rank", y = "Total Abundance")
  # labs(title = "Rank Abundance Curve (top 20)", x = "Rank", y = "Total Abundance")

## Top abundance - top 20 OTUs
top_otu_abundance <- rowSums(otu_table_data, na.rm = TRUE) %>%
  enframe(name = "OTU", value = "Total_Abundance") %>%
  arrange(desc(Total_Abundance)) %>%
  slice_head(n = 20)
# top_species_abundance <- species_abundance %>% 
#   slice_head(n = 20)
top_species_OTU <- otu_table_data[top_otu_abundance$OTU, ]
heatmap_data <- as.data.frame(as.matrix(top_species_abundance), stringsAsFactors = FALSE) %>% 
  mutate(Total_Abundance=as.numeric(Total_Abundance),
         Relative_Abundance=as.numeric(Relative_Abundance))

## Alpha Div - Shannon Index
alpha_diversity <- estimate_richness(fungal_physeq_filtered, measures = "Shannon")
alpha_diversity_df <- as.data.frame(alpha_diversity)
summary(alpha_diversity_df)

# Add participant IDs from sample data
sample_data_df <- as.data.frame(sample_data(fungal_physeq_filtered))
alpha_diversity_df$biome_id <- sample_data_df$biome_id

alpha_diversity_summary <- alpha_diversity_df %>%
  mutate(biome_id=as.numeric(biome_id)) %>% 
  filter(!is.na(biome_id)) %>% 
  # mutate(biome_id=as.factor(biome_id)) %>% 
  group_by(biome_id) %>%
  summarise(Mean_Shannon = mean(Shannon, na.rm = TRUE),
            SD_Shannon = sd(Shannon, na.rm = TRUE),
            Min_Shannon = min(Shannon, na.rm = TRUE),
            Max_Shannon = max(Shannon, na.rm = TRUE)) %>% 
  ungroup() 

alpha_diversity_summary <- alpha_diversity_summary %>%
  mutate(biome_id = reorder(biome_id, -Mean_Shannon))

summary(round(alpha_diversity_summary$Mean_Shannon, 2))

# Plot Shannon Index
ggplot(alpha_diversity_summary, aes(x = as.factor(biome_id), y = Mean_Shannon)) +
  geom_bar(stat = "identity", fill = "coral") +
  theme_minimal() +
  labs(title = " ", x = " ", y = "Shannon Index") +
  # labs(title = "Shannon Diversity Index by Participants", x = "Participant ID", y = "Shannon Index") +
  # geom_text(aes(label = round(Mean_Shannon, 2)), vjust = -0.5) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) 

## Beta Diversity - Bray-Curtis dissimilarity
bray_curtis <- vegdist(otu_table_data, method = "bray")
bray_curtis_df <- as.data.frame(as.matrix(bray_curtis))
biome_id <- sample_data_df$biome_id
participant_id_df <- data.frame(biome_id=biome_id)

# Plot Bray-Curtis dissimilarity
pheatmap(bray_curtis_df,
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# Random subset
subset_otu_table_data <- otu_table_data[1:50, ]
subset_sample_data <- sample_data_df[1:50, ]
bray_curtis_subset <- vegdist(subset_otu_table_data, method = "bray")
bray_curtis_subset_df <- as.data.frame(as.matrix(bray_curtis_subset))
biome_id_subset <- subset_sample_data$biome_id
participant_id_subset_df <- data.frame(biome_id=biome_id_subset)
pheatmap(bray_curtis_subset_df,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = participant_id_subset_df,
         annotation_colors = list(biome_id = RColorBrewer::brewer.pal(3, "Set1")))

# NMDS ordination
fungal_physeq_subset <- subset_samples(fungal_physeq_filtered, sampleType %in% c("vaginal"))
beta_div <- ordinate(fungal_physeq_subset, method = "NMDS", distance = "bray")

plot_ordination(fungal_physeq_subset, beta_div, color = "biome_id") + theme_minimal() +
  labs(title = "Beta Diversity (NMDS Plot) for Vaginal Samples by Participant ID")
