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
fungal.data <- readRDS("/Users/alicezhang/Desktop/microbiome_data/sequeced_data/Walther-Antonio_Project_022_ITS2.rds")
uminn_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_uminn_data.csv", header=TRUE)
samples.data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_samples.csv")

# bacterial.data <- readRDS("/Users/alicezhang/Desktop/microbiome_data/sequeced_data/Walther-Antonio_Project_022_16S.rds")
# bacterial_otu_table <- as.data.frame(otu_table(bacterial.data))
# View(bacterial_otu_table)

# uminn_data 
uminn_data <- uminn_data %>% 
  select(Sample.ID, Special.Notes, qr)

fungal.data

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
samples.data$qr <- sub("_.*", "", samples.data$qr)
samples.data <- samples.data %>% 
  rename(SampleID=qr)
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
# sample_sums(fungal_physeq_filtered)
# taxa_sums(fungal_physeq_filtered)   
ntaxa(fungal_physeq_filtered)
nsamples(fungal_physeq_filtered)

# OTU Confidence
# summary(as.numeric(taxa_table_fungal.data$confidence))

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


# physeq_phylum <- tax_glom(fungal_physeq_vag, taxrank = "Phylum")
# plot_bar(physeq_phylum, fill = "Phylum") + theme_minimal()

# Total OTU abundance per sample
# barplot(sample_sums(fungal_physeq_filtered))
# Total abundance per OTU across all samples
# taxa_sums(fungal_physeq_filtered)


# # Standardize sample names in OTU table
# rownames(otu_table_fungal.data) <- gsub(" ", "_", rownames(otu_table_fungal.data))
# rownames(otu_table_fungal.data) <- toupper(rownames(otu_table_fungal.data))
# rownames(otu_table_fungal.data) <- paste("sample_", rownames(otu_table_fungal.data), sep = "")
# # Add OTU col
# otu_table_fungal.data2 <- otu_table_fungal.data
# taxa_fungal.data2 <- taxa_fungal.data
# otu_table_fungal.data2$OTU <- rownames(otu_table_fungal.data2)
# taxa_fungal.data2$OTU <- rownames(taxa_fungal.data2)
# # Merge OTU and Taxa
# otu_taxa <- otu_table_fungal.data2 %>% 
#   left_join(taxa_fungal.data2, by = 'OTU') %>% 
#   # filtering for actual calssifications
#   filter(Kingdom != "") %>%
#   filter(confidence > 0) %>% 
#   select("OTU", everything())
# # Filter OTU table
# filtered_OTU_table <- otu_table_fungal.data[rownames(otu_table_fungal.data) %in% otu_taxa$OTU,]
# # Filter taxa table
# filtered_taxa_table <- taxa_fungal.data %>% 
#   filter(confidence > 0) %>% 
#   filter(Kingdom != "")


# Total Abundance
species_abundance <- otu_taxa %>% 
  group_by(Species) %>%
  summarise(total_abundance=sum(across(starts_with("sample_"), ~sum(.)))) %>% 
  filter(total_abundance != 0) %>% 
  filter(Species != "") %>% 
  arrange(desc(total_abundance))
dim(species_abundance)
# dim(otu_taxa)

# TOTAL ABUNDANCE
summary(species_abundance$total_abundance)

# Relative Abundance
relative_abundance <- species_abundance %>% 
  mutate(relative_abundance=((total_abundance/sum(total_abundance))*100))
dim(relative_abundance)

# RELATIVE ABUNDANCE
summary(relative_abundance$relative_abundance)

### Visualization

## Rank Abundance Curve
otu_rank_abundance <- filtered_OTU_table %>%
  rowSums(na.rm = TRUE) %>%
  enframe(name = "OTU", value = "Total_Abundance") %>%
  arrange(desc(Total_Abundance)) %>%
  mutate(Rank = row_number())
  # slice_head(n=20)

otu_rank_abundance_20 <- otu_rank_abundance %>% 
  slice_head(n=40)

ggplot(otu_rank_abundance_20, aes(x = Rank, y = Total_Abundance)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Rank Abundance Curve", x = "Rank", y = "Total Abundance")

# Get dominant Species
otu_rank_abundance <- otu_rank_abundance %>% 
  left_join(otu_taxa, by="OTU")
otu_rank_abundance <- otu_rank_abundance %>% 
  filter(Species != "") %>% 
  select(OTU, Total_Abundance, Rank, Species)

otu_rank_abundance_species <- otu_rank_abundance %>% 
  group_by(Species) %>% 
  summarise(Count=sum(Total_Abundance)) %>% 
  arrange(desc(Count))

head(otu_rank_abundance_species, 10)

# Subset samples
set.seed(360)
sample_subset <- sample(colnames(filtered_OTU_table), 60)

otu_sample_subset <- filtered_OTU_table[, sample_subset]

# Top 20 most abundant OTUs
top_otu_abundance <- otu_sample_subset %>%
  mutate(across(everything(), as.numeric)) %>% 
  rowSums(na.rm = TRUE) %>% 
  enframe(name = "OTU", value = "Total_Abundance") %>%
  arrange(desc(Total_Abundance)) %>% 
  slice_head(n = 20)

top_species_OTU <- otu_sample_subset %>%
  rownames_to_column(var = "OTU") %>%
  filter(OTU %in% top_otu_abundance$OTU)
dim(top_species_OTU)

heatmap_data <- top_species_OTU %>%
  pivot_longer(cols = starts_with("sample_"), names_to = "Sample", values_to = "Abundance")

ggplot(heatmap_data, aes(x = Sample, y = OTU, fill = Abundance)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(title = "Heatmap of Top 20 OTU Abundances for 60 Samples (random)", x = "Sample", y = "OTU")

## Diversity Metrics

# Shannon Index
shannon_diversity <- otu_sample_subset %>%
  summarise(across(starts_with("sample_"), ~ diversity(.x, index = "shannon"))) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Shannon_Index")
summary(shannon_diversity$Shannon_Index)

ggplot(shannon_diversity, aes(x = Sample, y = Shannon_Index)) +
  geom_bar(stat = "identity", fill = "coral") +
  theme_minimal() +
  labs(title = "Shannon Diversity Index Across 60 Samples (random)", x = "Sample", y = "Shannon Index") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Beta Diversity - Bray-Curtis dissimilarity matrix
bray_curtis <- vegdist(top_species_OTU[,-1], method = "bray", na.rm=TRUE)

bray_curtis_df <- as.data.frame(as.matrix(bray_curtis))
pheatmap(bray_curtis_df, cluster_rows = TRUE, cluster_cols = TRUE)

### Pseudo participants (cluster into 60)

top_otu_abundance <- filtered_OTU_table %>%
  mutate(across(everything(), as.numeric)) %>% 
  rowSums(na.rm = TRUE) %>% 
  enframe(name = "OTU", value = "Total_Abundance") %>%
  arrange(desc(Total_Abundance)) %>% 
  slice_head(n = 60)
otu_top_60 <- filtered_OTU_table[top_otu_abundance$OTU, ]
bray_curtis_dist <- vegdist(t(otu_top_60), method = "bray", na.rm=TRUE)

## k-means
set.seed(360)
kmeans_result <- kmeans(bray_curtis_dist, centers = 60, nstart = 25)
pseudo_participant_ids <- kmeans_result$cluster

# create mapping
sample_to_participant <- data.frame(
  sample_id = names(kmeans_result$cluster),
  participant_id = kmeans_result$cluster
)
otu_table_long <- otu_table_fungal.data %>%
  rownames_to_column(var = "OTU") %>%
  gather(key = "sample_id", value = "abundance", -OTU)
otu_table_with_participant <- otu_table_long %>%
  left_join(sample_to_participant, by = "sample_id")

otu_table_by_participant <- otu_table_with_participant %>%
  group_by(participant_id, OTU) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE))
 
# psuedo clustering
cluster_counts <- table(kmeans_result$cluster)
cluster_counts_df <- as.data.frame(cluster_counts)
colnames(cluster_counts_df) <- c("cluster_id", "sample_count")

ggplot(cluster_counts_df, aes(x = factor(cluster_id), y = sample_count)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  geom_text(aes(label = sample_count), vjust = -0.5) +
  labs(title = "Pseudo-Participant Groupings from K-means Clustering",
       x="Pseudo-Participant",
       y="Sample Count") +
  scale_x_discrete(breaks = seq(1, 60, by = 1)) +
  theme_minimal()

## OTU and Taxa
participant_id_vector <- 
  kmeans_result$cluster[match(colnames(filtered_OTU_table), 
                               names(kmeans_result$cluster))]

filtered_OTU_table["participant_id", ] <- participant_id_vector
filtered_OTU_table$OTU <- rownames(filtered_OTU_table)

filtered_taxa_table$OTU <- rownames(filtered_taxa_table)
otu_taxa <- filtered_OTU_table %>% 
  left_join(filtered_taxa_table, by = 'OTU') %>% 
  select("OTU", everything())
dim(otu_taxa)
# colnames(otu_taxa) <- tail(filtered_OTU_table, n=1)["participant_id",]
# otu_taxa <- otu_taxa[-9272,]

# participant_ids <- sub("\\..*", "", colnames(otu_taxa))

participant_ids <- t(otu_taxa[nrow(otu_taxa),])
participant_ids_name <- participant_ids[1, 1]
flat_participant_ids <- as.vector(participant_ids[-1, ])
participant_ids_df <- data.frame(item = flat_participant_ids)
colnames(participant_ids_df) <- c(participant_ids_name)

participant_ids <- as.data.frame(participant_ids["OTU"], participant_ids[["participant_id"]])
participant_ids$participant_id <- as.numeric(as.data.frame(participant_ids)$participant_id)
otu_taxa <- otu_taxa[-nrow(otu_taxa),]

colnames(participant_ids) <- "participant_id" 
otu_taxa2 <- cbind(otu_taxa, (participant_ids))

# Calculate total abundance by participant and species
total_abundance <- otu_taxa %>%
  gather(key = "sample_id", value = "abundance", -participant_id) %>%
  group_by(participant_id) %>%
  summarise(total_abundance = sum(abundance))

# relative abundance
relative_abundance <- otu_taxa_table %>%
  gather(key = "sample_id", value = "abundance", -participant_id) %>%
  group_by(participant_id) %>%
  mutate(relative_abundance = abundance / sum(abundance)) %>%
  summarise(total_relative_abundance = sum(relative_abundance))