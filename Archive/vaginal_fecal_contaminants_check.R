library(vegan)
library(pheatmap)
library(tidyverse)
library(Matrix)

library(cluster)
library("igraph")
library("markovchain")
library("RColorBrewer")
library("gridExtra")
library(viridis)

source("~/Microbiome Thesis/functions.R")


bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/vaginal_cleaned_max_taxa.rds")

###############################################################################################

bacteria_taxa_df <- tax_table(bacterial.data)
bacteria_metadata_df <- sample_data(bacterial.data)

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


###################################################################################################

###### Five Community State Type (CSTs) Clustering

# species to CST mapping
species_to_cst <- data.frame(
  Species = c("crispatus", "gasseri", "iners", "jensenii"),
  CST = c("I", "II", "III", "V") # IV is anaerobic/diverse cluster
)

# Get most abundant OTU per sample
max_taxa <- apply(relative_abundance_otu_t, 2, function(sample) {
  taxa_idx <- which.max(sample)
  taxa_names(vaginal_relative_abundances)[taxa_idx]
})

# Map most abundant OTU to the sample data
bacteria_metadata_df$max_taxa <- max_taxa

# Get max taxa names
# bacteria_metadata_df$OTU <- as.character(bacteria_taxa_df[bacteria_metadata_df$max_taxa, "Species_exact"]) # 2039 samples
bacteria_metadata_df$OTU <- as.character(bacteria_taxa_df[bacteria_metadata_df$max_taxa, "BLAST_species"])

# Set sample data in phyloseq obj
sample_data(bacterial.data) <- bacteria_metadata_df
# bacteria_taxa_df <- as.data.frame(bacteria_taxa_df)

# Assign CST col
bacteria_metadata_df <- as(bacteria_metadata_df, "data.frame")
class(bacteria_metadata_df)

bacteria_metadata_df <- bacteria_metadata_df %>%
  mutate(
    CST = case_when(
      str_detect(OTU, "gasseri") ~ "II",
      str_detect(OTU, "crispatus") ~ "I",
      str_detect(OTU, "iners") ~ "III",
      str_detect(OTU, "jensenii") ~ "V",
      TRUE ~ "IV" # Assign "IV" for diverse/anaerobic or unclassified species
    )
  )

table(bacteria_metadata_df$CST)
cst_summary <- bacteria_metadata_df %>% 
  count(CST) %>% 
  mutate(Percentage=100*(n/sum(n)))
print(cst_summary)

sample.cst <- bacteria_metadata_df$CST %>% 
  as.data.frame()
colnames(sample.cst) <- "CST"
rownames(sample.cst) <- rownames(bacteria_metadata_df)
sample.cst$SampleID <- rownames(sample.cst)



################################################################################

head(bacteria_metadata_df)
vaginal_relative_abundances <- transform_sample_counts(bacterial.data, function(x) x/sum(x))

## taxa frequency plot

# top 10 max taxa?

# # aggregate at Species level - with BLAST
# taxFreq <- tax_glom(vaginal_relative_abundances, "BLAST_species")
# taxFreq_df <- psmelt(taxFreq)
# taxFreq_df$biome_id <- sample_data(bacterial.data)$biome_id[match(taxFreq_df$Sample, sample_names(bacterial.data))]
# 
# # top 10
# taxFreq_df_top <- taxFreq_df %>%
#   group_by(biome_id) %>%
#   arrange(desc(Abundance)) %>%
#   mutate(rank = row_number()) %>%
#   mutate(BLAST_species = ifelse(rank > 5, "Other", BLAST_species)) %>%
#   ungroup() %>% 
#   filter(SampleID=="S1763_V3V5_S1176")
# 
# taxFreq_df_top <- taxFreq_df_top %>%
#   group_by(Sample, BLAST_species) %>%
#   summarise(Abundance = sum(Abundance)) %>%
#   ungroup()

# S1763_V3V5_S1176

# p <- ggplot(taxFreq_df_top, aes(x = Sample, y = Abundance, fill = BLAST_species)) +
#   geom_bar(stat = "identity") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(title = "Taxa Frequency Plot", y = "Relative Abundance", x = "Sample") #+
# # theme(legend.position = "none")
# p + theme(
#   legend.text = element_text(size = 8),
#   legend.key.size = unit(0.4, "cm"),
#   legend.title = element_text(size = 9)  # Optional: adjust title size
# )

# taxFreq_df_participant <- taxFreq_df %>%
#   filter(biome_id==1)
# group_by(biome_id, BLAST_species) %>%
# summarise(Abundance = mean(Abundance, na.rm = TRUE))

# taxFreq_df_participant_truncated <- taxFreq_df_participant %>%
#   mutate(Blast_species = ifelse(nchar(BLAST_species) > 10, 
#                                 paste0(substr(BLAST_species, 1, 10), "…"), 
#                                 BLAST_species),
#          Sample = ifelse(nchar(Sample) > 10, 
#                          paste0(substr(Sample, 1, 10), "…"), 
#                          Sample),
#   )


# ## reads per samples
# reads_per_sample_df <- data.frame(Sample = sample_names(bacterial.data),
#                                   Reads = sample_sums(bacterial.data))
# 
# ggplot(reads_per_sample_df, aes(x = Sample, y = Reads)) +
#   geom_bar(stat = "identity", fill = "steelblue") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(title = "Number of Reads per Sample", y = "Read Count", x = "Sample") +
#   theme(legend.position = "none")
# 
# ## Fecal contamination
# fecal_taxa <- c("bacteroides", "enterococcus", "escherichia", "clostridium")
