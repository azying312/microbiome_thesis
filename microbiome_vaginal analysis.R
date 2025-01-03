library(dplyr)
library(phyloseq)
library(vegan)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(Matrix)

# Use corrected data after
bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/vaginal_bacteria_cleaned.rds")

##########################################################################################

#### Exploratory Data Analysis

# Relative abundances
vaginal_relative_abundances <- transform_sample_counts(bacterial.data, function(x) x/sum(x))
# vaginal_relative_abundances <- transform_sample_counts(bacteria_subset_vaginal, function(x) x/sum(x))
relative_abundance_otu <- as.data.frame(otu_table(vaginal_relative_abundances))

relative_abundance_otu_t <- t(relative_abundance_otu) %>% as.data.frame()
relative_abundance_otu_t$SampleID <- rownames(relative_abundance_otu_t)

bacteria_metadata_df <- sample_data(bacterial.data)
bacteria_metadata_df <- bacteria_metadata_df[,-2]
otu_with_participant <- relative_abundance_otu_t %>%
  left_join(bacteria_metadata_df, by="SampleID")

# bray_curtis_dist <- vegdist(t(relative_abundance_otu), method = "bray")
specnumber.ra <- specnumber(relative_abundance_otu_t)

###### Five Community State Type (CSTs) Clustering

# species to CST mapping
species_to_cst <- data.frame(
  Species = c("Lactobacillus crispatus", "Lactobacillus gasseri", "Lactobacillus iners", "Lactobacillus jensenii"),
  CST = c("I", "II", "III", "V") # IV is anaerobic/diverse cluster
)

# Get most abundant OTU per sample
max_taxa <- apply(relative_abundance_otu, 2, function(sample) {
  taxa_idx <- which.max(sample)
  taxa_names(vaginal_relative_abundances)[taxa_idx]
})

# Map most abundant OTU to the sample data
bacteria_metadata_df$max_taxa <- max_taxa
bacteria_taxa_table <- tax_table(bacterial.data)

# Map max_taxa to name
bacteria_metadata_df$OTU <- bacteria_taxa_table[bacteria_metadata_df$max_taxa, "Species"]
# Set sample data in phyloseq obj
sample_data(bacterial.data) <- bacteria_metadata_df

# Add CST col - check what this doing
# bacteria_taxa_table$CST <- ifelse(
bacteria_taxa_table <- as.data.frame(bacteria_taxa_table)

# Assign CST col
bacteria_metadata_df
bacteria_metadata_df <- as.data.frame(bacteria_metadata_df)

bacteria_metadata_df$CST <- ifelse(
  bacteria_metadata_df$OTU %in% species_to_cst$Species,
  species_to_cst$CST[match(bacteria_metadata_df$OTU, species_to_cst$Species)],
  # Assign to CST IV if not in I, II, III, V
  "IV"  
)

table(bacteria_metadata_df$CST)
table(bacteria_metadata_df$OTU)

########################################################

# Aggregate data by participants by mean relative abundance for a given OTU
participant_otu <- tapply(sample_names(bacteria_metadata_df), 
                          sample_data(bacteria_metadata_df)$biome_id, 
                          function(samples) rowMeans(otu_table(vaginal_relative_abundances)[, samples, drop = FALSE]))
participant_otu <- do.call(cbind, participant_otu)
rownames(participant_otu) <- taxa_names(vaginal_relative_abundances)

## Alpha Div - Shannon Index
shannon.24 <- diversity(t(otu_table(bacterial.data)), index="shannon")

# Add participant IDs from sample data | Merge the calculated Shannon diversity values with metadata
bacteria_metadata_df <- as(bacteria_metadata_df, "data.frame")
shannon.df.24 <- data.frame("SampleID"=names(shannon.24), "shannon"=shannon.24)
shannon.qr.merged.24 <- merge(shannon.df.24, bacteria_metadata_df[,1:4], by="SampleID") %>% 
  mutate(biome_id=as.integer(biome_id))
# sample_data_df <- as.data.frame(sample_data(fungal_physeq_filtered))
# alpha_diversity_df$biome_id <- sample_data_df$biome_id

# Define unique user IDs and labels
unique.uid.24 <- sort(unique(shannon.qr.merged.24$biome_id))
# alph <- LETTERS[1:length(unique.uid.24)] 

par(mfrow = c(2, 2), mai = c(.6, .6, .3, .3))

# Loop through each user and plot diversity trends
i <- 1
for (i.uid in unique.uid.24[1:4]) {
  # Extract data for the current user
  current.div <- shannon.qr.merged.24$shannon[shannon.qr.merged.24$biome_id == i.uid]
  current.dates <- shannon.qr.merged.24$logDate[shannon.qr.merged.24$biome_id == i.uid]
  
  # Remove NA values from both dates and diversity scores
  valid_indices <- !is.na(current.dates) & !is.na(current.div)
  current.dates <- current.dates[valid_indices]
  current.div <- current.div[valid_indices]
  
  # Plot diversity trend
  plot(current.dates, current.div,
       xlim = range(as.Date(shannon.qr.merged.24$logDate)),
       ylim = c(0, max(shannon.qr.merged.24$shannon, na.rm = TRUE)),
       xlab = "Date",
       ylab = "Shannon Index",
       main = paste("User", i.uid),
       pch = 16, col = "blue")
  
  i <- i + 1
}

########################################################
## Corr with diet
diet.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/fully_merged_diet_data.csv", header = TRUE)

diet.data.collapsed <- diet.data %>% 
  group_by(study_id, Date) %>% 
  summarise(total_cal=sum(caloriesall)) %>% 
  rename(biome_id=study_id,
         logDate=Date)

ggplot(diet.data.collapsed, aes(x=total_cal, color=factor(biome_id))) +
  geom_histogram(binwidth=100)

# remove outliers with IQR
diet.data.collapsed.filtered <- diet.data.collapsed %>% 
  mutate(Q1=quantile(total_cal, 0.25),
         Q3=quantile(total_cal, 0.75),
         IQR=Q3-Q1) %>%
  filter(total_cal >= (Q1 - 1.5 * IQR) & total_cal <= (Q3 + 1.5 * IQR)) %>%
  select(-Q1, -Q3, -IQR)

ggplot(diet.data.collapsed.filtered, aes(x=total_cal, color=factor(biome_id))) +
  geom_histogram(binwidth=100)

shannon.diet <- shannon.qr.merged.24 %>% 
  left_join(diet.data.collapsed.filtered, by=c("biome_id", "logDate")) %>%
  filter(!is.na(total_cal)) # 2039 -> 1196 = 843 no diet data | before IQR remove

ggplot(shannon.diet, aes(x = total_cal, y = shannon, color=factor(biome_id))) +
  geom_point() +
  geom_smooth(method = "lm", col = "blue") +
  labs(
    title = "",
    x = "Total Calories",
    y = "Shannon Diversity Index"
  ) +
  theme_minimal()

cor.test(shannon.diet$shannon, shannon.diet$total_cal, method = "spearman")

# group_by biome_id
participant.spearman <- shannon.diet %>% 
  group_by(biome_id) %>% 
  summarise(spearman=cor(shannon, total_cal, method="spearman"))

hist(x=participant.spearman$spearman)

