library(dplyr)
library(phyloseq)
library(vegan)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(Matrix)

# Use corrected data after
bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/vaginal_bacteria_cleaned.rds")
# bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/old_data/vaginal_bacteria_cleaned.rds")
###############################################################################################

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
# bray_curtis_dist <- vegdist(t(relative_abundance_otu), method = "bray")
# specnumber.ra <- specnumber(relative_abundance_otu_t)

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


###################################################################################################

### Get Species variable cleaned in new data

bacteria_taxa_table <- tax_table(bacterial.data)
bacteria_taxa_df <- as.data.frame(bacteria_taxa_table)

dim(bacteria_taxa_df) # 92499    10

length(bacteria_taxa_df$Species) # 92499
sum(bacteria_taxa_df$Species=="") # 40397
sum(bacteria_taxa_df$Species_exact=="") # 90484

bacteria_taxa_df <- bacteria_taxa_df %>% 
  mutate(Species_exact=ifelse(Species_exact=="", Species, Species_exact))
sum(bacteria_taxa_df$Species_exact=="") # 39845
sum(bacteria_taxa_df$Species_exact!="") # 52654

# Subset taxa; Filter ASVs w/o Species_exact assignment
bacterial.data_subset <- subset_taxa(
  bacterial.data,
  Species_exact != ""
)

bacteria_taxa_table <- tax_table(bacterial.data_subset)
bacteria_taxa_df <- as.data.frame(bacteria_taxa_table)
dim(bacteria_taxa_df) # 2015 OTU <-- 52654 OTUs  10 cols

bacteria_metadata_df$OTU <- as.character(bacteria_taxa_df[bacteria_metadata_df$max_taxa, "Species_exact"]) # 2039 samples

# Set sample data in phyloseq obj
sample_data(bacterial.data) <- bacteria_metadata_df

bacteria_taxa_table <- as.data.frame(bacteria_taxa_table)

# Assign CST col
# dim(bacteria_metadata_df)
bacteria_metadata_df <- as(bacteria_metadata_df, "data.frame")
class(bacteria_metadata_df)

# Expand meta data for CST assignment
match_species <- species_to_cst$Species

bacteria_metadata_df_long <- bacteria_metadata_df %>% 
  mutate(separate_flag = grepl(paste(match_species, collapse = "|"), OTU))  %>%
  # Separate rows only for flagged rows
  separate_rows(OTU, sep = "/") %>%
  filter(separate_flag | (!separate_flag & !grepl("/", OTU))) %>%
  mutate(CST_species=OTU %in% match_species) %>% 
  group_by(across(-OTU)) %>%
  # Separate CST I, II, III, V species
  reframe(
    OTU = c(
      OTU[CST_species],                     
      paste(OTU[!CST_species], collapse = "/") # IV Species
    )) %>% 
  filter(!(OTU=="")) %>% 
  select(-c(separate_flag, CST_species))

bacteria_metadata_df_long$CST <- ifelse(
  bacteria_metadata_df_long$OTU %in% species_to_cst$Species,
  species_to_cst$CST[match(bacteria_metadata_df_long$OTU, species_to_cst$Species)],
  # Assign to CST IV if not in I, II, III, V
  "IV"  
)

table(bacteria_metadata_df_long$CST)
table(bacteria_metadata_df_long$OTU)

## Check repeats
freq_table <- table(bacteria_metadata_df_long$SampleID)
freq_df <- as.data.frame(freq_table)
colnames(freq_df) <- c("SampleID", "Freq")
summary(freq_df$Freq)

freq_table <- table(bacteria_metadata_df_long$max_taxa)
freq_df <- as.data.frame(freq_table)
colnames(freq_df) <- c("OTU", "Freq")
summary(freq_df$Freq)

freq_table <- table(bacteria_metadata_df$max_taxa)
freq_df <- as.data.frame(freq_table)
colnames(freq_df) <- c("OTU", "Freq")
summary(freq_df$Freq)

########################################################

# Aggregate data by participants by mean relative abundance for a given OTU
participant_otu <- tapply(sample_names(bacteria_metadata_df), 
                          sample_data(bacteria_metadata_df)$biome_id, 
                          function(samples) rowMeans(t(otu_table(vaginal_relative_abundances))[, samples, drop = FALSE]))
participant_otu <- do.call(cbind, participant_otu)
rownames(participant_otu) <- taxa_names(vaginal_relative_abundances)

## Alpha Div - Shannon Index
shannon.24 <- diversity((otu_table(bacterial.data)), index="shannon")

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

