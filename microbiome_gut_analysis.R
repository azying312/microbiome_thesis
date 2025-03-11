library(tidyverse)
library(vegan)
library(viridis)

source("~/Microbiome Thesis/functions.R")

# bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/fecal_bacteria_filteredv2.rds")
# RELABELED DATA
bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/relabeled_data/fecal_bacteria_cleanedv3.rds")

###############################################################################################

otu_table_df <- as(otu_table(bacterial.data), "matrix")
bacteria_taxa_df <- tax_table(bacterial.data)
bacteria_metadata_df <- sample_data(bacterial.data)

#### Exploratory Data Analysis

# Relative abundances
gut_relative_abundances <- transform_sample_counts(bacterial.data, function(x) x/sum(x))
relative_abundance_otu <- as.data.frame(otu_table(gut_relative_abundances)) # 1418 samples, 1015 OTU | 2015 OTU 1449 samples
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

# Get most abundant OTU per sample
# max_taxa <- apply(relative_abundance_otu_t, 2, function(sample) {
max_taxa <- apply(relative_abundance_otu_t, 1, function(sample) { # on relabeled data
  taxa_idx <- which.max(sample)
  taxa_names(gut_relative_abundances)[taxa_idx]
})

# Map most abundant OTU to the sample data
bacteria_metadata_df$max_taxa <- max_taxa

# Get max taxa names
bacteria_metadata_df$OTU <- as.character(bacteria_taxa_df[bacteria_metadata_df$max_taxa, "Species_exact"]) # 1418 samples

################################################################################

# Aggregate data by participants by mean relative abundance for a given OTU
# participant_otu <- tapply(sample_names(bacteria_metadata_df),
#                           sample_data(bacterial.data)$biome_id,
#                           function(samples) rowMeans(t(otu_table(gut_relative_abundances))[, samples, drop = FALSE]))
participant_otu <- tapply(sample_names(bacteria_metadata_df),
                          sample_data(bacterial.data)$biome_id,
                          function(samples) rowMeans(otu_table(gut_relative_abundances)[, samples, drop = FALSE]))

participant_otu <- do.call(cbind, participant_otu)
rownames(participant_otu) <- taxa_names(gut_relative_abundances)

## Alpha Div - Shannon Index
shannon.24 <- vegan::diversity(t(otu_table_df), "shannon")

# Add participant IDs from sample data | Merge the calculated Shannon diversity values with metadata
bacteria_metadata_df <- as(bacteria_metadata_df, "data.frame")
shannon.df.24 <- data.frame("SampleID"=names(shannon.24), "shannon"=shannon.24)
shannon.qr.merged.24 <- merge(shannon.df.24, bacteria_metadata_df, by="SampleID")
shannon.cst.qr.merged.24 <- shannon.qr.merged.24 %>% 
  mutate(biome_id=as.integer(biome_id)) %>% 
  filter(!is.na(biome_id)) %>% 
  # Filter the data to be within study days: 10-14 to 12-14
  filter(logDate > "2022-10-13" & logDate < "2022-12-15") # 2038 to 1971

dim(shannon.cst.qr.merged.24)

###############################################################################################

### Check UMinn Spreadsheet v. Sequenced Data
uminn_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_uminn_data.csv")
# samples.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_samples.csv")
samples.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_samplesv2.csv")

uminn_data <- uminn_data %>% 
  select(Sample.ID, Special.Notes) %>% 
  filter(Sample.ID != "BLANK") %>% 
  # errors in sample processing from UMinn
  filter(!str_detect(Special.Notes, "error")) %>% 
  filter(!str_detect(Special.Notes, "No swab in tube")) # 1418 -> 1369
uminn_data$qr <- sub("_.*", "", uminn_data$Sample.ID)

fecal_data <- samples.data %>% 
  filter(sampleType=="fecal")

uminn_data_fecal <- uminn_data %>% 
  filter(qr %in% fecal_data$qr)
uminn_data_fecal <- uminn_data_fecal %>% 
  left_join(fecal_data, by="qr")

# 345 - in sequenced data but not in the uminn returned swabs --> 340
length(setdiff(shannon.cst.qr.merged.24$qr, uminn_data_fecal$qr)) 
length(setdiff(uminn_data_fecal$qr, shannon.cst.qr.merged.24$qr)) # 19

diff.qr <- setdiff(shannon.cst.qr.merged.24$qr, uminn_data_fecal$qr)

setdiff(unique(shannon.cst.qr.merged.24$biome_id), unique(uminn_data_fecal$biome_id)) # 68
setdiff(unique(uminn_data_fecal$biome_id), unique(shannon.cst.qr.merged.24$biome_id))

write.csv(shannon.cst.qr.merged.24, "/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/relabeled_data/gut_shannon.cst.qr.merged.24.csv")

###############################################################################################

# Top 10 species per participant across samples for a participant
gut_relative_abundances_df <- psmelt(gut_relative_abundances)

gut_relative_abundances_df_subset <- gut_relative_abundances_df %>% 
  filter(biome_id %in% c(1,2,3))

top10_species <- gut_relative_abundances_df %>% 
  group_by(Species_exact) %>% 
  # group_by(biome_id, Species_exact) %>% 
  summarise(Total_Abundance = sum(Abundance), .groups = "drop") %>% 
  arrange(desc(Total_Abundance)) %>%  
  # arrange(biome_id, desc(Total_Abundance)) %>%  
  # group_by(biome_id) %>%
  slice_head(n = 10) %>% 
  pull(Species_exact)

top10_df <- gut_relative_abundances_df %>% 
  filter(Species_exact %in% top10_species) %>%
  group_by(biome_id, Species_exact) %>% 
  summarise(Total_Abundance = sum(Abundance), .groups = "drop")

top10_relative_df <- top10_df %>% 
  group_by(biome_id) %>%
  mutate(Relative_Abundance = Total_Abundance / sum(Total_Abundance) * 100)

ggplot(top10_relative_df, aes(x = factor(biome_id), y = Relative_Abundance, fill = Species_exact)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal() +
  scale_fill_viridis_d() +
  labs(title = "Top 10 Species per Participant", x = "Participant", y = "Relative Abundance") 
#+
  # theme(legend.position = "none")


###############################################################################################
# menses.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/imputed_menstruation_data_2_12.csv")
menses.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/relabeled_data/imputed_menstruation_data_3_11.csv")
menses.data <- menses.data %>% 
  rename_with(~gsub("X2022.", "2022.", .), starts_with("X2022.")) %>% 
  rename_with(~gsub("\\.", "-", .))

# Reshape
menses.data.long <- menses.data %>% 
  pivot_longer(cols=starts_with("2022-"), names_to="logDate", values_to="menses_status")

fecal.microbial.menses.24 <- shannon.cst.qr.merged.24 %>% 
  left_join(menses.data.long, by=c("biome_id", "logDate"))

fecal.microbial.menses.24 <- fecal.microbial.menses.24 %>% 
  mutate(menses_day = ifelse(menses_status %in% c(1,2,3,7,9,78), "menses", 
                             ifelse(menses_status %in% c(4,5,6,10), "not_menses", NA))) %>% 
  mutate(menses_day = ifelse(is.na(menses_day), "not_menses", "menses"))

write.csv(fecal.microbial.menses.24, file="/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/relabeled_data/gut.microbial.menses.24.csv")

menses.table.df <- fecal.microbial.menses.24 %>% 
  select(biome_id, logDate, shannon, qr, max_taxa, OTU, menses_status, menses_day)

##########################################################################################

# Plot Shannon diversity over time for each participant
participant_ids <- unique(fecal.microbial.menses.24$biome_id)
all_days <- seq.Date(as.Date("2022-10-13"), as.Date("2022-12-16"), by = "day")
all_days <- data.frame(logDate=as.character(all_days))
# file_path <- "/Volumes/T7/microbiome_data/graphics/Results from 2022 (compare to 2017-18)/gut_shannon_diversity_logDates/"
file_path <- "/Volumes/T7/microbiome_data/graphics/Results from 2022 (compare to 2017-18)/relabeled_data/gut_shannon_diversity_logDates/"
for(id in participant_ids) {
  # print(id)
  file_name_id <- paste0(file_path, id, "_id.png")
  fecal.microbial.menses.24.participant <- fecal.microbial.menses.24 %>% 
    filter(biome_id==id)
  fecal.microbial.menses.24.participant <- all_days %>%
    left_join(fecal.microbial.menses.24.participant, by = "logDate")
  fecal.microbial.menses.24.participant <- fecal.microbial.menses.24.participant %>% 
    mutate(menses_day=ifelse(!is.na(menses_day), menses_day, 
                             ifelse(!is.na(shannon), "not_recorded", NA)))
  dupe.day <- fecal.microbial.menses.24.participant %>% 
    group_by(logDate) %>% 
    filter(n() > 1) %>% 
    distinct(logDate)
  
  shannon_plt <- ggplot(fecal.microbial.menses.24.participant, 
                        aes(x=as.factor(logDate), y=shannon, color=as.factor(menses_day), 
                            shape=as.factor(menses_day))) +
    geom_point() +
    scale_color_manual(values = c("not_menses" = "black", "menses" = "red", "not_recorded" = "purple")) +  
    scale_shape_manual(values = c("not_menses" = 16, "menses" = 17, "not_recorded"=18)) +
    ylim(0,5)+
    # Add lines for the days with multiple samples
    geom_segment(data = dupe.day, 
                 aes(x = logDate, xend = logDate, y = -Inf, yend = Inf), 
                 linetype = "dashed", color = "gray50", inherit.aes = FALSE) +  
    theme_minimal() +
    labs(x = "Date", y = "Shannon Index", title=paste(id, "Shannon diversity days")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  ggsave(file_name_id, shannon_plt, width = 8, height = 6, dpi = 300)
  
  print(shannon_plt)
}

##########################################################################################

# vaginal.microbial.menses.24 <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/vaginal.microbial.menses.24.csv")
vaginal.microbial.menses.24 <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/relabeled_data/vaginal.microbial.menses.24.csv")

vaginal.microbial.menses.24 <- vaginal.microbial.menses.24 %>% 
  select(SampleID, shannon, qr, biome_id, logDate, sampleType, timestamp, max_taxa, OTU, survey_menstruate, menses_status, menses_day)
fecal.microbial.menses.24_subset <- fecal.microbial.menses.24 %>% 
  select(SampleID, shannon, qr, biome_id, logDate, sampleType, timestamp, max_taxa, OTU, survey_menstruate, menses_status, menses_day)

gut.vaginal.microbial <- rbind(fecal.microbial.menses.24_subset, vaginal.microbial.menses.24)
dim(gut.vaginal.microbial)

# Plot Shannon diversity over time for each participant
participant_ids <- unique(gut.vaginal.microbial$biome_id)
all_days <- seq.Date(as.Date("2022-10-13"), as.Date("2022-12-16"), by = "day")
all_days <- data.frame(logDate=as.character(all_days))
# file_path <- "/Volumes/T7/microbiome_data/graphics/Results from 2022 (compare to 2017-18)/vaginal_gut_shannon_diversity_logDates_10/"
file_path <- "/Volumes/T7/microbiome_data/graphics/Results from 2022 (compare to 2017-18)/relabeled_data/vaginal_gut_shannon_diversity_logDates_10/"

for(id in participant_ids) {
  # print(id)
  file_name_id <- paste0(file_path, id, "_id.png")
  gut.vaginal.microbial.participant <- gut.vaginal.microbial %>% 
    filter(biome_id==id)
  gut.vaginal.microbial.participant <- all_days %>%
    left_join(gut.vaginal.microbial.participant, by = "logDate")
  # gut.vaginal.microbial.participant <- gut.vaginal.microbial.participant %>% 
  #   mutate(menses_day=ifelse(!is.na(menses_day), menses_day, 
  #                            ifelse(!is.na(shannon), "not_recorded", NA)))
  dupe.day <- gut.vaginal.microbial.participant %>% 
    group_by(logDate) %>% 
    filter(n() > 1) %>% 
    distinct(logDate)
  
  shannon_plt <- ggplot(gut.vaginal.microbial.participant, 
                        aes(x=as.factor(logDate), y=shannon, color=as.factor(sampleType), 
                            shape=as.factor(sampleType))) +
    geom_point() +
    scale_color_manual(values = c("fecal" = "blue", "vaginal" = "red", "not_recorded" = "purple")) +  
    scale_shape_manual(values = c("fecal" = 16, "vaginal" = 17, "not_recorded"=18)) +
    ylim(0,10)+
    # Add lines for the days with multiple samples
    geom_segment(data = dupe.day, 
                 aes(x = logDate, xend = logDate, y = -Inf, yend = Inf), 
                 linetype = "dashed", color = "gray50", inherit.aes = FALSE) +  
    theme_minimal() +
    labs(x = "Date", y = "Shannon Index", title=paste(id, "Shannon diversity days")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  ggsave(file_name_id, shannon_plt, width = 8, height = 6, dpi = 300)
  
  print(shannon_plt)
}

####

fecal.microbial.menses.24.filtered <- fecal.microbial.menses.24 %>% 
  filter(!is.na(survey_menstruate))

ggplot(fecal.microbial.menses.24.filtered, aes(x = as.Date(logDate), y = shannon)) +
  geom_point(aes(color=as.factor(survey_menstruate))) +
  geom_smooth(method = "lm", se=FALSE, aes(color = as.factor(survey_menstruate))) +
  labs(
    x = "Days", 
    y = "Shannon Diversity Index",
    title = "",
    color="Menstruate v. No menstruate"
  ) +
  theme_minimal() +
  scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "1 day") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

