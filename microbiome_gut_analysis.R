########################
#
# Gut Microbiome Analysis with Lifestyle Factors
# Last updated: 04/03/2025
#
#########################

library(tidyverse)
library(vegan)
library(viridis)
library(phyloseq)

source("~/Microbiome Thesis/functions.R")

# bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/fecal_bacteria_filteredv2.rds")
# RELABELED DATA
bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/relabeled_data/fecal_cleaned_max_taxa.rds")

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
# bacteria_metadata_df$OTU <- as.character(bacteria_taxa_df[bacteria_metadata_df$max_taxa, "Species_exact"]) # 1418 samples
bacteria_metadata_df$OTU <- as.character(bacteria_taxa_df[bacteria_metadata_df$max_taxa, "BLAST_species"])

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
shannon.qr.merged.24 <- shannon.qr.merged.24 %>% 
  mutate(biome_id=as.integer(biome_id)) %>% 
  filter(!is.na(biome_id)) %>% 
  # Filter the data to be within study days: 10-14 to 12-14
  filter(logDate > "2022-10-13" & logDate < "2022-12-15") # 2038 to 1971

dim(shannon.qr.merged.24) # 1401   10

###############################################################################################

### Check UMinn Spreadsheet v. Sequenced Data
uminn_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_uminn_data.csv")
# samples.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_samples.csv")
samples.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/relabeled_data/cleaned_samplesv2.csv")

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

# 345 - in sequenced data but not in the uminn returned swabs --> 339
length(setdiff(shannon.qr.merged.24$qr, uminn_data_fecal$qr)) 
length(setdiff(uminn_data_fecal$qr, shannon.qr.merged.24$qr)) # 19

diff.qr <- setdiff(shannon.qr.merged.24$qr, uminn_data_fecal$qr)

setdiff(unique(shannon.qr.merged.24$biome_id), unique(uminn_data_fecal$biome_id)) # 68
setdiff(unique(uminn_data_fecal$biome_id), unique(shannon.qr.merged.24$biome_id))

# write.csv(shannon.qr.merged.24, "/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/relabeled_data/gut_shannon.qr.merged.24.csv")

###############################################################################################

# Top 10 species per participant across samples for a participant
gut_relative_abundances_df <- psmelt(gut_relative_abundances)

gut_relative_abundances_df_subset <- gut_relative_abundances_df %>% 
  filter(biome_id %in% c(1,2,3))

top10_species <- gut_relative_abundances_df %>% 
  group_by(BLAST_species) %>% 
  # group_by(biome_id, Species_exact) %>% 
  summarise(Total_Abundance = sum(Abundance), .groups = "drop") %>% 
  arrange(desc(Total_Abundance)) %>%  
  # arrange(biome_id, desc(Total_Abundance)) %>%  
  # group_by(biome_id) %>%
  slice_head(n = 10) %>% 
  pull(BLAST_species)

top10_df <- gut_relative_abundances_df %>% 
  filter(BLAST_species %in% top10_species) %>% # comment out for top 10 in each sample
  group_by(biome_id, BLAST_species) %>%
  summarise(Total_Abundance = sum(Abundance), .groups = "drop")

top10_relative_df <- top10_df %>% 
  group_by(biome_id) %>%
  mutate(Relative_Abundance = Total_Abundance / sum(Total_Abundance) * 100)

length(unique(top10_relative_df$biome_id))

# Gut: Top 10 Species Across Study in Gut Microbiome per Participant
ggplot(top10_relative_df, aes(x = factor(biome_id), y = Relative_Abundance, fill = BLAST_species)) +
  geom_bar(stat = "identity", position = "fill", show.legend=FALSE) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  theme_minimal() +
  scale_fill_viridis_d() +
  labs(title = " ", x = "Participant", y = "Relative Abundance",
       fill="Species") +
  theme(text=element_text(size=14))

unique(top10_relative_df$BLAST_species)

###############################################################################################
# menses.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/imputed_menstruation_data_2_12.csv")
menses.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/relabeled_data/imputed_menstruation_data_3_11.csv")
menses.data <- menses.data %>% 
  rename_with(~gsub("X2022.", "2022.", .), starts_with("X2022.")) %>% 
  rename_with(~gsub("\\.", "-", .))

# Reshape
menses.data.long <- menses.data %>% 
  pivot_longer(cols=starts_with("2022-"), names_to="logDate", values_to="menses_status")

fecal.microbial.menses.24 <- shannon.qr.merged.24 %>% 
  left_join(menses.data.long, by=c("biome_id", "logDate"))

fecal.microbial.menses.24 <- fecal.microbial.menses.24 %>% 
  mutate(menses_day = ifelse(menses_status %in% c(1,2,3,7,9,78), "menses", 
                             ifelse(menses_status %in% c(4,5,6,10), "not_menses", NA))) %>% 
  mutate(menses_day = ifelse(is.na(menses_day), "not_menses", menses_day))

table(fecal.microbial.menses.24$menses_day)

# 03/21 saved
# write.csv(fecal.microbial.menses.24, file="/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/relabeled_data/gut.microbial.menses.24.csv")

menses.table.df <- fecal.microbial.menses.24 %>% 
  select(biome_id, logDate, shannon, qr, max_taxa, OTU, menses_status, menses_day, regular_periods, sampleType)

# average shannon for each participant on menses vs. not on menses
menses.table.df.summary <- menses.table.df %>% 
  group_by(biome_id, menses_day) %>% 
  summarise(avg_shannon=mean(shannon), .groups="drop") %>% 
  pivot_wider(names_from = menses_day, values_from = avg_shannon, names_prefix = "menses_day_") %>% 
  filter(!is.na(menses_day_not_menses) & !is.na(menses_day_menses)) 
  # %>% select(!menses_day_NA)

### Wilcox test
wilcox.test(menses.table.df.summary$menses_day_not_menses,
            menses.table.df.summary$menses_day_menses, paired=TRUE)
t.test(menses.table.df.summary$menses_day_not_menses,
       menses.table.df.summary$menses_day_menses, paired=TRUE)

# Does Shannon diversity vary by menstruate/not menstruate?
# t.test(dass.participant.avg$avg_shannon ~ dass.participant.avg$person_menstruates)

length(unique(menses.table.df.summary$biome_id))
# Gut: Average Shannon Diversity for Participants on menses v. not
ggplot(
  menses.table.df.summary,
  aes(x = menses_day_not_menses, y = menses_day_menses)) +
  geom_point(color="orchid3", show.legend=FALSE) +
  xlim(3,4)+
  ylim(3,4)+
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  theme_minimal() +
  theme(text=element_text(size=18))+
  labs(x="Not Menses (Shannon)", y="Menses (Shannon)",
       color="Participant ID")

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
  fecal.microbial.menses.24.participant <- study_days(fecal.microbial.menses.24.participant)
  
  shannon_plt <- ggplot(fecal.microbial.menses.24.participant, 
                        aes(x=study_day, y=shannon, color=as.factor(menses_day), 
                            shape=as.factor(menses_day))) +
    geom_point() +
    scale_color_manual(name = "Menses Status",
                       values = c("not_menses" = "black", "menses" = "red", "not_recorded" = "purple"),
                       labels = c("not_menses" = "Not Menses", "menses" = "Menses", "not_recorded" = "Not Recorded")) +  
    scale_shape_manual(name = "Menses Status",
                       values = c("not_menses" = 16, "menses" = 17, "not_recorded"=18),
                       labels = c("not_menses" = "Not Menses", "menses" = "Menses", "not_recorded" = "Not Recorded")) +
    ylim(0,5)+
    # Add lines for the days with multiple samples
    # geom_segment(data = dupe.day, 
    #              aes(x = logDate, xend = logDate, y = -Inf, yend = Inf), 
    #              linetype = "dashed", color = "gray50", inherit.aes = FALSE) +  
    theme_minimal() +
    labs(x = "Study Day", y = "Shannon Diversity")+ #, title=paste(id, "Shannon diversity days")) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1),
          text=element_text(size=18)) +
    scale_x_continuous(breaks = seq(0, max(fecal.microbial.menses.24.participant$study_day), by = 5))
  
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
  
  dupe.day <- gut.vaginal.microbial.participant %>% 
    group_by(logDate) %>% 
    filter(n() > 1) %>% 
    distinct(logDate)
  dupe.day <- study_days(dupe.day)
  gut.vaginal.microbial.participant <- study_days(gut.vaginal.microbial.participant)
  
  shannon_plt <- ggplot(gut.vaginal.microbial.participant, 
                        aes(x=as.numeric(study_day), y=shannon, color=as.factor(sampleType), 
                            shape=as.factor(sampleType))) +
    geom_point() +
    scale_color_manual(name = "Microbiome",
                       values = c("fecal" = "blue", "vaginal" = "red", "not_recorded" = "purple"),
                       labels = c("fecal" = "Gut", "vaginal" = "Vaginal", "not_recorded" = "Not Recorded")) +  
    scale_shape_manual(name = "Microbiome",
                       values = c("fecal" = 16, "vaginal" = 17, "not_recorded"=18),
                       labels = c("fecal" = "Gut", "vaginal" = "Vaginal", "not_recorded" = "Not Recorded")) +
    ylim(0,10)+
    # Add lines for the days with multiple samples
    # geom_segment(data = dupe.day, 
    #              aes(x = study_day, xend = study_day, y = -Inf, yend = Inf), 
    #              linetype = "dashed", color = "gray50", inherit.aes = FALSE) +  
    # geom_segment(data = dupe.day,
    #              aes(x = logDate, xend = logDate, y = -Inf, yend = Inf),
    #              linetype = "dashed", color = "gray50", inherit.aes = FALSE) +
    theme_minimal() +
    labs(x = "Study Day", y = "Shannon Diversity")+ #, title=paste(id, "Shannon diversity days")) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1),
          text=element_text(size=16)) +
    scale_x_continuous(breaks = seq(0, max(gut.vaginal.microbial.participant$study_day), by = 5))
  
  ggsave(file_name_id, shannon_plt, width = 8, height = 6, dpi = 300)
  
  print(shannon_plt)
}

##########################################################################################
# Gut microbiome and menses

menses_participants <- fecal.microbial.menses.24 %>% 
  filter(menses_day=="menses")
menses_participants <- unique(menses_participants$biome_id)
survey_menses_participants <- fecal.microbial.menses.24 %>% 
  filter(survey_menstruate==1)
survey_menses_participants <- unique(survey_menses_participants$biome_id)
menses_participants <- unique(c(menses_participants, survey_menses_participants))

fecal.microbial.menses.24.filtered <- fecal.microbial.menses.24 %>%
  mutate(
    person_menstruates = ifelse(
      biome_id %in% menses_participants,
      "Menstruates",
      "No menstruate"
    )
  )

length(unique(fecal.microbial.menses.24.filtered$biome_id))

fecal.microbial.menses.24.filtered <- study_days(fecal.microbial.menses.24.filtered)
length(unique(fecal.microbial.menses.24.filtered$biome_id))

# Gut: scatterplot of gut microbiome shannon diversity and menses
# ggplot(fecal.microbial.menses.24.filtered, aes(x = as.Date(logDate), y = shannon)) +
ggplot(fecal.microbial.menses.24.filtered, aes(x = study_day, y = shannon)) +
  geom_point(aes(color=as.factor(menses_day)), alpha=0.5) +
  # geom_smooth(method = "lm", se=FALSE, aes(color = as.factor(menses_day))) +
  geom_smooth(method = "loess", se=FALSE, aes(color = as.factor(menses_day))) +
  labs(
    x = "Study Day", 
    y = "Shannon Diversity",
    title = "",
    color="Menstruating Status"
  ) +
  theme_minimal() +
  scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "1 day") +
  scale_color_discrete(labels = c("menses" = "Menstruating", "not_menses" = "Not menstruating")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),
        text=element_text(size=18)) +
  scale_x_continuous(breaks = seq(0, max(fecal.microbial.menses.24.filtered$study_day), by = 5))

##########################################################################################

## Changes in Genus level analysis in Menses

# relative abundance
bacterial.ra <- transform_sample_counts(bacterial.data, function(x) x / sum(x))  

# aggregate to genus level
genus.ra <- tax_glom(bacterial.ra, taxrank = "Genus")
genus.df <- psmelt(genus.ra) %>%
  rename(RelativeAbundance = Abundance)

dim(genus.df) # 563661 samples and otu combined

# merge with menses data
fecal.microbial.menses.24.subset <- fecal.microbial.menses.24.filtered %>% 
  select(SampleID, shannon, biome_id, logDate, sampleType, regular_periods, menses_status, menses_day)
genus.df.menses <- genus.df %>% 
  left_join(fecal.microbial.menses.24.subset)

# collapse data to be paired
genus.df.menses.paired <- genus.df.menses %>% 
  filter(!is.na(menses_day)) %>% 
  group_by(biome_id, menses_day, Genus) %>% 
  summarise(avg_shannon = sum(shannon) / n(),
            avg_ra = sum(RelativeAbundance) / n(),
            .groups="drop") %>% 
  group_by(biome_id, Genus) %>% 
  filter(n_distinct(menses_day) == 2) %>% 
  ungroup() %>%
  # reshape df
  pivot_wider(names_from = menses_day, values_from = c(avg_ra, avg_shannon),
              names_sep="_") %>%
  drop_na()

genus_stats <- genus.df.menses.paired %>%
  group_by(Genus) %>%
  summarise(
    p_value = t.test(avg_ra_menses, avg_ra_not_menses, paired = TRUE)$p.value
  ) %>%
  arrange(p_value)

# Menses: Log fold change - volcano plot
genus_stats <- genus.df.menses.paired %>%
  group_by(Genus) %>%
  summarise(
    # Log fold change
    fold_change = log2(mean(avg_ra_menses) / mean(avg_ra_not_menses)),  
    p_value = t.test(avg_ra_menses, avg_ra_not_menses, paired = TRUE)$p.value
  ) %>%
  # bonferonni correction
  mutate(
    p_adj = p.adjust(p_value, method = "bonferroni"),
    significance = ifelse(p_adj < 0.05, "Significant", "Not Significant")
  ) %>% 
  arrange(p_adj)

genus_stats.filtered <- genus_stats %>% 
  filter(!is.na(p_adj)) %>% 
  filter(fold_change != -Inf & fold_change != Inf) %>% 
  mutate(
    significance = ifelse(p_adj < 1, "Significant", "Not Significant"),
    color = case_when(
      significance == "Significant" & fold_change > 0 ~ "darkgray",    
      significance == "Significant" & fold_change < 0 ~ "blue",  
      TRUE ~ "gray" 
    )
  )

genus_stats.filtered <- genus_stats.filtered %>%
  mutate(hjust_text = ifelse(fold_change < 0, 1.1, -0.1))

# genus_stats <- genus.df.menses.paired %>%
#   group_by(Genus) %>%
#   summarise(
#     # Log fold change
#     fold_change = log2(mean(avg_ra_menses) / mean(avg_ra_not_menses)),  
#     p_value = t.test(avg_ra_menses, avg_ra_not_menses, paired = TRUE)$p.value
#   ) %>%
#   mutate(significance = ifelse(p_value < 0.05, "Significant", "Not Significant")) 

# genus_stats.filtered <- genus_stats %>% 
#   filter(!is.na(p_value)) %>% 
#   filter(fold_change != -Inf & fold_change != Inf) %>% 
#   mutate(
#     significance = ifelse(p_value < 0.05, "Significant", "Not Significant"),
#     color = case_when(
#       significance == "Significant" & fold_change > 0 ~ "red",    
#       significance == "Significant" & fold_change < 0 ~ "blue",  
#       TRUE ~ "gray" 
#     )
#   )

library(ggrepel)

# Gut: Log fold change plot in menses and not in menses
ggplot(genus_stats.filtered, aes(x = fold_change, y = -log10(p_adj), color = color)) +
  geom_point(alpha = 0.7) +
  # geom_segment(data = subset(genus_stats.filtered, significance == "Significant"), 
  #              aes(xend = fold_change, yend = -log10(p_adj)), 
  #              linetype = "dashed", color = "black", size = 0.5) + 
  # geom_text(data = subset(genus_stats.filtered, significance == "Significant"), 
  #           aes(label = Genus), vjust = -0.5, hjust = 1.2, size = 3) +  
  geom_text_repel(data = subset(genus_stats.filtered, significance == "Significant"), 
                  aes(label = Genus), size = 5, box.padding = 0.5, point.padding = 0.2, max.overlaps = 15) +  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (Menses vs. Not Menses)",
    y = "p-value (-log10)",
    title = " ",
    color = "Significance"
  ) +
  xlim(-10, 10) +
  scale_color_identity() +
  theme(text=element_text(size=18))

# top genus change
top_genera <- genus_stats %>%
  filter(p_value < 0.05) %>%
  arrange(p_value)


##########################################################################################

## Corr with DASS data/stress

### Imputed and cleaned DASS data
dass <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/DASS_0503_2024-final_df.csv")
dass <- dass %>%
  rename(biome_id=study_id)

# add stress severity
dass$stressseverity[dass$stress_score>=0 & dass$stress_score<=14] <- 0
dass$stressseverity[dass$stress_score>=15 & dass$stress_score<=18] <- 1
dass$stressseverity[dass$stress_score>=19 & dass$stress_score<=25] <- 2
dass$stressseverity[dass$stress_score>=26 & dass$stress_score<=33] <- 3
dass$stressseverity[dass$stress_score>=34] <- 4

dass2 <- dass %>% 
  mutate(Timestamp = as.Date(Timestamp)) %>% 
  rename(mood_date = Timestamp) %>% 
  select(biome_id, week, mood_date, stress_score, stressseverity, cisWoman, sport, probiotic, sexuallyActive)

shannon.dass <- fecal.microbial.menses.24.subset %>% # issue week variable shows up twice?
  left_join(dass2, by = c("biome_id")) %>%
  mutate(
    diff_days = abs(as.numeric(difftime(logDate, mood_date, units = "days")))
  ) %>%
  group_by(biome_id, logDate) %>%
  # keep closest survey
  filter(diff_days == min(diff_days)) %>%
  ungroup() %>%
  # only keep within 7 day survey
  mutate(stress_score = ifelse(diff_days > 7, NA, stress_score))
dim(shannon.dass)

# filter for dupes
dupes <- shannon.dass %>% 
  filter(duplicated(SampleID) | duplicated(SampleID, fromLast=TRUE))
dim(dupes) # 0 dupes, surveys are same days apart

unique_pairs <- dupes %>% 
  group_by(biome_id, SampleID) %>% 
  count()
dim(unique_pairs)

dupes$sampleWeek <- rep(NA, nrow(dupes))
dupes$sampleWeek[dupes$logDate >= 19276 & dupes$logDate <= 19280] <- 1
dupes$sampleWeek[dupes$logDate >= 19281 & dupes$logDate <= 19286] <- 2
dupes$sampleWeek[dupes$logDate >= 19290 & dupes$logDate <= 19293] <- 3
dupes$sampleWeek[dupes$logDate >= 19296 & dupes$logDate <= 19300] <- 4
dupes$sampleWeek[dupes$logDate >= 19302 & dupes$logDate <= 19308] <- 5
dupes$sampleWeek[dupes$logDate >= 19312 & dupes$logDate <= 19315] <- 6
dupes$sampleWeek[dupes$logDate >= 19319 & dupes$logDate <= 19321] <- 7
dupes$sampleWeek[dupes$logDate >= 19323 & dupes$logDate <= 19329] <- 8
dupes$sampleWeek[dupes$logDate >= 19330 & dupes$logDate <= 19336] <- 9
dupes$sampleWeek[dupes$logDate >= 19339] <- 10
table(dupes$sampleWeek)

dupes.filtered <- dupes %>%
  filter(week==sampleWeek) %>%
  ungroup()
dim(dupes.filtered)

# take out dupes
shannon.dass.filtered <- shannon.dass %>%
  anti_join(dupes)
dim(shannon.dass.filtered)
dim(shannon.dass)
names(shannon.dass.filtered)

# add back unique
shannon.dass.filtered <- shannon.dass.filtered %>% 
  left_join(dupes.filtered)
dim(shannon.dass.filtered) # 1321   35
names(shannon.dass.filtered)

# filter samples with no survey
shannon.dass.filtered2 <- shannon.dass.filtered %>% 
  filter(!is.na(stress_score))
length(unique(shannon.dass.filtered2$biome_id))

# average stress score
dass.avg <- shannon.dass.filtered %>% 
  group_by(biome_id) %>% 
  summarise(
    avg_stress_score=sum(stress_score, na.rm=TRUE)/n(),
    avg_shannon=sum(shannon, na.rm=TRUE)/n(),
  ) %>% 
  mutate(stress_severity = case_when(
    avg_stress_score <= 14 ~ "Normal",
    avg_stress_score >= 15 & avg_stress_score <= 18 ~ "Mild",
    avg_stress_score >= 19 & avg_stress_score <= 25 ~ "Moderate",
    avg_stress_score >= 26 & avg_stress_score <= 33 ~ "Severe",
    avg_stress_score >= 34 ~ "Extremely Severe",
    TRUE ~ NA_character_
  )) %>% 
  filter(!is.na(avg_stress_score))

# differences in avg shannon for stress severity category
lm.obj <- lm(avg_shannon~stress_severity, data=dass.avg)
summary(lm.obj)
TukeyHSD(aov(avg_shannon ~ stress_severity, data = dass.avg))

dass.avg$stress_severity <- factor(dass.avg$stress_severity, levels=c("Normal", "Mild", "Moderate", "Severe", "Extremely Severe"))
dim(dass.avg)

# Gut: avg shannon by avg stress severity category
dass.avg %>% 
  filter(!is.na(stress_severity)) %>% 
  ggplot(aes(x = factor(stress_severity), y = avg_shannon, 
                                   col=as.factor(biome_id))) +
  geom_boxplot(fill = "skyblue", color = "black", alpha=0.1) + 
  geom_point(position = position_jitter(width = 0.4, height = 0), alpha = 0.7) +
  labs(x = " ", y = "Average Shannon Diversity") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("0" = "Normal", "1" = "Mild", "2" = "Moderate", "3" = "Severe", "4" = "Extremely Severe")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        text=element_text(size=18))

# Gut: shannon by stress severity category
ggplot(shannon.dass.filtered2, aes(x = factor(stressseverity), y = shannon, 
                                   col=as.factor(biome_id))) +
  geom_boxplot(fill = "skyblue", color = "black", alpha=0.1) + 
  geom_point(position = position_jitter(width = 0.4, height = 0), alpha = 0.7) +
  labs(x = " ", y = "Shannon Diversity") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("0" = "Normal", "1" = "Mild", "2" = "Moderate", "3" = "Severe", "4" = "Extremely Severe")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        text=element_text(size=18))


# Gut: shannon by stress severity category
ggplot(shannon.dass.filtered2, aes(x = factor(stressseverity), y = shannon, 
                                   col=as.factor(biome_id))) +
  geom_boxplot(fill = "white", color = "black") + 
  geom_point(position = position_jitter(width = 0.4, height = 0), alpha = 0.7) +
  labs(x = "Stress Level", y = "Shannon Diversity") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("0" = "Normal", "1" = "Mild", "2" = "Moderate", "3" = "Severe", "4" = "Extremely Severe")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        text=element_text(size=18))

## Mixed effects models
library(lme4)
library(lmerTest)
library(performance)

# dim(shannon.dass.filtered2)
# shannon.dass.filtered2 <- study_days(shannon.dass.filtered2)
# 
# names(shannon.dass.filtered2)
# # levels=c("normal","mild","moderate","severe", "extremely severe"), 
# shannon.dass.filtered2$stressseverity <- factor(shannon.dass.filtered2$stressseverity, ordered = TRUE)
# contrasts(shannon.dass.filtered2$stressseverity) <- contr.poly(5)
# model <- lmer(shannon ~ stressseverity + (1|`biome_id`), data = shannon.dass.filtered2)
# 
# summary(model)
# library(emmeans)
# emmeans::emm(model, pairwise ~ stressseverity)

lmer.null <- lmer(shannon ~ (1|`biome_id`), data = shannon.dass.filtered2)
r2(lmer.null)
# lmer.shannon.stressscore <- lmer(shannon~stress_score + (1|`biome_id`), data = shannon.dass.filtered2)
# r2(lmer.shannon.stressscore)
# summary(lmer.shannon.stressscore)

# numerics
lmer.shannon.stress <- lmer(shannon~as.factor(stressseverity) + (1|`biome_id`), data = shannon.dass.filtered2)
r2(lmer.shannon.stress)
summary(lmer.shannon.stress)

anova(lmer.null, lmer.shannon.stress)

shannon.dass.filtered2 <- study_days(shannon.dass.filtered2)
lmer.shannon.stress.time <- lmer(shannon~study_day + as.factor(stressseverity) + (1|`biome_id`), data = shannon.dass.filtered2)
r2(lmer.shannon.stress.time)
summary(lmer.shannon.stress.time)

anova(lmer.shannon.stress, lmer.shannon.stress.time)

# quad time
lmer.shannon.stress.time2 <- lmer(shannon~study_day + I(study_day)^2 + as.factor(stressseverity) + (1|`biome_id`), data = shannon.dass.filtered2)
r2(lmer.shannon.stress.time2)
summary(lmer.shannon.stress.time2)

anova(lmer.shannon.stress, lmer.shannon.stress.time)

# categories
lmer.shannon.stress <- lmer(shannon~as.factor(stressseverity) + (1|`biome_id`), data = shannon.dass.filtered2)
r2(lmer.shannon.stress)
summary(lmer.shannon.stress)

lmer.shannon.stress.time <- lmer(shannon~study_day + as.factor(stressseverity) + (1|`biome_id`), data = shannon.dass.filtered2)
r2(lmer.shannon.stress.time)
summary(lmer.shannon.stress.time)

anova(lmer.shannon.stress.time, lmer.shannon.stress.time2)

lmer.shannon.stress.time2 <- lmer(shannon~study_day + I(study_day^2) + stressseverity + (1|`biome_id`), data = shannon.dass.filtered2)
r2(lmer.shannon.stress.time2)
summary(lmer.shannon.stress.time2)

anova(lmer.shannon.stress.time, lmer.shannon.stress.time2)

# random slopes
lmer.shannon.stress <- lmer(shannon~stressseverity + (stressseverity|`biome_id`), data = shannon.dass.filtered2)
r2(lmer.shannon.stress)
summary(lmer.shannon.stress)

### High and low stress comparisons
table(shannon.dass.filtered2$stressseverity)
shannon.dass.filtered2 <- shannon.dass.filtered2 %>% 
  mutate(stressbinary = ifelse(stressseverity %in% c(2,3,4), "High Stress", "Low Stress"))
table(shannon.dass.filtered2$stressbinary)

# alpha div change significantly in low v. high? - t-test

# average shannon for each participant for high and low stress
# shannon.dass.filtered2.summary <- shannon.dass.filtered2 %>% 
#   group_by(biome_id, stressbinary) %>% 
#   summarise(avg_shannon=mean(shannon), .groups="drop") %>% 
#   pivot_wider(names_from = stressbinary, values_from = avg_shannon) %>%  #, names_prefix = "stress_")
#   filter(!is.na(`High Stress`) & !is.na(`Low Stress`))
# length(unique(shannon.dass.filtered2.summary$biome_id))

### Wilcox test
# wilcox.test(shannon.dass.filtered2.summary$`High Stress`,
#             shannon.dass.filtered2.summary$`Low Stress`, paired=TRUE)
# t.test(shannon.dass.filtered2.summary$`High Stress`,
#        shannon.dass.filtered2.summary$`Low Stress`, paired=TRUE)

# Does Shannon diversity vary by stress 

# Gut: Average Shannon Diversity for Participants high and low stress
# ggplot(
#   shannon.dass.filtered2.summary,
#   aes(y = `High Stress`, x = `Low Stress`)) +
#   geom_point(color="orchid3", show.legend=FALSE) +
#   xlim(3,4)+
#   ylim(3,4)+
#   geom_abline(slope=1, intercept=0, linetype="dashed") +
#   theme_minimal() +
#   theme(text=element_text(size=18))+
#   labs(y="High Stress (Shannon)", x="Low Stress (Shannon)",
#        color="Participant ID")

# ## Log fold change - Changes in Genus level analysis in Menses
# 
# # relative abundance
# bacterial.ra <- transform_sample_counts(bacterial.data, function(x) x / sum(x))  
# 
# # aggregate to genus level
# genus.ra <- tax_glom(bacterial.ra, taxrank = "Genus")
# genus.df <- psmelt(genus.ra) %>%
#   rename(RelativeAbundance = Abundance)
# 
# dim(genus.df) # 563661 samples and otu combined
# 
# # merge with stress data
# genus.df.stress <- genus.df %>% 
#   left_join(shannon.dass.filtered2)
# 
# # collapse data to be paired
# genus.df.stress.paired <- genus.df.stress %>% 
#   filter(!is.na(stressbinary)) %>% 
#   group_by(biome_id, stressbinary, Genus) %>% 
#   summarise(avg_shannon = sum(shannon) / n(),
#             avg_ra = sum(RelativeAbundance) / n(),
#             .groups="drop") %>% 
#   group_by(biome_id, Genus) %>% 
#   filter(n_distinct(stressbinary) == 2) %>% 
#   ungroup() %>%
#   # reshape df
#   pivot_wider(names_from = stressbinary, values_from = c(avg_ra, avg_shannon),
#               names_sep="_") %>%
#   drop_na()
# 
# genus_stats <- genus.df.stress.paired %>%
#   group_by(Genus) %>%
#   summarise(
#     p_value = t.test(`avg_ra_High Stress`, `avg_ra_Low Stress`, paired = TRUE)$p.value
#   ) %>%
#   arrange(p_value)
# 
# # Stress: Log fold change - volcano plot
# genus_stats <- genus.df.stress.paired %>%
#   group_by(Genus) %>%
#   summarise(
#     # Log fold change
#     fold_change = log2(mean(`avg_ra_High Stress`) / mean(`avg_ra_Low Stress`)),  
#     p_value = t.test(`avg_ra_High Stress`, `avg_ra_Low Stress`, paired = TRUE)$p.value
#   ) %>%
#   mutate(significance = ifelse(p_value < 0.05, "Significant", "Not Significant")) 
# 
# genus_stats.filtered <- genus_stats %>% 
#   filter(!is.na(p_value)) %>% 
#   filter(fold_change != -Inf & fold_change != Inf) %>% 
#   mutate(
#     significance = ifelse(p_value < 0.05, "Significant", "Not Significant"),
#     color = case_when(
#       significance == "Significant" & fold_change > 0 ~ "red",    
#       significance == "Significant" & fold_change < 0 ~ "blue",  
#       TRUE ~ "gray" 
#     )
#   )

# library(ggrepel)
# # Gut: Log fold change plot in high and low stress
# ggplot(genus_stats.filtered, aes(x = fold_change, y = -log10(p_value), color = color)) +
#   geom_point(alpha = 0.7) +
#   geom_segment(data = subset(genus_stats.filtered, significance == "Significant"), 
#                aes(xend = fold_change, yend = -log10(p_value)), 
#                linetype = "dashed", color = "black", size = 0.5) +  
#   geom_text_repel(data = subset(genus_stats.filtered, significance == "Significant"), 
#                   aes(label = Genus), size = 5, box.padding = 0.5, point.padding = 0.2, max.overlaps = 15) +  
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + 
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
#   theme_minimal() +
#   labs(
#     x = "Log2 Fold Change (High Stress vs. Low Stress)",
#     y = "p-value (-log10)",
#     title = " ",
#     color = "Significance"
#   ) +
#   xlim(-10, 10) +
#   scale_color_identity() +
#   theme(text=element_text(size=18))
# 
# # top genus change
# top_genera <- genus_stats %>%
#   filter(p_value < 0.05) %>%
#   arrange(p_value)
# top_genera

##########################################################################################


### Birth Control

## Corr with volunteer history data
participant.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 9-Volunteer Medical History.csv", header = TRUE)

# filter non-hormonal & no samples
participant.data <- participant.data %>% 
  filter(birthControl!="Orilissa (Elagolix)")
# select birth control
birthControl.df <- participant.data %>% 
  select(biome_id, birthControl)

# Join dfs
shannon.dass.filtered.birthControl <- shannon.dass.filtered2 %>% 
  left_join(birthControl.df)

length(unique(shannon.dass.filtered.birthControl$biome_id))

# collapse birth control
shannon.dass.filtered.birthControl <- shannon.dass.filtered.birthControl %>% 
  mutate(birthControl_collapsed=ifelse(birthControl=="Systemic Combined (E&P)" | (birthControl=="Systemic P only"), "Systemic", 
                                       birthControl))

# relevel
shannon.dass.filtered.birthControl$birthControl_collapsed <- 
  factor(shannon.dass.filtered.birthControl$birthControl_collapsed, 
         levels=c("None", "Local P", "Systemic"))

shannon.dass.filtered.birthControl$birthControl <- 
  factor(shannon.dass.filtered.birthControl$birthControl, 
         levels=c("None", "Local P", "Systemic P only", "Systemic Combined (E&P)"))

shannon.dass.filtered.birthControl.avg <- shannon.dass.filtered.birthControl %>% 
  group_by(biome_id) %>% 
  summarise(avg_shannon=sum(shannon, na.rm=TRUE) / n(),
            birthControl = first(birthControl)) %>% 
  filter(!is.na(birthControl))

# Gut: avg shannon for birth controls
ggplot(shannon.dass.filtered.birthControl.avg, aes(x = factor(birthControl), y = avg_shannon, 
             col=as.factor(biome_id))) +
  geom_boxplot(fill = "skyblue", color = "black", alpha = 0.1) + 
  geom_point(position = position_jitter(width = 0.4, height = 0), alpha = 0.7) +
  labs(x = " ", y = "Average Shannon Diversity") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("0" = "Normal", "1" = "Mild", "2" = "Moderate", "3" = "Severe", "4" = "Extremely Severe")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        text=element_text(size=18))

# Gut: boxplot of birth control and shannon diversity
ggplot(shannon.dass.filtered.birthControl, aes(x = birthControl, y = shannon)) +
  geom_jitter(aes(color=as.factor(biome_id)), size=1, alpha=0.6) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
  labs(x = " ", y = "Shannon Diversity Index", title = "",
       color = "Biome ID") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),
        legend.position = "None", text = element_text(size = 18)) 

# add days
shannon.dass.filtered.birthControl <- study_days(shannon.dass.filtered.birthControl)

## Spline model - Non aggregate

spline.obj <- shannon.dass.filtered.birthControl %>%
  filter(!is.na(shannon)) %>% 
  filter(!is.na(stress_score)) %>% 
  group_by(birthControl) %>%
  # filter(n_distinct(stress_score) >= 4) %>%
  summarise(
    model = list(smooth.spline(x=stress_score, y=shannon, spar=0.7)),
    .groups = "drop"
  ) %>%
  mutate(
    x = map(model, ~ seq(min(.x$x), max(.x$x), length.out = 100)),
    y = map2(model, x, ~ predict(.x, .y)$y)
  ) %>%
  select(birthControl, x, y) %>%
  unnest(c(x, y))

# Gut: Scatterplot of shannon by stress with splines by birth control
ggplot(shannon.dass.filtered.birthControl, aes(x = stress_score, y = shannon)) +
  # geom_point(aes(color=as.factor(birthControl)), size=1, alpha=0.5) +
  geom_jitter(aes(color=as.factor(birthControl)), size=1, alpha=0.6) +
  geom_line(data = spline.obj, aes(x = x, y = y, color = birthControl), linewidth = 1) +
  # scale_color_viridis_d(option="D") +
  labs(x = "Stress Score", y = "Shannon Diversity Index", title = "",
       color = "Birth Control") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),
        text=element_text(size=18))

## add aggregate models (avg stress and avg shannon)

shannon.dass.filtered.birthControl.summary <- shannon.dass.filtered.birthControl %>% 
  group_by(biome_id) %>% 
  summarise(avg_shannon = sum(shannon)/n(),
            avg_stressscore = sum(stress_score)/n(),
            avg_stressseverity = sum(stressseverity)/n(),
            sport=first(sport),
            cisWoman=first(cisWoman),
            probiotic=first(probiotic),
            sexuallyActive=first(sexuallyActive),
            birthControl=first(birthControl),
            birthControl_collapsed = first(birthControl_collapsed))

# Regressing models
lm.obj <- lm(avg_shannon ~ avg_stressscore*birthControl, data=shannon.dass.filtered.birthControl.summary)
summary(lm.obj)
lm.obj2 <- lm(avg_shannon ~ avg_stressseverity*birthControl, data=shannon.dass.filtered.birthControl.summary)
summary(lm.obj2)

lm.obj3 <- lm(avg_shannon ~ avg_stressscore+birthControl, data=shannon.dass.filtered.birthControl.summary)
summary(lm.obj3)
lm.obj4 <- lm(avg_shannon ~ avg_stressseverity+birthControl, data=shannon.dass.filtered.birthControl.summary)
summary(lm.obj4)

anova(lm.obj,lm.obj3)
anova(lm.obj2,lm.obj4)

aov.obj <- aov(shannon.dass.filtered.birthControl.summary$avg_shannon ~ shannon.dass.filtered.birthControl.summary$birthControl)
summary(aov.obj)

## Collapsed birth control
table(shannon.dass.filtered.birthControl.summary$birthControl_collapsed)
lm.obj <- lm(avg_shannon ~ avg_stressscore*birthControl_collapsed, data=shannon.dass.filtered.birthControl.summary)
summary(lm.obj)
lm.obj2 <- lm(avg_shannon ~ avg_stressseverity*birthControl_collapsed, data=shannon.dass.filtered.birthControl.summary)
summary(lm.obj2)

lm.obj3 <- lm(avg_shannon ~ avg_stressscore+birthControl_collapsed, data=shannon.dass.filtered.birthControl.summary)
summary(lm.obj3)
lm.obj4 <- lm(avg_shannon ~ avg_stressseverity+birthControl_collapsed, data=shannon.dass.filtered.birthControl.summary)
summary(lm.obj4)

anova(lm.obj,lm.obj3)
anova(lm.obj2,lm.obj4)

aov.obj <- aov(shannon.dass.filtered.birthControl.summary$avg_shannon ~ shannon.dass.filtered.birthControl.summary$birthControl_collapsed)
summary(aov.obj)

# Mixed effects models
library(lme4)
library(lmerTest)
library(performance)

lmer.null <- lmer(shannon~ (1|`biome_id`), data = shannon.dass.filtered.birthControl)
r2(lmer.null)
summary(lmer.null)

lmer.shannon.bc <- lmer(shannon~as.factor(birthControl) + (1|`biome_id`), data = shannon.dass.filtered.birthControl)
r2(lmer.shannon.bc)
summary(lmer.shannon.bc)

anova(lmer.null, lmer.shannon.bc)

lmer.shannon.bc2 <- lmer(shannon~birthControl*stressseverity + (1|`biome_id`), data = shannon.dass.filtered.birthControl)
r2(lmer.shannon.bc2)
summary(lmer.shannon.bc2)

lmer.shannon.bc2 <- lmer(shannon~birthControl*stressseverity + (1|`biome_id`), data = shannon.dass.filtered.birthControl)
r2(lmer.shannon.bc2)
summary(lmer.shannon.bc2)

# collapsed birth control
lmer.shannon.bc_collapsed <- lmer(shannon~birthControl_collapsed + (1|`biome_id`), data = shannon.dass.filtered.birthControl)
r2(lmer.shannon.bc_collapsed)
summary(lmer.shannon.bc_collapsed)
lmer.shannon.bc2_collapsed <- lmer(shannon~birthControl_collapsed*stressseverity + (1|`biome_id`), data = shannon.dass.filtered.birthControl)
r2(lmer.shannon.bc2_collapsed)
summary(lmer.shannon.bc2_collapsed)

# add time
lmer.shannon.bc.time <- lmer(shannon~study_day + birthControl + (1|`biome_id`), data = shannon.dass.filtered.birthControl)
r2(lmer.shannon.bc.time)
summary(lmer.shannon.bc.time)

lmer.shannon.bc.time2 <- lmer(shannon~study_day + I(study_day^2) + birthControl + (1|`biome_id`), data = shannon.dass.filtered.birthControl)
r2(lmer.shannon.bc.time2)
summary(lmer.shannon.bc.time2)

anova(lmer.shannon.bc, lmer.shannon.bc.time)
anova(lmer.shannon.bc.time, lmer.shannon.bc.time2)

# random slopes
lmer.shannon.bc.slopes <- lmer(shannon~birthControl + (birthControl|`biome_id`), data = shannon.dass.filtered.birthControl)
r2(lmer.shannon.bc.slopes)
summary(lmer.shannon.bc.slopes)

lmer.shannon.bc.slopes <- lmer(shannon~birthControl_collapsed + (birthControl_collapsed|`biome_id`), data = shannon.dass.filtered.birthControl)
r2(lmer.shannon.bc.slopes)
summary(lmer.shannon.bc.slopes)

# Try birth control v. none
shannon.birthControl.binary <- shannon.birthControl %>% 
  mutate(bc_binary = ifelse(birthControl=="None", "None", "birtControl"))
lmer.shannon.bc.binary <- 
  lmer(shannon~bc_binary + (1|`biome_id`), data = shannon.birthControl.binary)
r2(lmer.shannon.bc.binary)
summary(lmer.shannon.bc.binary)

#### Make full df
full.data <- fecal.microbial.menses.24 %>%
  left_join(shannon.dass.filtered.birthControl)

# Figure: boxplot of collapsed systemic birth control and shannon diversity
# ggplot(shannon.birthControl.collapsed, aes(x = birthControl_collapsed, y = shannon)) +
#   geom_jitter(aes(color=as.factor(biome_id)), size=1, alpha=0.6) +
#   geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
#   scale_color_viridis_d(option="D") +
#   labs(x = " ", y = "Shannon Diversity Index", title = "",
#        color = "Biome ID") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),
#         text=element_text(size=18),
#         legend.position="None")

# Average Shannon diversity per person
shannon.birthControl.avg <- full.data %>%
  group_by(biome_id) %>% 
  summarise(avg_shannon = sum(shannon) / n(),
            birthControl = first(birthControl),
            birthControl_collapsed = first(birthControl_collapsed),
            avg_stressscore = sum(stress_score, na.rm=TRUE)/n(),
            avg_stressseverity = sum(stressseverity, na.rm=TRUE)/n())

# Gut: boxplot of birth control and avg shannon diversity
shannon.birthControl.avg %>% 
  filter(!is.na(birthControl)) %>% 
  ggplot(aes(x = birthControl, y = avg_shannon)) +
  geom_jitter(aes(color=as.factor(biome_id)), size=1, alpha=0.6) +
  # geom_point(aes(color=as.factor(biome_id)), size=1, alpha=0.7) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
  # geom_text(aes(label = as.factor(biome_id)), vjust = -0.5, hjust=1.5, size = 1.5) +
  scale_color_viridis_d(option="D") +
  labs(x = " ", y = "Average Shannon Diversity Index", title = "",
       color = "Biome ID") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),
        text=element_text(size=18),
        legend.position="None")
  # theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1))

# Figure: boxplot of collapsed Systemic birth control and avg shannon diversity
shannon.birthControl.avg %>% 
  filter(!is.na(birthControl)) %>% 
  ggplot(aes(x = birthControl_collapsed, y = avg_shannon)) +
  geom_jitter(aes(color=as.factor(biome_id)), size=1, alpha=0.6) +
  # geom_point(aes(color=as.factor(biome_id)), size=1, alpha=0.7) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
  # geom_text(aes(label = as.factor(biome_id)), vjust = -0.5, hjust=1.5, size = 1.5) +
  scale_color_viridis_d(option="D") +
  labs(x = " ", y = "Average Shannon Diversity Index", title = "",
       color = "Biome ID") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),
        text=element_text(size=18),
        legend.position="None")

### Alpha diversity over time, color maybe by birth control, menstruate, or both

# merge df
full.data$logDate <- as.Date(full.data$logDate)

full.data.filter <- full.data %>% 
  filter(!is.na(birthControl))
length(unique(full.data.filter$biome_id))


# Gut: scatterplot shannon over days by birth control
full.data %>% 
  filter(!is.na(birthControl)) %>% 
  ggplot(aes(x = logDate, y = shannon)) +
  geom_point(aes(color=birthControl), alpha=0.5) +
  geom_smooth(method = "lm", se=FALSE, aes(color = birthControl)) +
  labs(
    x = "Days", 
    y = "Shannon Diversity Index",
    title = "",
    color="Birth Control"
  ) +
  theme_minimal() +
  scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "1 day") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

##########################################################################################

## Diet Data
merged_diet_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/relabeled_data/nutrition_withVegPerc_data_3_20.csv")
past_two_days_diet_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/relabeled_data/past_two_days_diet_data_3_20.csv")

merged_diet_data <- merged_diet_data %>% 
  rename(logDate=Date) %>% 
  select(-X)
past_two_days_diet_data <- past_two_days_diet_data %>% 
  rename(logDate=Date) %>% 
  select(-X)

# filter data
merged_diet_data <- filter_days(merged_diet_data)
past_two_days_diet_data <- filter_days(past_two_days_diet_data)

dim(merged_diet_data) # 9332   42
dim(past_two_days_diet_data) # 1323   19

# collapse diet data
merged_diet_data_collapsed <- merged_diet_data %>% 
  group_by(biome_id, logDate) %>% 
  summarise(
    caloriesall = sum(caloriesall, na.rm=TRUE),
    cholesterolall = sum(cholesterolall, na.rm=TRUE),
    saturatedFatall = sum(saturatedFatall, na.rm=TRUE),
    sodiumall = sum(sodiumall, na.rm=TRUE),
    carbohydratesall = sum(carbohydratesall, na.rm=TRUE),
    dietaryFiberall = sum(dietaryFiberall, na.rm=TRUE),
    sugarsall = sum(sugarsall, na.rm=TRUE),
    proteinall = sum(proteinall, na.rm=TRUE),
    fatall = sum(fatall, na.rm=TRUE),
    caloriesFromFat = sum(caloriesFromFat, na.rm=TRUE),
    saturatedFat = sum(saturatedFat, na.rm=TRUE),
    caloriesFromSatFat = sum(caloriesFromSatFat, na.rm=TRUE),
    transFat = sum(transFat, na.rm=TRUE),
    addedSugarall = sum(addedSugarall, na.rm=TRUE),
    perc_veg = first(perc_veg),
    is_vegetarian = first(is_vegetarian)
  )

# fecal microbiota and specific nutrient intake

nutrient.diet.22 <-  merged_diet_data %>% #past_two_days_diet_data %>% #
  group_by(biome_id, logDate) %>% 
  summarise(caloriesall = sum(caloriesall, na.rm=TRUE),
            cholesterol_prop = sum(cholesterolall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            satFat_prop = sum(saturatedFatall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            sodium_prop = sum(sodiumall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            carb_prop = sum(carbohydratesall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            dietFib_prop = sum(dietaryFiberall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            sugar_prop = sum(sugarsall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            protein_prop = sum(proteinall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            fat_prop = sum(fatall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            # sat fat of all fat prop
            satfat_prop_fat = sum(saturatedFatall, na.rm=TRUE) / sum(fatall, na.rm=TRUE),
            fat_cal_prop = sum(caloriesFromFat, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            addedSugarall_prop = sum(addedSugarall, na.rm=TRUE) / sum(addedSugarall, na.rm=TRUE)) %>% 
  ungroup() %>% 
  mutate(logDate=as.character(logDate))
# add veg cols to the rolling window df
merged_diet_data_collapsed_subset <- merged_diet_data_collapsed %>%
  select(biome_id, logDate, perc_veg, is_vegetarian)
# past_two_days_diet_data <- past_two_days_diet_data %>% 
#   left_join(merged_diet_data_collapsed_subset, by=c("biome_id", "logDate"))

# join fecal data with diet data
full.data$logDate <- as.character(full.data$logDate)
full.data.subset <- full.data %>% 
  select(SampleID, shannon, biome_id, logDate, survey_menstruate,
         menses_status, menses_day, stress_score, stressseverity, cisWoman, sport, probiotic,
         sexuallyActive, 
         # stressbinary, 
         birthControl, birthControl_collapsed)
full.data.subset.collapsed <- full.data.subset %>% 
  group_by(biome_id, logDate) %>% 
  summarise(daily_shannon = sum(shannon, na.rm=TRUE) / n(),
            survey_menstruate = first(survey_menstruate),
            menses_status = first(menses_status), 
            menses_day = first(menses_day),
            stress_score = first(stress_score),
            stressseverity = first(stressseverity),
            sport = first(sport),
            sexuallyActive = first(sexuallyActive),
            birthControl = first(birthControl))
gut.diet.df <- nutrient.diet.22 %>% # past_two_days_diet_data %>% #
  left_join(full.data.subset.collapsed) %>% 
  filter(!is.na(daily_shannon))
gut.diet.df <- gut.diet.df %>% 
  left_join(merged_diet_data_collapsed_subset)

# aggregate analyses

nutrient.diet.22.shannon <- full.data.subset %>% 
  left_join(nutrient.diet.22)

nutrient_diet_avg_per_person <- nutrient.diet.22.shannon %>%
  group_by(biome_id) %>%
  summarise(across(
    .cols = c(shannon, caloriesall, cholesterol_prop, satFat_prop, sodium_prop, 
              carb_prop, dietFib_prop, sugar_prop, protein_prop, 
              fat_prop, satfat_prop_fat, fat_cal_prop, addedSugarall_prop),
    .fns = ~mean(.x, na.rm = TRUE)
  )) %>%
  ungroup()

lm.obj <- lm(shannon~., nutrient_diet_avg_per_person[,-1])
summary(lm.obj)

gut.diet.df <- study_days(gut.diet.df)
dim(gut.diet.df)
length(unique(gut.diet.df$biome_id))
head(gut.diet.df)

#### vegetarian analysis
gut.diet.summary.df <- gut.diet.df %>% 
  group_by(biome_id) %>% 
  summarise(avg_shannon = sum(daily_shannon, na.rm = TRUE) / n(),
            is_vegetarian = first(is_vegetarian),
            perc_veg = first(perc_veg))
table(gut.diet.summary.df$is_vegetarian)
summary(gut.diet.summary.df$avg_shannon)

# Gut: boxplot of average shannon diversity and vegetarian status
ggplot(gut.diet.summary.df, aes(x = is_vegetarian, y = avg_shannon)) +
  geom_boxplot(outlier.shape = NA, fill="skyblue", alpha=0.1) +
  geom_jitter(aes(color=as.factor(biome_id)), width = 0.2, alpha = 0.6, show.legend = FALSE) +
  labs(
    x = " ", 
    y = "Average Shannon Diversity", 
    title = ""
  ) +
  theme_minimal() +
  scale_x_discrete(labels = c("FALSE" = "Not vegetarian", "TRUE" = "Vegetarian")) +
  theme(text=element_text(size=18))

ggplot(gut.diet.df, aes(x = is_vegetarian, y = daily_shannon)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color="orchid") +
  labs(
    x = " ", 
    y = "Shannon Diversity", 
    title = ""
  ) +
  theme_minimal() +
  scale_x_discrete(labels = c("FALSE" = "Not vegetarian", "TRUE" = "Vegetarian")) +
  theme(text=element_text(size=18))

# testing
t.test(avg_shannon ~ is_vegetarian, data=gut.diet.summary.df)
wilcox.test(avg_shannon ~ is_vegetarian, data = gut.diet.summary.df)

### PERCENTAGE VEGETARIAN

# average shannon and percent veg
veg_perc_df_summary  <- gut.diet.df %>% 
  group_by(biome_id) %>% 
  summarise(avg_shannon = sum(daily_shannon, na.rm=TRUE) / n(),
            perc_veg = first(perc_veg))
summary(veg_perc_df_summary$avg_shannon)

lm.obj.summary <- lm(avg_shannon ~ perc_veg, data = veg_perc_df_summary)
summary(lm.obj.summary)

# Join max OTU
fecal.microbial.menses.24.OTU <- fecal.microbial.menses.24 %>% 
  select(biome_id, logDate, OTU)

veg_perc_df_summary2 <- veg_perc_df_summary %>%
  left_join(fecal.microbial.menses.24.OTU, by="biome_id")

# Gut: scatter plot of percent vegetarian on average shannon diversity
ggplot(veg_perc_df_summary2, aes(x = perc_veg, y = avg_shannon)) +
  geom_point(color="orchid") +
  labs(x = "Percent Vegetarian", y = "Average Shannon Diversity", 
       title = " ") +
  theme_minimal() + 
  scale_x_discrete(labels = c("FALSE" = "Vegetarian", "TRUE" = "Vegetarian")) +
  theme(text=element_text(size=18), legend.position = "None") #+
  # ylim(0, 4)

# Gut: scatter plot of percent vegetarian on shannon diversity colored by max taxa
ggplot(gut.diet.df, aes(x = perc_veg, y = daily_shannon)) +
  geom_point(color="orchid") +
  labs(x = "Percent Vegetarian", y = "Shannon Diversity", 
       title = " ") +
  theme_minimal() + 
  scale_x_discrete(labels = c("FALSE" = "Vegetarian", "TRUE" = "Vegetarian")) +
  theme(text=element_text(size=18), legend.position = "None") #+
# ylim(0, 4)

## Mixed Effects Model
library(lme4)
library(lmerTest)
library(performance)

## random intercept, fixed slopes

gut.diet.df.subset <- gut.diet.df %>%
  select(daily_shannon, biome_id, caloriesall, cholesterol_prop, satFat_prop,
         sodium_prop, carb_prop, dietFib_prop, sugar_prop, protein_prop,
         fat_prop, fat_cal_prop, addedSugarall_prop) %>%
  na.omit()

lmer.null <- lmer(daily_shannon~(1|`biome_id`), data = gut.diet.df.subset)
r2(lmer.null)

# Full model
lmer.full <- lmer(daily_shannon~caloriesall + cholesterol_prop + satFat_prop + sodium_prop +
                    carb_prop + dietFib_prop + sugar_prop + protein_prop + fat_prop +
                    fat_cal_prop + addedSugarall_prop +
                    (1|`biome_id`), data = gut.diet.df.subset)
r2(lmer.full)
summary(lmer.full)

anova(lmer.null, lmer.full)

# Time
lmer.time.full <- lmer(daily_shannon~study_day + caloriesall + cholesterol_prop + satFat_prop + sodium_prop +
                    carb_prop + dietFib_prop + sugar_prop + protein_prop + fat_prop +
                    fat_cal_prop + addedSugarall_prop +
                    (1|`biome_id`), data = gut.diet.df)
r2(lmer.time.full)
summary(lmer.time.full)

anova(lmer.full, lmer.time.full)

# Time + Time^2
lmer.time2.full <- lmer(daily_shannon~study_day + I(study_day^2) + caloriesall + cholesterol_prop + satFat_prop + sodium_prop +
                         carb_prop + dietFib_prop + sugar_prop + protein_prop + fat_prop +
                         fat_cal_prop + addedSugarall_prop +
                         (1|`biome_id`), data = gut.diet.df)
r2(lmer.time2.full)
summary(lmer.time2.full)

anova(lmer.time.full, lmer.time2.full)

# Check vegetarian variable
lmer.full.veg <- lmer(daily_shannon~caloriesall + cholesterol_prop + satFat_prop + sodium_prop +
                        carb_prop + dietFib_prop + sugar_prop + protein_prop + fat_prop +
                        perc_veg +
                        fat_cal_prop + addedSugarall_prop +
                        (1|`biome_id`), data = gut.diet.df)
r2(lmer.full.veg)
summary(lmer.full.veg)

lmer.full.veg <- lmer(daily_shannon~caloriesall + cholesterol_prop + satFat_prop + sodium_prop +
                        carb_prop + dietFib_prop + sugar_prop + protein_prop + fat_prop +
                        is_vegetarian +
                        fat_cal_prop + addedSugarall_prop +
                        (1|`biome_id`), data = gut.diet.df)
r2(lmer.full.veg)
summary(lmer.full.veg)

##########################################################################################

## Activity

activity.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 4-Physical Activity.csv")
activity.data <- activity.data %>% 
  filter(!is.na(as.numeric(biome_id)))
activity.data <- study_days(activity.data)

# gut microbiome data
full.data.subset <- full.data.subset %>% 
  left_join(veg_perc_df_summary)
dim(full.data.subset)

# join by date and participant
length(unique(full.data.subset$biome_id))

fecal.microbial.menses.24.activity <- full.data.subset %>% 
  left_join(activity.data) %>% 
  filter(!is.na(steps))
dim(fecal.microbial.menses.24.activity)
length(unique(fecal.microbial.menses.24.activity$biome_id)) # 46

### Summarize activity data
activity.data.summary <- activity.data %>% 
  group_by(biome_id) %>% 
  summarise(avg_cals_burned = sum(calories_burned, na.rm=TRUE) / n(),
            avg_steps = sum(steps, na.rm=TRUE) / n(),
            avg_distance = sum(distance, na.rm=TRUE) / n(),
            avg_minutes_sedentary = sum(minutes_sedentary, na.rm=TRUE) / n(),
            avg_minutes_lightly_active = sum(minutes_lightly_active, na.rm=TRUE) / n(),
            avg_minutes_fairly_active = sum(minutes_fairly_active, na.rm=TRUE) / n(),
            avg_minues_very_active = sum(minues_very_active, na.rm=TRUE) / n(),
            avg_activity_calories = sum(activity_calories, na.rm=TRUE) / n(),
            total_min_active = sum(minutes_lightly_active, minutes_fairly_active, minues_very_active, na.rm=TRUE) / n()
  )
fecal.microbial.menses.24.summary <- full.data.subset %>% 
  group_by(biome_id) %>% 
  summarise(avg_shannon=sum(shannon)/n())

activity.data.summary <- fecal.microbial.menses.24.summary %>% 
  left_join(activity.data.summary) %>% 
  # filter thorugh without activity data
  filter(!is.na(avg_steps))
length(unique(activity.data.summary$biome_id)) # 48

setdiff(unique(activity.data.summary$biome_id), unique(fecal.microbial.menses.24.activity$biome_id))

## Linear regression
library(car)

lm.shannon.activity <- lm(avg_shannon ~ avg_steps+avg_cals_burned+avg_activity_calories+
                            avg_distance+avg_minutes_sedentary+avg_minutes_lightly_active+
                            avg_minutes_fairly_active+avg_minues_very_active, data = activity.data.summary)
vif(lm.shannon.activity)

# lm.shannon.activity <- lm(avg_shannon ~ avg_steps+avg_distance+avg_minutes_sedentary+
#                             avg_minutes_fairly_active+avg_minues_very_active, data = activity.data.summary)
summary(lm.shannon.activity)

## Mixed Effects Model
library(lme4)
library(lmerTest)
library(performance)

fecal.microbial.menses.24.activity.filtered <- 
  fecal.microbial.menses.24.activity[complete.cases(fecal.microbial.menses.24.activity),]

lmer.null <- lmer(shannon~(1|`biome_id`), 
                data = fecal.microbial.menses.24.activity.filtered)
r2(lmer.null)

# full model - fixed slopes, rnd intercepts
lmer.full <- lmer(shannon~steps + distance +
                    calories_burned + minutes_sedentary + minutes_lightly_active +
                    minutes_fairly_active + minues_very_active + activity_calories +
                    (1|`biome_id`), 
                  data = fecal.microbial.menses.24.activity.filtered)
r2(lmer.full)
summary(lmer.full)

anova(lmer.null, lmer.full)

# full model - random slopes, rnd intercepts
rnd.slope.lmer.full <- lmer(shannon ~ steps + (steps||`biome_id`) + 
                              distance + (distance||`biome_id`) +
                              calories_burned + (calories_burned||`biome_id`) +
                              minutes_sedentary + (minutes_sedentary||`biome_id`) +
                              minutes_lightly_active + (minutes_lightly_active||`biome_id`) +
                              minutes_fairly_active + (minutes_fairly_active||`biome_id`) +
                              minues_very_active + (minues_very_active||`biome_id`) +
                              activity_calories + (activity_calories||`biome_id`),
                            data = fecal.microbial.menses.24.activity)
r2(rnd.slope.lmer.full)

anova(lmer.full, rnd.slope.lmer.full)

## Time
activity.time.model <- lmer(shannon ~ study_day + steps +
                              distance + calories_burned + minutes_sedentary +
                              minutes_lightly_active + minutes_fairly_active +
                              minues_very_active + activity_calories +
                              + (1 | biome_id),
                            data = fecal.microbial.menses.24.activity)
r2(activity.time.model)
summary(activity.time.model)
anova(lmer.full, activity.time.model)

## Time + Time^2
activity.time2.model <- lmer(shannon ~ study_day + I(study_day^2) + steps +
                              distance + calories_burned + minutes_sedentary +
                              minutes_lightly_active + minutes_fairly_active +
                              minues_very_active + activity_calories +
                              + (1 | biome_id),
                            data = fecal.microbial.menses.24.activity)
r2(activity.time2.model)
summary(activity.time2.model)
anova(activity.time.model, activity.time2.model)

## high and low
# table(fecal.microbial.menses.24.activity$minues_very_active)
# names(fecal.microbial.menses.24.activity)

## check if student athlete matters
participant.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 9-Volunteer Medical History.csv", header = TRUE)

participant.data <- participant.data %>% 
  mutate(sport_varsity = ifelse((sport=="In-Season" | sport == "Off-Season"), "Varsity", "Not Varsity"),
         sport_collapsed = ifelse((sport=="In-Season"), sport, "Not in-Season")) %>%
  # mutate(sport_collapsed = ifelse(sport=="In-Season" | sport == "Club", "Activly in Sport", "Not In-Season")) %>%
  select(biome_id, logDate, activity_level, sport, sport_collapsed, sport_varsity, field_hockey)

activity.sport.summary <- activity.data.summary %>% 
  left_join(participant.data, by="biome_id") %>% 
  filter(!is.na(sport))

length(unique(activity.sport.summary$biome_id))
table(activity.sport.summary$sport)

# Gut: boxplot of gut shannon diversity by sports type
ggplot(activity.sport.summary, aes(x = sport, y = avg_shannon)) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha=0.1) +
  geom_jitter(width = 0.2, alpha = 0.6, color="orchid") +
  labs(
    x = " ", 
    y = "Average Shannon Diversity", 
    title = ""
  ) +
  theme_minimal() +
  theme(text=element_text(size=18), legend.position = "None")

# testing
aov_sport <- aov(avg_shannon ~ sport, data = activity.sport.summary)
summary(aov_sport)

## Sports collapsed

length(unique(activity.sport.summary$biome_id))
table(activity.sport.summary$sport_collapsed)

# Gut: boxplot of shannon diversity by collapsed sports
ggplot(activity.sport.summary, aes(x = sport_collapsed, y = avg_shannon)) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha=0.1) +
  geom_jitter(width = 0.2, alpha = 0.6, color="orchid") +
  labs(
    x = " ", 
    y = "Average Shannon Diversity", 
    title = ""
  ) +
  theme_minimal() +
  theme(text=element_text(size=18), legend.position = "None")

t.test(avg_shannon ~ sport_collapsed, data=activity.sport.summary)
wilcox.test(avg_shannon ~ sport_collapsed, data = activity.sport.summary)

# Gut: boxplot of shannon diversity by varsity sports
ggplot(activity.sport.summary, aes(x = sport_varsity, y = avg_shannon)) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha=0.1) +
  geom_jitter(width = 0.2, alpha = 0.6, color="orchid") +
  labs(
    x = " ", 
    y = "Average Shannon Diversity", 
    title = ""
  ) +
  theme_minimal() +
  theme(text=element_text(size=18), legend.position = "None")

# testing
aov_sport <- aov(avg_shannon ~ sport_varsity, data = activity.sport.summary)
summary(aov_sport)

## Levels of activity
activity.data.summary.levels <- activity.data.summary %>% 
  mutate(exercise_level = case_when(
    total_min_active <= quantile(total_min_active, 0.25, na.rm = TRUE) ~ "Low",
    total_min_active <= quantile(total_min_active, 0.50, na.rm = TRUE) ~ "Moderate",
    total_min_active <= quantile(total_min_active, 0.75, na.rm = TRUE) ~ "High",
    TRUE ~ "Very High"
  ))

activity.data.summary.levels <- activity.data.summary.levels %>% 
  left_join(participant.data, by="biome_id") 

activity.data.summary.levels <- activity.data.summary.levels %>% 
  filter(!is.na(sport))

activity.data.summary.levels$exercise_level <- factor(activity.data.summary.levels$exercise_level, levels = c("Low", "Moderate", "High", "Very High"), ordered = TRUE)

table(activity.data.summary.levels$exercise_level)
table(activity.data.summary.levels$exercise_level, activity.data.summary.levels$sport)

lm.obj <- lm(avg_shannon ~ exercise_level, data = activity.data.summary.levels)
summary(lm.obj)


## Change activity graphic (in data streams) and add one for log date on the x-axis
# pick one activity variable (maybe high activity) show 67 ppl across time
#see if they are different over time
activity.data.levels <- fecal.microbial.menses.24.activity %>% 
  group_by(biome_id, logDate) %>% 
  mutate(total_min_active = sum(minutes_lightly_active, minutes_fairly_active, minues_very_active, na.rm=TRUE),
    exercise_level = case_when(
    total_min_active <= quantile(total_min_active, 0.25, na.rm = TRUE) ~ "Low",
    total_min_active <= quantile(total_min_active, 0.50, na.rm = TRUE) ~ "Moderate",
    total_min_active <= quantile(total_min_active, 0.75, na.rm = TRUE) ~ "High",
    TRUE ~ "Very High"
  ))

activity.data.levels <- activity.data.levels %>% 
  left_join(participant.data) 

activity.data.levels <- activity.data.levels %>% 
  filter(!is.na(sport))

activity.data.levels$exercise_level <- factor(activity.data.levels$exercise_level, levels = c("Low", "Moderate", "High", "Very High"), ordered = TRUE)
names(activity.data.levels)
length(unique(activity.data.levels$biome_id))

# add study days
activity.data.levels <- study_days(as.data.frame(activity.data.levels))

# Activity: scatterplot and splines of activity over time
activity.data.levels %>% 
  filter(steps != 0) %>% 
  ggplot(aes(x = study_day, y = steps, group=biome_id, color=as.factor(biome_id))) +
  geom_point(aes(color=as.factor(biome_id)), alpha=0.2, show.legend=FALSE) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5, show.legend=FALSE) +
  # geom_smooth(method = "lm", se=FALSE, aes(color = as.factor(biome_id)), size=0.5, show.legend=FALSE) +
  labs(
    x = "Days", 
    y = "Steps",
    title = ""
  ) +
  theme_minimal() +
  # scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "1 day") +
  # scale_color_discrete(labels = c("menses" = "Menstruating", "not_menses" = "Not menstruating")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),
        text=element_text(size=18)) +
  scale_x_continuous(breaks = seq(0, max(activity.data.levels$study_day), by = 5))

# Activity: scatterplot and splines of high levels activity over time
activity.data.levels %>% 
  filter(minues_very_active != 0) %>% 
  ggplot(aes(x = study_day, y = minues_very_active, group=biome_id, color=as.factor(biome_id))) +
  geom_point(aes(color=as.factor(biome_id)), alpha=0.2, show.legend=FALSE) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5, show.legend=FALSE) +
  # geom_smooth(method = "lm", se=FALSE, aes(color = as.factor(biome_id)), size=0.5, show.legend=FALSE) +
  labs(
    x = "Days", 
    y = "Minutes of High Activity",
    title = ""
  ) +
  theme_minimal() +
  scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "1 day") +
  # scale_color_discrete(labels = c("menses" = "Menstruating", "not_menses" = "Not menstruating")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),
        text=element_text(size=18)) +
  scale_x_continuous(breaks = seq(0, max(activity.data.levels$study_day), by = 5))

# Activity: scatterplot and splines of high levels activity over time
activity.data.levels %>% 
  filter(minues_very_active != 0) %>% 
  ggplot(aes(x = study_day, y = minues_very_active, group=biome_id, color=as.factor(biome_id))) +
  geom_point(aes(color=as.factor(biome_id)), alpha=0.2, show.legend=FALSE) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5, show.legend=FALSE) +
  # geom_smooth(method = "lm", se=FALSE, aes(color = as.factor(biome_id)), size=0.5, show.legend=FALSE) +
  labs(
    x = "Days", 
    y = "Minutes of High Activity",
    title = ""
  ) +
  theme_minimal() +
  scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "1 day") +
  # scale_color_discrete(labels = c("menses" = "Menstruating", "not_menses" = "Not menstruating")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),
        text=element_text(size=18)) +
  scale_x_continuous(breaks = seq(0, max(activity.data.levels$study_day), by = 5))

# Activity: scatterplot and splines of sedentary levels activity over time
activity.data.levels %>% 
  filter(minutes_sedentary != 0) %>% 
  ggplot(aes(x = study_day, y = minutes_sedentary, group=biome_id, color=as.factor(biome_id))) +
  geom_point(aes(color=as.factor(biome_id)), alpha=0.2, show.legend=FALSE) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5, show.legend=FALSE) +
  # geom_smooth(method = "lm", se=FALSE, aes(color = as.factor(biome_id)), size=0.5, show.legend=FALSE) +
  labs(
    x = "Days", 
    y = "Minutes of Sedentary Activity",
    title = ""
  ) +
  theme_minimal() +
  scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "1 day") +
  # scale_color_discrete(labels = c("menses" = "Menstruating", "not_menses" = "Not menstruating")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),
        text=element_text(size=18)) +
  scale_x_continuous(breaks = seq(0, max(activity.data.levels$study_day), by = 5))

# Activity: scatterplot and splines of light levels activity over time
activity.data.levels %>% 
  filter(minutes_lightly_active != 0) %>% 
  ggplot(aes(x = study_day, y = minutes_lightly_active, group=biome_id, color=as.factor(biome_id))) +
  geom_point(aes(color=as.factor(biome_id)), alpha=0.2, show.legend=FALSE) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5, show.legend=FALSE) +
  # geom_smooth(method = "lm", se=FALSE, aes(color = as.factor(biome_id)), size=0.5, show.legend=FALSE) +
  labs(
    x = "Days", 
    y = "Minutes of Light Activity",
    title = ""
  ) +
  theme_minimal() +
  scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "1 day") +
  # scale_color_discrete(labels = c("menses" = "Menstruating", "not_menses" = "Not menstruating")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),
        text=element_text(size=18)) +
  scale_x_continuous(breaks = seq(0, max(activity.data.levels$study_day), by = 5))

# Activity: scatterplot and splines of moderate levels activity over time
activity.data.levels %>% 
  filter(minutes_fairly_active != 0) %>% 
  ggplot(aes(x = study_day, y = minutes_fairly_active, group=biome_id, color=as.factor(biome_id))) +
  geom_point(aes(color=as.factor(biome_id)), alpha=0.2, show.legend=FALSE) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5, show.legend=FALSE) +
  # geom_smooth(method = "lm", se=FALSE, aes(color = as.factor(biome_id)), size=0.5, show.legend=FALSE) +
  labs(
    x = "Days", 
    y = "Minutes of Moderate Activity",
    title = ""
  ) +
  theme_minimal() +
  scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "1 day") +
  # scale_color_discrete(labels = c("menses" = "Menstruating", "not_menses" = "Not menstruating")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),
        text=element_text(size=18)) +
  scale_x_continuous(breaks = seq(0, max(activity.data.levels$study_day), by = 5))

# Gut: scatterplot and splines of moderate levels activity over time
activity.data.levels %>% 
  filter(minutes_fairly_active != 0) %>% 
  ggplot(aes(x = as.Date(logDate), y = minutes_fairly_active, group=biome_id, color=as.factor(biome_id))) +
  geom_point(aes(color=as.factor(biome_id)), alpha=0.2, show.legend=FALSE) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5, show.legend=FALSE) +
  # geom_smooth(method = "lm", se=FALSE, aes(color = as.factor(biome_id)), size=0.5, show.legend=FALSE) +
  labs(
    x = "Days", 
    y = "Minutes of High Activity",
    title = ""
  ) +
  theme_minimal() +
  scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "1 day") +
  # scale_color_discrete(labels = c("menses" = "Menstruating", "not_menses" = "Not menstruating")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
        text=element_text(size=18))

## Mixed effects model
## Levels of activity
fecal.microbial.menses.24.activity.levels <- activity.data.levels %>% 
  group_by(biome_id, logDate) %>% 
  # mutate(total_min_active = sum(minutes_lightly_active, minutes_fairly_active*2, minues_very_active*3, na.rm=TRUE)) %>% 
  ungroup() %>% 
  mutate(exercise_level = case_when(
    minues_very_active <= quantile(minues_very_active, 0.25, na.rm = TRUE) ~ "Low",
    minues_very_active <= quantile(minues_very_active, 0.50, na.rm = TRUE) ~ "Moderate",
    minues_very_active <= quantile(minues_very_active, 0.75, na.rm = TRUE) ~ "High",
    TRUE ~ "Very High"
  ))
table(fecal.microbial.menses.24.activity.levels$exercise_level)

# Self classification of participants against actual fitbit data
table(participant.data$activity_level)
participant.data.subset <- participant.data %>% 
  select(biome_id, activity_level)
fecal.microbial.menses.24.activity.levels <- fecal.microbial.menses.24.activity.levels %>% 
  left_join(participant.data.subset)
# relevel exercise level (fitbit)
fecal.microbial.menses.24.activity.levels$exercise_level <- 
  factor(fecal.microbial.menses.24.activity.levels$exercise_level,
         levels=c("Low", "Moderate", "High", "Very High"))

table(fecal.microbial.menses.24.activity.levels$activity_level, fecal.microbial.menses.24.activity.levels$exercise_level)

# maybe make it for the most common one for participant then compare?

lmer.obj <- lmer(shannon ~ exercise_level + (1|`biome_id`), data = fecal.microbial.menses.24.activity.levels)
r2(lmer.obj)
summary(lmer.obj)

################################################################################ 

## Merge dfs

## Join data sets
# fecal.microbial.menses.24.activity.levels.subset <- fecal.microbial.menses.24.activity.levels %>% 
#   left_join(participant.data)
# fecal.microbial.menses.24.activity.levels.subset <- fecal.microbial.menses.24.activity.levels.subset %>% 
#   select(!c(Date, study_day))
# 
# gut.lifestyle.data <- fecal.microbial.menses.24.activity.levels.subset %>% 
#   left_join(fecal.microbial.menses.24)
# 
# gut.lifestyle.data <- gut.lifestyle.data %>% 
#   left_join(gut.diet.df)

colSums(is.na(fecal.microbial.menses.24.filtered))
dim(fecal.microbial.menses.24.filtered)
fecal.microbial.menses.24.filtered.subset <- fecal.microbial.menses.24.filtered %>% 
  select(c(shannon, biome_id, logDate, OTU, SampleID)) %>% 
  mutate(logDate = as.Date(logDate)) %>% 
  group_by(biome_id, logDate) %>% 
  summarise(shannon = first(shannon), # get first sample for day
      OTU = first(OTU),
      SampleID = first(SampleID))
dim(fecal.microbial.menses.24.filtered.subset)

participant.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 9-Volunteer Medical History.csv", header = TRUE)
participant.data.subset <- participant.data %>% 
  select(biome_id, cisWoman, sport, probiotic, birthControl, study_menstruate,
         sexuallyActive)

gut.lifestyle.data <- fecal.microbial.menses.24.filtered.subset %>% 
  left_join(participant.data.subset)
colSums(is.na(gut.lifestyle.data))

nutrient.diet.subset <- nutrient.diet.22 %>% 
  select(c(biome_id, logDate, caloriesall, cholesterol_prop, satFat_prop, sodium_prop,
           carb_prop, dietFib_prop, sugar_prop, protein_prop, fat_prop, fat_cal_prop,
           addedSugarall_prop)) %>% 
  mutate(logDate=as.Date(logDate)) %>% 
  group_by(biome_id, logDate) %>% 
  summarise(caloriesall_avg = sum(caloriesall, na.rm=TRUE) / n(),
            cholesterol_prop = first(cholesterol_prop),
            satFat_prop = first(satFat_prop),
            sodium_prop = first(sodium_prop),
            carb_prop = first(carb_prop),
            dietFib_prop = first(dietFib_prop),
            sugar_prop = first(sugar_prop),
            protein_prop = first(protein_prop),
            fat_prop = first(fat_prop),
            fat_cal_prop = first(fat_cal_prop),
            addedSugarall_prop = first(addedSugarall_prop)
  )

gut.lifestyle.data <- gut.lifestyle.data %>% 
  left_join(nutrient.diet.subset)
colSums(is.na(gut.lifestyle.data))
dim(gut.lifestyle.data)

shannon.dass.filtered2.subset <- shannon.dass.filtered2 %>% 
  select(shannon, biome_id, logDate, menses_day, stress_score) %>% 
  mutate(logDate=as.Date(logDate)) %>% 
  group_by(biome_id, logDate) %>% 
  summarise(menses_day = first(menses_day),
            stress_score = first(stress_score))

gut.lifestyle.data <- gut.lifestyle.data %>% 
  left_join(shannon.dass.filtered2.subset)

activity.data.levels.subset <- activity.data.levels %>% 
  select(c(shannon, biome_id, logDate, calories_burned, steps,
           distance, floors, minutes_sedentary, minutes_lightly_active, minutes_fairly_active,
           minues_very_active, activity_calories)) %>% 
  mutate(logDate=as.Date(logDate)) %>% 
  group_by(biome_id, logDate) %>% 
  summarise(calories_burned = first(calories_burned),
            steps = first(steps),
            distance = first(distance),
            minutes_sedentary = first(minutes_sedentary),
            minutes_lightly_active = first(minutes_lightly_active),
            minutes_fairly_active = first(minutes_fairly_active),
            minues_very_active = first(minues_very_active),
            activity_calories = first(activity_calories))

gut.lifestyle.data <- gut.lifestyle.data %>% 
  left_join(activity.data.levels.subset)

colSums(is.na(gut.lifestyle.data))

gut.lifestyle.data <- study_days(gut.lifestyle.data)

# 04/13 uncomment to save - full gut df
# write.csv(gut.lifestyle.data, "/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/relabeled_data/gut.lifestyle.merged.csv")

colSums(is.na(gut.lifestyle.data))

lmer.obj <- lmer(shannon~.+(1|`biome_id`), data=gut.lifestyle.data[,-c(2,4,5)])
r2(lmer.obj)
summary(lmer.obj)

#### aggregate

# Aggregate model
gut.merged.df.summary <- gut.lifestyle.data %>% 
  group_by(biome_id) %>%
  summarise(avg_shannon=sum(shannon, na.rm=TRUE) / n(),
            avg_stress_score=sum(stress_score, na.rm=TRUE) / n(),
            cisWoman = first(cisWoman),
            sport = first(sport),
            probiotic = first(probiotic),
            birthControl = first(birthControl),
            sexuallyActive = first(sexuallyActive),
            caloriesall_avg = sum(caloriesall_avg, na.rm=TRUE) / n(),
            cholesterol_prop = sum(cholesterol_prop, na.rm=TRUE) / n(),
            satFat_prop = sum(satFat_prop, na.rm=TRUE) / n(),
            sodium_prop = sum(sodium_prop, na.rm=TRUE) / n(),
            carb_prop = sum(carb_prop, na.rm=TRUE) / n(),
            dietFib_prop = sum(dietFib_prop, na.rm=TRUE) / n(),
            sugar_prop = sum(sugar_prop, na.rm=TRUE) / n(),
            protein_prop = sum(protein_prop, na.rm=TRUE) / n(),
            fat_prop = sum(fat_prop, na.rm=TRUE) / n(),
            fat_cal_prop = sum(fat_cal_prop, na.rm=TRUE) / n(),
            addedSugarall_prop = sum(addedSugarall_prop, na.rm=TRUE) / n(),
            study_menstruate = first(study_menstruate),
            calories_burned = sum(calories_burned, na.rm=TRUE) / n(),
            steps = sum(steps, na.rm=TRUE) / n(),
            distance = sum(distance, na.rm=TRUE) / n(),
            minutes_sedentary = sum(minutes_sedentary, na.rm=TRUE) / n(),
            minutes_lightly_active = sum(minutes_lightly_active, na.rm=TRUE) / n(),
            minutes_fairly_active = sum(minutes_fairly_active, na.rm=TRUE) / n(),
            minues_very_active = sum(minues_very_active, na.rm=TRUE) / n(),
            activity_calories = sum(activity_calories, na.rm=TRUE) / n()) %>% 
  mutate(stress_severity = case_when(
    avg_stress_score <= 14 ~ "Normal",
    avg_stress_score >= 15 & avg_stress_score <= 18 ~ "Mild",
    avg_stress_score >= 19 & avg_stress_score <= 25 ~ "Moderate",
    avg_stress_score >= 26 & avg_stress_score <= 33 ~ "Severe",
    avg_stress_score >= 34 ~ "Extremely Severe",
    TRUE ~ NA_character_
  )) %>% 
  select(-avg_stress_score) 

colSums(is.na(gut.merged.df.summary))

gut.merged.df.summary.filtered <- 
  gut.merged.df.summary[complete.cases(gut.merged.df.summary),]


lm.obj.orig <- lm(avg_shannon~., data=gut.merged.df.summary.filtered[,-1])
summary(lm.obj.orig)

# menses
lm.obj <- lm(avg_shannon~., data=gut.merged.df.summary.filtered[,-c(1, 19)])
summary(lm.obj)
# 0.1359 
anova(lm.obj.orig, lm.obj)

# stress
lm.obj <- lm(avg_shannon~., data=gut.merged.df.summary.filtered[,-c(1, 28)])
summary(lm.obj)
# 0.1359 
anova(lm.obj.orig, lm.obj)

# HBC
lm.obj <- lm(avg_shannon~., data=gut.merged.df.summary.filtered[,-c(1, 6)])
summary(lm.obj)

anova(lm.obj.orig, lm.obj)

# nutrition
lm.obj <- lm(avg_shannon~., data=gut.merged.df.summary.filtered[,-c(1, 8:18)])
summary(lm.obj)
# 0.1359 

anova(lm.obj.orig, lm.obj)

# activity
lm.obj <- lm(avg_shannon ~ ., 
             data=gut.merged.df.summary.filtered[,-c(1, 20:27)])
summary(lm.obj)
# 0.1359

anova(lm.obj.orig, lm.obj)

# sport
lm.obj <- lm(avg_shannon~., data=gut.merged.df.summary.filtered[,-c(1, 4)])
summary(lm.obj)
# 0.1359
anova(lm.obj.orig, lm.obj)

# probiotic
lm.obj <- lm(avg_shannon~., data=gut.merged.df.summary.filtered[,-c(1, 5)])
summary(lm.obj)
# 0.1359
anova(lm.obj.orig, lm.obj)

# sex act
lm.obj <- lm(avg_shannon~., data=gut.merged.df.summary.filtered[,-c(1, 7)])
summary(lm.obj)
# 0.1359
anova(lm.obj.orig, lm.obj)

# gender - cisWoman
lm.obj <- lm(avg_shannon~., data=gut.merged.df.summary.filtered[,-c(1, 3)])
summary(lm.obj)
# 0.1359
anova(lm.obj.orig, lm.obj)
