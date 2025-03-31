########################
#
# Gut Microbiome Analysis with Lifestyle Factors
# Last updated: 03/21/2025
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
  filter(BLAST_species %in% top10_species) %>%
  group_by(biome_id, BLAST_species) %>%
  summarise(Total_Abundance = sum(Abundance), .groups = "drop")

top10_relative_df <- top10_df %>% 
  group_by(biome_id) %>%
  mutate(Relative_Abundance = Total_Abundance / sum(Total_Abundance) * 100)

# Figure: Top 10 Species Across Study in Gut Microbiome per Participant
ggplot(top10_relative_df, aes(x = factor(biome_id), y = Relative_Abundance, fill = BLAST_species)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal() +
  scale_fill_viridis_d() +
  labs(title = "Top 10 Species per Participant", x = "Participant", y = "Relative Abundance")

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
  select(biome_id, logDate, shannon, qr, max_taxa, OTU, menses_status, menses_day)


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
t.test(dass.participant.avg$avg_shannon ~ dass.participant.avg$person_menstruates)

# Figure: Average Shannon Diversity for Participants on menses v. not
ggplot(
  menses.table.df.summary,
  aes(x = menses_day_not_menses, y = menses_day_menses,
      color=as.factor(biome_id))) +
  geom_point() +
  xlim(3,4)+
  ylim(3,4)+
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  theme_minimal() +
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

# Figure: scatterplot of gut microbiome shannon diversity and menses
ggplot(fecal.microbial.menses.24.filtered, aes(x = as.Date(logDate), y = shannon)) +
  geom_point(aes(color=as.factor(person_menstruates))) +
  geom_smooth(method = "lm", se=FALSE, aes(color = as.factor(person_menstruates))) +
  labs(
    x = "Days", 
    y = "Shannon Diversity Index",
    title = "",
    color="Menstruates v. Doesn't menstruate"
  ) +
  theme_minimal() +
  scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "1 day") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

### Birth Control

## Corr with volunteer history data
participant.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 9-Volunteer Medical History.csv", header = TRUE)

# filter non-hormonal & no samples
participant.data <- participant.data %>% 
  filter(birthControl!="Orilissa (Elagolix)")

# select birth control
birthControl.df <- participant.data %>% 
  select(biome_id, birthControl)
birthControl.collapsed <- birthControl.df %>% 
  count(birthControl, name="frequency")

# merge df for shannon index with birth control
shannon.birthControl <- menses.table.df %>% 
  left_join(birthControl.df, by="biome_id") %>% 
  filter(!is.na(birthControl))

# add days
shannon.birthControl <- study_days(shannon.birthControl)

# Figure: boxplot of birth control and shannon diversity
ggplot(shannon.birthControl, aes(x = birthControl, y = shannon)) +
  geom_jitter(aes(color=as.factor(biome_id)), size=1, alpha=0.6) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
  scale_color_viridis_d(option="D") +
  labs(x = "Contraceptive", y = "Shannon Diversity Index", title = "",
       color = "Biome ID") +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1))

# Mixed effects models
library(lme4)
library(lmerTest)
library(performance)

lmer.shannon.bc <- lmer(shannon~birthControl + (1|`biome_id`), data = shannon.birthControl)
r2(lmer.shannon.bc)
summary(lmer.shannon.bc)

lmer.shannon.bc.time <- lmer(shannon~study_day + birthControl + (1|`biome_id`), data = shannon.birthControl)
r2(lmer.shannon.bc.time)
summary(lmer.shannon.bc.time)

lmer.shannon.bc.time2 <- lmer(shannon~study_day + I(study_day^2) + birthControl + (1|`biome_id`), data = shannon.birthControl)
r2(lmer.shannon.bc.time2)
summary(lmer.shannon.bc.time2)

# Try birth control v. none
shannon.birthControl.binary <- shannon.birthControl %>% 
  mutate(bc_binary = ifelse(birthControl=="None", "None", "birtControl"))
lmer.shannon.bc.binary <- 
  lmer(shannon~bc_binary + (1|`biome_id`), data = shannon.birthControl.binary)
r2(lmer.shannon.bc.binary)
summary(lmer.shannon.bc.binary)

# Collapse birth control variable
shannon.birthControl.collapsed <- shannon.birthControl %>% 
  mutate(birthControl_collapsed=ifelse(birthControl=="Systemic Combined (E&P)" | (birthControl=="Systemic P only"), "Systemic", 
                                       birthControl)) %>% 
  mutate(birthControl=as.factor(birthControl))

shannon.birthControl.collapsed$birthControl_collapsed <- factor(shannon.birthControl.collapsed$birthControl_collapsed, levels=c("None", "Local P", "Systemic"))

## Make full df
shannon.birthControl.collapsed.subset <- shannon.birthControl.collapsed %>% 
  select(biome_id, logDate, birthControl, birthControl_collapsed) %>% 
  distinct()
fecal.microbial.menses.24.birthControl <- fecal.microbial.menses.24 %>% 
  left_join(shannon.birthControl.collapsed.subset, by=c("biome_id", "logDate"))

## Mixed effects models

lmer.shannon.bc.collapsed <- lmer(shannon~birthControl_collapsed + (1|`biome_id`), data = shannon.birthControl.collapsed)
r2(lmer.shannon.bc.collapsed)
summary(lmer.shannon.bc.collapsed)

# Figure: boxplot of collapsed systemic birth control and shannon diversity
ggplot(shannon.birthControl.collapsed, aes(x = birthControl_collapsed, y = shannon)) +
  geom_jitter(aes(color=as.factor(biome_id)), size=1, alpha=0.6) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
  scale_color_viridis_d(option="D") +
  labs(x = "Contraceptive", y = "Shannon Diversity Index", title = "",
       color = "Biome ID") +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1))

# Average Shannon diversity per person
shannon.birthControl.avg <- shannon.birthControl %>%
  group_by(biome_id) %>% 
  summarise(avg_shannon = sum(shannon) / n(),
            birthControl = first(birthControl))

# Figure: boxplot of birth control and avg shannon diversity
ggplot(shannon.birthControl.avg, aes(x = birthControl, y = avg_shannon)) +
  geom_jitter(aes(color=as.factor(biome_id)), size=1, alpha=0.6) +
  geom_point(aes(color=as.factor(biome_id)), size=1, alpha=0.7) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
  # geom_text(aes(label = as.factor(biome_id)), vjust = -0.5, hjust=1.5, size = 1.5) +
  scale_color_viridis_d(option="D") +
  labs(x = "Contraceptive", y = "Average Shannon Diversity Index", title = "",
       color = "Biome ID") +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1))

lm.avgshannon.bc <- lm(avg_shannon~birthControl, data = shannon.birthControl.avg)
summary(lm.avgshannon.bc)

# Average Shannon diversity per person
shannon.birthControl.collapsed.avg <- shannon.birthControl.collapsed %>%
  group_by(biome_id) %>% 
  summarise(avg_shannon = sum(shannon) / n(),
            birthControl_collapsed = first(birthControl_collapsed))

# Figure: boxplot of collapsed Systemic birth control and avg shannon diversity
ggplot(shannon.birthControl.collapsed.avg, aes(x = birthControl_collapsed, y = avg_shannon)) +
  geom_jitter(aes(color=as.factor(biome_id)), size=1, alpha=0.6) +
  geom_point(aes(color=as.factor(biome_id)), size=1, alpha=0.7) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
  # geom_text(aes(label = as.factor(biome_id)), vjust = -0.5, hjust=1.5, size = 1.5) +
  scale_color_viridis_d(option="D") +
  labs(x = "Contraceptive", y = "Average Shannon Diversity Index", title = "",
       color = "Biome ID") +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1))

lm.avgshannon.bc.collapsed <- lm(avg_shannon~birthControl_collapsed, data = shannon.birthControl.collapsed.avg)
summary(lm.avgshannon.bc.collapsed)

##########################################################################################

# Menses
# full_menses_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/relabeled_data/cleaned_menstruation_data.csv", header=TRUE)
# 
# # Q: how many days with fecal swabs do I have menses info for (from vaginal)
# table(shannon.birthControl$menses_day)

##########################################################################################

## Corr with DASS data/stress
dass <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_dass.csv")

## DASS Survey completion - aggregate by week

all_days <- seq.Date(as.Date(min(dass$Timestamp, na.rm=TRUE)), 
                     as.Date(max(dass$Timestamp, na.rm=TRUE)), by = "day")
# all_days <- as.character(all_days)
# complete_grid <- expand.grid(biome_id = unique(dass$biome_id),
#                              Timestamp = all_days)
complete_grid <- expand.grid(
  biome_id = unique(dass$biome_id),
  Timestamp = all_days
) %>%
  mutate(week = cut(Timestamp, breaks = "week", start.on.monday = TRUE))
# 
# heatmap_data <- complete_grid %>%
#   left_join(dass, by = c("biome_id", "Timestamp"))
heatmap_data <- complete_grid %>%
  left_join(dass %>% mutate(week = cut(Timestamp, breaks = "week")), 
            by = c("biome_id", "week"))

heatmap_data$Value <- ifelse(is.na(heatmap_data$stress_score),
                             "No Survey", "Survey")

heatmap_data <- heatmap_data %>% 
  group_by(biome_id) %>%
  mutate(survey_count = sum(Value == "Survey", na.rm = TRUE)) %>%
  arrange(desc(survey_count), biome_id) %>% 
  ungroup()

ggplot(heatmap_data, aes(x = as.character(week), y = reorder(factor(biome_id), survey_count), fill = Value)) +
  geom_tile(color = "gray25") +
  scale_fill_manual(values = c("white", "darkblue"), na.value = "white", 
                    labels = c("No survey", "Survey")) +
  labs(title = "DASS Completion", 
       x = "Week", 
       y = "Biome ID", 
       fill = "Submission") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

### Imputed and cleaned DASS data
dass <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/DASS_0503_2024-final_df.csv")
dass <- dass %>% 
  rename(biome_id=study_id)

## Aggregate

# Average stress score
dass.avg <- dass %>% 
  group_by(biome_id) %>% 
  summarise(
    # avg_depr=sum(depression_score, na.rm=TRUE)/n(),
    avg_anx=sum(anxiety_score, na.rm=TRUE)/n(),
    avg_stress=sum(stress_score, na.rm=TRUE)/n()
  )

# average shannon diversity
shannon.birthControl.avg <- shannon.birthControl %>% 
  group_by(biome_id) %>% 
  summarise(avg_shannon = sum(shannon) / n(),
            birthControl = first(birthControl))
dass.participant.avg <- merge(shannon.birthControl.avg, dass.avg, by="biome_id")

spline.obj <- smooth.spline(x=dass.participant.avg$avg_stress, y=dass.participant.avg$avg_shannon)
spline.df <- data.frame(avg_stress=spline.obj$x,
                        avg_shannon=spline.obj$y)

# Figure: scatter plot of avg shannon idx by avg stress score, color by birth control
ggplot(dass.participant.avg, aes(x = avg_stress, y = avg_shannon)) +
  geom_point(aes(color=as.factor(birthControl)), size=1, alpha=1) +
  geom_line(data=spline.df, aes(x=avg_stress, y=avg_shannon), color="orchid", linewidth=1, alpha=0.7) +
  scale_color_viridis_d(option="D") +
  labs(x = "Average Stress Score", y = "Average Shannon Index", title = "",
       color = "Birth Control") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

# Average stress v. avg shannon diversity index, birth control
lm.results <- dass.participant.avg %>%
  group_by(birthControl) %>%
  summarise(
    model = list(lm(avg_shannon ~ avg_stress, data = cur_data())),
    .groups = "drop"
  ) %>%
  mutate(
    slope = sapply(model, function(m) coef(m)["avg_stress"]),
    p_value = sapply(model, function(m) summary(m)$coefficients["avg_stress", "Pr(>|t|)"])
  ) %>%
  select(birthControl, slope, p_value)

spline.obj <- dass.participant.avg %>%
  group_by(birthControl) %>%
  filter(n_distinct(avg_stress) >= 4) %>%
  summarise(
    model = list(smooth.spline(x=avg_stress, y=avg_shannon, spar=1)),
    .groups = "drop"
  ) %>%
  mutate(
    x = map(model, ~ seq(min(.x$x), max(.x$x), length.out = 100)),
    y = map2(model, x, ~ predict(.x, .y)$y)
  ) %>%
  select(birthControl, x, y) %>%
  unnest(c(x, y))

# Figure: Scatterplot of avg shannon by avg stress with splines
ggplot(dass.participant.avg, aes(x = avg_stress, y = avg_shannon)) +
  geom_point(aes(color=as.factor(birthControl)), size=1, alpha=1) +
  geom_line(data = spline.obj, aes(x = x, y = y, color = birthControl), linewidth = 1) +
  scale_color_viridis_d(option="D") +
  labs(x = "Average Stress Score", y = "Average Shannon Diversity Index", title = "",
       color = "Birth Control") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

### Alpha diversity over time, color maybe by birth control, menstruate, or both

# merge df
dass.participant <- merge(shannon.birthControl, dass.avg, by="biome_id")
dass.participant$logDate <- as.Date(dass.participant$logDate)

# Figure: scatterplot shannon over days by birth control
ggplot(dass.participant, aes(x = logDate, y = shannon)) +
  geom_point(aes(color=birthControl)) +
  geom_smooth(method = "lm", se=FALSE, aes(color = birthControl)) +
  labs(
    x = "Days", 
    y = "Shannon Diversity Index",
    title = ""
  ) +
  theme_minimal() +
  scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "1 day") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
  scale_color_viridis_d(option="D")

# Figure: scatterplot shannon over days by menses day
ggplot(dass.participant, aes(x = logDate, y = shannon)) +
  geom_point(aes(color=menses_day)) +
  geom_smooth(method = "lm", se=FALSE, aes(color = menses_day)) +
  labs(
    x = "Days", 
    y = "Shannon Diversity Index",
    title = ""
  ) +
  theme_minimal() +
  scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "1 day") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
  scale_color_viridis_d(option="D")

fecal.microbial.menses.24.menstruates <- fecal.microbial.menses.24.filtered %>% 
  select(biome_id, person_menstruates) %>% 
  distinct()
dass.participant <- dass.participant %>% 
  left_join(fecal.microbial.menses.24.menstruates, by="biome_id")

dass.participant.avg <- dass.participant %>% 
  group_by(biome_id) %>% 
  summarise(avg_shannon = sum(shannon) / n(),
            person_menstruates = first(person_menstruates))

# Figure: Scatterplot of avg shannon over time by if someone menstruates
ggplot(dass.participant, aes(x = logDate, y = shannon)) +
  geom_point(aes(color=person_menstruates)) +
  geom_smooth(method = "lm", se=FALSE, aes(color = person_menstruates)) +
  labs(
    x = "Days", 
    y = "Shannon Diversity Index",
    title = ""
  ) +
  theme_minimal() +
  scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "1 day") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
  scale_color_viridis_d(option="D")

## Regress shannon diversity on stress + stress*birth control, if interaction terms (2) is significant

# collapse birth control
dass.participant.collapsed <- merge(shannon.birthControl.collapsed, dass.avg, by="biome_id")
dass.participant.collapsed$birthControl_collapsed <- factor(
  dass.participant.collapsed$birthControl_collapsed,
  levels = c("None", "Local P", "Systemic")
)

dass.participant.collapsed.avg <- dass.participant.collapsed %>% 
  group_by(biome_id) %>% 
  summarise(avg_shannon = sum(shannon)/n(),
            birthControl = first(birthControl),
            birthControl_collapsed = first(birthControl_collapsed),
            # avg_depr = first(avg_depr),
            avg_anx = first(avg_anx),
            avg_stress = first(avg_stress))

# Regressing models
lm.obj <- lm(avg_shannon ~ avg_stress*birthControl, data=dass.participant.collapsed.avg)
summary(lm.obj)
lm.obj2 <- lm(avg_shannon ~ avg_stress+birthControl, data=dass.participant.collapsed.avg)
summary(lm.obj2)
anova(lm.obj,lm.obj2)

# Collapsed birth control
lm.obj <- lm(avg_shannon ~ avg_stress*birthControl_collapsed, data=dass.participant.collapsed.avg)
summary(lm.obj)
lm.obj2 <- lm(avg_shannon ~ avg_stress+birthControl_collapsed, data=dass.participant.collapsed.avg)
summary(lm.obj2)
anova(lm.obj,lm.obj2)


########################################################

## Individual DASS surveys

shannon.birthControl$logDate <- as.Date(shannon.birthControl$logDate)
dass$Timestamp <- as.Date(dass$Timestamp)
dass.filtered <- dass %>% 
  select(!birthControl)

# matching: select closest survey
shannon.birthControl.dass <- shannon.birthControl %>%
  left_join(dass, by = "biome_id") %>%
  mutate(DayDifference = abs(as.numeric(logDate - Timestamp))) %>%
  group_by(biome_id, logDate) %>% 
  slice_min(DayDifference, with_ties = FALSE) %>% 
  ungroup() %>% 
  select("biome_id", "logDate", "Timestamp", "DayDifference", everything())

# Mixed effects models
library(lme4)
library(lmerTest)
library(performance)

# mixed effects with stress
lmer.stress <-lmer(shannon~stress_score + 
                     (1|`biome_id`), 
                   data = shannon.birthControl.dass)
r2(lmer.stress)
summary(lmer.stress)

shannon.birthControl.dass <- study_days(shannon.birthControl.dass)

# mixed effects with stress and time
lmer.stress.days <-lmer(shannon~stress_score + study_day +
                     (1|`biome_id`), 
                   data = shannon.birthControl.dass)
r2(lmer.stress.days)
summary(lmer.stress.days)

anova(lmer.stress, lmer.stress.days)

# mixed effects with stress and time, time^2
lmer.stress.days2 <-lmer(shannon~stress_score + study_day +
                          I(study_day^2) +
                          (1|`biome_id`), 
                        data = shannon.birthControl.dass)
r2(lmer.stress.days2)
summary(lmer.stress.days2)

anova(lmer.stress.days, lmer.stress.days2)
anova(lmer.stress, lmer.stress.days2)

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
merged_diet_data <- study_days(merged_diet_data)
past_two_days_diet_data <- study_days(past_two_days_diet_data)
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
    is_vegetarian = first(is_vegetarian),
    study_day = first(study_day)
  )

# add veg cols to the rolling window df
merged_diet_data_collapsed_subset <- merged_diet_data_collapsed %>% 
  select(biome_id, logDate, perc_veg, is_vegetarian)
past_two_days_diet_data <- past_two_days_diet_data %>% 
  left_join(merged_diet_data_collapsed_subset, by=c("biome_id", "logDate"))

# join fecal data with diet data
gut.diet.df <- merged_diet_data_collapsed %>% # past_two_days_diet_data %>% #
  left_join(fecal.microbial.menses.24, by=c("biome_id", "logDate")) %>% 
  filter(!is.na(shannon))

# Figure: Macronutrients Plot with Shannon Diversity Index
par(mfrow = c(4,3))

plot(col = gut.diet.df$biome_id, y = gut.diet.df$shannon , x = gut.diet.df$caloriesall)
plot(col = gut.diet.df$biome_id, y = gut.diet.df$shannon , x = gut.diet.df$cholesterolall)
plot(col = gut.diet.df$biome_id, y = gut.diet.df$shannon , x = gut.diet.df$saturatedFatall)
plot(col = gut.diet.df$biome_id, y = gut.diet.df$shannon , x = gut.diet.df$sodiumall)
plot(col = gut.diet.df$biome_id, y = gut.diet.df$shannon , x = gut.diet.df$carbohydratesall)
plot(col = gut.diet.df$biome_id, y = gut.diet.df$shannon , x = log(gut.diet.df$dietaryFiberall))
plot(col = gut.diet.df$biome_id, y = gut.diet.df$shannon , x = gut.diet.df$sugarsall)
plot(col = gut.diet.df$biome_id, y = gut.diet.df$shannon , x = gut.diet.df$proteinall)
plot(col = gut.diet.df$biome_id, y = gut.diet.df$shannon , x = gut.diet.df$fatall)
plot(col = gut.diet.df$biome_id, y = gut.diet.df$shannon , x = gut.diet.df$saturatedFat)
plot(col = gut.diet.df$biome_id, y = gut.diet.df$shannon , x = log(gut.diet.df$transFat))
plot(col = gut.diet.df$biome_id, y = gut.diet.df$shannon , x = gut.diet.df$addedSugarall)

par(mfrow = c(1,1))

## vegetarian analysis
gut.diet.summary.df <- gut.diet.df %>% 
  group_by(biome_id) %>% 
  summarise(avg_shannon = sum(shannon, na.rm = TRUE) / n(),
            is_vegetarian = first(is_vegetarian),
            perc_veg = first(perc_veg))
table(gut.diet.summary.df$is_vegetarian)
summary(gut.diet.summary.df$avg_shannon)

# Figure: boxplot of average shannon diversity and vegetarian
ggplot(gut.diet.summary.df, aes(x = is_vegetarian, y = avg_shannon)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color="orchid") +
  labs(
    x = " ", 
    y = "Average Shannon Diversity", 
    title = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none") 

ggplot(gut.diet.df, aes(x = is_vegetarian, y = shannon)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color="orchid") +
  labs(
    x = " ", 
    y = "Shannon Diversity", 
    title = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none") 

# testing
t.test(avg_shannon ~ is_vegetarian, data=gut.diet.summary.df)
wilcox.test(avg_shannon ~ is_vegetarian, data = gut.diet.summary.df)

### PERCENTAGE VEGETARIAN

# average shannon and percent veg
veg_perc_df_summary  <- gut.diet.df %>% 
  group_by(biome_id) %>% 
  summarise(avg_shannon = sum(shannon, na.rm=TRUE) / n(),
            perc_veg = first(perc_veg))
summary(veg_perc_df_summary$avg_shannon)

lm.obj.summary <- lm(avg_shannon ~ perc_veg, data = veg_perc_df_summary)
summary(lm.obj.summary)

# lm.obj <- lm(shannon ~ perc_veg, data = gut.diet.df)
# summary(lm.obj)

# Join max OTU
fecal.microbial.menses.24.OTU <- fecal.microbial.menses.24 %>% 
  select(biome_id, logDate, OTU)

veg_perc_df_summary2 <- veg_perc_df_summary %>%
  left_join(fecal.microbial.menses.24.OTU, by="biome_id")

# Figure: scatter plot of percent vegetarian on average shannon diversity colored by max taxa
ggplot(veg_perc_df_summary2, aes(x = perc_veg, y = avg_shannon, color = as.factor(OTU))) +
  geom_point() +
  labs(x = "Percent Vegetarian", y = "Shannon Diversity", 
       title = "Scatter Plot of Average Shannon Diversity v. Percent Vegetarian") +
  theme_minimal() + 
  theme(legend.position = "None") +
  ylim(0, 5)

# Figure: scatter plot of percent vegetarian on shannon diversity colored by max taxa
ggplot(gut.diet.df, aes(x = perc_veg, y = shannon, color = as.factor(OTU))) +
  geom_point() +
  labs(x = "Percent Vegetarian", y = "Shannon Diversity", 
       title = " ") +
  theme_minimal() + 
  theme(legend.position = "None") +
  ylim(0, 5)

## Mixed Effects Model
library(lme4)
library(lmerTest)
library(performance)

nutrient.diet.22.day <- gut.diet.df %>% 
  group_by(biome_id, logDate) %>% 
  summarise(daily_shannon = sum(shannon, na.rm=TRUE) / n(), # average the samples for the day
            caloriesall_avg = sum(caloriesall, na.rm=TRUE) / n(), # get avg cals per day
            cholesterol_prop = sum(cholesterolall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            satFat_prop = sum(saturatedFatall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            sodium_prop = sum(sodiumall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            carb_prop = sum(carbohydratesall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            dietFib_prop = sum(dietaryFiberall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            sugar_prop = sum(sugarsall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            protein_prop = sum(proteinall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            fat_prop = sum(fatall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            fat_cal_prop = sum(caloriesFromFat, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            addedSugarall_prop = sum(addedSugarall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            perc_veg = first(perc_veg),
            is_vegetarian = first(is_vegetarian)) %>% 
  ungroup() %>% 
  mutate(logDate = as.character(logDate))
nutrient.diet.22.day <- study_days(nutrient.diet.22.day)

## random intercept, fixed slopes

# Full model
lmer.full <- lmer(daily_shannon~caloriesall_avg + cholesterol_prop + satFat_prop + sodium_prop +
                    carb_prop + dietFib_prop + sugar_prop + protein_prop + fat_prop +
                    fat_cal_prop + addedSugarall_prop +
                    (1|`biome_id`), data = nutrient.diet.22.day)
r2(lmer.full)
summary(lmer.full)

# Time
lmer.time.full <- lmer(daily_shannon~study_day + caloriesall_avg + cholesterol_prop + satFat_prop + sodium_prop +
                    carb_prop + dietFib_prop + sugar_prop + protein_prop + fat_prop +
                    fat_cal_prop + addedSugarall_prop +
                    (1|`biome_id`), data = nutrient.diet.22.day)
r2(lmer.time.full)
summary(lmer.time.full)

anova(lmer.full, lmer.time.full)

# Time + Time^2
lmer.time2.full <- lmer(daily_shannon~study_day + I(study_day^2) + caloriesall_avg + cholesterol_prop + satFat_prop + sodium_prop +
                         carb_prop + dietFib_prop + sugar_prop + protein_prop + fat_prop +
                         fat_cal_prop + addedSugarall_prop +
                         (1|`biome_id`), data = nutrient.diet.22.day)
r2(lmer.time2.full)
summary(lmer.time2.full)

anova(lmer.time.full, lmer.time2.full)

# Check vegetarian variable
lmer.full.veg <- lmer(daily_shannon~caloriesall_avg + cholesterol_prop + satFat_prop + sodium_prop +
                        carb_prop + dietFib_prop + sugar_prop + protein_prop + fat_prop +
                        perc_veg +
                        fat_cal_prop + addedSugarall_prop +
                        (1|`biome_id`), data = nutrient.diet.22.day)
r2(lmer.full.veg)
summary(lmer.full.veg)

lmer.full.veg <- lmer(daily_shannon~caloriesall_avg + cholesterol_prop + satFat_prop + sodium_prop +
                        carb_prop + dietFib_prop + sugar_prop + protein_prop + fat_prop +
                        is_vegetarian +
                        fat_cal_prop + addedSugarall_prop +
                        (1|`biome_id`), data = nutrient.diet.22.day)
r2(lmer.full.veg)
summary(lmer.full.veg)

##########################################################################################
library(ggrepel)

bacterial.taxa <- tax_table(bacterial.data)

# aggregate at genus level & relative abundances
genus.taxa <- tax_glom(bacterial.data, taxrank = "Genus")
genus.ra <- transform_sample_counts(genus.taxa, function(x) x / sum(x)) 
bactieral.meta <- as(sample_data(genus.ra), "data.frame")

# make meta data obj
gut.diet.df.subset <- gut.diet.df %>%
  select(SampleID, shannon, biome_id, logDate, sampleType, perc_veg, is_vegetarian)
# genus.ra.merged <- bactieral.meta %>% 
#   left_join(veg_perc_df.subset, by=c("SampleID", "biome_id", "logDate", "sampleType"))
# nutrient.diet.22.day.microbiome.filtered.subset <- nutrient.diet.22.day.microbiome.filtered %>% 
#   select(!c(qr, status, timestamp, sampleType))
genus.ra.nutrient.merged <- gut.diet.df.subset %>% 
  left_join(nutrient.diet.22.day, by=c("biome_id", "logDate", "perc_veg", "is_vegetarian")) %>% 
  as.data.frame()
rownames(genus.ra.nutrient.merged) <- genus.ra.nutrient.merged$SampleID
sample_data(genus.ra) <- sample_data(genus.ra.nutrient.merged)

# ordination
genus.ordination <- ordinate(genus.ra, method = "PCoA", distance = "bray")

# Figure: cluster by vegetarian
plot_ordination(genus.ra, genus.ordination, color = "is_vegetarian") + 
  geom_point(size=3) +
  geom_text_repel(aes(label = biome_id), size = 3) +
  theme_minimal()

# Figure: cluster by percent vegetarian
plot_ordination(genus.ra, genus.ordination, color = "perc_veg") + 
  geom_point(size=3) +
  geom_text_repel(aes(label = biome_id), size = 3) +
  theme_minimal()

# Figure: cluster by Cholesterol Intake
plot_ordination(genus.ra, genus.ordination, color = "cholesterol_prop") + 
  geom_point(size=3) +
  geom_text_repel(aes(label = biome_id), size = 3) +
  theme_minimal()

# Figure: cluster by Study Day
plot_ordination(genus.ra, genus.ordination, color = "study_day") + 
  geom_point(size=3) +
  geom_text_repel(aes(label = biome_id), size = 3) +
  theme_minimal()

# Figure: cluster by satFat_prop
plot_ordination(genus.ra, genus.ordination, color = "satFat_prop") + 
  geom_point(size=3) +
  geom_text_repel(aes(label = biome_id), size = 3) +
  theme_minimal()

# Figure: cluster by carb_prop
plot_ordination(genus.ra, genus.ordination, color = "carb_prop") + 
  geom_point(size=3) +
  geom_text_repel(aes(label = biome_id), size = 3) +
  theme_minimal()

# Figure: cluster by dietFib_prop
plot_ordination(genus.ra, genus.ordination, color = "dietFib_prop") + 
  geom_point(size=3) +
  geom_text_repel(aes(label = biome_id), size = 3) +
  theme_minimal()

# clean df
meta_df <- as(sample_data(genus.ra), "data.frame")
meta_df_clean <- na.omit(meta_df)
genus.ra.clean <- prune_samples(rownames(meta_df_clean), genus.ra)

# analyze if microbiome composition (Bray-Curtis distances) is significantly associated with the dietary factors
vegan::adonis2(distance(genus.ra.clean, method="bray") ~ cholesterol_prop + satFat_prop, #+ 
               # sodium_prop + carb_prop + dietFib_prop + sugar_prop +
               # protein_prop + fat_prop + fat_cal_prop + addedSugarall_prop, 
               data = meta_df_clean)

# time and dietary category
genus.df <- psmelt(genus.ra)
genus.df.time <- study_days(genus.df)

lmer_genus <- lmer(Abundance ~ study_day + Abundance + (1 | biome_id), data = genus.df.time)
summary(lmer_genus)

##########################################################################################

## Activity

activity.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 4-Physical Activity.csv")
activity.data <- activity.data %>% 
  filter(!is.na(as.numeric(biome_id)))
activity.data <- study_days(activity.data)

# gut microbiome data
dim(fecal.microbial.menses.24.filtered)

# join by date and participant
fecal.microbial.menses.24.activity <- fecal.microbial.menses.24.filtered %>% 
  left_join(activity.data, by=c("biome_id", "logDate"))

### Summarize activity data
activity.data.summary <- activity.data %>% 
  group_by(biome_id) %>% 
  summarise(avg_cals_burned = sum(calories_burned) / n(),
            avg_steps = sum(steps) / n(),
            avg_distance = sum(distance) / n(),
            avg_minutes_sedentary = sum(minutes_sedentary) / n(),
            avg_minutes_lightly_active = sum(minutes_lightly_active) / n(),
            avg_minutes_fairly_active = sum(minutes_fairly_active) / n(),
            avg_minues_very_active = sum(minues_very_active) / n(),
            avg_activity_calories = sum(activity_calories) / n(),
            total_min_active = sum(minutes_lightly_active, minutes_fairly_active, minues_very_active) / n()
  )
fecal.microbial.menses.24.summary <- fecal.microbial.menses.24.filtered %>% 
  group_by(biome_id) %>% 
  summarise(avg_shannon=sum(shannon)/n())

activity.data.summary <- fecal.microbial.menses.24.summary %>% 
  left_join(activity.data.summary, by="biome_id")


# steps
lm.shannon.activity <- lm(avg_shannon ~ avg_steps, data = activity.data.summary)
summary(lm.shannon.activity)
# total cals burned - significant
lm.shannon.activity <- lm(avg_shannon ~ avg_cals_burned, data = activity.data.summary)
summary(lm.shannon.activity)
# activity cals burned 
lm.shannon.activity <- lm(avg_shannon ~ avg_activity_calories, data = activity.data.summary)
summary(lm.shannon.activity)
# distance burned 
lm.shannon.activity <- lm(avg_shannon ~ avg_distance, data = activity.data.summary)
summary(lm.shannon.activity)
# minutes sedentary
lm.shannon.activity <- lm(avg_shannon ~ avg_minutes_sedentary, data = activity.data.summary)
summary(lm.shannon.activity)
# minutes lightly active
lm.shannon.activity <- lm(avg_shannon ~ avg_minutes_lightly_active, data = activity.data.summary)
summary(lm.shannon.activity)
# minutes fairly active
lm.shannon.activity <- lm(avg_shannon ~ avg_minutes_fairly_active, data = activity.data.summary)
summary(lm.shannon.activity)
# minutes very active
lm.shannon.activity <- lm(avg_shannon ~ avg_minues_very_active, data = activity.data.summary)
summary(lm.shannon.activity)

## Mixed Effects Model
library(lme4)
library(lmerTest)

## Different baseline, same slope, different intercept

source("~/Microbiome Thesis/functions.R")

fixed_list <- c("steps", "distance", "calories_burned", "minutes_sedentary", "minutes_lightly_active",
                "minutes_fairly_active", "minues_very_active", "activity_calories")

lmer.models <- mixed_effects_fixed_slope(fecal.microbial.menses.24.activity, 10, "shannon", fixed_list)

r2(lmer.models$steps)
r2(lmer.models$distance)
r2(lmer.models$calories_burned)
r2(lmer.models$minutes_sedentary)
r2(lmer.models$minutes_lightly_active)
r2(lmer.models$minutes_fairly_active)
r2(lmer.models$minues_very_active)
r2(lmer.models$activity_calories)

summary(lmer.models$steps)
summary(lmer.models$distance)
summary(lmer.models$calories_burned)
summary(lmer.models$minutes_sedentary)
summary(lmer.models$minutes_lightly_active)
summary(lmer.models$minutes_fairly_active)
summary(lmer.models$minues_very_active)
summary(lmer.models$activity_calories)

## random intercepts, random slops
rnd.slope.lmer.models <- mixed_effects_rnd_slope(fecal.microbial.menses.24.activity, 10, "shannon", fixed_list)

r2(rnd.slope.lmer.models$steps)
r2(rnd.slope.lmer.models$distance)
r2(rnd.slope.lmer.models$calories_burned)
r2(rnd.slope.lmer.models$minutes_sedentary)
r2(rnd.slope.lmer.models$minutes_lightly_active)
r2(rnd.slope.lmer.models$minutes_fairly_active)
r2(rnd.slope.lmer.models$minues_very_active)
r2(rnd.slope.lmer.models$activity_calories)

summary(rnd.slope.lmer.models$steps)
summary(rnd.slope.lmer.models$distance)
summary(rnd.slope.lmer.models$calories_burned)
summary(rnd.slope.lmer.models$minutes_sedentary)
summary(rnd.slope.lmer.models$minutes_lightly_active)
summary(rnd.slope.lmer.models$minutes_fairly_active)
summary(rnd.slope.lmer.models$minues_very_active)
summary(rnd.slope.lmer.models$activity_calories)

# full model - fixed slopes, rnd intercepts
lmer.full <- lmer(shannon~steps + distance +
                    calories_burned + minutes_sedentary + minutes_lightly_active +
                    minutes_fairly_active + minues_very_active + activity_calories +
                    (1|`biome_id`), 
                  data = fecal.microbial.menses.24.activity)
r2(lmer.full)
summary(lmer.full)

# full model - random slopes, rnd intercepts
rnd.slope.lmer.full <- lmer(shannon~ steps + (steps||`biome_id`) + 
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

## check if student athlete matters
participant.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 9-Volunteer Medical History.csv", header = TRUE)

participant.data <- participant.data %>% 
  mutate(sport_collapsed = ifelse(sport=="In-Season", sport, "Not In-Season")) %>%
  # mutate(sport_collapsed = ifelse(sport=="In-Season" | sport == "Club", "Activly in Sport", "Not In-Season")) %>%
  select(biome_id, logDate, activity_level, sport, sport_collapsed, field_hockey)

activity.sport.summary <- activity.data.summary %>% 
  left_join(participant.data, by="biome_id") %>% 
  filter(!is.na(sport))

# Figure: boxplot of sannon diversity by sports type
ggplot(activity.sport.summary, aes(x = sport, y = avg_shannon)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color="orchid") +
  labs(
    x = " ", 
    y = "Average Shannon Diversity", 
    title = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none") 

# testing
aov_sport <- aov(avg_shannon ~ sport, data = activity.sport.summary)
summary(aov_sport)

## Sports collapsed

# Figure: boxplot of sannon diversity by collapsed sports
ggplot(activity.sport.summary, aes(x = sport_collapsed, y = avg_shannon)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color="orchid") +
  labs(
    x = " ", 
    y = "Average Shannon Diversity", 
    title = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none") 

t.test(avg_shannon ~ sport_collapsed, data=activity.sport.summary)
wilcox.test(avg_shannon ~ sport_collapsed, data = activity.sport.summary)

### Minutes active 
lm.shannon.activity <- lm(avg_shannon ~ total_min_active, data = activity.data.summary)
summary(lm.shannon.activity)

## field hockey team
activity.sport.summary2 <- activity.sport.summary %>% 
  mutate(field_hockey=ifelse(field_hockey==TRUE, "fieldHockey", 
                             ifelse(sport_collapsed=="In-Season",
                                    "In-Season", "Off-Season")))

# Figure: boxplot of sannon diversity by field hockey
ggplot(activity.sport.summary2, aes(x = field_hockey, y = avg_shannon)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color="orchid") +
  labs(
    x = " ", 
    y = "Average Shannon Diversity", 
    title = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none") 

t.test(avg_shannon ~ field_hockey, data=activity.sport.summary)
wilcox.test(avg_shannon ~ field_hockey, data = activity.sport.summary)


## Levels of activity
activity.data.summary.levels <- activity.data.summary %>% 
  mutate(exercise_level = case_when(
    total_min_active <= quantile(total_min_active, 0.25, na.rm = TRUE) ~ "Low",
    total_min_active <= quantile(total_min_active, 0.50, na.rm = TRUE) ~ "Moderate",
    total_min_active <= quantile(total_min_active, 0.75, na.rm = TRUE) ~ "High",
    TRUE ~ "Very High"
  ))

activity.data.summary.levels <- activity.data.summary.levels %>% 
  left_join(participant.data, by="biome_id") %>% 
  filter(!is.na(sport))
activity.data.summary.levels$exercise_level <- factor(activity.data.summary.levels$exercise_level, levels = c("Low", "Moderate", "High", "Very High"), ordered = TRUE)

table(activity.data.summary.levels$exercise_level)
table(activity.data.summary.levels$exercise_level, activity.data.summary.levels$sport)

## Sports collapsed

# Figure: boxplot of sannon diversity by activity levels
ggplot(activity.data.summary.levels, aes(x = exercise_level, y = avg_shannon)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color="orchid") +
  labs(
    x = " ", 
    y = "Average Shannon Diversity", 
    title = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none") 

# testing
aov_sport <- aov(avg_shannon ~ exercise_level, data = activity.data.summary.levels)
summary(aov_sport)






