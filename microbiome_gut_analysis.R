########################
#
# Gut Microbiome Analysis with Lifestyle Factors
# Last updated: 03/21/2025
#
#########################

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
shannon.qr.merged.24 <- shannon.qr.merged.24 %>% 
  mutate(biome_id=as.integer(biome_id)) %>% 
  filter(!is.na(biome_id)) %>% 
  # Filter the data to be within study days: 10-14 to 12-14
  filter(logDate > "2022-10-13" & logDate < "2022-12-15") # 2038 to 1971

dim(shannon.qr.merged.24)

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

# write.csv(fecal.microbial.menses.24, file="/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/relabeled_data/gut.microbial.menses.24.csv")

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


##########################################################################################
# Gut microbiome and menses

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

# join fecal data with diet data
gut.diet.df <- merged_diet_data_collapsed %>% 
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

lm.obj <- lm(shannon ~ perc_veg, data = gut.diet.df)
summary(lm.obj)

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
  summarise(avg_shannnon = sum(shannon, na.rm=TRUE) / n(),
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

# Full model (no avg cal)
lmer.full <- lmer(avg_shannnon~caloriesall_avg + cholesterol_prop + satFat_prop + sodium_prop +
                    carb_prop + dietFib_prop + sugar_prop + protein_prop + fat_prop +
                    fat_cal_prop + addedSugarall_prop +
                    (1|`biome_id`), data = nutrient.diet.22.day)
r2(lmer.full)
summary(lmer.full)

# Time
lmer.time.full <- lmer(avg_shannnon~study_day + caloriesall_avg + cholesterol_prop + satFat_prop + sodium_prop +
                    carb_prop + dietFib_prop + sugar_prop + protein_prop + fat_prop +
                    fat_cal_prop + addedSugarall_prop +
                    (1|`biome_id`), data = nutrient.diet.22.day)
r2(lmer.time.full)
summary(lmer.time.full)

anova(lmer.full, lmer.time.full)

# Time + Time^2
lmer.time2.full <- lmer(avg_shannnon~study_day + I(study_day^2) + caloriesall_avg + cholesterol_prop + satFat_prop + sodium_prop +
                         carb_prop + dietFib_prop + sugar_prop + protein_prop + fat_prop +
                         fat_cal_prop + addedSugarall_prop +
                         (1|`biome_id`), data = nutrient.diet.22.day)
r2(lmer.time2.full)
summary(lmer.time2.full)

anova(lmer.time.full, lmer.time2.full)

# Check vegetarian variable
lmer.full.veg <- lmer(avg_shannnon~caloriesall_avg + cholesterol_prop + satFat_prop + sodium_prop +
                        carb_prop + dietFib_prop + sugar_prop + protein_prop + fat_prop +
                        perc_veg +
                        fat_cal_prop + addedSugarall_prop +
                        (1|`biome_id`), data = nutrient.diet.22.day)
r2(lmer.full.veg)
summary(lmer.full.veg)

lmer.full.veg <- lmer(avg_shannnon~caloriesall_avg + cholesterol_prop + satFat_prop + sodium_prop +
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




