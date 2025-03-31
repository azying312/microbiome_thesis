library(tidyverse)

sleep_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_sleep.csv")
dass_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/DASS_0503_2024-final_df.csv") # read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_dass.csv")
merged_diet_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/fully_merged_diet_data.csv")
fitbit_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 4-Physical Activity.csv")
med_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_med_data.csv")
menses_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/relabeled_data/cleaned_menstruation_data.csv", header=TRUE) # read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_menstruation_data.csv", header=TRUE)
survey_data_full <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 9-Volunteer Medical History.csv", header=TRUE)
sex_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 3-Sexual Activity.csv")

missing_survey_ids <- c(70, 71, 69, 73, 68)
missing_survey <- survey_data_full %>% filter(biome_id %in% missing_survey_ids)
unique(missing_survey$biome_id)

missing_sleep <- sleep_data %>% filter(biome_id %in% missing_survey_ids)
unique(missing_sleep$biome_id)
missing_dass <- dass_data %>% filter(biome_id %in% missing_survey_ids)
unique(missing_dass$biome_id)
missing_diet <- merged_diet_data %>% filter(biome_id %in% missing_survey_ids)
unique(missing_diet$biome_id)

missing_activity <- fitbit_data %>% filter(biome_id %in% missing_survey_ids)
dim(missing_activity)
unique(missing_activity$biome_id)

miss_p <- missing_activity %>% 
  filter(biome_id==70)
dim(miss_p)
unique(miss_p$activity_date)
length(unique(miss_p$Date))

missing_med <- med_data %>% filter(biome_id %in% missing_survey_ids)
unique(missing_med$biome_id)
missing_menses <- menses_data %>% filter(biome_id %in% missing_survey_ids)
unique(missing_menses$biome_id)
# 69 70 71 73
miss_p <- missing_menses %>% 
  filter(biome_id==73)
dim(miss_p)
unique(miss_p$logDate)
length(unique(miss_p$logDate))

missing_sex <- sex_data %>% filter(biome_id %in% missing_survey_ids)
unique(missing_sex$biome_id)

