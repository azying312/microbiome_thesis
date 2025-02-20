########################
#
# Clean Medication Data
# v1: 8 October 2024
#
#########################

library(tidyverse)

med_data <- read.csv("/Volumes/T7/microbiome_data/original_data/Report 2-Medications.csv")
id_mapping <- read.csv("/Volumes/T7/microbiome_data/original_data/Original Study Mapping - Sheet3.csv", header = TRUE)

# Data Prep
med_data <- med_data %>%
  rename("biome_id" = "uid")
med_data <- study_mapping(med_data, id_mapping)

### Save final data output
write.csv(med_data,
          file = "/Volumes/T7/microbiome_data/cleaned_data/prehandcoded_med_data.csv",
          row.names = FALSE)

# Handcoded: "/Volumes/T7/microbiome_data/cleaned_data/cleaned_med_data.csv"
