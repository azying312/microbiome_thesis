########################
#
# Clean Sleep Data
# v1: 8 October 2024
#
#########################

source("~/Microbiome Thesis/functions.R")

library(tidyverse)

sleep_data <- read.csv("/Volumes/T7/microbiome_data/original_data/Report 5-Sleep.csv")
id_mapping <- read.csv("/Volumes/T7/microbiome_data/original_data/Original Study Mapping - Sheet3.csv", header = TRUE)

# Data Prep
sleep_data <- sleep_data %>%
  rename("biome_id" = "uid")
sleep_data <- study_mapping(sleep_data, id_mapping)

### Save final data output
write.csv(sleep_data,
          file = "/Volumes/T7/microbiome_data/cleaned_data/cleaned_sleep.csv",
          row.names = FALSE)

