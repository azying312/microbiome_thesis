########################
#
# Clean Samples Data
# v1: 8 October 2024
#
#########################

library(tidyverse)

source("~/Microbiome Thesis/functions.R")

samples_data <- read.csv("/Volumes/T7/microbiome_data/original_data/Report 8-Sample.csv")
id_mapping <- read.csv("/Volumes/T7/microbiome_data/Original Study Mapping - Sheet3.csv", header = TRUE)
uminn_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_uminn_data.csv", header=TRUE)

### Map uid to study id
samples_data <- samples_data %>% 
  rename("biome_id"="uid")
samples_data <- study_mapping(samples_data, id_mapping)

# UMinn Subset
uminn_data_subset <- uminn_data %>% 
  select(Sample.ID, Well, Special.Notes) 
uminn_data_subset$qr <- sub("_.*", "", uminn_data_subset$Sample.ID)
uminn_data_subset$inUminn <- TRUE

### Duplicates and Errors in Vaginal Swabs
vaginal_data <- samples_data %>% 
  filter(sampleType=="vaginal") %>% 
  filter(!is.na(biome_id)) %>%  
  filter(logDate != "0000-00-00")
dim(vaginal_data) # 4024
vaginal_data$vaginal_swab <- TRUE
vaginal_data <- uminn_data_subset %>% 
  left_join(vaginal_data, by="qr") # 2894 actual samples

## Filter out any errors
vaginal_data <- vaginal_data %>% 
  filter(logDate!="0000-00-00") %>% # no corresponding logDate
  mutate(uMinn_menstruation=ifelse(str_detect(Special.Notes, "Blood")==TRUE, 1, 0)) %>%  # set to menstruation true
  select(biome_id, logDate, qr, timestamp, uMinn_menstruation, inUminn, vaginal_swab, Special.Notes) %>% 
  mutate(biome_id=as.numeric(biome_id)) %>% 
  # errors in sample processing from UMinn
  filter(!str_detect(Special.Notes, "error")) %>% 
  filter(!str_detect(Special.Notes, "No swab in tube")) %>% 
  filter(!is.na(biome_id))
dim(vaginal_data) # 1536    8

vaginal_data <- vaginal_data %>% 
  select(biome_id, logDate, qr, timestamp) %>% 
  mutate(sampleType="vaginal")

## Join back cleaned vaginal data
samples_data <- samples_data %>% 
  filter(sampleType!="vaginal")
samples_data <- samples_data %>% 
  full_join(vaginal_data)

# Fix log dates
samples_data <- samples_data %>% 
  mutate(logDate = as.Date(timestamp),
         timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S"))

dim(samples_data) # 4754 6

write.csv(samples_data,
          file = "/Volumes/T7/microbiome_data/cleaned_data/cleaned_samples.csv",
          row.names = FALSE)

