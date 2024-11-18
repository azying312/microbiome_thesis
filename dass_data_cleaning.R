########################
#
# Clean DASS Data
# v1: 8 October 2024
#
#########################

library(tidyverse)

dass_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Report 10-DASS-21.csv")
id_mapping <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Original Study Mapping - Sheet3.csv", header = TRUE)

### Map uid to study id

# Get unique study and biome health pairings
study_and_u_id <- unique(id_mapping %>% 
                           select(STUDY.ID, Biome.Health.App.ID))
# Match and join columns
study_and_u_id <- study_and_u_id %>% 
  rename("study_id" = "STUDY.ID") %>% 
  rename("biome_id" = "Biome.Health.App.ID")
dass_data <- dass_data %>% 
  rename("biome_id" = "Your.Biome.Health.App.ID")
study_and_u_id$study_id <- as.character(study_and_u_id$study_id)

# Map ids
dass_data <- dass_data %>%
  left_join(study_and_u_id, by = "biome_id") %>%
  mutate(biome_id = coalesce(study_id, biome_id)) %>%
  select(-study_id)

# Check missing ids
missing_list <- dass_data %>%
  filter(is.na(as.numeric(biome_id)))
print(unique(missing_list$biome_id))

dim(dass_data)

