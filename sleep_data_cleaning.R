########################
#
# Clean Sleep Data
# v1: 8 October 2024
#
#########################

library(tidyverse)

sleep_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Report 5-Sleep.csv")
id_mapping <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Original Study Mapping - Sheet3.csv", header = TRUE)

### Map uid to study id

# Get unique study and biome health pairings
study_and_u_id <- unique(id_mapping %>% 
                           select(STUDY.ID, Biome.Health.App.ID))
# Match and join columns
study_and_u_id <- study_and_u_id %>% 
  rename("study_id" = "STUDY.ID") %>% 
  rename("biome_id" = "Biome.Health.App.ID")
sleep_data <- sleep_data %>% 
  rename("biome_id" = "uid")
study_and_u_id$study_id <- as.character(study_and_u_id$study_id)

# Map ids
sleep_data <- sleep_data %>%
  left_join(study_and_u_id, by = "biome_id") %>%
  mutate(biome_id = coalesce(study_id, biome_id)) %>%
  select(-study_id)

# Check missing ids
missing_list <- sleep_data %>%
  filter(is.na(as.numeric(biome_id)))
print(unique(missing_list$biome_id))

dim(sleep_data)

### Save final data output
write.csv(sleep_data,
          file = "/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_sleep.csv",
          row.names = FALSE)

