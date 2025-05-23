########################
#
# Merge the Handcoded Diet Data
# Last updated: 1 October 2024
#
#########################

# Packages
library(tidyverse)

merged_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/diet_intermediary/no_dupes_diet_data.csv", header = TRUE)
# merged_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/diet_intermediary/manual_merged_diet_data.csv", header = TRUE)
other_food_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/diet_intermediary/check_dupes_diet_data - handcoded.csv", header = TRUE)
type_other_food_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/diet_intermediary/type_other_diet_data - handcoded.csv", header = TRUE)
no_nutrition_map <-read.csv("/Volumes/T7/microbiome_data/cleaned_data/diet_intermediary/no_nutrition - handcoded.csv", header = TRUE)

merged_data <- merged_data %>% 
  mutate(study_id=as.numeric(study_id),
         servings=as.numeric(servings))

merged_data <- merged_data %>% 
  filter(!is.na(study_id)) %>% 
  # Take out other_food == TRUE
  filter(food_other==FALSE) %>%
  # Take out type other
  filter(type !="other")

### Clean hand_code

# Food Other; filter dupes
dim(other_food_data) # 2701   73
other_food_data <- other_food_data %>% 
  select(!c(Assigner, Comments.Column)) %>%
  mutate(study_id=as.numeric(study_id),
         servings=as.numeric(servings)) %>% 
  # get the food other
  # filter(food_other==TRUE) %>% 
  filter(type != "other") %>%
  filter(IS_DUPE==FALSE)
dim(other_food_data) # 2642   73 | 59 dupes

# what to do about foods that are not in food other

other_food_data <- other_food_data %>% # other_food_data rows
  filter(!(is.na(study_id))) %>% 
  mutate(servings=as.numeric(servings))
dim(other_food_data) # 2629   73 | filter out data without study id

# Type other food | 10
type_other_food_data <- type_other_food_data %>% 
  mutate(
    study_id=as.numeric(type_other_food_data$study_id),
    id=as.integer(id),
    servings=as.numeric(servings)
  )

# glimpse(other_food_data)
other_food_data$study_id <- as.numeric(other_food_data$study_id)
other_food_data$caloriesall <- as.numeric(other_food_data$caloriesall)
other_food_data$sodiumall <- as.numeric(other_food_data$sodiumall)
other_food_data$proteinall <- as.numeric(other_food_data$proteinall)
other_food_data$caloriesFromFat <- as.numeric(other_food_data$caloriesFromFat)
other_food_data$fatall <- as.numeric(other_food_data$fatall)
other_food_data$addedSugarall <- as.numeric(other_food_data$addedSugarall)
other_food_data$carbohydratesall <- as.numeric(other_food_data$carbohydratesall)
# other_food_data <- other_food_data %>% 
#   select(-c(X, Source, Notes, Do.we.think.this.is.a.duplicate.)) %>% 
#   rename(name = name..highlighted.items...either.couldn.t.find.nutritional.information.for.the.exact.item..or.the.entry.seemed.ambiguous.could.mean.multiple.things..)

# Merge other_food
dim(merged_data) #  8093   68
fully_merged_data <- merged_data %>% 
  bind_rows(other_food_data) %>%
  bind_rows(type_other_food_data)
dim(fully_merged_data) # 10732    74

## Add nutrition data
no_nutrition_map <- no_nutrition_map %>% 
  select(c(type, servings, id, name, caloriesall, cholesterolall, saturatedFatall, sodiumall, carbohydratesall, 
           dietaryFiberall, sugarsall, proteinall, fatall)) %>%
  mutate(id=as.integer(id)) 
dim(no_nutrition_map) # 36 13 

d_merged <- fully_merged_data %>% 
  rows_update(no_nutrition_map, by="name")

dim(d_merged) # 11321; NEW: 11262 x 73 | 10732 x 74
# names(fully_merged_data)

# Remove all cols that are all NA
m_clean <- d_merged %>%
  select(where(~ !all(is.na(.)))) # 11262 x 68
dim(m_clean) # 10732

m_filtered <- m_clean %>% 
  mutate(Date=as.Date(Date)) %>% 
  filter(Date < as.Date("2022-12-17")) 
dim(m_filtered) # 10644    37

table(m_filtered$type)
dim(m_filtered) # 11013 x 68 | 10644    37

m_filtered <- m_filtered %>% 
  filter(!is.na(caloriesall)) # 10176
dim(m_filtered) # 9807   37

m_filtered <- m_filtered %>% 
  mutate(servings=ifelse(is.na(servings), 1, servings))

m_filtered <- m_filtered %>% 
  rename(biome_id=study_id)

dim(m_filtered)
length(unique(m_filtered$biome_id))

### Save final data output
write.csv(m_filtered,
          file = "/Volumes/T7/microbiome_data/cleaned_data/fully_merged_diet_data.csv",
          row.names = FALSE)
