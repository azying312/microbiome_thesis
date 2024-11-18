########################
#
# Merge the Handcoded Diet Data
# v1: 1 October 2024
#
#########################

# Packages
library(tidyverse)

merged_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/premanual_merged_diet_data.csv", header = TRUE)
other_food_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/other_food_handcoded.csv", header = TRUE)

# other_food_data_NA <- other_food_data %>% 
#   filter(is.na(food_other))
# View(other_food_data_NA)
# Set all food_other to TRUE

# glimpse(other_food_data)
# other_food_data$study_id <- as.numeric(other_food_data$study_id)
other_food_data$caloriesall <- as.numeric(other_food_data$caloriesall)
other_food_data$sodiumall <- as.numeric(other_food_data$sodiumall)
other_food_data$proteinall <- as.numeric(other_food_data$proteinall)
other_food_data$caloriesFromFat <- as.numeric(other_food_data$caloriesFromFat)
other_food_data$fatall <- as.numeric(other_food_data$fatall)
other_food_data$addedSugarall <- as.numeric(other_food_data$addedSugarall)
other_food_data$carbohydratesall <- as.numeric(other_food_data$carbohydratesall)

other_food_data <- other_food_data %>% 
  select(-c(X, Source, Notes, Do.we.think.this.is.a.duplicate.)) %>% 
  rename(name = name..highlighted.items...either.couldn.t.find.nutritional.information.for.the.exact.item..or.the.entry.seemed.ambiguous.could.mean.multiple.things..)

### Clean hand_code
other_food_data$food_other <- TRUE
other_food_data <- other_food_data %>% 
  filter(!(is.na(study_id)))
sum(is.na(other_food_data$study_id))

# Take out other_food == TRUE
no_other_food <- merged_data %>% 
  filter(food_other != TRUE)

# Merge other_food
fully_merged_data <- no_other_food %>% 
  bind_rows(other_food_data)

dim(fully_merged_data) # 11321
names(fully_merged_data)

# Remove all cols that are all NA



### Save final data output
write.csv(fully_merged_data,
          file = "/Users/alicezhang/Desktop/microbiome_data/alice_fully_merged_diet_data.csv",
          row.names = FALSE)
