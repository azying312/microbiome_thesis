# Packages
library(tidyverse)

merged_first <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Fully Merged Diet (Dining Hall, Bitesnap, Food_Other, Beverages, Desserts) - nutritiondata_DHBSFOHC.csv")
premanual_merged_alice <- read.csv("/Users/alicezhang/Desktop/microbiome_data/premanual_merged_diet_data.csv")
old_not_food_other <- read.csv("/Users/alicezhang/Desktop/microbiome_data/old_not_food_other.csv")
fully_merged_alice <- read.csv("/Users/alicezhang/Desktop/microbiome_data/alice_fully_merged_diet_data.csv")
fully_merged_old <- read.csv("/Users/alicezhang/Desktop/microbiome_data/old_merged.csv")

# Change type
fully_merged_old$study_id <- as.numeric(fully_merged_old$study_id)
fully_merged_old$id <- as.numeric(fully_merged_old$id)
fully_merged_old$caloriesall <- as.numeric(fully_merged_old$caloriesall)
fully_merged_old$sodiumall <- as.numeric(fully_merged_old$sodiumall)
fully_merged_old$proteinall <- as.numeric(fully_merged_old$proteinall)
fully_merged_old$caloriesFromFat <- as.numeric(fully_merged_old$caloriesFromFat)
fully_merged_old$fatall <- as.numeric(fully_merged_old$fatall)
fully_merged_old$addedSugarall <- as.numeric(fully_merged_old$addedSugarall)

## Difference between fully merged - I have 31 more rows
diff_df1 <- anti_join(fully_merged_alice, fully_merged_old, by = "id")
dim(diff_df1)
View(diff_df1)

diff_df2 <- anti_join(old_not_food_other, merged_alice, by = "id") # I have all of their rows
dim(diff_df2)
View(diff_df2)

## Look at food_other
food_other <- fully_merged_alice %>% 
  filter(food_other==TRUE)
dim(food_other)

fully_merged_alice %>% 
  filter(study_id==45 & Date=='2022-10-20' & type=='breakfast')

# Investigate
fully_merged_alice %>% 
  filter(study_id == 65 & Date == '2022-10-14' & type == 'lunch')
fully_merged_alice %>% 
  filter(study_id == 7 & Date == '2022-10-17' & type == 'dinner')

fully_merged_alice %>% 
  filter(name == "Spaghetti and Meatballs")
fully_merged_alice %>% 
  filter(name == "Pudding")

## Explanation: this meal was somehow dropped
fully_merged_alice %>% 
  filter(study_id == 76 & Date == '2022-12-13' & type == 'dinner')
fully_merged_alice %>% # this has others that are in the old fully merged
  filter(name == "Apple")
fully_merged_alice %>% 
  filter(name == "Honey") # this has others that are in the old fully merged
fully_merged_alice %>% 
  filter(name == "Turkey Sub Sandwich")
fully_merged_alice %>% # this has others that are in the old fully merged
  filter(name == "Cheese Pizza")

