########################
#
# Clean Diet Data
# v1: 17 September 2024
#
#########################

# Packages
library(tidyverse)
library(jsonlite)
library(utils)

original_dh_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Report 6-DH Food.csv", header = TRUE)
avi_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/AVI Foods (correct) - avi_foods.csv", header = TRUE)
id_mapping <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Original Study Mapping - Sheet3.csv", header = TRUE)
original_bitesnap_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Bitesnap-filtered - Bitesnap-filtered.csv", header = TRUE)
avi_bd_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/BDNutrition - Sheet1.csv", header = TRUE)
bd_mapping <- read.csv("/Users/alicezhang/Desktop/microbiome_data/BDNutrition - Sheet3.csv", header = TRUE)
food_other <- read.csv("/Users/alicezhang/Desktop/microbiome_data/food_other.csv", header = TRUE)

## Merge avi_data and avi_bd_data
mapped_avi_bd_data <- avi_bd_data %>% 
  rename(id=ID) %>% 
  select(-c(Notes, Source))

avi_data$id <- as.character(avi_data$id)
bd_mapping$new_id <- as.character(bd_mapping$new_id)

merged_avi_data <- avi_data %>% 
  bind_rows(mapped_avi_bd_data)

### Recode uid to numbers 1-75

# Get unique study and biome health pairings
study_and_u_id <- unique(id_mapping %>% 
  select(STUDY.ID, Biome.Health.App.ID))

# Match and join columns
study_and_u_id <- study_and_u_id %>% 
  rename("study_id" = "STUDY.ID") %>% 
  rename("biome_id" = "Biome.Health.App.ID")

original_dh_data <- original_dh_data %>% 
  rename("biome_id" = "uid")

study_and_u_id$study_id <- as.character(study_and_u_id$study_id)

studyID_dh_data <- original_dh_data %>%
  left_join(study_and_u_id, by = "biome_id") %>%
  mutate(biome_id = coalesce(study_id, biome_id)) %>%
  select(-study_id)

### Some IDs are missing
missing_list <- studyID_dh_data %>%
  filter(is.na(as.numeric(biome_id)))
print(unique(missing_list$biome_id))

### Splitting days into individual item entries
# Onerow generates warnings if there are IDs that are not numbers (beverages & desserts)
# Note: Desserts and Beverages are dropped here && prints "uneven number of entries" if serving size is NULL

# Initialize an empty data frame to store dropped IDs
# dropped_ids <- data.frame(biome_id = character(), logDate = character(), id = character(), servings = character(), stringsAsFactors = FALSE)
# dropped_ids_list <- list()

# 25 - missing serving
# food_other is filtered out
onerow<-function(dh, row){ 
  #print(row) #only important for de-bugging
  temp<-dh[row,"foods"]
  temp2<-strsplit(temp, ",") #splits row where comma
  temp3<-gsub("[[:punct:]]","", temp2) #takes out punctuation
  
  # Replace all occurrences of "null" with "1" in the string
  temp3 <- gsub("null", "1", temp3)
  
  temp4<-unlist(strsplit(temp3, " ")) #separate string when space
  isnumber<-grepl("[0-9]", temp4) #pull out numbers 
  temp5<-temp4[isnumber]
  temp5 <- gsub("\\n", "", temp5)
  
  if (length(temp5)==0) return (data.frame())
  
  ids <- gsub("foodsid","", temp5[seq(1,length(temp5)-1, 2)])
  ids <- gsub("id","",ids)
  
  servings<-gsub("servings","",temp5[seq(2,length(temp5),2)])
  
  #running for loop because some ids appear multiple times in one day, and 
  ## some food ids aren't in the avi dataset
  noserving <- c()
  
  food <- data.frame()
  
  for(item in 1:length(ids)) {
    newfood <- merged_avi_data[merged_avi_data$id%in%ids[item],]
    # newfood <- avi_data[avi_data$id%in%ids[item],]
    food<-rbind(food, newfood)
    if(nrow(newfood)==0) noserving<-c(noserving,item)
  }
  if (length(noserving)>0) servings<-servings[-noserving] 
  
  
  if(nrow(food)==0) {
    return(data.frame())
  } else {
    return (data.frame(dh$biome_id[row], dh$logDate[row], servings, food))
  }
}

# Include type (ie, breakfast) --------------------------------------------
parsed_foods <- lapply(studyID_dh_data$foods, fromJSON) #study id and date gone 
big <- data.frame()

for(ii in 1:nrow(studyID_dh_data)) {
  # print(ii)
  temp1 <- onerow(studyID_dh_data, ii)
  if (nrow(temp1) > 0) { #if temp1 is NOT empty
    temp2 <- parsed_foods[[ii]]
    temp1$type <- NA
    
    for(i in 1:nrow(temp1)){
      id <- temp1$id[i]
      temp_row <- grep(id, temp2$foods)
      temp1$type[i] <- temp2$type[temp_row]
    }
    
    big <- rbind(big, temp1)
  }
}

# Change names
colnames(big)[1] <- "study_id"
colnames(big)[2] <- "Date"

## Use BDNutrition Numerical IDs
big <- big %>%
  left_join(bd_mapping, by = "id") %>%
  mutate(id = coalesce(new_id, id)) %>%
  select(-new_id)
# Change IDs back to numeric
big$id <- as.numeric(big$id)

# Multiply all vars by # of servings --------------------
big$caloriesall<-big$calories*as.numeric(big$servings)
big$cholesterolall<-big$cholesterol*as.numeric(big$servings)
big$caloriesFromFatall<-big$caloriesFromFat*as.numeric(big$servings)
big$saturatedFatall<-big$saturatedFat*as.numeric(big$servings)
big$caloriesFromSatFatall<-as.numeric(big$caloriesFromSatFat)*as.numeric(big$servings)
big$transFatall<-big$transFat*as.numeric(big$servings)
big$sodiumall<-big$sodium*as.numeric(big$servings)
big$carbohydratesall<-big$carbohydrates*as.numeric(big$servings)
big$dietaryFiberall<-big$dietaryFiber*as.numeric(big$servings)
big$sugarsall<-big$sugars*as.numeric(big$servings)
big$addedSugarall<-big$addedSugar*as.numeric(big$servings)
big$proteinall<-big$protein*as.numeric(big$servings)
big$fatall<-big$fat*as.numeric(big$servings)

### Add food other ------------------

## Clean up food_other df - hand cleaned
food_other_df <- food_other %>% 
  filter(IS_DUPE==FALSE) %>% 
  filter(food_other==TRUE) %>% 
  select(!c(Assigner, Comments.Column, IS_DUPE, servings)) %>% 
  # set servings
  rename(servings=serving.size) %>% 
  mutate(servings=ifelse(!grepl("^[0-9]+$", servings), 1, servings)) %>% 
  # clean df
  filter(!is.na(study_id)) %>% 
  filter(!is.na(name)) %>% 
  filter(!is.na(Date)) %>% 
  # get rid of all NA cols
  select(where(~ !all(is.na(.)))) 

### Append big and food_other
food_other_df_cleaned <- food_other_df %>%
  mutate(food_other = TRUE)
big <- big %>%
  mutate(food_other = FALSE) %>% 
  mutate(caloriesFromSatFat=as.numeric(caloriesFromSatFat))
big <- bind_rows(big, food_other_df_cleaned)

# Clean up big df
big_cleaned <- big %>% 
  select(c("study_id","Date","type","servings","id","name","caloriesall","cholesterolall","saturatedFatall","sodiumall",
           "carbohydratesall","dietaryFiberall","sugarsall","proteinall","fatall","caloriesFromFat","saturatedFat",
           "caloriesFromSatFat","transFat","addedSugarall","food_other","place"))

# Empty Type -- remove
sum(is.na(big_cleaned$type))
big_cleaned %>% 
  filter(is.na(big_cleaned$type))
big_cleaned <- big_cleaned %>% 
  filter(!is.na(type))

#### Bitesnap Data
# map uid to study_id
bitesnap_data <- original_bitesnap_data %>% 
  rename("biome_id" = "uid")

studyID_bitesnap_data <- bitesnap_data %>%
  left_join(study_and_u_id, by = "biome_id") %>%
  mutate(biome_id = coalesce(study_id, biome_id)) %>%
  select(-study_id) %>% 
  rename("study_id" = "biome_id")

#use times to make a meal column 
#how many meals per person? 
posix_date_time <- as.POSIXct(studyID_bitesnap_data$mod_timestamp, format = "%Y-%m-%d %H:%M:%S")
studyID_bitesnap_data$times <- format(posix_date_time, format = "%H:%M:%OS")
studyID_bitesnap_data$Date <- format(posix_date_time, format = "%Y-%m-%d")
studyID_bitesnap_data$food_other <- FALSE

#classifying by meal
#dhall hours
#if between 7 and 10, breakfast
#if between 1130 and 2, lunch
#if between 5 and 7, dinner

## our accommodations 
#if between 5 and 11, breakfast
#if between 11 and 3, lunch
#if between 3 and 12, dinner
#### CHECK TIMES ON THE THINGS THAT WERE NOT ASSIGNED A CATEGORY

studyID_bitesnap_data$type <- rep("other", nrow(studyID_bitesnap_data)) 
studyID_bitesnap_data$type[studyID_bitesnap_data$times >= "05:00:00" & studyID_bitesnap_data$times <= "10:59:59"] <- "breakfast" #5 to 11
studyID_bitesnap_data$type[studyID_bitesnap_data$times >= "11:00:00" & studyID_bitesnap_data$times <= "14:59:59"] <- "lunch" #11 to 3
studyID_bitesnap_data$type[studyID_bitesnap_data$times >= "15:00:00" & studyID_bitesnap_data$times <= "23:59:59"] <- "dinner" #3 to 12
table(studyID_bitesnap_data$type)

# other_type <- studyID_bitesnap_data %>%
#   filter(type=="other")
# other_type$times
# other_type$study_id
# Rename bitesnap_data cols
studyID_bitesnap_data <- studyID_bitesnap_data %>% 
  rename("name" = "foodItemName") %>% 
  rename("servings" = "serving") %>% 
  rename("caloriesall" = "calories") %>% 
  rename("cholesterolall" = "cholesterol") %>% 
  rename("sodiumall" = "sodium") %>% 
  rename("carbohydratesall" = "totalCarb") %>% 
  rename("dietaryFiberall" = "dietaryFiber") %>% 
  rename("sugarsall" = "sugars") %>% 
  rename("proteinall" = "protein") %>%
  rename("fatall" = "totalFat") %>% 
  rename("saturatedFatall" = "saturatedFat") %>%
  rename("addedSugarall" = "addedSugars")

m <- bind_rows(studyID_bitesnap_data, big_cleaned) 

# Clean up merged df
m_cleaned <- m %>% 
  select(c("study_id","Date","type","servings","id","name","caloriesall","cholesterolall","saturatedFatall","sodiumall",
           "carbohydratesall","dietaryFiberall","sugarsall","proteinall","fatall","caloriesFromFat","saturatedFat",
           "caloriesFromSatFat","transFat","addedSugarall","food_other","place", everything()))

m_cleaned$caloriesFromSatFat <- as.numeric(m_cleaned$caloriesFromSatFat)

m_cleaned <- m_cleaned %>% 
  select(where(~ !all(is.na(.)))) %>% 
  filter(m_cleaned$name!="")

### Save final data output
write.csv(m_cleaned,
          file = "/Users/alicezhang/Desktop/microbiome_data/manual_merged_diet_data.csv",
          row.names = FALSE)

################################################################################################