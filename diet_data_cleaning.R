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
bd_mapping <- read.csv("/Users/alicezhang/Desktop/microbiome_data/BDNutrition - Mapping.csv", header = TRUE)
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
  
  # # This code assumes that uneven number of entries is because serving size was null
  # if(length(temp5)%%2!=0) { #len() of temp 5 needed to be even
  #   print("uneven number of entries")
  #   temp3<-gsub("null","1", temp3)
  #   temp4<-unlist(strsplit(temp3, " "))
  #   isnumber<-grepl("[0-9]", temp4)
  #   temp5<-temp4[isnumber]
  #   if(length(temp5)%%2!=0) print("still uneven number")
  # }
  # 
  # for (jj in seq(1, length(temp5), by = 2)) {
  #   id <- gsub("foodsid|id", "", temp5[jj])  # Clean ID
  #   serving <- ifelse(jj + 1 <= length(temp5), gsub("servings", "", temp5[jj + 1]), NA)  # Clean serving
  # }
  
  if (length(temp5)==0) return (data.frame())
  
  ids <- gsub("foodsid","", temp5[seq(1,length(temp5)-1, 2)])
  ids <- gsub("id","",ids)
  
  servings<-gsub("servings","",temp5[seq(2,length(temp5),2)])
  
  # Get the dessert and beverage ids
  # bd_vector <- ids[ids %in% avi_bd_data$ID]
  # if (length(bd_vector) == 0) {
  #   # Create a data frame with NA values if no matching IDs
  #   bd_df <- data.frame(
  #     biome_id = as.character(dh$biome_id[row]),  # Ensure character type
  #     logDate = as.character(dh$logDate[row]),    # Ensure character type
  #     id = NA,                                    # NA for id
  #     servings = NA,                              # NA for servings
  #     stringsAsFactors = FALSE
  #   )
  # } else {
  #   # Create bd_df with matched IDs
  #   bd_df <- data.frame(
  #     biome_id = as.character(dh$biome_id[row]),  # Ensure character type
  #     logDate = as.character(dh$logDate[row]),    # Ensure character type
  #     id = as.character(bd_vector),                # Ensure character type
  #     servings = as.character(servings[ids %in% avi_bd_data$ID]),  # Ensure character type
  #     stringsAsFactors = FALSE
  #   )
  # }
  # 
  # dropped_ids_list <<- append(dropped_ids_list, list(bd_df))
  # 
  # ids <- as.numeric(gsub("id","",ids)) # at this step beverages & desserts become NA

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

# head(big)
# dim(big) # 8928

# Dropped IDs
# dropped_ids <- do.call(rbind, dropped_ids_list) # 4489

# Add Nutrition info from AVI for beverage and dessert
# dropped_ids_merged <- dropped_ids %>%
#   left_join(avi_bd_data, by = c("id" = "ID"))

# Map to numeric IDs
# bd_ids_merged <- big %>% 
#   left_join(bd_mapping, by = "id")
# dropped_ids_merged <- dropped_ids_merged %>% 
#   mutate(id=new_id) %>% 
#   select(-new_id)

## Use BDNutrition Numerical IDs
big <- big %>%
  left_join(bd_mapping, by = "id") %>%
  mutate(id = coalesce(new_id, id)) %>%
  select(-new_id)
# Change IDs back to numeric
big$id <- as.numeric(big$id)

### Add food other ------------------
# food_other
# Parse food_other
studyID_dh_data$food_other_2 <- lapply(seq_len(nrow(original_dh_data)), function(i) {
  # Skip rows with no "food_other:
  if (original_dh_data$food_other[i] == "[]") {
    return(NULL)  
  }
  # Parse the JSON data
  food_other_parsed <- as.data.frame(fromJSON(original_dh_data$food_other[i]))
  # Create a new column for study_id and logDate
  food_other_parsed$study_id <- studyID_dh_data$biome_id[i]
  food_other_parsed$logDate <- studyID_dh_data$logDate[i]
  
  return(food_other_parsed)
})

# Make food_other DF
food_other_df <- data.frame(
  study_id = numeric(),
  logDate = character(),
  type = character(),
  place = character(),
  foods = character(),
  stringsAsFactors = FALSE
)

# Iterate through each row of dh
for (i in 1:nrow(studyID_dh_data)) {
  x <- studyID_dh_data$food_other_2[i][[1]]
  
  # Check if there are foods
  if (!all(x$foods == ""))  {
    # Create a df for the current row
    # print(i)
    df <- data.frame(
      study_id = studyID_dh_data$biome_id[i],
      logDate = studyID_dh_data$logDate[i],
      type = c(x$type),
      place = ifelse(is.null(x$place), NA, c(x$place)),
      foods = c(x$foods)
    )
    
    # append the df to food_other_df
    food_other_df <- rbind(food_other_df, df)
  }
}

## Clean up food_other df
colnames(food_other_df)[1] <- "study_id"
colnames(food_other_df)[2] <- "Date"
colnames(food_other_df)[5] <- "name"
food_other_df$place[food_other_df$place == ""] <- as.character(NA)

# Filter out dessert and drinks for now
# food_other_bd <- food_other_df %>% 
#   filter(type %in% c("beverage", "dessert"))
# food_other_df_cleaned <- food_other_df %>%
#   filter(name != "") %>% 
#   filter(!(type %in% c("beverage", "dessert")))
food_other_df_cleaned <- food_other_df

### Append big and food_other
food_other_df_cleaned <- food_other_df_cleaned %>%
  mutate(food_other = TRUE)
big <- big %>%
  mutate(food_other = FALSE)
big <- bind_rows(big, food_other_df_cleaned)
# big$id <- as.character(big$id)

# Bind the drinks and desserts
# dropped_ids_merged <- dropped_ids_merged %>% 
#   # select(-c(Source, Notes)) %>% 
#   rename("study_id" = "biome_id") %>% 
#   rename(Date=logDate) %>% 
#   mutate(place = NA) %>%
#   mutate(food_other = FALSE)
# Note: type for the dropped_ids was not kept

# big <- bind_rows(big, dropped_ids_merged)

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

# common <- intersect(colnames(big_cleaned), colnames(studyID_bitesnap_data))
# final_columns <- c(common, c("caloriesFromFat", "saturatedFat", "caloriesFromSatFat", "place"))

# studyID_bitesnap_data$id <- as.character(studyID_bitesnap_data$id)
m <- bind_rows(studyID_bitesnap_data, big_cleaned) 
# %>%
#   select(all_of(final_columns))

# Clean up merged df
m_cleaned <- m %>% 
  select(c("study_id","Date","type","servings","id","name","caloriesall","cholesterolall","saturatedFatall","sodiumall",
           "carbohydratesall","dietaryFiberall","sugarsall","proteinall","fatall","caloriesFromFat","saturatedFat",
           "caloriesFromSatFat","transFat","addedSugarall","food_other","place", everything()))
# 
# sum(is.na(m_cleaned$type))
# length(m_cleaned$type)
# dim(m_cleaned)
# 
# bitesnap_names <- names(studyID_bitesnap_data)
# avi_names <- names(big_cleaned)
# merged_names <- names(m_cleaned)
# 
# # Elements in list1 but not in list2
# merged_names[!merged_names %in% bitesnap_names]
# # m_cleaned$study_id <- as.numeric(m_cleaned$study_id)
# # sum(is.na(m_cleaned$study_id))

# Remove rows with NA in the "study_id" column
# m_cleaned <- m_cleaned %>% filter(!is.na(study_id))

m_cleaned$caloriesFromSatFat <- as.numeric(m_cleaned$caloriesFromSatFat)

### Save final data output
write.csv(m_cleaned,
          file = "/Users/alicezhang/Desktop/microbiome_data/premanual_merged_diet_data.csv",
          row.names = FALSE)

################################################################################################