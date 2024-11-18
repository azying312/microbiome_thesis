########################
#
# Exploration of the Merged Diet Data
# 11 September 2024
#
#########################

## Packages
packs <- c("tidyverse")
lapply(packs, require, character.only = TRUE)

## Load Data
merged_diet_data <- read_csv("/Users/alicezhang/Desktop/microbiome_data/alice_fully_merged_diet_data.csv")

sum(is.na(merged_diet_data$serving.size))
dim(merged_diet_data)
sum(is.na(merged_diet_data$servings))

dim(merged_diet_data)
# glimpse(merged_diet_data)
names(merged_diet_data)

head(merged_diet_data)

#### Vegetarian
common_meats <- c("beef", "fish", "pork", "chicken", "steak", "turkey", "shrimp", "sushi", "burger king hamburger", 
                    "salmon", "tonkotsu", "mcrib", "bacon", "meatballs", "taco de birria", "meat", "sausage", "burger",
                  "tilapia", "cod", "kielbasa", "clam", "meatloaf", "chop suey", "crab", "ham", "dumplings", "lgblt",
                  "enchilada", "hawaiian pizza", "hot dog", "italian wedding soup", "jimmy johns sandwich", "lamb",
                  "leaky beaker sandwich", "london broil", "mto burrito bar", "mto baked potato bar", "mto deli bar",
                  "mto gyro bar", "pot roast", "meatball calzone", "mcdonald's hamburger", "mulligatawny", "bolognese",
                  "pepperoni", "pho", "chorizo", "ribs", "prosciutto", "scallops", "southwest breakfast wrap",
                  "salami", "taco pizza", "tuna", "wonton", "calamari", "chick fil a sandwich",
                  "duck", "el table sandwich", "wonton soup")

exclude_foods <- c("goldfish", "plant-based chicken nuggets", "goldfish, original", "plant based", "vegan", "vegetarian", 
                   "vegetarian sausage patty", "vegetarian sausage patties", "beyond burger", "no meat", "vegetable sushi")
# Create regex patterns for meat and exclusion lists
meat_pattern <- paste0("\\b(", paste(common_meats, collapse = "|"), ")\\b")
exclude_pattern <- paste0("\\b(", paste(exclude_foods, collapse = "|"), ")\\b")

# Check the common meats
merged_diet_data <- merged_diet_data %>%
  mutate(name_lower = str_to_lower(name),  # Convert names to lowercase
         vegetarian=TRUE,
         # Check for common meats first
         vegetarian = ifelse(str_detect(name_lower, meat_pattern), FALSE, vegetarian)) 

# Look at exclusions
merged_diet_data <- merged_diet_data %>%
  mutate(vegetarian = ifelse(str_detect(name_lower, exclude_pattern), TRUE, vegetarian)) %>% 
  select(-name_lower)

unique_vegetarian_pairs <- merged_diet_data %>%
  select(name, vegetarian) %>%  # Select the relevant columns
  distinct() %>%                   # Get unique combinations
  arrange(name)                    # Optional: arrange by name for easier viewing

## Check if participant vegetarian
participant_vegetarian_status <- merged_diet_data %>%
  group_by(study_id) %>%
  # Check if all foods for a participant are vegetarian
  summarise(is_vegetarian = all(vegetarian))  
participant_vegetarian_status$study_id <- as.numeric(participant_vegetarian_status$study_id)
participant_vegetarian_status <- participant_vegetarian_status %>%
  drop_na()

table(participant_vegetarian_status$is_vegetarian) # 7 vegetarian

# Add to diet data
participant_vegetarian_status$study_id <- as.character(participant_vegetarian_status$study_id)
merged_diet_data <- merged_diet_data %>%
  left_join(participant_vegetarian_status, by = "study_id")

### Save final data output
# write.csv(merged_diet_data,
#           file = "/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_diet.csv",
#           row.names = FALSE)

###############

## Creating dataset with participant, date, item, and nutritional info
big <- merged_diet_data

### Collapse by Day
small <- big %>%
  group_by(study_id, Date) %>%
  summarise(across(
    c(caloriesall,cholesterolall,carbohydratesall,sodiumall,
      proteinall,transFat,caloriesFromSatFat,caloriesFromFat,
      sugarsall,fatall,addedSugarall,dietaryFiberall), sum
  ))



# Week Variable
small$Timestamp <- as.Date(small$Date, format="%m/%d/%y", tz="UTC")
time_values <- sort(unique(small$Timestamp))
num_time_values <- as.numeric(time_values)

## Figuring out numbers and dates
as.numeric(as.Date("2022-10-12", format="%Y-%m-%d"))


small$week <- rep(NA, nrow(small))
# Week 1: October 12-October 18
small$week[small$Timestamp >= 19272 & small$Timestamp <= 19280] <- 1
small$week[small$Timestamp >= 19281 & small$Timestamp <= 19286] <- 2
small$week[small$Timestamp >= 19287 & small$Timestamp <= 19293] <- 3
small$week[small$Timestamp >= 19294 & small$Timestamp <= 19300] <- 4
small$week[small$Timestamp >= 19301 & small$Timestamp <= 19308] <- 5
small$week[small$Timestamp >= 19309 & small$Timestamp <= 19315] <- 6
small$week[small$Timestamp >= 19316 & small$Timestamp <= 19321] <- 7
small$week[small$Timestamp >= 19322 & small$Timestamp <= 19329] <- 8
small$week[small$Timestamp >= 19330 & small$Timestamp <= 19336] <- 9
small$week[small$Timestamp >= 19337] <- 10

table(small$week)
sum(is.na(small$week)) # 10

small_naweek10 <- small %>% 
  filter(is.na(week))
dim(small_naweek10)
#  1   2   3   4   5   6   7   8   9  10 
# 150 239 216 184 177 120  49 111  89  57 

# Creating a variable with weeks
small$week <- rep(NA, nrow(small))
small$week[small$Date >= "2022-10-08" & small$Date <= "2022-10-14"] <- 1
small$week[small$Date >= "2022-10-15" & small$Date <= "2022-10-21"] <- 2
small$week[small$Date >= "2022-10-22" & small$Date <= "2022-10-27"] <- 3
small$week[small$Date >= "2022-10-28" & small$Date <= "2022-11-03"] <- 4
small$week[small$Date >= "2022-11-04" & small$Date <= "2022-11-10"] <- 5
small$week[small$Date >= "2022-11-11" & small$Date <= "2022-11-17"] <- 6
small$week[small$Date >= "2022-11-18" & small$Date <= "2022-11-24"] <- 7
small$week[small$Date >= "2022-11-25" & small$Date <= "2022-12-01"] <- 8
small$week[small$Date >= "2022-12-02" & small$Date <= "2022-12-08"] <- 9
small$week[small$Date >= "2022-12-09" & small$Date <= "2022-12-15"] <- 10
small$week[small$Date >= "2022-12-15" & small$Date <= "2022-12-21"] <- 11

# 1   2   3   4   5   6   7   8   9  10  11 
# 109 280 190 184 173 120  76  89  91  73   4 
# 13 NA
sum(is.na(small$week)) # 13
small_naweek <- small %>% 
  filter(is.na(week))
small_week11<- small %>% 
  filter(week==11)
dim(small_week11)

#########################

# Collapse by Meal
small2<-big%>%
  group_by(study_id, Date,type)%>%
  summarise(across(c(caloriesall, cholesterolall, carbohydratesall,
                     sodiumall, proteinall, transFat, caloriesFromSatFat,
                     caloriesFromFat, sugarsall, fatall, addedSugarall, dietaryFiberall), sum))


# Contigency table for Study Participant and Meal Type frequency
table(small2$study_id, small2$type)

#### TO DO: beverage and dessert need to be added to the collapse by meal

# Look at beverages
table(merged_diet_data$type)

beverage_data <- merged_diet_data %>% 
  filter(type == "beverage")
dim(beverage_data)
head(beverage_data)

table(beverage_data$name)

# Look at other
other_data <- merged_diet_data %>% 
  filter(type == "other")
dim(other_data)
head(other_data)

table(other_data$name)
# These are drinks? Tequila, 100% Honeycrisp Apple Cider, Apple Cider, Beer, Boba, Brewed Tea, Hot chocolate, latte, Nesquik Chocolate Powder, Pina Colada, Rose wine, sangria, shirley temple, tea, vodka


