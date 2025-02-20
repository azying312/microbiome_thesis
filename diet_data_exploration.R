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
merged_diet_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/fully_merged_diet_data.csv")

##########################################################################################
# Visualize diet completeness map

merged_diet_data <- merged_diet_data %>% 
  mutate(Date=as.Date(Date)) %>% 
  filter(Date < as.Date("2022-12-17")) 

all_days <- seq.Date(as.Date(min(merged_diet_data$Date, na.rm=TRUE)), 
                     as.Date(max(merged_diet_data$Date, na.rm=TRUE)), by = "day")

# Expand to include all days
expanded_data <- expand.grid(biome_id = unique(merged_diet_data$biome_id), Date = all_days)
expanded_data <- left_join(expanded_data, merged_diet_data, by = c("biome_id", "Date"))

diet_data_subset <- merged_diet_data %>% 
  select(biome_id, Date, type, name, caloriesall, proteinall, sugarsall, fatall, food_other)

# get meals per day (only 3 meals max)
meals_per_day <- diet_data_subset %>%
  filter(type %in% c("breakfast", "lunch", "dinner")) %>% 
  group_by(biome_id, Date) %>%
  summarize(
    MealsPerDay = n_distinct(type),
    .groups = "drop"
  ) %>% 
  mutate(study_id=as.numeric(biome_id)) %>% 
  filter(!is.na(biome_id))

complete_grid <- expand.grid(biome_id = unique(meals_per_day$biome_id),
                             Date = all_days)
heatmap_data <- complete_grid %>%
  left_join(meals_per_day, by = c("biome_id", "Date"))

## All days plot
# heatmap_data$Value <- ifelse((is.na(heatmap_data$type)), "Not menstruating", "Menstruating")
# heatmap_data$Value <- factor(heatmap_data$Value, levels = c("Not menstruating", "Menstruating"))

heatmap_data_plot <- heatmap_data %>% 
  group_by(biome_id) %>%
  mutate(meals_count = sum(MealsPerDay, na.rm = TRUE)) %>%
  mutate(MealsPerDay=ifelse(is.na(MealsPerDay), 0, MealsPerDay)) %>% 
  mutate(MealsPerDay=as.factor(MealsPerDay)) %>% 
  arrange(desc(meals_count), biome_id) %>% 
  ungroup()

# all dates showing
ggplot(heatmap_data_plot, aes(x = as.factor(Date), y = reorder(factor(biome_id), meals_count), fill = MealsPerDay)) +
  geom_tile(color = "gray25") +
  scale_fill_manual(values = c("white", "lightblue", "blue", "darkblue"), 
                    na.value = "white")+
  labs(title = " ", 
       x = "Date", 
       y = "Biome ID", 
       fill = "Number of Meals") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

##########################################################################################
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
                   "vegetarian sausage patty", "vegetarian sausage patties", "beyond burger", "no meat", "vegetable sushi",
                   "Chipotle Black Bean Burger")
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
  select(name, vegetarian) %>%  
  distinct() %>%                
  arrange(name) 

## Check if participant vegetarian
participant_vegetarian_status <- merged_diet_data %>%
  group_by(biome_id) %>%
  # Check if all foods for a participant are vegetarian
  summarise(is_vegetarian = all(vegetarian))  
participant_vegetarian_status$biome_id <- as.numeric(participant_vegetarian_status$biome_id)
participant_vegetarian_status <- participant_vegetarian_status %>%
  drop_na()

table(participant_vegetarian_status$is_vegetarian) # 9 vegetarian

# Add to diet data
# participant_vegetarian_status$study_id <- as.character(participant_vegetarian_status$study_id)
merged_diet_data <- merged_diet_data %>%
  left_join(participant_vegetarian_status, by = "biome_id")

##########################################################################################
merged_diet_data_vegetarian <- merged_diet_data %>% 
  select(biome_id, is_vegetarian) %>% 
  distinct(biome_id, is_vegetarian)

# Correlate with vegetarian status and vaginal microbiome
vaginal.microbial.menses.24 <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifetyle/vaginal.microbial.menses.24.csv")
vaginal.microbial.menses.24 <- vaginal.microbial.menses.24[,-1]

vaginal.microbial.menses.24.veg <- vaginal.microbial.menses.24 %>% 
  left_join(merged_diet_data_vegetarian, by="biome_id")

vaginal.microbial.menses.24.summary <- vaginal.microbial.menses.24 %>% 
  group_by(biome_id) %>% 
  summarise(avg_shannon=sum(shannon)/n())

vaginal.microbial.menses.24.summary <- vaginal.microbial.menses.24.summary %>% 
  left_join(merged_diet_data_vegetarian, by="biome_id")
vaginal.microbial.menses.24.summary %>% 
  filter(is.na(is_vegetarian))

vaginal.microbial.menses.24.summary <- vaginal.microbial.menses.24.summary %>% 
  filter(!is.na(is_vegetarian)) %>% 
  mutate(is_vegetarian_bin=ifelse(is_vegetarian==TRUE, "Vegetarian", "Non-vegetarian")) 

table(vaginal.microbial.menses.24.summary$is_vegetarian_bin)

ggplot(vaginal.microbial.menses.24.summary, aes(x = is_vegetarian_bin, y = avg_shannon)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color="orchid") +
  labs(
    x = " ", 
    y = "Average Shannon Diversity", 
    title = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none") 

# testing
t.test(avg_shannon ~ is_vegetarian, data=vaginal.microbial.menses.24.summary)
wilcox.test(avg_shannon ~ is_vegetarian, data = vaginal.microbial.menses.24.summary)

##########################################################################################

# For a participant on a given day, get cals that are vegetarian v. not vegetarian
# collapse data into foods for that participant by day and the percent of food that is vegetarian

# vegetarian_percentage <- merged_diet_data %>% 
#   mutate(study_id=as.numeric(study_id)) %>% 
#   filter(!is.na(study_id)) %>% 
#   group_by(study_id, Date) %>% 
#   summarise(full_day_cal=sum(caloriesall, na.rm = TRUE),
#             vegetarian_perc=sum(caloriesall[is_vegetarian==TRUE], na.rm=TRUE)) %>% 
#   ungroup()
# 
# View(vegetarian_percentage)

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


