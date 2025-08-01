library(tidyverse)

diet_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/fully_merged_diet_data.csv")
names(diet_data)
dim(diet_data)

### (1) how many days do we have for each person out of all days
day_per_person <- diet_data %>%
  group_by(biome_id) %>%
  summarise(days_recorded = n_distinct(Date))

ggplot(day_per_person, aes(x=factor(biome_id), y=days_recorded, fill="orchid")) +
  geom_bar(stat="identity") +
  theme_minimal() +
  labs(x = "Participant ID", y = "Days Recorded", title = "Distribution of Days Recorded") +
  theme(legend.position = "none") 

# meal missingness
100 * sum(complete_days$meals_per_day) / (length(unique(complete_days$Date))*length(unique(complete_days$biome_id))*3)

# day missingness
100*(dim(complete_days)[1]) / (length(unique(complete_days$Date))*length(unique(complete_days$biome_id)))

### (2) how complete is the days (meals per day, calories per day, look at proportion of calories)
complete_days <- diet_data %>%
  group_by(biome_id, Date) %>%
  # TODO: meals per day & calories per day
  summarise(
    meals_per_day = n_distinct(type),
    total_calories = sum(caloriesall, na.rm = TRUE),
    .groups = 'drop'
  ) %>% 
  # TODO: look at proportion of calories
  mutate(expected_daily_cal = 2200, # online source said 2,000–2,400 calories
         proportion_daily_cal = total_calories / expected_daily_cal)

dim(complete_days) # 1383 complete days
# View(complete_days)

### (3) Completeness by person

complete_days_by_person <- complete_days %>%
  group_by(biome_id) %>%
  summarise(
    avg_meals_per_day = mean(meals_per_day, na.rm = TRUE),
    avg_calories_per_day = mean(total_calories, na.rm = TRUE),
    avg_proportion_of_calories = mean(proportion_daily_cal, na.rm = TRUE) 
  )

dim(complete_days_by_person) # 72 people
# View(complete_days_by_person)

### (4) Missingness for each col
sum(is.na(diet_data$caloriesall))
colSums(is.na(diet_data))

diet_data_notype <- diet_data %>% 
  filter(is.na(type))
# View(diet_data_notype)
