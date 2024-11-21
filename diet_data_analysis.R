library(tidyverse)
library(ggplot2)

diet_data <-read.csv("/Users/alicezhang/Desktop/microbiome_data/cleaned_data/manual_merged_diet_data.csv")

# heatmap for the days - same as menses one
diet_data <- diet_data %>% 
  mutate(Date=as.Date(Date)) %>% 
  filter(Date < as.Date("2022-12-17")) 

all_days <- seq.Date(as.Date(min(diet_data$Date, na.rm=TRUE)), 
                     as.Date(max(diet_data$Date, na.rm=TRUE)), by = "day")
# all_days <- as.character(all_days)

# Expand to include all days
expanded_data <- expand.grid(study_id = unique(diet_data$study_id), Date = all_days)
expanded_data <- left_join(expanded_data, diet_data, by = c("study_id", "Date"))

diet_data_subset <- diet_data %>% 
  select(study_id, Date, type, name, caloriesall, proteinall, sugarsall, fatall, food_other)

# get meals per day (only 3 meals max)
meals_per_day <- diet_data_subset %>%
  filter(type %in% c("breakfast", "lunch", "dinner")) %>% 
  group_by(study_id, Date) %>%
  summarize(
    MealsPerDay = n_distinct(type),
    .groups = "drop"
  ) %>% 
  mutate(study_id=as.numeric(study_id)) %>% 
  filter(!is.na(study_id))

complete_grid <- expand.grid(study_id = unique(meals_per_day$study_id),
                             Date = all_days)
heatmap_data <- complete_grid %>%
  left_join(meals_per_day, by = c("study_id", "Date"))

## All days plot
# heatmap_data$Value <- ifelse((is.na(heatmap_data$type)), "Not menstruating", "Menstruating")
# heatmap_data$Value <- factor(heatmap_data$Value, levels = c("Not menstruating", "Menstruating"))

heatmap_data_plot <- heatmap_data %>% 
  group_by(study_id) %>%
  mutate(meals_count = sum(MealsPerDay, na.rm = TRUE)) %>%
  mutate(MealsPerDay=ifelse(is.na(MealsPerDay), 0, MealsPerDay)) %>% 
  mutate(MealsPerDay=as.factor(MealsPerDay)) %>% 
  arrange(desc(meals_count), study_id) %>% 
  ungroup()

# all dates showing
ggplot(heatmap_data_plot, aes(x = as.factor(Date), y = reorder(factor(study_id), meals_count), fill = MealsPerDay)) +
  geom_tile(color = "gray25") +
  scale_fill_manual(values = c("white", "lightblue", "blue", "darkblue"), 
                    na.value = "white")+
  labs(title = "Number of Meals per person", 
       x = "Date", 
       y = "Biome ID", 
       fill = "Number of Meals") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

# some dates showing
ggplot(heatmap_data_plot, aes(x = (Date), y = reorder(factor(study_id), meals_count), fill = MealsPerDay)) +
  geom_tile(color = "gray25") +
  scale_fill_manual(values = c("white", "lightblue", "blue", "darkblue"), 
                    na.value = "white")+
  labs(title = "Number of Meals per person", 
       x = "Date", 
       y = "Biome ID", 
       fill = "Number of Meals") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

