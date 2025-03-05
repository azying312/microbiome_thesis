library(tidyverse)

activity.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 4-Physical Activity.csv")

activity.data <- activity.data %>% 
  filter(!is.na(as.numeric(biome_id)))

# distribution of steps
ggplot(activity.data, aes(x=factor(biome_id), y=steps, fill="orchid")) +
  geom_bar(stat="identity") +
  theme_minimal() +
  labs(x = "Participant ID", y = "Steps", title = "Distribution of Steps") +
  theme(legend.position = "none") 

#### Correlate with shannon data
# vaginal.microbial.menses.24 <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifetyle/vaginal.microbial.menses.24.csv")
vaginal.microbial.menses.24 <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/relabeled_data/vaginal.microbial.menses.24.csv")
vaginal.microbial.menses.24 <- vaginal.microbial.menses.24[,-1]

# join by date and participant
activity.data.shannon <- vaginal.microbial.menses.24 %>% 
  left_join(activity.data, by=c("biome_id", "logDate"))

# steps
# lm.shannon.activity <- lm(shannon ~ steps, data = activity.data.shannon)
# summary(lm.shannon.activity)
# # total cals burned - significant
# lm.shannon.activity <- lm(shannon ~ calories_burned, data = activity.data.shannon)
# summary(lm.shannon.activity)
# # activity cals burned 
# lm.shannon.activity <- lm(shannon ~ activity_calories, data = activity.data.shannon)
# summary(lm.shannon.activity)
# # distance burned 
# lm.shannon.activity <- lm(shannon ~ distance, data = activity.data.shannon)
# summary(lm.shannon.activity)
# # minutes sedentary
# lm.shannon.activity <- lm(shannon ~ minutes_sedentary, data = activity.data.shannon)
# summary(lm.shannon.activity)
# # minutes lightly active
# lm.shannon.activity <- lm(shannon ~ minutes_lightly_active, data = activity.data.shannon)
# summary(lm.shannon.activity)
# # minutes fairly active
# lm.shannon.activity <- lm(shannon ~ minutes_fairly_active, data = activity.data.shannon)
# summary(lm.shannon.activity)
# # minutes very active
# lm.shannon.activity <- lm(shannon ~ minues_very_active, data = activity.data.shannon)
# summary(lm.shannon.activity)

### Summary of activity data
activity.data.summary <- activity.data %>% 
  group_by(biome_id) %>% 
  summarise(avg_cals_burned = sum(calories_burned) / n(),
            avg_steps = sum(steps) / n(),
            avg_distance = sum(distance) / n(),
            avg_minutes_sedentary = sum(minutes_sedentary) / n(),
            avg_minutes_lightly_active = sum(minutes_lightly_active) / n(),
            avg_minutes_fairly_active = sum(minutes_fairly_active) / n(),
            avg_minues_very_active = sum(minues_very_active) / n(),
            avg_activity_calories = sum(activity_calories) / n(),
            total_min_active = sum(minutes_lightly_active, minutes_fairly_active, minues_very_active) / n()
            )
vaginal.microbial.menses.24.summary <- vaginal.microbial.menses.24 %>% 
  group_by(biome_id) %>% 
  summarise(avg_shannon=sum(shannon)/n())

activity.data.summary <- vaginal.microbial.menses.24.summary %>% 
  left_join(activity.data.summary, by="biome_id")

# steps
lm.shannon.activity <- lm(avg_shannon ~ avg_steps, data = activity.data.summary)
summary(lm.shannon.activity)
# total cals burned - significant
lm.shannon.activity <- lm(avg_shannon ~ avg_cals_burned, data = activity.data.summary)
summary(lm.shannon.activity)
# activity cals burned 
lm.shannon.activity <- lm(avg_shannon ~ avg_activity_calories, data = activity.data.summary)
summary(lm.shannon.activity)
# distance burned 
lm.shannon.activity <- lm(avg_shannon ~ avg_distance, data = activity.data.summary)
summary(lm.shannon.activity)
# minutes sedentary
lm.shannon.activity <- lm(avg_shannon ~ avg_minutes_sedentary, data = activity.data.summary)
summary(lm.shannon.activity)
# minutes lightly active
lm.shannon.activity <- lm(avg_shannon ~ avg_minutes_lightly_active, data = activity.data.summary)
summary(lm.shannon.activity)
# minutes fairly active
lm.shannon.activity <- lm(avg_shannon ~ avg_minutes_fairly_active, data = activity.data.summary)
summary(lm.shannon.activity)
# minutes very active
lm.shannon.activity <- lm(avg_shannon ~ avg_minues_very_active, data = activity.data.summary)
summary(lm.shannon.activity)


##########################################################################################

## check if student athlete matters
participant.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 9-Volunteer Medical History.csv", header = TRUE)

participant.data <- participant.data %>% 
  mutate(sport_collapsed = ifelse(sport=="In-Season", sport, "Not In-Season")) %>%
  # mutate(sport_collapsed = ifelse(sport=="In-Season" | sport == "Club", "Activly in Sport", "Not In-Season")) %>%
  select(biome_id, logDate, activity_level, sport, sport_collapsed)

activity.sport.summary <- activity.data.summary %>% 
  left_join(participant.data, by="biome_id") %>% 
  filter(!is.na(sport))

ggplot(activity.sport.summary, aes(x = sport, y = avg_shannon)) +
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
aov_sport <- aov(avg_shannon ~ sport, data = activity.sport.summary)
summary(aov_sport)

## Sports collapsed
ggplot(activity.sport.summary, aes(x = sport_collapsed, y = avg_shannon)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color="orchid") +
  labs(
    x = " ", 
    y = "Average Shannon Diversity", 
    title = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none") 

t.test(avg_shannon ~ sport_collapsed, data=activity.sport.summary)
wilcox.test(avg_shannon ~ sport_collapsed, data = activity.sport.summary)

### Minutes active 
lm.shannon.activity <- lm(avg_shannon ~ total_min_active, data = activity.data.summary)
summary(lm.shannon.activity)

##########################################################################################

# levels of activity

activity.data.summary.levels <- activity.data.summary %>% 
  mutate(exercise_level = case_when(
    total_min_active <= quantile(total_min_active, 0.25, na.rm = TRUE) ~ "Low",
    total_min_active <= quantile(total_min_active, 0.50, na.rm = TRUE) ~ "Moderate",
    total_min_active <= quantile(total_min_active, 0.75, na.rm = TRUE) ~ "High",
    TRUE ~ "Very High"
  ))

activity.data.summary.levels <- activity.data.summary.levels %>% 
  left_join(participant.data, by="biome_id") %>% 
  filter(!is.na(sport))
activity.data.summary.levels$exercise_level <- factor(activity.data.summary.levels$exercise_level, levels = c("Low", "Moderate", "High", "Very High"), ordered = TRUE)

table(activity.data.summary.levels$exercise_level)
table(activity.data.summary.levels$exercise_level, activity.data.summary.levels$sport)

## Sports collapsed
ggplot(activity.data.summary.levels, aes(x = exercise_level, y = avg_shannon)) +
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
aov_sport <- aov(avg_shannon ~ exercise_level, data = activity.data.summary.levels)
summary(aov_sport)

##########################################################################################

participant.data.subset <- participant.data %>% 
  filter(activity_level > 2)

# activity.data.summary
table(participant.data.subset$activity_level)

activity.data.summary.subset <- activity.data.summary %>% 
  filter(biome_id %in% unique(participant.data.subset$biome_id)) %>% 
  left_join(participant.data.subset, by="biome_id")

ggplot(activity.data.summary.subset, aes(x = factor(activity_level), y = avg_shannon)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color="orchid") +
  labs(
    x = " ", 
    y = "Average Shannon Diversity", 
    title = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none") 


### Minutes active 
lm.shannon.activity <- lm(avg_shannon ~ total_min_active, data = activity.data.summary.subset)
summary(lm.shannon.activity)

