########################
#
# Exploration of the Activity Data
# Last updated: 03/18/2025
#
#########################

library(tidyverse)

source("functions.R")

activity.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 4-Physical Activity.csv")
activity.data <- activity.data %>% 
  filter(!is.na(as.numeric(biome_id)))
activity.data <- study_days(activity.data)

##########################################################################################

length(unique(activity.data$biome_id))

activity.data.collapsed <- activity.data %>% 
  group_by(biome_id) %>% 
  summarise(avg_steps = sum(steps) / n()) 

head(activity.data.collapsed, 7)

# Activity: distribution of average steps throughout study
ggplot(activity.data.collapsed, aes(x=factor(biome_id), y=avg_steps)) +
  geom_bar(stat="identity", fill ="orchid", color="black") +
  theme_minimal() +
  labs(x = "Participant ID", y = "Average Steps", title = " ") +
  theme(legend.position = "none") 

summary(activity.data.collapsed$avg_steps)

dim(activity.data) # 2532   13
dim(activity.data %>% 
      group_by(biome_id, Date) %>% 
      distinct())

# filter for study days
activity.data <- study_days(activity.data)
dim(activity.data)

# filter out 0 step count data
summary(activity.data$steps)
activity.data <- activity.data %>%
  filter(steps != 0)

#### Correlate with shannon data
# vaginal.microbial.menses.24 <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifetyle/vaginal.microbial.menses.24.csv")
vaginal.microbial.menses.24 <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/relabeled_data/vaginal.microbial.menses.24.csv")
vaginal.microbial.menses.24 <- vaginal.microbial.menses.24[,-1]
dim(vaginal.microbial.menses.24)

# join by date and participant
activity.data.shannon <- vaginal.microbial.menses.24 %>% 
  left_join(activity.data, by=c("biome_id", "logDate"))
dim(activity.data.shannon)
sum(!is.na(activity.data.shannon$study_day))
sum(is.na(activity.data.shannon$study_day))

activity.data.shannon <- activity.data.shannon %>% 
  filter(!is.na(study_day))
length(unique(activity.data.shannon$biome_id))

# 449 activity days not matched with swabs before filter 0 step count days
# 1122 activity days do match before filtering

# after filtering for 0
# 476 activity not matched with swabs
# 1095 match with swabs

activity.data.collapsed <- activity.data %>% 
  group_by(biome_id) %>% 
  summarise(avg_steps = sum(steps) / n()) 

# Activity: filtered distribution of average steps throughout study
ggplot(activity.data.collapsed, aes(x=factor(biome_id), y=avg_steps)) +
  geom_bar(stat="identity", fill ="orchid", color="black") +
  theme_minimal() +
  labs(x = "Participant ID", y = "Average Steps", title = " ") +
  theme(legend.position = "none") 

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

### Summarize activity data
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

## plot relationship between preds and avg shannon
plot(activity.data.summary$avg_cals_burned, activity.data.summary$avg_shannon)
plot(activity.data.summary$avg_steps, activity.data.summary$avg_shannon)
plot(activity.data.summary$avg_distance, activity.data.summary$avg_shannon)
plot(activity.data.summary$avg_minutes_sedentary, activity.data.summary$avg_shannon)
plot(activity.data.summary$avg_minutes_lightly_active, activity.data.summary$avg_shannon)
plot(activity.data.summary$avg_minutes_fairly_active, activity.data.summary$avg_shannon)
plot(activity.data.summary$avg_minues_very_active, activity.data.summary$avg_shannon)
plot(activity.data.summary$avg_activity_calories, activity.data.summary$avg_shannon)
plot(activity.data.summary$total_min_active, activity.data.summary$avg_shannon)



## simple linear regression
library(car)

lm.shannon.activity <- lm(avg_shannon ~ avg_steps+#avg_cals_burned+
                            #avg_distance+ # alias for avg_steps
                            avg_minutes_sedentary+
                            avg_minutes_lightly_active+
                            avg_minutes_fairly_active+
                            avg_minues_very_active #+avg_activity_calories+
                            #total_min_active
                          , data = activity.data.summary)
#pairs(activity.data.summary[,-c(1,2)])
# cor(activity.data.summary[,-c(1,2)])
summary(lm.shannon.activity)


vif(lm.shannon.activity)

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

## Mixed Effects Model
library(lme4)
library(lmerTest)
library(performance)
# activity.data.shannon.filtered <- activity.data.shannon %>% 
#   group_by(biome_id) %>% 
#   filter(n() > 10)
# table(activity.data.shannon.filtered$biome_id)

## Different baseline, same slope, different intercept

source("~/Microbiome Thesis/functions.R")

fixed_list <- c("steps", "distance", "calories_burned", "minutes_sedentary", "minutes_lightly_active",
                "minutes_fairly_active", "minues_very_active", "activity_calories")

lmer.models <- mixed_effects_fixed_slope(activity.data.shannon, 10, "shannon", fixed_list)

r2(lmer.models$steps)
r2(lmer.models$distance)
r2(lmer.models$calories_burned)
r2(lmer.models$minutes_sedentary)
r2(lmer.models$minutes_lightly_active)
r2(lmer.models$minutes_fairly_active)
r2(lmer.models$minues_very_active)
r2(lmer.models$activity_calories)

summary(lmer.models$steps)
summary(lmer.models$distance)
summary(lmer.models$calories_burned)
summary(lmer.models$minutes_sedentary)
summary(lmer.models$minutes_lightly_active)
summary(lmer.models$minutes_fairly_active)
summary(lmer.models$minues_very_active)
summary(lmer.models$activity_calories)

## random intercepts, random slopes
rnd.slope.lmer.models <- mixed_effects_rnd_slope(activity.data.shannon, 10, "shannon", fixed_list)

# r2(rnd.slope.lmer.models$steps)
r2(rnd.slope.lmer.models$distance)
r2(rnd.slope.lmer.models$calories_burned)
r2(rnd.slope.lmer.models$minutes_sedentary)
r2(rnd.slope.lmer.models$minutes_lightly_active)
r2(rnd.slope.lmer.models$minutes_fairly_active)
r2(rnd.slope.lmer.models$minues_very_active)
r2(rnd.slope.lmer.models$activity_calories)

summary(rnd.slope.lmer.models$steps)
summary(rnd.slope.lmer.models$distance)
summary(rnd.slope.lmer.models$calories_burned)
summary(rnd.slope.lmer.models$minutes_sedentary)
summary(rnd.slope.lmer.models$minutes_lightly_active)
summary(rnd.slope.lmer.models$minutes_fairly_active)
summary(rnd.slope.lmer.models$minues_very_active)
summary(rnd.slope.lmer.models$activity_calories)

# full model - fixed slopes, rnd intercepts
lmer.full <- lmer(shannon~steps + # distance + calories_burned + 
                  minutes_sedentary + minutes_lightly_active +
                    minutes_fairly_active + minues_very_active + #activity_calories +
                    (1|`biome_id`), 
                  data = activity.data.shannon)
r2(lmer.full)
summary(lmer.full)

# full model - random slopes, rnd intercepts
rnd.slope.lmer.full <- lmer(shannon~ steps + (steps||`biome_id`) + 
                              #distance + (distance||`biome_id`) +
                    #calories_burned + (calories_burned||`biome_id`) +
                      # activity_calories + (activity_calories||`biome_id`)
                      minutes_sedentary + (minutes_sedentary||`biome_id`) +
                      minutes_lightly_active + (minutes_lightly_active||`biome_id`) +
                    minutes_fairly_active + (minutes_fairly_active||`biome_id`) +
                      minues_very_active + (minues_very_active||`biome_id`),
                  data = activity.data.shannon)
r2(rnd.slope.lmer.full)

anova(lmer.full, rnd.slope.lmer.full)

##########################################################################################

## Time
activity.time.model <- lmer(shannon ~ study_day + steps + distance + activity_calories +
                              #calories_burned + 
                              minutes_sedentary +
                               minutes_lightly_active + minutes_fairly_active +
                               minues_very_active + 
                               + (1 | biome_id),
                             data = activity.data.shannon)
r2(activity.time.model)
summary(activity.time.model)
anova(lmer.full, activity.time.model)

## Time and Time^2
activity.time2.model <- lmer(shannon ~ study_day + I(study_day^2) + steps +
                              #distance + calories_burned + activity_calories +
                               minutes_sedentary +
                              minutes_lightly_active + minutes_fairly_active +
                              minues_very_active + 
                              + (1 | biome_id),
                            data = activity.data.shannon)
r2(activity.time2.model)
anova(activity.time.model, activity.time2.model)
anova(lmer.full, activity.time2.model)

##########################################################################################

## check if student athlete matters
participant.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 9-Volunteer Medical History.csv", header = TRUE)

participant.data <- participant.data %>% 
  mutate(sport_collapsed = ifelse(sport=="In-Season", sport, "Not In-Season")) %>%
  # mutate(sport_collapsed = ifelse(sport=="In-Season" | sport == "Club", "Activly in Sport", "Not In-Season")) %>%
  select(biome_id, logDate, activity_level, sport, sport_collapsed, field_hockey)

activity.sport.summary <- activity.data.summary %>% 
  left_join(participant.data, by="biome_id") %>% 
  filter(!is.na(sport))

table(activity.sport.summary$sport)
sum(is.na(activity.sport.summary$sport))
length(activity.sport.summary$sport)

# Activity: avg shannon and type of sport
ggplot(activity.sport.summary, aes(x = sport, y = avg_shannon)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.9, color="orchid") +
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

table(activity.sport.summary$sport_collapsed)

# Activity: avg shannon on Sports collapsed
ggplot(activity.sport.summary, aes(x = sport_collapsed, y = avg_shannon)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.9, color="orchid") +
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

## field hockey team
activity.sport.summary2 <- activity.sport.summary %>% 
  mutate(field_hockey=ifelse(field_hockey==TRUE, "fieldHockey", 
                             ifelse(sport_collapsed=="In-Season",
                                    "In-Season", "Off-Season")))

table(activity.sport.summary2$field_hockey)

ggplot(activity.sport.summary2, aes(x = field_hockey, y = avg_shannon)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color="orchid") +
  labs(
    x = " ", 
    y = "Average Shannon Diversity", 
    title = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none") 

t.test(avg_shannon ~ field_hockey, data=activity.sport.summary)
wilcox.test(avg_shannon ~ field_hockey, data = activity.sport.summary)

##########################################################################################

# levels of activity - skipped in writing because more arbitrary analysis

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
# Investigate intense exercise (greater than 2)
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

##########################################################################################

## mixed effects models
activity.sport.participant <- activity.data.shannon %>% 
  left_join(participant.data, by="biome_id") %>% 
  filter(!is.na(sport))

fixed_list <- c("steps*field_hockey", "distance*field_hockey", "calories_burned*field_hockey", 
                "minutes_sedentary*field_hockey", "minutes_lightly_active*field_hockey",
                "minutes_fairly_active*field_hockey", "minues_very_active*field_hockey", "activity_calories*field_hockey")

lmer.models <- mixed_effects_fixed_slope(activity.sport.participant, 10, "shannon", fixed_list)

summary(lmer.models$`steps*field_hockey`)
summary(lmer.models$`distance*field_hockey`)
summary(lmer.models$`calories_burned*field_hockey`)
summary(lmer.models$`minutes_sedentary*field_hockey`)
summary(lmer.models$`minutes_lightly_active*field_hockey`)
summary(lmer.models$`minutes_fairly_active*field_hockey`)
summary(lmer.models$`minues_very_active*field_hockey`)
summary(lmer.models$`activity_calories*field_hockey`)

##########################################################################################

head(activity.data.shannon)

