########################
#
# Clean Volunteer Medical Survey Data
# v1: 8 October 2024
#
#########################

library(tidyverse)
library(lubridate) # date transformation

history_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Report 9-Volunteer Medical History.csv")
id_mapping <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Original Study Mapping - Sheet3.csv", header = TRUE)

### Map uid to study id

# Get unique study and biome health pairings
study_and_u_id <- unique(id_mapping %>% 
                           select(STUDY.ID, Biome.Health.App.ID))
# Match and join columns
study_and_u_id <- study_and_u_id %>% 
  rename("study_id" = "STUDY.ID") %>% 
  rename("biome_id" = "Biome.Health.App.ID")
history_data <- history_data %>% 
  rename("biome_id" = "Your.Biome.Health.App.ID")
study_and_u_id$study_id <- as.character(study_and_u_id$study_id)

# Map ids
history_data <- history_data %>%
  left_join(study_and_u_id, by = "biome_id") %>%
  mutate(biome_id = coalesce(study_id, biome_id)) %>%
  select(-study_id)

# Check missing ids
missing_list <- history_data %>%
  filter(is.na(as.numeric(biome_id)))
print(unique(missing_list$biome_id))

#### Cleaning

## Repeat surveys
dim(history_data)
length(unique(history_data$biome_id))
which(history_data$biome_id==50)
history_data <- history_data[-68, ] # remove dupe

## Recode
activity_scale <- c("Sedentary (little to no exercise/physical activity)",
                    "Lightly active (1-3 hours of physical exercise per week)",
                    "Moderately active (4-6 hours of physical exercise per week)",
                    "Very active (7+ hours of physical exercise per week)")



history_data_recode <- history_data %>% 
  # Clean menstruate variable
  rename(menstruate=Do.you.menstruate...Choose.one.) %>% 
  # Date
  rename(logDate=Today.s.Date) %>% 
  mutate(logDate = mdy(logDate)) %>%
  # Ethnicity
  rename(ethnicity=How.would.you.describe.your.racial.ethnic.background...Indicate.all.that.apply.) %>% 
  # Sexuality
  rename(sexuality=How.would.you.describe.your.sexual.identity...Indicate.all.that.apply.) %>% 
  # Activity Level
  rename(activity_level=How.would.you.describe.your.general.activity.levels.on.average.) %>% 
  mutate(across(where(~ any(. %in% activity_scale)), ~ case_when(
    . == "Sedentary (little to no exercise/physical activity)" ~ 1,
    . == "Lightly active (1-3 hours of physical exercise per week)" ~ 2,
    . == "Moderately active (4-6 hours of physical exercise per week)" ~ 3,
    . == "Very active (7+ hours of physical exercise per week)" ~ 4,
    TRUE ~ NA_real_
  ))) %>% 
  # Student Athlete
  rename(student_athelete=Are.you.a.student.athlete.or.involved.in.club.sports..If.yes..please.indicate.what.activity.) %>% 
  mutate(student_athelete=ifelse(student_athelete=="", NA, student_athelete)) %>% # not working 
  # Taken Antibiotics, 0 Yes 1 NO
  rename(taken_antibiotics=Have.you.taken.any.antibiotics.within.the.last.month...Circle.one.) %>% 
  mutate(taken_antibiotics=ifelse(taken_antibiotics=="Yes", 0, 1)) %>% # not working 
  # Vegetarian
  rename("dietary_habits"="How.would.you.describe.your.dietary.habits...Choose.one.") %>% 
  mutate(dietary_habits = case_when(
    dietary_habits
      %in% c("omnivore but i don't eat red meat", "no red meat", "No red meat", 
             "Omnivore but don't eat red meat, beef (the only meat i eat is turkey, seafood, chicken)", 
             "omnivore but i don't eat red meat  ") ~ "No red meat",
    dietary_habits
      %in% c("Omnivore (Meat and Plant)", "Omnivore (Meat and Plant)") ~ "Omnivore",
    dietary_habits
      %in% c("Vegetarian and Gluten Free", "Vegetarian") ~ "Vegetarian",
    TRUE ~ dietary_habits
  )) %>% 
  # meds
  rename(meds=Are.you.currently.taking.any.hormonal.medications.and.or.birth.control...Indicate.all.that.apply.) %>%
  # regular periods
  rename(regular_periods=Do.you.have.regular.periods..every.21.to.35.days.and.last.2.to.7.days..or.irregular.periods...Choose.one.) %>% 
  mutate(regular_periods=ifelse(regular_periods=="Regular", 1, 
                                ifelse(regular_periods=="Irregular", 0, NA))) %>%
  # vagina bacterial infection
  rename(vag_infection=Have.you.ever.been.diagnosed.with.vaginal.bacterial.infection..vaginosis..) %>%
  mutate(vag_infection=ifelse(vag_infection=="Yes", 1,
                                ifelse(vag_infection=="60", 0, NA))) %>%
  # menstrual products
  rename(menstrual_prod=Which.menstrual.products.do.you.typically.use...Select.all.that.apply.) %>%
  mutate(
    menstrual_cup=ifelse(str_detect(menstrual_prod, "Menstrual Cup"), 1, 0),
    tampon=ifelse(str_detect(menstrual_prod, "Tampons"), 1, 0),
    pad=ifelse(str_detect(menstrual_prod, "Pads"), 1, 0),
    no_menstrual_product=ifelse(str_detect(menstrual_prod, "none"), 1, 0)
  ) %>% 
  # notes about period
  rename(vaginal_notes=If.you.would.like.to.report.any.other.condition.or.event.related.to.vaginal.health.and.or.menstruation.that.you.think.may.be.relevant.to.the.study..please.do.so.in.the.space.below.)

table(history_data_recode$dietary_habits)
history_data_recode$logDate

table(history_data_recode$dietary_habits)
length(history_data_recode$menstruate)
### Save final data output
# write.csv(history_data_recode,
#           file = "/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_Report 9-Volunteer Medical History.csv",
#           row.names = FALSE)

table(history_data_recode$ethnicity)
