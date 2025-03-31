########################
#
# Clean Volunteer Medical Survey Data
# v1: 8 October 2024
#
#########################

library(tidyverse)
library(lubridate) # date transformation

source("~/Microbiome Thesis/functions.R")

# Load data
history_data <- read.csv("/Volumes/T7/microbiome_data/original_data/Report 9-Volunteer Medical History.csv")
id_mapping <- read.csv("/Volumes/T7/microbiome_data/original_data/Original Study Mapping - Sheet3.csv", header = TRUE)

# Data Prep
history_data <- history_data %>% 
  rename("biome_id" = "Your.Biome.Health.App.ID")
history_data <- study_mapping(history_data, id_mapping)

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
  # Make logDate
  mutate(logDate=format(mdy_hms(Timestamp), "%Y-%m-%d")) %>% 
  # Clean menstruate variable
  rename(survey_menstruate=Do.you.menstruate...Choose.one.) %>% 
  mutate(survey_menstruate=ifelse(survey_menstruate=="Yes", 1, 0)) %>% 
  # Remove this col: Date
  select(!Today.s.Date) %>%
  # mutate(logDate = mdy(logDate)) %>%
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
  rename(sport=Are.you.a.student.athlete.or.involved.in.club.sports..If.yes..please.indicate.what.activity.) %>% 
  mutate(sport=tolower(ifelse(sport=="", "no", sport))) %>%
  mutate(athele=ifelse(sport=="no", 0, 1)) %>% 
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
  # birthControl
  rename(birthControl=Are.you.currently.taking.any.hormonal.medications.and.or.birth.control...Indicate.all.that.apply.) %>%
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
  rename(vaginal_notes=If.you.would.like.to.report.any.other.condition.or.event.related.to.vaginal.health.and.or.menstruation.that.you.think.may.be.relevant.to.the.study..please.do.so.in.the.space.below.) %>% 
  # Irregular period notes If.irregular..please.elaborate.on.the.frequency.
  rename(irreg_period_notes=If.irregular..please.elaborate.on.the.frequency.) %>% 
  rename(period_len=What.is.the.number.of.days.in.your.monthly.cycle..that.is..how.many.days.are.there.from.the.first.day.of.one.period.to.the.first.day.of.your.next.period...) %>% 
  rename(vag_discomfort=Are.you.currently.experiencing.vaginal.pain.or.discomfort...Choose.one.) %>% 
  rename(period_last_day=What.was.the.date.of.the.first.day.of.your.last.menstrual.period.) %>% 
  rename(cisWoman=How.would.you.describe.your.gender.identity...Indicate.all.that.apply.) %>% 
  mutate(cisWoman=ifelse(cisWoman == "Cisgender woman", 1, 0)) %>% 
  rename(probiotic=Are.you.currently.taking.a.probiotic.) %>% 
  mutate(probiotic=ifelse(probiotic=="Yes", 1, 0)) %>% 
  rename(sexuallyActive=Are.you.sexually.active...Choose.one.) %>% 
  mutate(sexuallyActive=ifelse(sexuallyActive=="Yes", 1, 0))

# Subset to cols of interest
history_data_subset <- history_data_recode %>% 
  select(biome_id, logDate, ethnicity, sexuality, cisWoman,
         activity_level, sport, taken_antibiotics, probiotic, birthControl,
         dietary_habits, survey_menstruate, period_last_day, regular_periods, irreg_period_notes, period_len,
         menstrual_prod, vaginal_notes, sexuallyActive, menstrual_cup, tampon, pad, no_menstrual_product)

### Add menstruation data from UMinn for menstruation variable
menses_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_menstruation_data.csv", header=TRUE)
# menstruation status 1,2,3,7,9
menses_data_collapsed <- menses_data %>% 
  group_by(biome_id) %>% 
  summarise(
    minn_menstruate = as.integer(any(menstruation_status %in% c(1, 2, 3, 7, 9)))
  )
history_data_subset <- history_data_subset %>% 
  left_join(menses_data_collapsed, by=as.character("biome_id"))
history_data_subset <- history_data_subset %>% 
  mutate(study_menstruate=ifelse(is.na(minn_menstruate), survey_menstruate, minn_menstruate)) %>% 
  select(!minn_menstruate)

## Dietary Habits Cleaning
diet_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/fully_merged_diet_data.csv")

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
                   "chipotle black bean burger")
# Create regex patterns for meat and exclusion lists
meat_pattern <- paste0("\\b(", paste(common_meats, collapse = "|"), ")\\b")
exclude_pattern <- paste0("\\b(", paste(exclude_foods, collapse = "|"), ")\\b")

# Check the common meats
diet_data <- diet_data %>%
  mutate(name_lower = str_to_lower(name),  # Convert names to lowercase
         vegetarian=TRUE,
         # Check for common meats first
         vegetarian = ifelse(str_detect(name_lower, meat_pattern), FALSE, vegetarian)) 
# Look at exclusions
diet_data <- diet_data %>%
  mutate(vegetarian = ifelse(str_detect(name_lower, exclude_pattern), TRUE, vegetarian)) %>% 
  select(-name_lower)
unique_vegetarian_pairs <- diet_data %>%
  select(name, vegetarian) %>%  
  distinct() %>%                
  arrange(name)

## Check if participant vegetarian
participant_vegetarian_status <- diet_data %>%
  group_by(biome_id) %>%
  # Check if all foods for a participant are vegetarian
  summarise(is_vegetarian = all(vegetarian)) %>% 
  # rename(biome_id=study_id) %>% 
  mutate(is_vegetarian=ifelse(is_vegetarian==TRUE, 1, 0))
participant_vegetarian_status$biome_id <- as.numeric(participant_vegetarian_status$biome_id)
participant_vegetarian_status <- participant_vegetarian_status %>%
  drop_na()
history_data_subset <- history_data_subset %>% 
  left_join(participant_vegetarian_status, by=as.character("biome_id")) %>% 
  mutate(dietary_habits=ifelse(is.na(is_vegetarian), dietary_habits, is_vegetarian)) %>% 
  select(!is_vegetarian) %>% 
  mutate(dietary_habits=ifelse(dietary_habits=="Vegetarian", 1, ifelse(dietary_habits=="Omnivore", 0, dietary_habits)))

## Birth Control
history_data_subset$birthControl[history_data_subset$birthControl=="" | 
                                   history_data_subset$birthControl=="None"] <- "None"
history_data_subset$birthControl[history_data_subset$birthControl=="Combination birth control pill"] <- "Systemic Combined (E&P)"
history_data_subset$birthControl[history_data_subset$birthControl=="Contraceptive implant (Nexplanon)" | 
                                   history_data_subset$birthControl=="Progesterone (Aygestin, provera, depo-provera, mini progestin-only birth control pill or norethindrone)"] <- "Systemic P only"
history_data_subset$birthControl[history_data_subset$birthControl=="Hormonal intrauterine device (Liletta, Skyla, Kyleena, or Mirena IUD)"
                     | history_data_subset$birthControl =="Vaginal ring"] <- "Local P"

unique(history_data_subset$birthControl) # Orilissa (Elagolix) is nonhormonal
table(history_data_subset$birthControl)

## Organize sport to be "off-season", "in-season", or "no sport" 

table(history_data_subset$activity_level)

history_data_subset <- history_data_subset %>% 
  mutate(field_hockey = ifelse(grepl("field hockey", history_data_subset$sport, ignore.case = TRUE), TRUE, FALSE))

history_data_subset$sport[grepl("frisbee", history_data_subset$sport, ignore.case = TRUE)] <- "Ultimate Frisbee"
history_data_subset$sport[grepl("climbing", history_data_subset$sport, ignore.case = TRUE) |
                            grepl("squash", history_data_subset$sport, ignore.case = TRUE) | 
                            grepl("Ultimate Frisbee", history_data_subset$sport, ignore.case = TRUE) |
                            grepl("rugby", history_data_subset$sport, ignore.case = TRUE) |
                            grepl("masters swim team at mit", history_data_subset$sport, ignore.case = TRUE) |
                            grepl("equestrian", history_data_subset$sport, ignore.case = TRUE) |
                            grepl("nordic ski", history_data_subset$sport, ignore.case = TRUE)] <- "Club"
history_data_subset$sport[grepl("tennis", history_data_subset$sport, ignore.case = TRUE) |
                            grepl("crew", history_data_subset$sport, ignore.case = TRUE) | 
                            grepl("rowing", history_data_subset$sport, ignore.case = TRUE) | 
                            grepl("basketball", history_data_subset$sport, ignore.case = TRUE) |
                            grepl("lacrosse", history_data_subset$sport, ignore.case = TRUE) |
                            # check if track was in season?
                            grepl("track", history_data_subset$sport, ignore.case = TRUE)] <- "Off-Season"
history_data_subset$sport[grepl("no", history_data_subset$sport, ignore.case = TRUE) |
                            grepl("pe class", history_data_subset$sport, ignore.case = TRUE)] <- "None"
history_data_subset$sport[grepl("field hockey", history_data_subset$sport, ignore.case = TRUE) |
                            grepl("cross country", history_data_subset$sport, ignore.case = TRUE) |
                            grepl("track", history_data_subset$sport, ignore.case = TRUE)] <- "In-Season"
# 9, 46
table(history_data_subset$sport)

### Save final data output
write.csv(history_data_subset,
          file = "/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 9-Volunteer Medical History.csv",
          row.names = FALSE)
# write.csv(history_data_recode,
#           file = "/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 9-Volunteer Medical History.csv",
#           row.names = FALSE)
