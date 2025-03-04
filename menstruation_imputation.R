library(tidyverse)
library(purrr)
library(zoo)

source("~/Microbiome Thesis/functions.R")

# RELABELED DATA
full_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/relabeled_data/cleaned_menstruation_data.csv", header=TRUE)

# full_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_menstruation_data.csv", header=TRUE)
survey_data_full <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 9-Volunteer Medical History.csv", header=TRUE)

# filter data
full_data <- filter_days(full_data)

survey_data <- survey_data_full %>% 
  select(biome_id,
         # activity_level, student_athelete, taken_antibiotics, meds, dietary_habits, vag_infection, vag_discomfort,
         survey_menstruate, period_last_day, regular_periods, irreg_period_notes, period_len,
         menstrual_prod, vaginal_notes, menstrual_cup, tampon, pad, no_menstrual_product) %>% 
  filter(!is.na(biome_id))

## Inconsistency levels (imputation pt. 1)
full_data <- full_data %>% 
  mutate(menstruation_status=ifelse((biome_id==56 & logDate=="2022-10-19"), 9, menstruation_status))

### Menses or no col
full_data <- full_data %>% 
  mutate(menses_decision=ifelse(menstruation_status %in% c(1,2,3,7,9), 1, 
                                ifelse(menstruation_status %in% c(4,5,10), 0, 
                                       # no blood, vaginal sample
                                       ifelse(menstruation_status==6, 99, NA))))

### Widen data and add cols for each day (0 no menses, 1 menses)
all_days <- seq.Date(as.Date(min(full_data$logDate, na.rm=TRUE)),
                     as.Date(max(full_data$logDate, na.rm=TRUE)), by = "day")
all_days <- as.character(all_days)

collapsed_days <- full_data %>%
  mutate(logdate = as.character(logDate)) %>%
  pivot_wider(
    id_cols = biome_id,
    names_from = logDate,
    values_from = menses_decision,
    values_fill = NA
  )

# Add back missing days (no data for them)
missing_days <- setdiff(all_days, colnames(collapsed_days))
if (length(missing_days) > 0) {
  collapsed_data[missing_days] <- NA
}

# Menstruation Data
menses_data <- survey_data %>% 
  left_join(collapsed_days, by="biome_id")

## Plot Initial Data
# heatmap_plot(menses_data)
heatmap_plot1(menses_data)

####################################################################
# Mismatch IDs - those who say they don't menstruate

# Make menstruation variable
full_data <- full_data %>% 
  mutate(person_menstruating=ifelse((menstruation_status==1 | menstruation_status==2 | menstruation_status==3 
                                     | menstruation_status==7), 1, 
                                    ifelse(menstruation_status==8, -1, 0)))
menstruate_df <- full_data %>% 
  left_join(survey_data, by="biome_id") %>%
  filter(person_menstruating==1) 

mismatch_menstruate_df <- menstruate_df %>%
  filter(((menstruation==1 | uMinn_menstruation==1) & survey_menstruate == 0))
mismatch_ids <- unique(mismatch_menstruate_df$biome_id)
mismatch_menstruate_df <- menses_data %>%
  filter(biome_id %in% mismatch_ids) %>%
  filter(!is.na(biome_id))
length(unique(mismatch_menstruate_df$biome_id))

# Fig 2: Heatmap of Biome ID by Date - self-report don't menstruate
heatmap_plot1(mismatch_menstruate_df)

# Check Volunteer Surveys
survey_data_mismatch <- survey_data %>% 
  filter(biome_id %in% mismatch_ids)

## Everyone who says they don't menstruate
person_says_no_menses <- survey_data %>% 
  filter(survey_menstruate==0) %>% 
  left_join(collapsed_days, by="biome_id")

# Fig 3: Heatmap of Biome ID by Date - self-report don't menstruate
heatmap_plot1(person_says_no_menses)

## Everyone who says they do menstruate
person_says_menses <- survey_data %>% 
  filter(survey_menstruate==1) %>% 
  left_join(collapsed_days, by="biome_id")

person_says_menses_data <- person_says_menses %>% 
  # filter if they have data
  filter(!(if_all(starts_with("2022-"), ~is.na(.)))) %>% 
  # filter for if they have menses data
  filter(if_any(starts_with("2022-"), ~. == 1) )
length(unique(person_says_menses_data$biome_id))

# Fig 4: Heatmap of Biome ID by Date - self-report do menstruate
heatmap_plot1(person_says_menses_data)

## Participants to Impute
participants_to_impute <- c(unique(person_says_menses_data$biome_id), unique(mismatch_menstruate_df$biome_id))
participants_to_impute_df <- menses_data %>% 
  filter(biome_id %in% participants_to_impute)

# Fig 5: Heatmap of Biome ID by Date - participants who have menses data to work with
heatmap_plot1(participants_to_impute_df)

####################################################################
# Imputation pt. 2: Check Volunteer surveys for cross-check menses (hard code)
survey_data_imputation <- survey_data %>% 
  filter(biome_id %in% participants_to_impute)

### Widen data and add cols for each day
all_days <- seq.Date(as.Date(min(full_data$logDate, na.rm=TRUE)),
                     as.Date(max(full_data$logDate, na.rm=TRUE)), by = "day")
all_days <- as.character(all_days)

collapsed_days <- full_data %>%
  mutate(logdate = as.character(logDate)) %>%
  pivot_wider(
    id_cols = biome_id,
    names_from = logDate,
    values_from = menstruation_status,
    values_fill = NA
  )

# Menstruation Data
imputation_data <- survey_data %>% 
  left_join(collapsed_days, by="biome_id") #%>% 
  # get only ppl that menstruated throughout study
  # filter(biome_id %in% participants_to_impute) #%>%
  # mutate(regular_periods=ifelse(is.na(regular_periods), 0, regular_periods))

heatmap_plot2(imputation_data)

####################################################################
# Imputation pt. 3: regular cycle participants

regular_cycle_df <- imputation_data %>% 
  filter(regular_periods==1)

irregular_cycle_df <- imputation_data %>% 
  filter(regular_periods!=1)

noinfo_cycle_df <- imputation_data %>% 
  filter(is.na(regular_periods))

dim(regular_cycle_df)
dim(irregular_cycle_df)
dim(noinfo_cycle_df)

heatmap_plot2(regular_cycle_df)

heatmap_plot2(irregular_cycle_df)

heatmap_plot2(noinfo_cycle_df)

####################################################################
# Imputation pt. 4: add fecal samples

fecal.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_fecal_samples_data.csv")

# If there is a fecal sample and a number already, no change
# no number and fecal sample, 50
imputation_data_fecal <- imputation_data %>% 
  rowwise() %>% 
  mutate(across(starts_with("2022-"),
    ~ ifelse(is.na(.) & biome_id %in% fecal.data$biome_id & 
                cur_column() %in% fecal.data$logDate[fecal.data$biome_id == biome_id & fecal.data$submitted],
             50, .)
  ))

regular_cycle_df <- imputation_data_fecal %>% 
  filter(regular_periods==1)

irregular_cycle_df <- imputation_data_fecal %>% 
  filter(regular_periods!=1)

ininfo_cycle_df <- imputation_data_fecal %>% 
  filter(is.na(regular_periods))

heatmap_plot_fecal(regular_cycle_df)

heatmap_plot_fecal(irregular_cycle_df)

heatmap_plot_fecal(ininfo_cycle_df)

####################################################################
# Imputation pt. 5: Add imputed dates (fecal samples not represented here)

# volunteer.df_subset <- volunteer.df %>% 
#   select(biome_id, birthControl, survey_menstruate, regular_periods, irreg_period_notes, period_len, study_menstruate, sexuallyActive)
# View(volunteer.df_subset %>% filter(biome_id==48))

# Regular (self-reported); person 66, 16, 14, 58, 49, 55, 51, 29, 15, 52, 37, 9 can't really impute
# person 50, 17, 4 no impute needed
person65_days <- c("2022-10-15", "2022-10-16", "2022-10-18", "2022-10-19", "2022-11-11", "2022-11-14", "2022-11-14", "2022-11-16", "2022-11-17", "2022-11-18")
person64_days <- c("2022-10-21", "2022-10-23", "2022-11-17", "2022-11-18")
person60_days <- c("2022-10-20", "2022-10-21", "2022-10-22", "2022-10-23", "2022-11-15", "2022-11-16", "2022-11-17")
person63_days <- c("2022-10-25")
person66_days <- c("2022-10-21", "2022-10-22", "2022-10-23")
person43_days <- c("2022-12-01", "2022-12-02", "2022-12-04", "2022-12-05", "2022-12-06")
person38_days <- c("2022-11-02", "2022-11-03", "2022-11-05", "2022-11-06")
person12_days <- c("2022-10-22")
person35_days <- c("2022-11-16", "2022-11-18")
person25_days <- c("2022-11-01", "2022-11-03")
person41_days <- c("2022-10-28", "2022-11-24", "2022-11-25", "2022-11-26", "2022-11-27", "2022-11-28") # actually imputed, not just fill in cycles
person39_days <- c("2022-11-05")
person23_days <- c("2022-10-21", "2022-10-22", "2022-10-23", "2022-10-24", "2022-10-25", "2022-10-26") # actually imputed, not just fill in cycles
person45_days <- c("2022-10-27") # impute since next to a self-report day
person27_days <- c("2022-10-21")
person26_days <- c("2022-10-24") # impute since next to a self-report day

# Irregular (self-reported); person 33, 2, 54, 59, 47, 11, 56, 44, 42, 3, 40, 36, 13 can't really impute
# person 53, 10 no impute needed
person48_days <- c("2022-10-25", "2022-10-27", "2022-10-28", "2022-11-28", "2022-12-01", "2022-11-29", "2022-12-02", "2022-12-03", "2022-12-04", "2022-12-06", "2022-12-02")
person1_days <- c("2022-10-24", "2022-10-25", "2022-10-28")
person7_days <- c("2022-10-31")
person61_days <- c("2022-10-31")
person34_days <- c("2022-10-18", "2022-10-20", "2022-10-20", "2022-10-21", "2022-10-22", "2022-10-23", "2022-10-25")
person32_days <- c("2022-11-12", "2022-11-13")

imputation_data_v2 <- imputation_data %>% 
  mutate(across(all_of(person65_days),
                ~ ifelse(biome_id==65, 78, .) )) %>% 
  mutate(across(all_of(person64_days),
                ~ ifelse(biome_id==64, 78, .) )) %>% 
  mutate(across(all_of(person66_days),
                ~ ifelse(biome_id==66, 78, .) )) %>% 
  mutate(across(all_of(person60_days),
                ~ ifelse(biome_id==60, 78, .) )) %>% 
  mutate(across(all_of(person63_days),
                ~ ifelse(biome_id==63, 78, .) )) %>% 
  mutate(across(all_of(person43_days),
                ~ ifelse(biome_id==43, 78, .) )) %>% 
  mutate(across(all_of(person38_days),
                ~ ifelse(biome_id==38, 78, .) )) %>% 
  mutate(across(all_of(person12_days),
                ~ ifelse(biome_id==12, 78, .) )) %>% 
  mutate(across(all_of(person35_days),
                ~ ifelse(biome_id==35, 78, .) )) %>% 
  mutate(across(all_of(person25_days),
                ~ ifelse(biome_id==25, 78, .) )) %>% 
  mutate(across(all_of(person41_days),
                ~ ifelse(biome_id==41, 78, .) )) %>% 
  mutate(across(all_of(person39_days),
                ~ ifelse(biome_id==39, 78, .) )) %>% 
  mutate(across(all_of(person23_days),
                ~ ifelse(biome_id==23, 78, .) )) %>% 
  mutate(across(all_of(person45_days),
                ~ ifelse(biome_id==45, 78, .) )) %>% 
  mutate(across(all_of(person27_days),
                ~ ifelse(biome_id==27, 78, .) )) %>% 
  mutate(across(all_of(person26_days),
                ~ ifelse(biome_id==26, 78, .) )) %>% 
  # irregular cycles
  mutate(across(all_of(person48_days),
                ~ ifelse(biome_id==48, 78, .) )) %>% 
  mutate(across(all_of(person1_days),
                ~ ifelse(biome_id==1, 78, .) )) %>% 
  mutate(across(all_of(person7_days),
                ~ ifelse(biome_id==7, 78, .) )) %>% 
  mutate(across(all_of(person61_days),
                ~ ifelse(biome_id==61, 78, .) )) %>% 
  mutate(across(all_of(person34_days),
                ~ ifelse(biome_id==34, 78, .) )) %>% 
  mutate(across(all_of(person32_days),
                ~ ifelse(biome_id==32, 78, .) ))

heatmap_plot_imputation(imputation_data_v2)

### Save final data output
# write.csv(imputation_data_v2,
#           file = "/Volumes/T7/microbiome_data/cleaned_data/imputed_menstruation_data_2_12.csv",
#           row.names = FALSE)

# Relabeled data
write.csv(imputation_data_v2,
          file = "/Volumes/T7/microbiome_data/cleaned_data/relabeled_data/imputed_menstruation_data_2_12.csv",
          row.names = FALSE)
