library(tidyverse)
library(purrr)
library(zoo)

source("~/Microbiome Thesis/functions.R")

full_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_menstruation_data.csv", header=TRUE)
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
  left_join(collapsed_days, by="biome_id") %>% 
  filter(biome_id %in% participants_to_impute) %>% 
  mutate(regular_periods=ifelse(is.na(regular_periods), 0, regular_periods))

heatmap_plot2(imputation_data)

####################################################################
# Imputation pt. 3: regular cycle participants

regular_cycle_df <- imputation_data %>% 
  filter(regular_periods==1)

irregular_cycle_df <- imputation_data %>% 
  filter(regular_periods!=1)
dim(regular_cycle_df)

heatmap_plot2(regular_cycle_df)

# 1. Get D1 of p0
logDates <- names(regular_cycle_df) %>% 
  grep("^2022-", ., value = TRUE)

imputation_regular_cycle_df <- regular_cycle_df %>%
  mutate(day1 = apply(select(., all_of(logDates)), 1, function(row) {
    matched_idx <- which(row %in% c(1, 2, 3, 7, 9))
    if (length(matched_idx) > 0) {
      return(logDates[matched_idx[1]])
    } else {
      return(NA)
    }
  })) %>% 
  mutate(period_length=sapply(str_extract_all(regular_cycle_df$period_len, "\\d+"), function(nums) max(as.numeric(nums))))

imputation_regular_cycle_df <- imputation_regular_cycle_df %>% 
  mutate(menstruation_length=ifelse((period_length < 10 & period_length > 0), period_length, NA),
         full_cycle_length=ifelse(period_length >= 10, period_length, NA))

# Reshape to long format
expanded_data_long <- imputation_regular_cycle_df %>%
  pivot_longer(cols = starts_with("2022-"),
               names_to = "logDate", 
               values_to = "person_menstruating") %>%
  arrange(biome_id, logDate)

all_days <- seq.Date(as.Date(min(full_data$logDate, na.rm=TRUE)),
                     as.Date(max(full_data$logDate, na.rm=TRUE)), by = "day")
all_days <- as.character(all_days)

# Expand to include all days
expanded_data <- expand.grid(biome_id = unique(expanded_data_long$biome_id), logDate = all_days)
expanded_data <- left_join(expanded_data, expanded_data_long, by = c("biome_id", "logDate"))

# Fill in cycles
imputed_filtered_data <- expanded_data %>%
  group_by(biome_id) %>%
  mutate(
    person_menstruating = ifelse(is.na(person_menstruating), 0, person_menstruating),
    imputation_menses = person_menstruating,
    # prev menstruation date
    prev_menstruating_date = ifelse(person_menstruating == 1, lag(logDate), NA),
    prev_menstruating_date = zoo::na.locf(prev_menstruating_date, na.rm = FALSE),
    
    next_menstruating_date = sapply(1:n(), function(i) {
      # Check if there's a `1` later in the vector
      next_index <- which(person_menstruating[i:n()] == 1)[1] + (i - 1)  # + (i - 1) to correct the offset from the slice
      if (length(next_index) > 0) logDate[next_index] else NA
    }),
    next_menstruating_date = zoo::na.locf(next_menstruating_date, na.rm = FALSE),
    next_menstruating_date = ifelse(next_menstruating_date > logDate, next_menstruating_date, NA),
    
    # Calculate gaps in days
    gap_between_menstruations = as.numeric(difftime(next_menstruating_date, prev_menstruating_date, units = "days")),
    day_gap = round(as.numeric(difftime(logDate, prev_menstruating_date, units = "days")), 4),
    
    # Impute menstruation if gap between menses is less than 10 days, person_menstruating == 0, and next_menstruating_date is not NA
    imputation_menses = ifelse(
      gap_between_menstruations > 0 & gap_between_menstruations < 10 & person_menstruating == 0 & !is.na(prev_menstruating_date) & !is.na(next_menstruating_date),
      1, imputation_menses
    )
  ) %>%
  ungroup()

# Pivot wide
# imputed_filtered_data_wide <- imputed_filtered_data %>%
#   pivot_wider(names_from = logDate, values_from = imputation_menses)

# Plot of imputed data colored with imputation
full_data_subset <- imputed_filtered_data %>% 
  select(biome_id, logDate, imputation_menses, uMinn_menstruation, survey_menstruate, person_menstruating)

complete_grid <- expand.grid(biome_id = unique(full_data_subset$biome_id),
                             logDate = all_days)
heatmap_data <- complete_grid %>%
  left_join(full_data_subset, by = c("biome_id", "logDate")) %>% 
  mutate(imputation_identifier=case_when(
    (imputation_menses==1 & person_menstruating==0) ~ 1, # imputed data point
    # person_menstruating==1 ~ 2,
    (uMinn_menstruation==1 & menstruation==1) ~ 2,
    (uMinn_menstruation==1 & (menstruation==0 | is.na(menstruation))) ~ 3,
    ((uMinn_menstruation==0 | is.na(uMinn_menstruation)) & menstruation==1) ~ 4,
    TRUE ~ NA_real_
  ))

heatmap_data$Value <- case_when(
  heatmap_data$imputation_identifier == 1 ~ "Imputed",
  heatmap_data$imputation_identifier == 2 ~ "Both sources agree",
  heatmap_data$imputation_identifier == 3 ~ "UMinn sample",
  heatmap_data$imputation_identifier == 4 ~ "Self-report",
  TRUE ~ "Missing"
)
heatmap_data$Value <- factor(heatmap_data$Value, levels = c("Missing", "Imputed", "Both sources agree", "UMinn sample", "Self-report"))
heatmap_data <- heatmap_data %>% 
  group_by(biome_id) %>%
  mutate(menstruating_days_count = sum((Value == "Imputed" | Value == "Both sources agree" | Value == "Self-report"), na.rm = TRUE)) %>%
  arrange(desc(menstruating_days_count), biome_id) %>% 
  ungroup()
ggplot(heatmap_data, aes(x = logDate, y = reorder(factor(biome_id), menstruating_days_count), fill = Value)) +
  geom_tile(color = "gray25") +
  scale_fill_manual(values = c("No Menstruation" = "white", 
                               "Imputed" = "lightpink", 
                               "Both sources agree" = "darkblue", 
                               "UMinn sample" = "lightblue", 
                               "Self-report" = "purple"),
                    na.value = "white") +
  labs(title = "Menstruation Heatmap by Identifier", 
       x = "Date", 
       y = "Biome ID", 
       fill = "Menstruation Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
