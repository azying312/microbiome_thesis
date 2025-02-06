########################
#
# Menstruation Data Analysis
# 9 November 2024
#
#########################

library(tidyverse)

full_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_menstruation_data.csv", header=TRUE)
survey_data_full <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 9-Volunteer Medical History.csv", header=TRUE)

survey_data <- survey_data_full %>% 
  select(biome_id, survey_menstruate) #%>%
  # rename(survey_menstruate=menstruate) # 14 no, 54 yes

## Q: What proportion of self-reported days had blood; days that are inUminn & in self report
selfReport_blood <- full_data %>% 
  filter(inSelfReport==TRUE & uMinn_menstruation==1)
dim(selfReport_blood) # 42 days
sum(full_data$inSelfReport==TRUE)
sum((full_data$inSelfReport==TRUE & full_data$menstruation==1), na.rm=TRUE) # 108

## Q: What proportion of self reported have corresponding samples? - includes no menses ppl too
selfReport_blood <- full_data %>% 
  filter(inUminn==TRUE & inSelfReport==TRUE)
dim(selfReport_blood)
sum((full_data$inSelfReport==TRUE), na.rm=TRUE) # 128

### Agreement we have btw UMinn data & Self Reported (put in a slide)
## Q: Do samples for the same day agree on blood/no-blood?
temp_full_data <- full_data %>% 
  mutate(menstruation=ifelse(menstruation==1, 1, 0))
table(temp_full_data$uMinn_menstruation, temp_full_data$menstruation)

## Q: What proportion of self reported no menstruation had samples - we don't have samples for days reported no menses...
dim(full_data %>% filter(inSelfReport==TRUE & menstruation==0 & inUminn==TRUE))

# Q: What number of days self reported menses had samples
dim(full_data %>% filter(inSelfReport==TRUE & uMinn_menstruation==1 & inUminn==TRUE))

## Q: For days we have self report & samples, do they correspond
selfReport_samples <- full_data %>%
  filter(inSelfReport==TRUE & inUminn==TRUE) 
dim(selfReport_samples) # 72 with both self report & samples
table(selfReport_samples$menstruation_status)
table(selfReport_samples$uMinn_menstruation, selfReport_samples$menstruation)

dim(full_data %>% filter(inSelfReport==TRUE)) # 128 in self report
dim(full_data %>% filter(inUminn==TRUE)) # 1422 in UMinn

## Q: Days in UMinn with Blood
dim(full_data %>% filter(inUminn==TRUE & uMinn_menstruation==TRUE))
# sum(full_data$sampleType=="vaginal", na.rm=TRUE)

# Make menstruation variable
full_data <- full_data %>% 
  mutate(person_menstruating=ifelse((menstruation_status==1 | menstruation_status==2 | menstruation_status==3 
                                     | menstruation_status==7), 1, 
                                    ifelse(menstruation_status==8, -1, 0)))

#### Contingency Tables -- smthg going wrong here

# filter for menstruating and add survey data
menstruate_df <- full_data %>% 
  filter(person_menstruating==1) %>% 
  left_join(survey_data, by="biome_id")
survey_data <- survey_data %>%
  mutate(menses_data_entry = ifelse(biome_id %in% menstruate_df$biome_id, 1, 0)) #%>%
  # mutate(survey_menstruate=ifelse(survey_menstruate=="Yes", 1, 0))

# find mismatched entries - if they say menstruate 0 menses_data_entry 1
mismatched_entries <- survey_data %>%
  filter(survey_menstruate == 0 & menses_data_entry == 1)
dim(mismatched_entries)
length(unique(mismatched_entries$biome_id))

mismatch_menstruate_df <- menstruate_df %>%
  mutate(survey_menstruate=ifelse(survey_menstruate=="Yes", 1, 0)) %>% 
  filter(((menstruation==1 | uMinn_menstruation==1) & survey_menstruate == 0))

# All ppl who say they don't menstruate
person_says_no_menses <- full_data %>% 
  left_join(survey_data, by="biome_id") %>% 
  filter(survey_menstruate==0)

# Check people who said they do menstruate with no entries (missingness)
missing_entries <- survey_data %>%
  filter(survey_menstruate == 1 & menses_data_entry == 0)
dim(missing_entries) # 3
length(unique(missing_entries$biome_id))

# Number of people who don't menstruate (say they don't & no menstruation reports for them)
no_menses <- survey_data %>%
  filter(survey_menstruate == 0 & menses_data_entry == 0)
dim(no_menses)
unique(no_menses$biome_id)
length(unique(no_menses$biome_id))

# Do menstruate with menses entries
yes_menses <- survey_data %>%
  filter(survey_menstruate == 1 & menses_data_entry == 1)
dim(yes_menses)
(unique(yes_menses$biome_id))
length((unique(yes_menses$biome_id)))

## Contingency table
menses_cont_table <- table(survey_data$menses_data_entry, survey_data$survey_menstruate) # menstruate is from the volunteer survey
dimnames(menses_cont_table) <- list("Menses in Study" = c("No", "Yes"), 
                                    "Volunteer Survey Data" = c("No", "Yes"))
menses_cont_table

### Plot of menses

# filter to study days
filtered_data <- full_data %>%
  filter(as.Date(logDate) < as.Date("2022-12-15")) %>% 
  filter(!is.na(biome_id))

# no menstruation
menstruate_df <- filtered_data %>% 
  group_by(biome_id) %>% 
  summarise(count = sum(person_menstruating == 1, na.rm = TRUE))
menstruate_df[menstruate_df$count==0,]$biome_id

##################### Heatmap Construction - Blue & White

# Plot 1: all participants
all_days <- seq.Date(as.Date(min(filtered_data$logDate, na.rm=TRUE)), 
                     as.Date(max(filtered_data$logDate, na.rm=TRUE)), by = "day")
all_days <- as.character(all_days)
full_data_subset <- filtered_data %>% 
  select(biome_id, logDate, person_menstruating)
complete_grid <- expand.grid(biome_id = unique(full_data_subset$biome_id),
                             logDate = all_days)
heatmap_data <- complete_grid %>%
  left_join(full_data_subset, by = c("biome_id", "logDate"))

heatmap_data$Value <- ifelse((is.na(heatmap_data$person_menstruating) | heatmap_data$person_menstruating==0), "Not menstruating", "Menstruating")
heatmap_data$Value <- factor(heatmap_data$Value, levels = c("Not menstruating", "Menstruating"))

heatmap_data <- heatmap_data %>% 
  group_by(biome_id) %>%
  mutate(menstruating_days_count = sum(Value == "Menstruating", na.rm = TRUE)) %>%
  arrange(desc(menstruating_days_count), biome_id) %>% 
  ungroup()

ggplot(heatmap_data, aes(x = logDate, y = reorder(factor(biome_id), menstruating_days_count), fill = Value)) +
  geom_tile(color = "gray25") +
  scale_fill_manual(values = c("white", "darkblue"), na.value = "white", 
                    labels = c("Not menstruating", "Menstruating")) +
  labs(title = "Menstruation Heatmap", 
       x = "Date", 
       y = "Biome ID", 
       fill = "Menstruation Scale") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

# Plot 2: Mismatch participants
mismatch_ids <- unique(mismatch_menstruate_df$biome_id)
mismatch_menstruate_df <- full_data %>%
  filter(biome_id %in% mismatch_ids) %>% 
  filter(as.Date(logDate) < as.Date("2022-12-15")) %>% 
  filter(!is.na(biome_id))
all_days <- seq.Date(as.Date(min(mismatch_menstruate_df$logDate, na.rm=TRUE)), 
                     as.Date(max(mismatch_menstruate_df$logDate, na.rm=TRUE)), by = "day")
all_days <- as.character(all_days)
full_data_subset <- mismatch_menstruate_df %>% 
  select(biome_id, logDate, person_menstruating, menstruation_status, menstruation, uMinn_menstruation)
complete_grid <- expand.grid(biome_id = unique(full_data_subset$biome_id),
                             logDate = all_days)
heatmap_data <- complete_grid %>%
  left_join(full_data_subset, by = c("biome_id", "logDate"))

heatmap_data$Value <- ifelse((is.na(heatmap_data$person_menstruating) | heatmap_data$person_menstruating==0), 
                             "Not menstruating", "Menstruating")
heatmap_data$Value <- factor(heatmap_data$Value, levels = c("Not menstruating", "Menstruating"))

heatmap_data <- heatmap_data %>% 
  group_by(biome_id) %>%
  mutate(menstruating_days_count = sum(Value == "Menstruating", na.rm = TRUE)) %>%
  arrange(desc(menstruating_days_count), biome_id) %>% 
  ungroup()

ggplot(heatmap_data, aes(x = logDate, y = reorder(factor(biome_id), menstruating_days_count), fill = Value)) +
  geom_tile(color = "gray25") +
  scale_fill_manual(values = c("white", "darkblue"), na.value = "white", 
                    labels = c("Not menstruating", "Menstruating")) +
  labs(title = "Menstruation Heatmap", 
       x = "Days from First Menstruating Day", 
       y = "Biome ID", 
       fill = "Menstruation Scale") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

# Plot 2.5: Ppl who say no menses
person_says_no_menses_ids <- unique(person_says_no_menses$biome_id)

no_menses_data <- person_says_no_menses %>% 
  group_by(biome_id) %>% 
  summarise(
    has_menses = if_else(any(menstruation == 1 | uMinn_menstruation == 1), 1, 0)
  ) %>% 
  ungroup()

## Survey data
survey_data_subset <- survey_data_full %>% 
  filter(biome_id %in% person_says_no_menses_ids) %>% 
  select(biome_id, activity_level, sport, taken_antibiotics, birthControl, survey_menstruate, regular_periods, irreg_period_notes,
         period_len,
         menstrual_prod, vaginal_notes, menstrual_cup, tampon, pad, no_menstrual_product) %>% 
  full_join(no_menses_data) %>% 
  filter(has_menses==1) %>% 
  select(has_menses, everything())

person_says_no_menses <- full_data %>%
  filter(biome_id %in% person_says_no_menses_ids) %>% 
  filter(as.Date(logDate) < as.Date("2022-12-15")) %>% 
  filter(!is.na(biome_id))
all_days <- seq.Date(as.Date(min(person_says_no_menses$logDate, na.rm=TRUE)), 
                     as.Date(max(person_says_no_menses$logDate, na.rm=TRUE)), by = "day")
all_days <- as.character(all_days)
full_data_subset <- person_says_no_menses %>% 
  select(biome_id, logDate, person_menstruating, menstruation_status, menstruation, uMinn_menstruation)
complete_grid <- expand.grid(biome_id = unique(full_data_subset$biome_id),
                             logDate = all_days)
heatmap_data <- complete_grid %>%
  left_join(full_data_subset, by = c("biome_id", "logDate"))

heatmap_data$Value <- ifelse((is.na(heatmap_data$person_menstruating) | heatmap_data$person_menstruating==0), 
                             "Not menstruating", "Menstruating")
heatmap_data$Value <- factor(heatmap_data$Value, levels = c("Not menstruating", "Menstruating"))

heatmap_data <- heatmap_data %>% 
  group_by(biome_id) %>%
  mutate(menstruating_days_count = sum(Value == "Menstruating", na.rm = TRUE)) %>%
  arrange(desc(menstruating_days_count), biome_id) %>% 
  ungroup()

ggplot(heatmap_data, aes(x = logDate, y = reorder(factor(biome_id), menstruating_days_count), fill = Value)) +
  geom_tile(color = "gray25") +
  scale_fill_manual(values = c("white", "darkblue"), na.value = "white", 
                    labels = c("Not menstruating", "Menstruating")) +
  labs(title = "Menstruation Heatmap", 
       x = "Date", 
       y = "Biome ID", 
       fill = "Menstruation Scale") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

### Reconstruct cycles

# Plot 3: all days, pushed to relative day 1 (first menses data)
all_days <- seq.Date(as.Date(min(filtered_data$logDate, na.rm=TRUE)), 
                     as.Date(max(filtered_data$logDate, na.rm=TRUE)), by = "day")
all_days <- as.character(all_days)
full_data_subset <- filtered_data %>% 
  select(biome_id, logDate, person_menstruating)
complete_grid <- expand.grid(biome_id = unique(full_data_subset$biome_id),
                             logDate = all_days)
heatmap_data <- complete_grid %>%
  left_join(full_data_subset, by = c("biome_id", "logDate"))

heatmap_data$Value <- ifelse((is.na(heatmap_data$person_menstruating) | heatmap_data$person_menstruating==0), "Not menstruating", "Menstruating")
heatmap_data$Value <- factor(heatmap_data$Value, levels = c("Not menstruating", "Menstruating"))

# Calculate the relative day for each biome_id starting from their first "Menstruating" day
relative_day_heatmap_data <- heatmap_data %>%
  mutate(logDate = as.Date(logDate)) %>% 
  group_by(biome_id) %>%
  arrange(logDate) %>%
  mutate(RelativeDay = row_number() - min(row_number()[Value == "Menstruating"], na.rm = TRUE)) %>%
  filter(RelativeDay >= 0) %>% 
  ungroup() %>% 
  group_by(biome_id) %>%
  mutate(menstruating_days_count = sum(Value == "Menstruating", na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(menstruating_days_count), biome_id)

ggplot(relative_day_heatmap_data, aes(x = RelativeDay, y = reorder(factor(biome_id), menstruating_days_count), fill = Value)) +
  geom_tile(color = "gray25") +
  scale_fill_manual(values = c("white", "darkblue"), na.value = "white", 
                    labels = c("Not menstruating", "Menstruating")) +
  labs(title = "Menstruation Heatmap (Aligned to First Menstruating Day)", 
       x = "Days from First Menstruating Day", 
       y = "Biome ID", 
       fill = "Menstruation Scale") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

# Plot 4: Color by what filled in the day

full_data_subset <- filtered_data %>% 
  select(biome_id, logDate, menstruation_status, uMinn_menstruation, menstruation, person_menstruating)
complete_grid <- expand.grid(biome_id = unique(full_data_subset$biome_id),
                             logDate = all_days)

heatmap_data <- complete_grid %>%
  left_join(full_data_subset, by = c("biome_id", "logDate")) %>%
  mutate(menstruation_identifier=case_when(
    # both UMinn Sample has blood & Self-report menses
    (uMinn_menstruation==1 & menstruation==1) ~ 1,
    # UMinn sample with blood & no self report
    (uMinn_menstruation==1 & (menstruation==0 | is.na(menstruation))) ~ 2,
    ((uMinn_menstruation==0 | is.na(uMinn_menstruation)) & menstruation==1) ~ 3,
    TRUE ~ NA_real_
  ))

heatmap_data$Value <- ifelse((is.na(heatmap_data$menstruation_identifier) | heatmap_data$menstruation_identifier==0), "Not menstruating", "Menstruating")
heatmap_data$Value <- factor(heatmap_data$Value, levels = c("Not menstruating", "Menstruating", "No report"))

heatmap_data <- heatmap_data %>% 
  group_by(biome_id) %>%
  mutate(menstruating_days_count = sum(Value == "Menstruating", na.rm = TRUE)) %>%
  arrange(desc(menstruating_days_count), biome_id) %>% 
  ungroup()

ggplot(heatmap_data, aes(x = logDate, y = reorder(factor(biome_id), menstruating_days_count), fill = factor(menstruation_identifier))) +
  geom_tile(color = "gray25") +
  scale_fill_manual(values = c("1" = "darkblue", "2" = "lightblue", "3" = "purple", 
                               "NA" = "white"),  # Specify "NA" as a string in the values argument
                    labels = c("1" = "Both sources agree", 
                               "2" = "UMinn Sample", 
                               "3" = "User Input", 
                               "NA" = "Missing"), 
                    na.value = "white") +
  labs(title = "Menstruation Heatmap by Identifier", 
       x = "Date", 
       y = "Biome ID", 
       fill = "Menstruation Identifier") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

# Plot 5: Calculate the relative day for each biome_id starting from their first "Menstruating" day
relative_day_heatmap_data <- heatmap_data %>%
  mutate(logDate = as.Date(logDate)) %>% 
  group_by(biome_id) %>%
  arrange(logDate) %>%
  mutate(RelativeDay = row_number() - min(row_number()[Value == "Menstruating"], na.rm = TRUE)) %>%
  filter(RelativeDay >= 0) %>% 
  ungroup() %>% 
  group_by(biome_id) %>%
  mutate(menstruating_days_count = sum(Value == "Menstruating", na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(menstruating_days_count), biome_id)

ggplot(relative_day_heatmap_data, aes(x = RelativeDay, y = reorder(factor(biome_id), menstruating_days_count), fill = factor(menstruation_identifier))) +
  geom_tile(color = "gray25") +
  scale_fill_manual(values = c("1" = "darkblue", "2" = "lightblue", "3" = "purple", 
                               "NA" = "white"),  # Specify "NA" as a string in the values argument
                    labels = c("1" = "Both sources agree", 
                               "2" = "UMinn Sample", 
                               "3" = "User Input", 
                               "NA" = "Missing"), 
                    na.value = "white") +
  labs(title = "Menstruation Heatmap by Identifier", 
       x = "Date", 
       y = "Biome ID", 
       fill = "Menstruation Identifier") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))


################################################################# Imputation try 1

all_days <- seq.Date(as.Date(min(filtered_data$logDate, na.rm=TRUE)), 
                     as.Date(max(filtered_data$logDate, na.rm=TRUE)), by = "day")
all_days <- as.character(all_days)
# Expand to include all days
expanded_data <- expand.grid(biome_id = unique(filtered_data$biome_id), logDate = all_days)
expanded_data <- left_join(expanded_data, filtered_data, by = c("biome_id", "logDate"))

imputed_filtered_data <- expanded_data %>%  
  group_by(biome_id) %>%
  arrange(logDate) %>%
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
  select(biome_id, logDate, person_menstruating, imputation_menses, prev_menstruating_date, 
         next_menstruating_date, gap_between_menstruations, day_gap, 
         # timestamp, 
         inUminn, inSelfReport, uMinn_menstruation, menstruation) %>%
  ungroup()

sum(imputed_filtered_data$imputation_menses==1, na.rm=TRUE)
sum(filtered_data$person_menstruating==1)

# Plot of imputed data - blue & white
all_days <- seq.Date(as.Date(min(imputed_filtered_data$logDate, na.rm=TRUE)), 
                     as.Date(max(imputed_filtered_data$logDate, na.rm=TRUE)), by = "day")
all_days <- as.character(all_days)
full_data_subset <- imputed_filtered_data %>% 
  select(biome_id, logDate, imputation_menses)
complete_grid <- expand.grid(biome_id = unique(full_data_subset$biome_id),
                             logDate = all_days)
heatmap_data <- complete_grid %>%
  left_join(full_data_subset, by = c("biome_id", "logDate"))

heatmap_data$Value <- ifelse((is.na(heatmap_data$imputation_menses) | heatmap_data$imputation_menses==0), "Not menstruating", "Menstruating")
heatmap_data$Value <- factor(heatmap_data$Value, levels = c("Not menstruating", "Menstruating"))

heatmap_data <- heatmap_data %>% 
  group_by(biome_id) %>%
  mutate(menstruating_days_count = sum(Value == "Menstruating", na.rm = TRUE)) %>%
  arrange(desc(menstruating_days_count), biome_id) %>% 
  ungroup()

ggplot(heatmap_data, aes(x = logDate, y = reorder(factor(biome_id), menstruating_days_count), fill = Value)) +
  geom_tile(color = "gray25") +
  scale_fill_manual(values = c("white", "darkblue"), na.value = "white", 
                    labels = c("Not menstruating", "Menstruating")) +
  labs(title = "Menstruation Heatmap", 
       x = "Date", 
       y = "Biome ID", 
       fill = "Menstruation Scale") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

# Plot of imputed data colored with imputation
full_data_subset <- imputed_filtered_data %>% 
  select(biome_id, logDate, imputation_menses, uMinn_menstruation, menstruation, person_menstruating)

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

##################### Heatmap Construction - fill in type

full_data_subset <- filtered_data %>% 
  select(biome_id, logDate, menstruation_status, uMinn_menstruation, menstruation, person_menstruating)
complete_grid <- expand.grid(biome_id = unique(full_data_subset$biome_id),
                             logDate = all_days)

# Plot: Map to with all ppl
heatmap_data <- complete_grid %>%
  left_join(full_data_subset, by = c("biome_id", "logDate")) %>%
  mutate(
    menstruation_identifier = case_when(
      # Assigning menstruation_identifier based on menstruation_status
      menstruation_status == 1 ~ "Menstruating (Sample with Blood)",
      menstruation_status %in% c(2, 7) ~ "Menstruating (Self-Reported)",
      menstruation_status == 3 ~ "Menstruating (Both Sources)",
      menstruation_status %in% c(5, 6) ~ "Not Menstruating (Sample with No Blood)",
      menstruation_status == 4 ~ "Not Menstruating (Self-Reported)",
      menstruation_status == 8 ~ "No Report",
      TRUE ~ "No Data (can impute)"
    ),
    menstruation_identifier = factor(menstruation_identifier, levels = c(
      "Not Menstruating (Sample with No Blood)",
      "Not Menstruating (Self-Reported)",
      "Menstruating (Sample with Blood)",
      "Menstruating (Self-Reported)",
      "Menstruating (Both Sources)",
      "No Report",
      "No Data (can impute)"
    ))
  ) %>%
  # Sorting data for plotting
  group_by(biome_id) %>%
  mutate(
    menstruating_days_count = sum(menstruation_identifier %in% c(
      "Menstruating (Sample with Blood)",
      "Menstruating (Self-Reported)",
      "Menstruating (Both Sources)"
    ), na.rm = TRUE)
  ) %>%
  arrange(desc(menstruating_days_count), biome_id) %>%
  ungroup()

ggplot(heatmap_data, aes(
  x = logDate,
  y = reorder(factor(biome_id), menstruating_days_count),
  fill = menstruation_identifier
)) +
  geom_tile(color = "gray25") +
  scale_fill_manual(
    values = c(
      "Not Menstruating (Sample with No Blood)" = "blue",
      "Not Menstruating (Self-Reported)" = "lightblue",
      "Menstruating (Sample with Blood)" = "red",
      "Menstruating (Self-Reported)" = "pink",
      "Menstruating (Both Sources)" = "darkred",
      "No Report" = "gray",
      "No Data (can impute)" = "white"
    ),
    na.value = "white" # Handle missing data
  ) +
  labs(
    title = "Menstruation Heatmap by Identifier",
    x = "Date",
    y = "Biome ID",
    fill = "Menstruation Status"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)
  )

### Plot: say no menses

# All ppl who say they don't menstruate
full_data_subset <- person_says_no_menses %>% 
  select(biome_id, logDate, menstruation_status, uMinn_menstruation, menstruation, person_menstruating)
complete_grid <- expand.grid(biome_id = unique(full_data_subset$biome_id),
                             logDate = all_days)

# Plot: Map to with all ppl
heatmap_data <- complete_grid %>%
  left_join(full_data_subset, by = c("biome_id", "logDate")) %>%
  mutate(
    menstruation_identifier = case_when(
      # Assigning menstruation_identifier based on menstruation_status
      menstruation_status == 1 ~ "Menstruating (Sample with Blood)",
      menstruation_status %in% c(2, 7) ~ "Menstruating (Self-Reported)",
      menstruation_status == 3 ~ "Menstruating (Both Sources)",
      menstruation_status %in% c(5, 6) ~ "Not Menstruating (Sample with No Blood)",
      menstruation_status == 4 ~ "Not Menstruating (Self-Reported)",
      menstruation_status == 8 ~ "No Report",
      TRUE ~ "No Data (can impute)"
    ),
    menstruation_identifier = factor(menstruation_identifier, levels = c(
      "Not Menstruating (Sample with No Blood)",
      "Not Menstruating (Self-Reported)",
      "Menstruating (Sample with Blood)",
      "Menstruating (Self-Reported)",
      "Menstruating (Both Sources)",
      "No Report",
      "No Data (can impute)"
    ))
  ) %>%
  # Sorting data for plotting
  group_by(biome_id) %>%
  mutate(
    menstruating_days_count = sum(menstruation_identifier %in% c(
      "Menstruating (Sample with Blood)",
      "Menstruating (Self-Reported)",
      "Menstruating (Both Sources)"
    ), na.rm = TRUE)
  ) %>%
  arrange(desc(menstruating_days_count), biome_id) %>%
  ungroup()

ggplot(heatmap_data, aes(
  x = logDate,
  y = reorder(factor(biome_id), menstruating_days_count),
  fill = menstruation_identifier
)) +
  geom_tile(color = "gray25") +
  scale_fill_manual(
    values = c(
      "Not Menstruating (Sample with No Blood)" = "blue",
      "Not Menstruating (Self-Reported)" = "lightblue",
      "Menstruating (Sample with Blood)" = "red",
      "Menstruating (Self-Reported)" = "pink",
      "Menstruating (Both Sources)" = "darkred",
      "No Report" = "gray",
      "No Data (can impute)" = "white"
    ),
    na.value = "white" # Handle missing data
  ) +
  labs(
    title = "Menstruation Heatmap by Identifier \n Self-report don't menstruate",
    x = "Date",
    y = "Biome ID",
    fill = "Menstruation Status"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)
  )

### Plot: Map to with menses relative day 1
full_data_subset <- filtered_data %>% 
  select(biome_id, logDate, menstruation_status, uMinn_menstruation, menstruation, person_menstruating)
complete_grid <- expand.grid(biome_id = unique(full_data_subset$biome_id),
                             logDate = all_days)
heatmap_data <- complete_grid %>%
  left_join(full_data_subset, by = c("biome_id", "logDate"))
heatmap_data$Value <- ifelse((heatmap_data$menstruation_status==8), "No report", 
                             ifelse(heatmap_data$menstruation_status==1 | heatmap_data$menstruation_status==2 | heatmap_data$menstruation_status==3 | heatmap_data$menstruation_status==7, "Menstruating", 
                                    ifelse(heatmap_data$menstruation_status==4 | heatmap_data$menstruation_status==5 | heatmap_data$menstruation_status==6, "Not menstruating", NA)))
heatmap_data$Value <- factor(heatmap_data$Value, levels = c("Not menstruating", "Menstruating", "No report"))

relative_day_heatmap_data <- heatmap_data %>%
  mutate(logDate = as.Date(logDate)) %>% 
  group_by(biome_id) %>%
  arrange(logDate) %>%
  mutate(RelativeDay = row_number() - min(row_number()[Value == "Menstruating"], na.rm = TRUE)) %>%
  filter(RelativeDay >= 0) %>% 
  ungroup() %>% 
  group_by(biome_id) %>%
  mutate(menstruating_days_count = sum(Value == "Menstruating", na.rm = TRUE),
         sample_count=sum(menstruation_status %in% c(1,3,5,6,7), na.rm=TRUE)) %>%
  ungroup() %>%
  arrange(desc(menstruating_days_count), biome_id)

ggplot(relative_day_heatmap_data, 
       aes(x = RelativeDay,
           y=reorder(factor(biome_id), menstruating_days_count + sample_count / max(sample_count)),
           # y = reorder(factor(biome_id), menstruating_days_count), 
           fill = factor(menstruation_status)))+
  geom_tile(color = "gray25") +
  scale_fill_manual(values = c("1" = "red", "2" = "maroon", "3" = "darkred", 
                               "4"="lightblue", "5"="darkblue","6"="blue",
                               "7"="pink","8"="gray",
                               "NA" = "white"),  # NA - technically can impute these?
                    labels = c("1" = "Blood UMinn Sample", 
                               "2" = "User Input Menses", 
                               "3" = "Both Sources Menses", 
                               "4"="User Input no Menses, no Sample",
                               "5"="Both Sources no Menses",
                               "6"="No Blood Sample and no User Input",
                               "7"="User Input Menses, no Blood Sample",
                               "8"="No data",
                               "NA" = "Missing"), 
                    na.value = "white") +
  labs(title = "Menstruation Heatmap by Identifier", 
       x = "Relative Day 1 to Menstruation Day 1", 
       y = "Biome ID", 
       fill = "Menstruation Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

### Plot: Map to with menses relative day 1 - order by samples
ggplot(relative_day_heatmap_data, 
       aes(x = RelativeDay,
           y=reorder(factor(biome_id), sample_count),
           # y = reorder(factor(biome_id), menstruating_days_count), 
           fill = factor(menstruation_status)))+
  geom_tile(color = "gray25") +
  scale_fill_manual(values = c("1" = "red", "2" = "maroon", "3" = "darkred", 
                               "4"="lightblue", "5"="darkblue","6"="blue",
                               "7"="pink","8"="gray",
                               "NA" = "white"),  # NA - technically can impute these?
                    labels = c("1" = "Blood UMinn Sample", 
                               "2" = "User Input Menses", 
                               "3" = "Both Sources Menses", 
                               "4"="User Input no Menses, no Sample",
                               "5"="Both Sources no Menses",
                               "6"="No Blood Sample and no User Input",
                               "7"="User Input Menses, no Blood Sample",
                               "8"="No data",
                               "NA" = "Missing"), 
                    na.value = "white") +
  labs(title = "Menstruation Heatmap by Identifier", 
       x = "Relative Day 1 to Menstruation Day 1", 
       y = "Biome ID", 
       fill = "Menstruation Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

### Plot: Map to impute from Uminn sample & menses relative day 1

relative_day_heatmap_data <- heatmap_data %>%
  mutate(logDate = as.Date(logDate)) %>% 
  group_by(biome_id) %>%
  arrange(logDate) %>%
  # filter out ppl with no menses at all
  filter(any(menstruation_status %in% c(1, 2, 3, 7))) %>% 
  # Identify the earliest logDate where we have sample
  mutate(first_sample_day = min(logDate[menstruation_status %in% c(1, 3, 5, 6, 7)], na.rm = TRUE),
         first_menstruating_day = min(logDate[Value == "Menstruating"], na.rm = TRUE))%>%
  # filter if no data before first menses report
  filter(logDate >= first_sample_day | 
           (logDate < first_sample_day & !is.na(menstruation_status))) %>%
  # Calculate RelativeDay with potential negative values
  mutate(RelativeDay = as.numeric(logDate - first_menstruating_day)) %>%
  ungroup() %>% 
  group_by(biome_id) %>%
  # Count menstruating days for sorting
  mutate(menstruating_days_count = sum(Value == "Menstruating", na.rm = TRUE),
         sample_count=sum(menstruation_status %in% c(1,3,5,6,7), na.rm=TRUE)) %>%
  ungroup() %>%
  arrange(desc(menstruating_days_count), desc(sample_count), biome_id)

ggplot(relative_day_heatmap_data, 
       aes(x = RelativeDay,
           y=reorder(factor(biome_id), menstruating_days_count + sample_count / max(sample_count)),
           # y = reorder(factor(biome_id), menstruating_days_count), 
           fill = factor(menstruation_status)))+
  geom_tile(color = "gray25") +
  scale_fill_manual(values = c("1" = "red", "2" = "maroon", "3" = "darkred", 
                               "4"="lightblue", "5"="darkblue","6"="blue",
                               "7"="pink","8"="gray",
                               "NA" = "white"),  # NA - technically can impute these?
                    labels = c("1" = "Blood UMinn Sample", 
                               "2" = "User Input Menses", 
                               "3" = "Both Sources Menses", 
                               "4"="User Input no Menses, no Sample",
                               "5"="Both Sources no Menses",
                               "6"="No Blood Sample and no User Input",
                               "7"="User Input Menses, no Blood Sample",
                               "8"="No data",
                               "NA" = "Missing"), 
                    na.value = "white") +
  labs(title = "Menstruation Heatmap by Identifier", 
       x = "Relative Day from Menstruation Day 1", 
       y = "Biome ID", 
       fill = "Menstruation Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

## Plot: Map to impute from Uminn sample & menses relative day 1 - order by num samples

ggplot(relative_day_heatmap_data, 
       aes(x = RelativeDay,
           y=reorder(factor(biome_id), sample_count),
           fill = factor(menstruation_status)))+
  geom_tile(color = "gray25") + 
  scale_fill_manual(values = c("1" = "red", "2" = "maroon", "3" = "darkred", 
                               "4"="lightblue", "5"="darkblue","6"="blue",
                               "7"="pink","8"="gray",
                               "NA" = "white"),  # NA - technically can impute these?
                    labels = c("1" = "Blood UMinn Sample", 
                               "2" = "User Input Menses", 
                               "3" = "Both Sources Menses", 
                               "4"="User Input no Menses, no Sample",
                               "5"="Both Sources no Menses",
                               "6"="No Blood Sample and no User Input",
                               "7"="User Input Menses, no Blood Sample",
                               "8"="No data",
                               "NA" = "Missing"), 
                    na.value = "white") +
  labs(title = "Menstruation Heatmap by Identifier", 
       x = "Relative Day from Menstruation Day 1", 
       y = "Biome ID", 
       fill = "Menstruation Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

## Cut participants 
# from samples_analysis.R
# (ids15)
# unique(no_menses_data$biome_id)
# ids15
# 
# relative_day_heatmap_data <- heatmap_data %>%
#   mutate(logDate = as.Date(logDate)) %>% 
#   group_by(biome_id) %>%
#   arrange(logDate) %>%
#   # filter out ppl with no menses at all
#   # filter(any(menstruation_status %in% c(1, 2, 3, 7))) %>% 
#   # filter ppl
#   filter(biome_id %in% ids15) %>% 
#   # Identify the earliest logDate where we have sample
#   mutate(first_sample_day = min(logDate[menstruation_status %in% c(1, 3, 5, 6, 7)], na.rm = TRUE),
#          first_menstruating_day = min(logDate[Value == "Menstruating"], na.rm = TRUE))%>%
#   # set first menstruating day to first data day if no menses
#   mutate(first_menstruating_day = ifelse((first_menstruating_day==Inf), first_sample_day, first_menstruating_day)) %>%
#   # filter if no data before first menses report
#   filter(logDate >= first_sample_day | 
#            (logDate < first_sample_day & !is.na(menstruation_status))) %>%
#   # Calculate RelativeDay with potential negative values
#   mutate(RelativeDay = as.numeric(logDate - first_menstruating_day)) %>%
#   ungroup() %>% 
#   group_by(biome_id) %>%
#   # Count menstruating days for sorting
#   mutate(menstruating_days_count = sum(Value == "Menstruating", na.rm = TRUE),
#          sample_count=sum(menstruation_status %in% c(1,3,5,6,7), na.rm=TRUE)) %>%
#   ungroup() %>%
#   arrange(desc(menstruating_days_count), desc(sample_count), biome_id)
# 
# ggplot(relative_day_heatmap_data, 
#        aes(x = RelativeDay,
#            y=reorder(factor(biome_id), sample_count),
#            fill = factor(menstruation_status)))+
#   geom_tile(color = "gray25") + 
#   scale_fill_manual(values = c("1" = "red", "2" = "maroon", "3" = "darkred", 
#                                "4"="lightblue", "5"="darkblue","6"="blue",
#                                "7"="pink","8"="gray",
#                                "NA" = "white"),  # NA - technically can impute these?
#                     labels = c("1" = "Blood UMinn Sample", 
#                                "2" = "User Input Menses", 
#                                "3" = "Both Sources Menses", 
#                                "4"="User Input no Menses, no Sample",
#                                "5"="Both Sources no Menses",
#                                "6"="No Blood Sample and no User Input",
#                                "7"="User Input Menses, no Blood Sample",
#                                "8"="No data",
#                                "NA" = "Missing"), 
#                     na.value = "white") +
#   labs(title = "Menstruation Heatmap by Identifier", 
#        x = "Relative Day from Menstruation Day 1", 
#        y = "Biome ID", 
#        fill = "Menstruation Status") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
# 

### Plot: days to impute

# add variable for days to impute - mutate menstruation status
relative_day_heatmap_data_impute <- relative_day_heatmap_data %>% 
  mutate(impute_status=case_when(
    # menses
    (menstruation_status %in% c(1,2,3,7)) ~ "Menses",
    # no menses
    (menstruation_status %in% c(4,5,6)) ~ "No menses",
    # impute
    (menstruation_status==8 | is.na(menstruation_status)) ~ "Impute",
    TRUE ~ NA_character_
  ))



# make the days that should have imputation done marked in gray
ggplot(relative_day_heatmap_data_impute, 
       aes(x = RelativeDay,
           y=reorder(factor(biome_id), menstruating_days_count + sample_count / max(sample_count)),
           # y = reorder(factor(biome_id), menstruating_days_count), 
           fill = factor(impute_status)))+
  geom_tile(color = "gray25") +
  scale_fill_manual(values = c("Menses" = "red", 
                               "No menses"= "black",
                               "Impute"="white",
                               "NA" = "gray"),
                    labels = c("Menses" = "Confirmed menses", 
                               "No menses" = "Confirmed no menses", 
                               "Impute" = "Can impute",
                               "NA" = "Missing"), 
                    na.value = "white") +
  labs(title = "Menstruation Heatmap by Identifier (prep for imputation)", 
       x = "Date", 
       y = "Biome ID", 
       fill = "Menstruation Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

################################################################# Imputation try 2

# Work from here, imputation doesn't work
# idea: fill in cycles first, count how many cycles per person i have
# impute cycles

# get how many "cycles" we have for a person

table(relative_day_heatmap_data_impute$impute_status)

imputed_filtered_data <- relative_day_heatmap_data_impute %>%  
  group_by(biome_id) %>%
  arrange(logDate) %>%
  mutate(
    person_menstruating = ifelse(is.na(person_menstruating), 0, person_menstruating),
    imputation_menses = person_menstruating,
    
    # Identify the last menstruation date
    prev_menstruating_date = ifelse(person_menstruating == 1, lag(logDate), NA),
    prev_menstruating_date = zoo::na.locf(prev_menstruating_date, na.rm = FALSE),
    
    # Identify the next menstruation date
    next_menstruating_date = sapply(1:n(), function(i) {
      next_index <- which(person_menstruating[i:n()] == 1)[1] + (i - 1)
      if (length(next_index) > 0) logDate[next_index] else NA
    }),
    next_menstruating_date = zoo::na.locf(next_menstruating_date, na.rm = FALSE),
    next_menstruating_date = ifelse(next_menstruating_date > logDate, next_menstruating_date, NA),
    
    # Calculate gaps in days
    gap_between_menstruations = as.numeric(difftime(next_menstruating_date, prev_menstruating_date, units = "days")),
    day_gap = round(as.numeric(difftime(logDate, prev_menstruating_date, units = "days")), 4),
    
    # Create an indicator for whether imputation should proceed
    valid_imputation = cumsum(impute_status == "No menses"),
    valid_imputation = ifelse(impute_status == "Impute" & lag(valid_imputation, default = 0) == 0, 1, 0),
    valid_imputation = zoo::na.locf(valid_imputation, na.rm = FALSE),
    
    # Impute menstruation if all conditions are met
    imputation_menses = ifelse(
      valid_imputation == 1 & 
        gap_between_menstruations > 0 & 
        gap_between_menstruations < 10 & 
        person_menstruating == 0 & 
        impute_status == "Impute" & 
        !is.na(prev_menstruating_date) & 
        !is.na(next_menstruating_date),
      1, imputation_menses
    )
  ) %>% 
  select(
    biome_id, logDate, person_menstruating, imputation_menses, 
    menstruation_status,
    prev_menstruating_date, next_menstruating_date, gap_between_menstruations, 
    day_gap, valid_imputation, impute_status
  ) %>%
  ungroup()

sum(imputed_filtered_data$imputation_menses==1, na.rm=TRUE)
sum(filtered_data$person_menstruating==1)

### plotting
all_days <- seq.Date(as.Date(min(imputed_filtered_data$logDate, na.rm=TRUE)), 
                     as.Date(max(imputed_filtered_data$logDate, na.rm=TRUE)), by = "day")
all_days <- as.character(all_days)
# full_data_subset <- imputed_filtered_data %>% 
#   select(biome_id, logDate, imputation_menses)
complete_grid <- expand.grid(biome_id = unique(imputed_filtered_data$biome_id),
                             logDate = all_days) %>% 
  mutate(logDate = as.Date(logDate))
# heatmap_data <- complete_grid %>%
#   left_join(imputed_filtered_data, by = c("biome_id", "logDate"))

heatmap_data <- complete_grid %>%
  left_join(imputed_filtered_data, by = c("biome_id", "logDate")) %>% 
  mutate(imputation_identifier=case_when(
    (imputation_menses==1 & person_menstruating==0) ~ 1, # imputed data point
    (menstruation_status %in% c(1,2,3,7)) ~ 2, # menstruate
    (menstruation_status %in% c(4,5,6)) ~ 3, # no menstruate
    (menstruation_status==8) ~ 4, # no conclusion
    TRUE ~ NA_real_
  ))

heatmap_data$Value <- case_when(
  heatmap_data$imputation_identifier == 1 ~ "Imputed",
  heatmap_data$imputation_identifier == 2 ~ "Menstruate",
  heatmap_data$imputation_identifier == 3 ~ "No Menstruate",
  heatmap_data$imputation_identifier == 4 ~ "Empty",
  TRUE ~ "Missing"
)
heatmap_data$Value <- factor(heatmap_data$Value, levels = c("Imputed", "Menstruate", "No Menstruate", "Empty"))
heatmap_data <- heatmap_data %>% 
  group_by(biome_id) %>%
  mutate(menstruating_days_count = sum((Value == "Imputed" | Value == "Menstruate"), na.rm = TRUE)) %>%
  arrange(desc(menstruating_days_count), biome_id) %>% 
  ungroup()
ggplot(heatmap_data, aes(x = logDate, y = reorder(factor(biome_id), menstruating_days_count), fill = Value)) +
  geom_tile(color = "gray25") +
  scale_fill_manual(values = c("No Menstruate" = "black", 
                               "Imputed" = "lightpink", 
                               "Menstruate" = "red",
                               "Empty"="white"),
                    na.value = "gray") +
  labs(title = "Menstruation Heatmap by Identifier", 
       x = "Date", 
       y = "Biome ID", 
       fill = "Menstruation Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

