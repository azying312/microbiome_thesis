library(tidyverse)
library(viridis)

volunteer.survey <- read.csv("/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_Report 9-Volunteer Medical History.csv")
menses.data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_Report 1-Menstruation.csv")
med.data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_med_data.csv")
sexact.data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_Report 3-Sexual Activity.csv")
physical.data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_Report 4-Physical Activity.csv")
sleep.data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_sleep.csv")
diet.data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_diet.csv")
samples.data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_samples.csv")
DASS.data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/cleaned_data/DASS.csv")

## Create missingness grid
missingness_grid <- data.frame(matrix(NA, nrow = nrow(volunteer.survey), ncol = 11))

## 09 - Volunteer Survey
volunteer.survey <- volunteer.survey %>%
  mutate(across(where(is.character), ~ na_if(., "")))
volunteer.survey$Missingness <- round(rowMeans(!is.na(volunteer.survey)), 4) * 100
missingness_grid[, 1] <- volunteer.survey$Missingness
colnames(missingness_grid)[1] <- "Volunteer Survey"
rownames(missingness_grid) <- volunteer.survey$ParticipantID

## 01 - Menses Data
all_days <- seq.Date(as.Date(min(menses.data$logDate)), as.Date(max(menses.data$logDate)), by = "day")
length(all_days)
menses.data <- menses.data %>%
  group_by(biome_id) %>%
  summarise(Entries = n(), .groups = "drop") %>%
  mutate(Percentage = (Entries / 81) * 100)
# Get menses data
menses_percentages <- menses.data$Percentage
names(menses_percentages) <- menses.data$biome_id
# Join to missingness grid
missingness_grid[, 2] <- menses_percentages[rownames(missingness_grid)]
colnames(missingness_grid)[2] <- "Menses Data"

## 02 - Medications Data
med.data <- med.data %>%
  mutate(across(where(is.character), ~ na_if(., "")))
med.data <- med.data %>% 
  select(biome_id, other_medicine) %>% 
  mutate(med_recorded = ifelse(is.na(other_medicine), 0, 1))
# Get med data
med_recorded <- med.data$med_recorded
names(med_recorded) <- med.data$biome_id
# Join to missingness grid
missingness_grid[, 3] <- med_recorded[rownames(missingness_grid)]
colnames(missingness_grid)[3] <- "Med Data (binary)"

## 03 - Sexual Activity
cleaned_sexact.data <- sexact.data %>%
  filter(!is.na(logDate))

cleaned_sexact.data$logDate <- as.Date(cleaned_sexact.data$logDate)
all_days <- seq.Date(as.Date(min(cleaned_sexact.data$logDate, na.rm = TRUE)), 
                     as.Date(max(cleaned_sexact.data$logDate, na.rm = TRUE)), by = "day")
length(all_days)
cleaned_sexact.data <- cleaned_sexact.data %>%
  mutate(across(where(is.character), ~ na_if(., "")))

cleaned_sexact.data <- cleaned_sexact.data %>%
  group_by(biome_id) %>%
  summarise(Entries = n(), .groups = "drop") %>%
  mutate(Percentage = round((Entries / length(all_days)), 4) * 100)
# Get sex active data
cleaned_sexact_percentage <- cleaned_sexact.data$Percentage
names(cleaned_sexact_percentage) <- cleaned_sexact.data$biome_id
# Join to missingness grid
missingness_grid[, 4] <- cleaned_sexact_percentage[rownames(missingness_grid)]
colnames(missingness_grid)[4] <- "Sexual Activity (%)"

# Sex active binary
sexact.data.subset <- sexact.data %>% 
  select(biome_id, contraceptions) %>% 
  mutate(contraception_recorded = ifelse(is.na(contraceptions), 0, 1))
# Get sex act data
contraception_recorded <- sexact.data.subset$contraceptions
names(contraception_recorded) <- sexact.data.subset$biome_id
contraception_recorded[] <- 1
contraception_recorded <- as.numeric(contraception_recorded)
names(contraception_recorded) <- sexact.data.subset$biome_id

# Join to missingness grid
# missingness_grid[, 5] <- contraception_recorded[rownames(missingness_grid)]
# colnames(missingness_grid)[5] <- "Sexual Activity (binary)"

## 04 - Physical Activity
cleaned_physical.data <- physical.data %>%
  filter(!is.na(activity_date))

cleaned_physical.data$logDate <- as.Date(cleaned_physical.data$activity_date)
all_days <- seq.Date(as.Date(min(cleaned_physical.data$logDate)), as.Date(max(cleaned_physical.data$logDate)), by = "day")
length(all_days)
cleaned_physical.data <- cleaned_physical.data %>%
  mutate(across(where(is.character), ~ na_if(., "")))

# Percentage phys activity days
cleaned_physical.data <- cleaned_physical.data %>%
  group_by(biome_id) %>%
  summarise(Entries = n(), .groups = "drop") %>%
  mutate(Percentage = round((Entries / length(all_days)), 4) * 100)
# Get sex active data
cleaned_physical_percentage <- cleaned_physical.data$Percentage
names(cleaned_physical_percentage) <- cleaned_physical.data$biome_id
# Join to missingness grid
missingness_grid[, 6] <- cleaned_physical_percentage[rownames(missingness_grid)]
colnames(missingness_grid)[6] <- "Physical Activity (%)"

## 05 - Sleep
cleaned_sleep.data <- sleep.data %>%
  filter(!is.na(start_time))

cleaned_sleep.data$logDate <- as.Date(cleaned_sleep.data$start_time)
all_days <- seq.Date(as.Date(min(cleaned_sleep.data$logDate)), as.Date(max(cleaned_sleep.data$logDate)), by = "day")
length(all_days)
cleaned_sleep.data <- cleaned_sleep.data %>%
  mutate(across(where(is.character), ~ na_if(., "")))

# Percentage sleep days
cleaned_sleep.data <- cleaned_sleep.data %>%
  group_by(biome_id) %>%
  summarise(Entries = n(), .groups = "drop") %>%
  mutate(Percentage = round((Entries / length(all_days)), 4) * 100)
# Get sleep data
cleaned_sleep_percentage <- cleaned_sleep.data$Percentage
names(cleaned_sleep_percentage) <- cleaned_sleep.data$biome_id
# Join to missingness grid
missingness_grid[, 7] <- cleaned_sleep_percentage[rownames(missingness_grid)]
colnames(missingness_grid)[7] <- "Sleep Data (%)"

## 07 - Food

# If there is food data for the person - col 8
diet.data$logDate <- as.Date(diet.data$Date)

# Filter for dates before 12/17
# over_diet.data <- diet.data %>% 
#   filter(logDate > as.Date("2022-12-17"))
filtered_diet.data <- diet.data %>% 
  filter(logDate <= as.Date("2022-12-17"))

all_days <- seq.Date(as.Date(min(filtered_diet.data$logDate, na.rm=TRUE)),
                     as.Date(max(filtered_diet.data$logDate, na.rm=TRUE)), by = "day")
length(all_days)
filtered_diet.data <- filtered_diet.data %>%
  mutate(across(where(is.character), ~ na_if(., "")))
participants.diet <- unique(as.numeric(filtered_diet.data$study_id))

# % days participant submitted
diet_submission <- filtered_diet.data %>%
  # get distinct days & participants
  distinct(study_id, logDate) %>% 
  group_by(study_id) %>%
  summarise(Entries = n(), .groups = "drop") %>%
  mutate(Percentage = round((Entries / length(all_days)), 4) * 100)

# Get diet data
percent_diet_submission <- diet_submission$Percentage
names(percent_diet_submission) <- diet_submission$study_id
# Join to missingness grid
missingness_grid[, 8] <- percent_diet_submission[rownames(missingness_grid)]
colnames(missingness_grid)[8] <- "Diet Data (% days)"

# Col 05 - % meals submitted
# assume 3 meals a day
all_meals <- length(all_days)*3

# % meals participant submitted
diet_submission <- filtered_diet.data %>%
  filter(type %in% c("breakfast", "lunch", "dinner")) %>% 
  # get distinct days, meals & participants
  distinct(study_id, type, logDate) %>% 
  group_by(study_id, type) %>%
  summarise(Entries = n(), .groups = "drop") %>%
  mutate(Percentage = round((Entries / all_meals), 4) * 100)

# Get diet data
percent_diet_submission <- diet_submission$Percentage
names(percent_diet_submission) <- diet_submission$study_id
# Join to missingness grid
missingness_grid[, 5] <- percent_diet_submission[rownames(missingness_grid)]
colnames(missingness_grid)[5] <- "Diet Data (% meals)"


## Percentage Vegetarian - col 9 (of all cals, how many for each person are vegetarian)
# veg.diet.data <- diet.data %>%
#   group_by(study_id) %>% 
#   summarise(
#     totalcals = sum(caloriesall, na.rm = TRUE),
#     total_veg_cal = sum(caloriesall[vegetarian == TRUE], na.rm = TRUE),
#     perc_veg_cal = round((total_veg_cal / totalcals), 4) * 100
#   )
# 
# # Get list
# veg.diet.perc <- veg.diet.data$perc_veg_cal
# names(veg.diet.perc) <- veg.diet.data$study_id

# Join to missingness grid
# missingness_grid[, 9] <- veg.diet.perc[rownames(missingness_grid)]
# colnames(missingness_grid)[9] <- "Plant-based (%)"
# 

## 08 - Sample
cleaned_samples.data <- samples.data %>%
  filter(!is.na(logDate)) %>% 
  filter(logDate!="0000-00-00")

cleaned_samples.data$logDate <- as.Date(cleaned_samples.data$logDate)
all_days <- seq.Date(as.Date(min(cleaned_samples.data$logDate)), as.Date(max(cleaned_samples.data$logDate)), by = "day")
length(all_days)
cleaned_samples.data <- cleaned_samples.data %>%
  mutate(across(where(is.character), ~ na_if(., "")))

# % Vaginal - col 9
vaginal.samples <- cleaned_samples.data %>% 
  filter(sampleType=="vaginal")
all_days <- seq.Date(as.Date(min(vaginal.samples$logDate)), as.Date(max(vaginal.samples$logDate)), by = "day")
length(all_days)

# Percentage vaginal swab submissions
vaginal.samples <- vaginal.samples %>%
  group_by(biome_id) %>%
  # get distinct days
  distinct(biome_id, logDate) %>% 
  summarise(Entries = n(), .groups = "drop") %>%
  mutate(Percentage = round((Entries / length(all_days)), 4) * 100)

# Get vaginal sample data
cleaned_vag_percentage <- vaginal.samples$Percentage
names(cleaned_vag_percentage) <- vaginal.samples$biome_id
# Join to missingness grid
missingness_grid[,9] <- cleaned_vag_percentage[rownames(missingness_grid)]
colnames(missingness_grid)[9] <- "Vaginal (%)"

# % Gut - col 10
gut.samples <- cleaned_samples.data %>% 
  filter(sampleType=="fecal") %>% 
  select(biome_id, logDate, sampleType)
all_days <- seq.Date(as.Date(min(gut.samples$logDate)), as.Date(max(gut.samples$logDate)), by = "day")
length(all_days)

# Percentage fecal swab submissions
length(unique(as.numeric(gut.samples$biome_id)))

gut.samples <- gut.samples %>%
  group_by(biome_id) %>%
  # get distinct days
  distinct(biome_id, logDate) %>% 
  summarise(Entries = n(), .groups = "drop") %>%
  mutate(Percentage = round((Entries / length(all_days)), 4) * 100)
# Get fecal sample data
cleaned_gut_percentage <- gut.samples$Percentage
names(cleaned_gut_percentage) <- gut.samples$biome_id
# Join to missingness grid
missingness_grid[, 10] <- cleaned_gut_percentage[rownames(missingness_grid)]
colnames(missingness_grid)[10] <- "Gut (%)"

## 10 - DASS - col 12
DASS.data$biome_id <- DASS.data$study_id
posix_date_time <- as.POSIXct(DASS.data$Timestamp, format = "%m/%d/%y %H:%M")
DASS.data$Converted_Timestamp <- posix_date_time
DASS.data$logDate <- format(posix_date_time, format = "%Y-%m-%d")
cleaned_DASS.data <- DASS.data %>%
  filter(!is.na(logDate))

cleaned_DASS.data$logDate <- as.Date(cleaned_DASS.data$logDate)
all_days <- seq.Date(as.Date(min(cleaned_DASS.data$logDate)), as.Date("2022-12-15"), by = "day")

weeks<-10 #length(all_days)/7

cleaned_DASS.data <- cleaned_DASS.data %>%
  filter(logDate >= min(cleaned_DASS.data$logDate) & logDate <= as.Date("2022-12-15"))

cleaned_DASS.data <- cleaned_DASS.data %>%
  mutate(across(where(is.character), ~ na_if(., "")))

cleaned_DASS.data <- cleaned_DASS.data %>%
  group_by(biome_id) %>%
  summarise(Entries = n(), .groups = "drop") %>%
  mutate(Percentage = round((Entries / weeks), 4) * 100)

# Get DASS data
cleaned_DASS_percentage <- cleaned_DASS.data$Percentage
names(cleaned_DASS_percentage) <- cleaned_DASS.data$biome_id
# Join to missingness grid
missingness_grid[, 11] <- cleaned_DASS_percentage[rownames(missingness_grid)]
colnames(missingness_grid)[11] <- "DASS (%)"

### Missingness Plot

# Set NA to 0
missingness_grid[is.na(missingness_grid)] <- 0

# Convert the missingness matrix to a tidy format
missingness_df <- as.data.frame(missingness_grid)

missingness_df$biome_ids <- rownames(missingness_grid)
# Filter out binary rows 
missingness_df_cont <- missingness_df[,-c(3)]
length(names(missingness_df_cont))

# Reorder
missingness_df_cont <- missingness_df_cont %>% 
  select("Diet Data (% days)", "Diet Data (% meals)", # "Plant-based (%)",
         "Menses Data", "Vaginal (%)", "Gut (%)",
         "Physical Activity (%)", "Sleep Data (%)",
         "DASS (%)",
         "Sexual Activity (%)", "Volunteer Survey", everything())

# Reshape to long format
missingness_long <- pivot_longer(missingness_df_cont,
                                 # cols=everything(),
                                 cols = -biome_ids,
                                 names_to = "Variable", 
                                 values_to = "Missing")

## Reorder
missingness_long_count <- missingness_long %>% 
  group_by(biome_ids) %>%
  mutate(total_sum = sum(across(where(is.numeric)), na.rm = TRUE)) %>%
  # arrange rows
  arrange(desc(total_sum), biome_ids) %>% 
  ungroup() #%>%
  # arrange cols
  # select(-total_sum) %>%
  # group_by(Variable) %>%
  # mutate(columnSum = sum(Missing, na.rm = TRUE)) %>% 
  # arrange(desc(columnSum), biome_ids) %>% 
  # ungroup()

colsum_df <- missingness_long_count %>%
  group_by(Variable) %>%
  mutate(columnSum = sum(Missing, na.rm = TRUE))

missingness_long_count <- missingness_long_count %>%
  mutate(
    Variable = factor(Variable, levels = c("Volunteer Survey", "Physical Activity (%)", "Sleep Data (%)",
                                           "Diet Data (% days)", "Diet Data (% meals)",
                                           "Gut (%)", "Vaginal (%)",
                                           "Menses Data", "Sexual Activity (%)",
                                           "DASS (%)")),
    biome_ids = factor(biome_ids, levels = rev(unique(biome_ids)))
  )

# missingness_long_count <- missingness_long_count %>%
#   mutate(
#     Variable = factor(Variable, levels = c("Volunteer Survey", "DASS (%)",
#                                            "Physical Activity (%)", "Sleep Data (%)",
#                                            "Vaginal (%)", "Gut (%)", "Diet Data (% days)", "Diet Data (% meals)",
#                                            "Menses Data", "Sexual Activity (%)"
#                                            )),
#     biome_ids = factor(biome_ids, levels = rev(unique(biome_ids)))
#   )

# Missingness Plot
ggplot(missingness_long_count, aes(x = Variable, y = biome_ids, fill = Missing)) +
  geom_tile(color = "gray25") +
  scale_fill_viridis_c(option = "viridis", direction=-1,
                       na.value = "grey50") +
  labs(title = "Missingness Heatmap", 
       x = "\n Data Type", 
       y = "Study ID \n", 
       fill = "Percentage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) 

# Samples # Missingness Plot
missingness_long_samples<- missingness_long %>% 
  filter(Variable=="Gut (%)" | Variable=="Vaginal (%)")
ggplot(missingness_long_samples, aes(x = Variable, y = factor(biome_ids), fill = Missing)) +
  geom_tile(color = "gray25") +
  scale_fill_viridis_c(option = "viridis", direction=-1,
                       na.value = "grey50") +
  labs(title = "Samples Missingness Heatmap", 
       x = "Data Type", 
       y = "Study ID", 
       fill = "Percentage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) 

# DASS # Missingness Plot
# missingness_long_DASS<- missingness_long %>% 
#   filter(Variable=="DASS (%)")
# ggplot(missingness_long_DASS, aes(x = Variable, y = factor(biome_ids), fill = Missing)) +
#   geom_tile(color = "gray25") +
#   scale_fill_viridis_c(option = "viridis", direction=-1,
#                        na.value = "grey50") +
#   labs(title = "DASS Missingness Heatmap", 
#        x = "Data Type", 
#        y = "Study ID", 
#        fill = "Percentage") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) 

# # Missingness Plot - Binary
# ggplot(missingness_long_bin, aes(x = Variable, y = factor(biome_ids), fill = factor(Missing))) +
#   geom_tile(color = "black") +
#   scale_fill_manual(values = c("white", "black"), 
#                     na.value = "white") +  # Specify colors for binary values
#   labs(title = "Missingness Heatmap", 
#        x = "Variable", 
#        y = "Biome ID", 
#        fill = "Missingness") +  # Adjust labels as needed
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
