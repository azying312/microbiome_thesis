# Merge the samples data set with each of the data sets that has one observation for each person (just the medical history?). 
# For each of the other data sets, check whether there are any IDs that appear >1 time per day. 
#If yes (eg, diet), ignore it. If each ID appears max 1 time per day, merge it in.
# Do not include the interpolated DASS values, only original data.
# Make a list of which data sets are included in your new data set and which are excluded.
# In every merge, hold on to all sample rows, even if there is no match in the new data set
#(all.x=TRUE, if the sample data set is first and a new data set is second), 
#but don't hold on to any rows that don't match the ID and date of a sample (all.y=FALSE).

library(dplyr)

### Report 8: Samples
samples_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Report 8-Sample.csv")
id_mapping <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Original Study Mapping - Sheet3.csv", header = TRUE)

### Map uid to study id

# Get unique study and biome health pairings
study_and_u_id <- unique(id_mapping %>% 
                           select(STUDY.ID, Biome.Health.App.ID))
# Match and join columns
study_and_u_id <- study_and_u_id %>% 
  rename("study_id" = "STUDY.ID") %>% 
  rename("biome_id" = "Biome.Health.App.ID")
samples_data <- samples_data %>% 
  rename("biome_id" = "uid")
study_and_u_id$study_id <- as.character(study_and_u_id$study_id)

# Map ids
samples_data <- samples_data %>%
  left_join(study_and_u_id, by = "biome_id") %>%
  mutate(biome_id = coalesce(study_id, biome_id)) %>%
  select(-study_id)

# Filter out saliva samples
samples_data <- samples_data %>% 
  filter(sampleType!="saliva")

# Check missing ids
missing_list <- samples_data %>%
  filter(is.na(as.numeric(biome_id)))
print(unique(missing_list$biome_id))

### Report 1: Menses
menses_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Report 1-Menstruation - Report 1-Menstruation.csv")
# id_mapping <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Original Study Mapping - Sheet3.csv", header = TRUE)

# View(menses_data)

### Map uid to study id

# Get unique study and biome health pairings
study_and_u_id <- unique(id_mapping %>% 
                           select(STUDY.ID, Biome.Health.App.ID))
# Match and join columns
study_and_u_id <- study_and_u_id %>% 
  rename("study_id" = "STUDY.ID") %>% 
  rename("biome_id" = "Biome.Health.App.ID")
menses_data <- menses_data %>% 
  rename("biome_id" = "uid")
study_and_u_id$study_id <- as.character(study_and_u_id$study_id)

# Map ids
menses_data <- menses_data %>%
  left_join(study_and_u_id, by = "biome_id") %>%
  mutate(biome_id = coalesce(study_id, biome_id)) %>%
  select(-study_id)

## Check if more than 1 a day for each participant
menses_data %>%
  group_by(biome_id, logDate) %>%           
  summarise(count = n()) %>%
  filter(count > 1) 

### Merge Report 8 and 1
full_data <- samples_data %>%
  left_join(menses_data, by=c('biome_id', 'logDate'))
dim(full_data)
# Dropped from med_data
dropped_rows <- menses_data %>%
  anti_join(full_data, by = c('biome_id', 'logDate'))
dim(dropped_rows)
# View(dropped_rows)

### Report 2 - Medications
med_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Report 2-Medications.csv")
# id_mapping <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Original Study Mapping - Sheet3.csv", header = TRUE)

### Map uid to study id

# Get unique study and biome health pairings
study_and_u_id <- unique(id_mapping %>% 
                           select(STUDY.ID, Biome.Health.App.ID))
# Match and join columns
study_and_u_id <- study_and_u_id %>% 
  rename("study_id" = "STUDY.ID") %>% 
  rename("biome_id" = "Biome.Health.App.ID")
med_data <- med_data %>% 
  rename("biome_id" = "uid")
study_and_u_id$study_id <- as.character(study_and_u_id$study_id)

# Map ids
med_data <- med_data %>%
  left_join(study_and_u_id, by = "biome_id") %>%
  mutate(biome_id = coalesce(study_id, biome_id)) %>%
  select(-study_id)

# Check missing ids
missing_list <- med_data %>%
  filter(is.na(as.numeric(biome_id)))
print(unique(missing_list$biome_id))

## Check if more than 1 a day for each participant
med_data %>%
  group_by(biome_id, logDate) %>%           
  summarise(count = n()) %>%
  filter(count > 1) 
names(med_data)

### Merge Full and Medications
full_data <- full_data %>%
  left_join(med_data, by=c('biome_id', 'logDate'))
dim(full_data)
# Dropped from med_data
dropped_rows <- med_data %>%
  anti_join(full_data, by = c('biome_id', 'logDate'))
dim(dropped_rows)
# View(dropped_rows)

### Report 3 - Sexual Activity
sex_act_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Report 3-Sexual Activity.csv")
# id_mapping <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Original Study Mapping - Sheet3.csv", header = TRUE)

### Map uid to study id

# Get unique study and biome health pairings
study_and_u_id <- unique(id_mapping %>% 
                           select(STUDY.ID, Biome.Health.App.ID))
# Match and join columns
study_and_u_id <- study_and_u_id %>% 
  rename("study_id" = "STUDY.ID") %>% 
  rename("biome_id" = "Biome.Health.App.ID")
sex_act_data <- sex_act_data %>% 
  rename("biome_id" = "uid")
study_and_u_id$study_id <- as.character(study_and_u_id$study_id)

# Map ids
sex_act_data <- sex_act_data %>%
  left_join(study_and_u_id, by = "biome_id") %>%
  mutate(biome_id = coalesce(study_id, biome_id)) %>%
  select(-study_id)

# Check missing ids
missing_list <- sex_act_data %>%
  filter(is.na(as.numeric(biome_id)))
print(unique(missing_list$biome_id))

dim(sex_act_data)

## Check if more than 1 a day for each participant
sex_act_data %>%
  group_by(biome_id, logDate) %>%           
  summarise(count = n()) %>%
  filter(count > 1) 
names(sex_act_data)

### Merge Full and Medications
full_data <- full_data %>%
  left_join(sex_act_data, by=c('biome_id', 'logDate'))
dim(full_data)
# Dropped from sex_act_data
dropped_rows <- sex_act_data %>%
  anti_join(full_data, by = c('biome_id', 'logDate'))
dim(dropped_rows)
# View(dropped_rows)

### Report 4-Physical Activity
activity.data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Report 4-Physical Activity.csv")
# id_mapping <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Original Study Mapping - Sheet3.csv", header = TRUE)

### Mapping IDs
study_and_u_id <- unique(id_mapping %>% 
                           select(STUDY.ID, Biome.Health.App.ID))
study_and_u_id <- study_and_u_id %>% 
  rename("study_id" = "STUDY.ID") %>% 
  rename("biome_id" = "Biome.Health.App.ID")

activity.data <- activity.data %>% 
  rename("biome_id" = "uid")

study_and_u_id$study_id <- as.character(study_and_u_id$study_id)

studyID_activity_data <- activity.data %>%
  left_join(study_and_u_id, by = "biome_id") %>%
  mutate(biome_id = coalesce(study_id, biome_id)) %>%
  select(-study_id)

### Some IDs are missing
missing_list <- studyID_activity_data %>%
  filter(is.na(as.numeric(biome_id)))
print(unique(missing_list$biome_id))

## Check if more than 1 a day for each participant
studyID_activity_data <- studyID_activity_data %>% 
  rename(logDate=activity_date)
studyID_activity_data %>%
  group_by(biome_id, logDate) %>%           
  summarise(count = n()) %>%
  filter(count > 1) 
names(studyID_activity_data)

### Merge Full and Medications
full_data <- full_data %>%
  left_join(studyID_activity_data, by=c('biome_id', 'logDate'))
dim(full_data)
# Dropped from sex_act_data
dropped_rows <- studyID_activity_data %>%
  anti_join(full_data, by = c('biome_id', 'logDate'))
dim(dropped_rows)
# View(dropped_rows)

### Report 5 - Sleep
sleep_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Report 5-Sleep.csv")
# id_mapping <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Original Study Mapping - Sheet3.csv", header = TRUE)

### Map uid to study id

# Get unique study and biome health pairings
study_and_u_id <- unique(id_mapping %>% 
                           select(STUDY.ID, Biome.Health.App.ID))
# Match and join columns
study_and_u_id <- study_and_u_id %>% 
  rename("study_id" = "STUDY.ID") %>% 
  rename("biome_id" = "Biome.Health.App.ID")
sleep_data <- sleep_data %>% 
  rename("biome_id" = "uid")
study_and_u_id$study_id <- as.character(study_and_u_id$study_id)

# Map ids
sleep_data <- sleep_data %>%
  left_join(study_and_u_id, by = "biome_id") %>%
  mutate(biome_id = coalesce(study_id, biome_id)) %>%
  select(-study_id)

# Check missing ids
missing_list <- sleep_data %>%
  filter(is.na(as.numeric(biome_id)))
print(unique(missing_list$biome_id))

## Check if more than 1 a day for each participant
names(sleep_data)
sleep_data <- sleep_data %>% 
  mutate(logDate=as.Date(start_time))
sleep_data %>%
  group_by(biome_id, logDate) %>%           
  summarise(count = n()) %>%
  filter(count > 1) 

multi_sleep_data <- sleep_data %>%
  group_by(biome_id, logDate) %>% 
  filter(n() > 1) %>%     
  ungroup()
# View(sleep_data)


### Report 9 - Volunteer Medical History
history_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Report 9-Volunteer Medical History.csv")
# id_mapping <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Original Study Mapping - Sheet3.csv", header = TRUE)

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

## Check if more than 1 a day for each participant
names(history_data)
history_data <- history_data %>% 
  mutate(logDate=as.character(as.Date(Timestamp, format = "%m/%d/%Y %H:%M:%S")))
history_data %>%
  group_by(biome_id, logDate) %>%           
  summarise(count = n()) %>%
  filter(count > 1) 
history_data <- history_data %>% 
  rename(surveyDate=logDate)

### Merge Full and Medications
full_data <- full_data %>%
  left_join(history_data, by=c('biome_id'))
dim(full_data)
# Dropped from history_data
dropped_rows <- history_data %>%
  anti_join(full_data, by = c('biome_id'))
dim(dropped_rows)
# View(dropped_rows)

### Report 10 - DASS Data
dass <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Report 10-DASS-21.csv")

### Map uid to study id
# Get unique study and biome health pairings
study_and_u_id <- unique(id_mapping %>% 
                           select(STUDY.ID, Biome.Health.App.ID))
# Match and join columns
study_and_u_id <- study_and_u_id %>% 
  rename("study_id" = "STUDY.ID") %>% 
  rename("biome_id" = "Biome.Health.App.ID")
dass <- dass %>% 
  rename("biome_id" = "Your.Biome.Health.App.ID")
study_and_u_id$study_id <- as.character(study_and_u_id$study_id)

# Map ids
dass <- dass %>%
  left_join(study_and_u_id, by = "biome_id") %>%
  mutate(biome_id = coalesce(study_id, biome_id)) %>%
  select(-study_id)

# Check missing ids
missing_list <- dass %>%
  filter(is.na(as.numeric(biome_id)))
print(unique(missing_list$biome_id))

## Check if more than 1 a day for each participant
names(dass)
dass <- dass %>% 
  mutate(logDate=as.character(as.Date(Timestamp, format = "%m/%d/%Y %H:%M:%S")))
dass %>%
  group_by(biome_id, logDate) %>%           
  summarise(count = n()) %>%
  filter(count > 1) 

multi_dass <- dass %>%
  group_by(biome_id, logDate) %>% 
  filter(n() > 1) %>%     
  ungroup()

# View(multi_dass)

dim(full_data)
## replace all blank with NA

is.na(full_data) <- full_data == ""

### Save final data output
write.csv(full_data,
          file = "/Users/alicezhang/Desktop/microbiome_data/meta_data_mayo.csv",
          row.names = FALSE)


