########################
#
# Clean FitBit Data
# v1: 3 October 2024
#
#########################

library(tidyverse)

activity.data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Report 4-Physical Activity.csv")
id_mapping <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Original Study Mapping - Sheet3.csv", header = TRUE)

head(activity.data)
head(id_mapping)

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

# Only keep numeric IDs
studyID_activity_data$biome_id <- as.numeric(studyID_activity_data$biome_id)
studyID_activity_data <- studyID_activity_data %>% 
  filter(!is.na(biome_id)) %>% 
  select(-id)

### Cleaning
studyID_activity_data$activity_date <- as.Date(studyID_activity_data$activity_date, format='%Y-%m-%d')
studyID_activity_data$calories_burned <- as.numeric(gsub(",", "", studyID_activity_data$calories_burned))
studyID_activity_data$steps <- as.numeric(gsub(",", "", studyID_activity_data$steps))
studyID_activity_data$minutes_sedentary <- as.numeric(gsub(",", "", studyID_activity_data$minutes_sedentary))
studyID_activity_data$activity_calories <- as.numeric(gsub(",", "", studyID_activity_data$activity_calories))

### Save final data output
write.csv(studyID_activity_data,
          file = "/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_Report 4-Physical Activity.csv",
          row.names = FALSE)

# Week Variable
# studyID_activity_data$Timestamp <- as.Date(studyID_activity_data$activity_date, format="%m/%d/%y", tz="UTC")
# time_values <- sort(unique(studyID_activity_data$Timestamp))
# num_time_values <- as.numeric(time_values)
# 
# studyID_activity_data$Timestamp[1] <= 19280
# 
# studyID_activity_data$week <- rep(NA, nrow(studyID_activity_data))
# studyID_activity_data$week[studyID_activity_data$Timestamp >= 19272 & studyID_activity_data$Timestamp <= 19280] <- 1
# studyID_activity_data$week[studyID_activity_data$Timestamp >= 19281 & studyID_activity_data$Timestamp <= 19286] <- 2
# studyID_activity_data$week[studyID_activity_data$Timestamp >= 19287 & studyID_activity_data$Timestamp <= 19293] <- 3
# studyID_activity_data$week[studyID_activity_data$Timestamp >= 19294 & studyID_activity_data$Timestamp <= 19300] <- 4
# studyID_activity_data$week[studyID_activity_data$Timestamp >= 19301 & studyID_activity_data$Timestamp <= 19308] <- 5
# studyID_activity_data$week[studyID_activity_data$Timestamp >= 19309 & studyID_activity_data$Timestamp <= 19315] <- 6
# studyID_activity_data$week[studyID_activity_data$Timestamp >= 19316 & studyID_activity_data$Timestamp <= 19321] <- 7
# studyID_activity_data$week[studyID_activity_data$Timestamp >= 19322 & studyID_activity_data$Timestamp <= 19329] <- 8
# studyID_activity_data$week[studyID_activity_data$Timestamp >= 19330 & studyID_activity_data$Timestamp <= 19336] <- 9
# studyID_activity_data$week[studyID_activity_data$Timestamp >= 19337] <- 10

# NA week?
# table(studyID_activity_data[is.na(studyID_activity_data$week),]$biome_id)
# Table of weeks
# table(studyID_activity_data$week)

# Reorder values by date
studyID_activity_data <- studyID_activity_data[order(studyID_activity_data$activity_date),]

id_values <- unique(studyID_activity_data$biome_id)

dim(studyID_activity_data)
View(studyID_activity_data)

### Data Analysis
# Plot calories burned
par(cex.main=1.5) 
plot(studyID_activity_data$activity_date[studyID_activity_data$biome_id==1], studyID_activity_data$calories_burned[studyID_activity_data$biome_id==1], type="l",
     main = "Calories Burned over Time", ylab = "Calories", xlab = "Date", xlim = c(19272.00, 19351.00), 
     ylim = c(154, 4963), 
     col="white")

for (id in id_values){
  points(studyID_activity_data$activity_date[studyID_activity_data$biome_id==id], studyID_activity_data$calories_burned[studyID_activity_data$biome_id==id], 
         type="l", lwd=0.3, col="darkgrey")}

sp <- smooth.spline(y = studyID_activity_data$calories_burned, x = studyID_activity_data$activity_date)
points(sp$x, sp$y, col = "deeppink4", type="o")


#PLOT ON STEPS
plot(d$activity_date[d$study.id==1], d$steps[d$study.id==1], type="l",
     main = "Steps over Time", ylab= "Steps Taken", xlab = "Date", xlim = c(19272.00, 19351.00), 
     ylim = c(0, 54675), col="white")

for (id in id_values){
  points(d$activity_date[d$study.id==id], d$steps[d$study.id==id], 
         type="l", lwd=0.3, col="darkgrey")}

sp1 <- smooth.spline(y = d$steps, x = d$activity_date)
points(sp1$x, sp1$y, col = "deeppink4", type="o")