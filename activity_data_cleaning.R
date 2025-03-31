########################
#
# Clean FitBit Data
# v1: 3 October 2024
#
#########################

source("~/Microbiome Thesis/functions.R")

library(tidyverse)

activity.data <- read.csv("/Volumes/T7/microbiome_data/original_data/Report 4-Physical Activity.csv")
id_mapping <- read.csv("/Volumes/T7/microbiome_data/original_data/Original Study Mapping - Sheet3.csv", header = TRUE)


# Data Prep
activity.data <- activity.data %>%
  rename("biome_id" = "uid")
activity.data <- study_mapping(activity.data, id_mapping)

length(unique(activity.data$biome_id))

### Cleaning
activity.data <- activity.data %>% 
  rename(logDate=activity_date) %>% 
  select(!c(id, mod_timestamp))
activity.data$calories_burned <- as.numeric(gsub(",", "", activity.data$calories_burned))
activity.data$steps <- as.numeric(gsub(",", "", activity.data$steps))
activity.data$minutes_sedentary <- as.numeric(gsub(",", "", activity.data$minutes_sedentary))
activity.data$activity_calories <- as.numeric(gsub(",", "", activity.data$activity_calories))

activity.data$minutes_lightly_active <- as.numeric(activity.data$minutes_lightly_active)

### Save final data output
write.csv(activity.data,
          file = "/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 4-Physical Activity.csv",
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
activity.data <- activity.data[order(activity.data$logDate),]

id_values <- unique(activity.data$biome_id)

dim(activity.data)
# View(studyID_activity_data)

### Data Analysis
# Plot calories burned
par(cex.main=1.5) 
activity.data$logDate <- as.Date(activity.data$logDate)

plot(activity.data$logDate[activity.data$biome_id == 1], 
     activity.data$calories_burned[activity.data$biome_id == 1], 
     type = "l", 
     main = "Calories Burned over Time", 
     ylab = "Calories", 
     xlab = "Date", 
     col = "blue")

plot(activity.data$logDate[activity.data$biome_id==1], 
     activity.data$calories_burned[activity.data$biome_id==1], type="l",
     main = "Calories Burned over Time", ylab = "Calories", xlab = "Date",
     col="white")
for (id in id_values){
  points(activity.data$logDate[activity.data$biome_id==id], 
         activity.data$calories_burned[activity.data$biome_id==id], 
         type="l", lwd=0.3, col="darkgrey")}

sp <- smooth.spline(y = activity.data$calories_burned, x = activity.data$logDate)
points(sp$x, sp$y, col = "deeppink4", type="o")

plot(activity.data$logDate, activity.data$calories_burned, type = "n",
     main = "Calories Burned over Time by Participant", 
     ylab = "Calories", xlab = "Date")
participant_ids <- unique(activity.data$biome_id)

source("~/Microbiome Thesis/functions.R")
activity.data.filtered <- filter_days(activity.data)
dim(activity.data.filtered)
activity.data.filtered$logDate_numeric <- as.numeric(activity.data.filtered$logDate)
activity.data.filtered <- activity.data.filtered %>%
  filter(!is.na(logDate_numeric) & !is.na(calories_burned))

# figure shows that cals burned over time doesn't change rlly
ggplot(activity.data.filtered, aes(x = logDate, y = calories_burned, group = biome_id)) +
  #geom_line(color = "darkgrey", size = 0.3) +
  geom_smooth(aes(color = as.factor(biome_id)), method = "loess", se = FALSE, size = 0.8) +
  labs(title = "Calories Burned over Time by Participant",
       x = "Date", y = "Calories Burned") +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_date(breaks = unique(activity.data.filtered$logDate), date_labels = "%Y-%m-%d")

# could plot smthg else since the figure above isn't that interesting

#PLOT ON STEPS
plot(activity.data$logDate[activity.data$biome_id==1], activity.data$steps[activity.data$biome_id==1], type="l",
     main = "Steps over Time", ylab= "Steps Taken", xlab = "Date", xlim = c(19272.00, 19351.00), 
     ylim = c(0, 54675), col="white")

for (id in id_values){
  points(activity.data$logDate[activity.data$biome_id==id], activity.data$steps[activity.data$biome_id==id], 
         type="l", lwd=0.3, col="darkgrey")}

sp1 <- smooth.spline(y = activity.data$steps, x = activity.data$logDate)
points(sp1$x, sp1$y, col = "deeppink4", type="o")
