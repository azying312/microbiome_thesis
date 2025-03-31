library(dplyr)
rm(list=ls()) # Clearing all saved objects in R workspace

dass <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Report 10-DASS-21.csv")

sample.data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Report 8-Sample.csv")
# View(sample.data)
# dass$Timestamp
# posix_date_time <- as.POSIXct(dass$Timestamp, format = "%m-%d-%Y %H:%M:%S")
# studyID_bitesnap_data$times <- format(posix_date_time, format = "%H:%M:%OS")
# studyID_bitesnap_data$Date <- format(posix_date_time, format = "%Y-%m-%d")
# studyID_bitesnap_data$food_other <- FALSE

# doesn't process correct
dass$Timestamp <- as.Date(dass$Timestamp, format="%m/%d/%y", tz="UTC")
id_values <- unique(dass$study_id) 
time_values <- sort(unique(dass$Timestamp))
num_time_values <- as.numeric(time_values)
dass <- dass[order(dass$Timestamp),]

#Remove any values that occurred after 12/16 (end of semester)
dass <- dass[dass$Timestamp<= 19342, ]

#remove people with unidentifiable study ID
dass <- dass[dass$study_id!="#N/A", ] 
dass$study_id <- as.factor(dass$study_id)

#Below we created a "week" variable to organize entries into 1 of 10 weeks.
dass$week <- rep(NA, nrow(dass))
dass$week[dass$Timestamp >= 19276 & dass$Timestamp <= 19280] <- 1
dass$week[dass$Timestamp >= 19281 & dass$Timestamp <= 19286] <- 2
dass$week[dass$Timestamp >= 19290 & dass$Timestamp <= 19293] <- 3
dass$week[dass$Timestamp >= 19296 & dass$Timestamp <= 19300] <- 4
dass$week[dass$Timestamp >= 19302 & dass$Timestamp <= 19308] <- 5
dass$week[dass$Timestamp >= 19312 & dass$Timestamp <= 19315] <- 6
dass$week[dass$Timestamp >= 19319 & dass$Timestamp <= 19321] <- 7
dass$week[dass$Timestamp >= 19323 & dass$Timestamp <= 19329] <- 8
dass$week[dass$Timestamp >= 19330 & dass$Timestamp <= 19336] <- 9
dass$week[dass$Timestamp >= 19339] <- 10
