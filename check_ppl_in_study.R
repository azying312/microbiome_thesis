library(tidyverse)

samples_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Report 8-Sample.csv")
history_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Report 9-Volunteer Medical History.csv")

# In samples, no volunteer survey
setdiff(samples_data$uid, history_data$Your.Biome.Health.App.ID)
# [1] "ef08154dacea08250a1f944befd52a26" "3dee14d32d23c3c1b1911fa236b0c874" "ae7cdadf86dad503089cdd9b623b1447"
# [4] "6411edfb62f428083dd35681f0862698" "4b4c694a838af7ccc89056289bfb3522" "d57c7a6a18fcc23860d951edb9a67a55"
# [7] "f4e1e88ac8f0db5923f0f72782e04d14" "fe199a72bf18ff0881a02c35a84f6a9a"

# In history data, not in sample IDs
setdiff(history_data$Your.Biome.Health.App.ID, samples_data$uid)
# [1] "3b540366ee825a64aa92fff4f98571f7" "dadd3e52822f17796ccce484ff487305" "5683d28b8ea12d4344a109885dfb19c2"

cleaned_history_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_Report 9-Volunteer Medical History.csv")
cleaned_samples_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_samples.csv")

# has survey, no samples
setdiff(cleaned_history_data$biome_id, cleaned_samples_data$biome_id)
# have samples, but no survey
setdiff(cleaned_samples_data$biome_id, cleaned_history_data$biome_id)

### Find samples for ppl w/o surveys

samples_data[which(samples_data$uid=="ef08154dacea08250a1f944befd52a26"),]
samples_data[which(samples_data$uid=="3dee14d32d23c3c1b1911fa236b0c874"),]
samples_data[which(samples_data$uid=="ae7cdadf86dad503089cdd9b623b1447"),]
samples_data[which(samples_data$uid=="6411edfb62f428083dd35681f0862698"),]
samples_data[which(samples_data$uid=="4b4c694a838af7ccc89056289bfb3522"),]
samples_data[which(samples_data$uid=="d57c7a6a18fcc23860d951edb9a67a55"),]
samples_data[which(samples_data$uid=="f4e1e88ac8f0db5923f0f72782e04d14"),]
samples_data[which(samples_data$uid=="fe199a72bf18ff0881a02c35a84f6a9a"),]

### Check timestamps on multiple sample submissions
dim(samples_data)

vaginal_samples <- samples_data %>% 
  filter(sampleType=="vaginal")
dim(vaginal_samples) # 4179
fecal_samples <- samples_data %>% 
  filter(sampleType=="fecal")

# remove ppl to drop
drop_participants <- c("ef08154dacea08250a1f944befd52a26", "ae7cdadf86dad503089cdd9b623b1447",
                       "4b4c694a838af7ccc89056289bfb3522", "fe199a72bf18ff0881a02c35a84f6a9a")
samples_data_subset <- samples_data %>% 
  filter(!uid %in% drop_participants)

# get logDate from timeStamp
samples_data_subset$logDate <- as.character(as.Date(samples_data_subset$timestamp, format = "%Y-%m-%d %H:%M:%S"))

###### vaginal samples

# check multiple submissions per day
multiple_vaginal_submission <- samples_data_subset %>% 
  group_by(uid, logDate) %>% 
  filter(sampleType=="vaginal") %>% 
  filter(logDate != "0000-00-00") %>% 
  summarise(Entries = n(), .groups = "drop") %>% 
  filter(Entries > 1)

multiple_vaginal_uids <- unique(multiple_vaginal_submission$uid)

multiple_vaginal_submission <- samples_data_subset %>% 
  group_by(uid, logDate) %>% 
  filter(sampleType=="vaginal") %>% 
  filter(logDate != "0000-00-00") %>% 
  filter(uid %in% multiple_vaginal_uids)

uminn_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Swabs with blood - Sheet1.csv", header=TRUE)
uminn_data_qr <- sapply(strsplit(uminn_data$Sample.ID, "_"), '[', 1)
vaginal_data_qr <- multiple_vaginal_submission$qr

# in vaginal qr, not in uminn
vaginal_notUminn <- setdiff(vaginal_data_qr, uminn_data_qr)
length(vaginal_notUminn)

multiple_vaginal_submission <- samples_data_subset %>% 
  group_by(uid, logDate) %>% 
  filter(sampleType=="vaginal") %>% 
  filter(logDate != "0000-00-00") %>% 
  filter(qr %in% vaginal_notUminn)

# Uminn having the multiple samples
multiple_vaginal_submission <- samples_data_subset %>% 
  group_by(uid, logDate) %>% 
  filter(sampleType=="vaginal") %>% 
  filter(logDate != "0000-00-00") %>% 
  filter(!qr %in% vaginal_notUminn) %>% 
  mutate(Entries = n()) %>% 
  filter(Entries > 1)

View(multiple_vaginal_submission)

# Uminn doesn't have the multiple samples
multiple_vaginal_submission <- samples_data_subset %>% 
  group_by(uid, logDate) %>% 
  filter(sampleType=="vaginal") %>% 
  filter(logDate != "0000-00-00") %>% 
  filter(qr %in% vaginal_notUminn) %>% 
  mutate(Entries = n()) %>% 
  filter(Entries > 1)

# View(multiple_vaginal_submission)

multiple_vaginal_submission_notUminn <- multiple_vaginal_submission %>% 
  group_by(uid, logDate) %>% 
  summarise(Entries = n(), .groups = "drop") %>% 
  filter(Entries > 1)

# View(multiple_vaginal_submission)

###### fecal samples

# check multiple submissions per day
multiple_fecal_submission <- samples_data_subset %>% 
  group_by(uid, logDate) %>% 
  filter(sampleType=="fecal") %>% 
  filter(logDate != "0000-00-00") %>% 
  summarise(Entries = n(), .groups = "drop") %>% 
  filter(Entries > 1)

multiple_fecal_uids <- unique(multiple_fecal_submission$uid)

multiple_fecal_submission <- samples_data_subset %>% 
  group_by(uid, logDate) %>% 
  filter(sampleType=="fecal") %>% 
  filter(logDate != "0000-00-00") %>% 
  filter(uid %in% multiple_fecal_uids)

fecal_data_qr <- multiple_fecal_submission$qr

# in fecal qr, not in uminn
fecal_notUminn <- setdiff(fecal_data_qr, uminn_data_qr)
length(fecal_notUminn)

multiple_fecal_submission <- samples_data_subset %>% 
  group_by(uid, logDate) %>% 
  filter(sampleType=="fecal") %>% 
  filter(logDate != "0000-00-00") %>% 
  filter(qr %in% fecal_notUminn)

# View(multiple_fecal_submission)

multiple_fecal_submission <- samples_data_subset %>% 
  group_by(uid, logDate) %>% 
  filter(sampleType=="fecal") %>% 
  filter(logDate != "0000-00-00") %>% 
  filter(!qr %in% fecal_notUminn)


