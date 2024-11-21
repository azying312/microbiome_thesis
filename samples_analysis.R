library(tidyverse)

samples.data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_samples.csv")
uminn_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Swabs with blood - Sheet1.csv", header=TRUE)
survey_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/alice_cleaned_survey_data.csv", header=TRUE)

table(samples.data$sampleType)

uminn_data_subset <- uminn_data %>% 
  select(Sample.ID, Well, Special.Notes) 
uminn_data_subset$qr <- sub("_.*", "", uminn_data_subset$Sample.ID)
uminn_data_subset$inUminn <- TRUE

## vaginal samples
vaginal_data <- samples.data %>% 
  filter(sampleType=="vaginal") %>% 
  filter(!is.na(biome_id)) %>%  
  filter(logDate != "0000-00-00")
dim(vaginal_data) # 4179 -> 4026

vaginal_data$vaginal_swab <- TRUE

vaginal_data <- uminn_data_subset %>% 
  left_join(vaginal_data, by="qr") # 1558 actual samples
 
## Filter out any errors
vaginal_data <- vaginal_data %>% 
  filter(logDate!="0000-00-00") %>% # no corresponding logDate
  mutate(uMinn_menstruation=ifelse(str_detect(Special.Notes, "Blood")==TRUE, 1, 0)) %>%  # set to menstruation true
  select(biome_id, logDate, qr, timestamp, uMinn_menstruation, inUminn, Special.Notes) %>% 
  mutate(biome_id=as.numeric(biome_id)) %>% 
  # errors in sample processing from UMinn
  filter(!str_detect(Special.Notes, "error")) %>% 
  filter(!str_detect(Special.Notes, "No swab in tube")) %>% 
  filter(!is.na(biome_id))

## Get duplicates

vaginal_data_diff <- vaginal_data %>%
  mutate(logDate = as.Date(timestamp),
         timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S")) %>% 
  group_by(biome_id, logDate) %>%
  mutate(count=n()) %>% 
  arrange(timestamp, .by_group = TRUE) %>%
  mutate(
    time_diff = as.numeric(difftime(timestamp, lag(timestamp), units = "mins"))
  ) %>%
  ungroup()

identical_vaginal_data <- vaginal_data_diff %>% 
  group_by(biome_id, logDate) %>%
  filter(n() == 2) %>%
  mutate(UMinn_mismatch = uMinn_menstruation[1] != uMinn_menstruation[2]) %>%
  filter(UMinn_mismatch == TRUE) %>% 
  mutate(
    time_diff = as.numeric(difftime(timestamp, lag(timestamp), units = "mins"))
  ) %>%
  # filter(n() > 1) %>%
  ungroup() 

# View(identical_vaginal_data)

# Filtering
count_vaginal_data_diff <- vaginal_data_diff %>% 
  group_by(biome_id) %>% 
  mutate(count=n()) %>% 
  filter(count > 11) %>% 
  ungroup()
dim(count_vaginal_data_diff)
length(unique(vaginal_data_diff$biome_id)) # 67
length(unique(count_vaginal_data_diff$biome_id)) # 49

## Add surveys
survey_data_ids <- unique(survey_data$biome_id)

vaginal_data_diff_wSurvey <- count_vaginal_data_diff %>% 
  filter(biome_id %in% survey_data_ids)

length(unique(vaginal_data_diff_wSurvey$biome_id))
dim(vaginal_data_diff_wSurvey)

## fecal samples
fecal_data <- samples.data %>% 
  filter(sampleType=="fecal") %>% 
  filter(!is.na(biome_id)) %>%  
  filter(logDate != "0000-00-00")
dim(fecal_data) # 2911 -> 2842

fecal_data$fecal_swab <- TRUE

fecal_data <- uminn_data_subset %>% 
  left_join(fecal_data, by="qr") 

## Filter out any errors
fecal_data <- fecal_data %>% 
  filter(logDate!="0000-00-00") %>% # no corresponding logDate
  select(biome_id, logDate, qr, timestamp, inUminn, Special.Notes) %>% 
  mutate(biome_id=as.numeric(biome_id)) %>% 
  # errors in sample processing from UMinn
  filter(!str_detect(Special.Notes, "error")) %>% 
  filter(!str_detect(Special.Notes, "No swab in tube")) %>% 
  filter(!is.na(biome_id)) # 1094 actual samples

## Get duplicates

fecal_data_diff <- fecal_data %>%
  mutate(logDate = as.Date(timestamp),
         timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S")) %>% 
  group_by(biome_id, logDate) %>%
  mutate(count=n()) %>% 
  arrange(timestamp, .by_group = TRUE) %>%
  mutate(
    time_diff = as.numeric(difftime(timestamp, lag(timestamp), units = "mins"))
  ) %>%
  ungroup()


#### Both vaginal and fecal
vaginal_data_unique <- vaginal_data %>% 
  group_by(biome_id, logDate) %>% 
  mutate(count=n()) %>% 
  filter(
    # keep menses swab
    (uMinn_menstruation == 1 & row_number() == 1) |
      # otherwise keep first swab
      (uMinn_menstruation != 1 & row_number() == 1)
  ) %>%
  ungroup() %>% 
  select(biome_id, logDate, uMinn_menstruation, vaginal_swab)
  
fecal_data_unique <- fecal_data %>% 
  group_by(biome_id, logDate) %>% 
  mutate(count=n()) %>% 
  filter(row_number() == 1) %>%
  ungroup() %>% 
  select(biome_id, logDate, fecal_swab)

vag_and_fec <- vaginal_data_unique %>% 
  left_join(fecal_data_unique) %>% 
  group_by(biome_id, logDate) %>% 
  filter(vaginal_swab==TRUE & fecal_swab==TRUE) %>% 
  ungroup()

length(unique(vag_and_fec$biome_id))
dim(vag_and_fec)

vag_and_fec <- vag_and_fec %>% 
  group_by(biome_id) %>% 
  mutate(count=n()) %>% 
  filter(count>12)
length(unique(vag_and_fec$biome_id))
dim(vag_and_fec)

## Add surveys
survey_data_ids <- unique(survey_data$biome_id)

vag_and_fec <- vag_and_fec %>% 
  filter(biome_id %in% survey_data_ids)

length(unique(vag_and_fec$biome_id))
dim(vag_and_fec)
