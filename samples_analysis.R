library(tidyverse)
library(ggplot2)

samples.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_samples.csv")
# uminn_data <- read.csv("/Volumes/T7/microbiome_data/Swabs with blood - Sheet1.csv", header=TRUE)
uminn_data <- read.csv("/Volumes/T7/microbiome_data/original_data/UMinn Samples - Sheet1.csv")
survey_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 9-Volunteer Medical History.csv", header = TRUE)

saliva <- samples.data %>%  
  filter(sampleType=="saliva")
table(samples.data$sampleType)

uminn_data_subset <- uminn_data %>% 
  select(Sample.ID, Well, Special.Notes) 
uminn_data_subset$qr <- sub("_.*", "", uminn_data_subset$Sample.ID)
uminn_data_subset$inUminn <- TRUE

## vaginal samples
vaginal_data <- samples.data %>% 
  filter(sampleType=="vaginal")
vaginal_data$vaginal_swab <- TRUE # 1536

# note that uminn_data_subset has 67 ppl and doesn't have participant 68 who is in the SILVA pipeline data
# issue is missing in the uminn data
vaginal_data <- uminn_data_subset %>% 
  left_join(vaginal_data, by="qr") %>%  # 1536 actual samples
  filter(vaginal_swab==TRUE) # there might be dupe QR codes

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

# Get time difference between samples on the same day for a given participant
identical_vaginal_data <- vaginal_data_diff %>% 
  group_by(biome_id, logDate) %>%
  filter(n() == 2) %>% 
  mutate(
    time_diff = as.numeric(difftime(timestamp, lag(timestamp), units = "mins"))
  ) %>%
  ungroup() 

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
  select(biome_id, logDate, qr, timestamp, inUminn, fecal_swab, Special.Notes) %>% 
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

## Distribution of vaginal and fecal samples
vag_and_fec <- vaginal_data_unique %>% 
  full_join(fecal_data_unique) %>% 
  mutate(has_both=ifelse((vaginal_swab==TRUE & fecal_swab==TRUE), 1, 0))

vag_and_fec_counts <- vag_and_fec %>%
  group_by(biome_id) %>%
  summarize(
    total_samples = sum(vaginal_swab, na.rm = TRUE) + sum(fecal_swab, na.rm = TRUE),
    total_vaginal = sum(vaginal_swab, na.rm = TRUE),
    total_fecal = sum(fecal_swab, na.rm = TRUE),
    has_both=sum(has_both,na.rm=TRUE),
    .groups = "drop"
  )

## potential cutoff: 15
all_sample_15 <- vag_and_fec_counts %>% 
  filter(has_both >= 15)
ids15 <- unique(all_sample_15$biome_id)

# vag_and_fec_46 <- vag_and_fec %>% 
  # filter(biome_id==46)

cutoff_levels <- c(0:20)

# Check # samples that meet cut off
cutoff_summary <- lapply(cutoff_levels, function(cutoff) {
  vag_and_fec_counts %>%
    summarize(
      cutoff = cutoff,
      fecal_meets_cutoff = sum(total_fecal >= cutoff),
      vaginal_meets_cutoff = sum(total_vaginal >= cutoff),
      has_both_meets_cutoff=sum(has_both >= cutoff),
      total_meets_cutoff = sum(total_samples >= cutoff)
    )
  }) %>%
  bind_rows()

cutoff_long <- cutoff_summary %>%
  pivot_longer(
    cols = c(fecal_meets_cutoff, vaginal_meets_cutoff, has_both_meets_cutoff, total_meets_cutoff),
    names_to = "sampleType",
    values_to = "count"
  ) %>% 
  mutate(sampleType=as.factor(sampleType))

cutoff_long$sampleType <- factor(
  cutoff_long$sampleType,
  levels = c("fecal_meets_cutoff", "vaginal_meets_cutoff", "has_both_meets_cutoff", "total_meets_cutoff")
)

# Plot
###LEGEND???
ggplot(cutoff_long, aes(x = factor(cutoff), y = count, fill = sampleType)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label=count), position = position_dodge(width = 0.9),
            vjust = -0.2) +
  labs(
    title = " ",
    x = " ",
    y = " ",
    fill = "Swab Type"
  ) +
  scale_fill_manual(
    values = c("fecal_meets_cutoff" = "blue", 
               "vaginal_meets_cutoff" = "red", 
               "has_both_meets_cutoff" = "purple",
               "total_meets_cutoff" = "lightgreen"),
    labels = c("Fecal", "Vaginal", "Both", "Total")
  ) +
  scale_y_continuous(
    breaks = seq(0, 70, by = 5),
    limits=c(0,75)
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

