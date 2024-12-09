########################
#
# Clean Menstruation Data
# 8 October 2024
#
# Dataset Outputs: cleaned_menstruation_data.csv & cleaned_vaginal_samples_data.csv
#
#########################

library(tidyverse)
library(reshape2)
library(stringr)

# Load data
menses_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Report 1-Menstruation - Report 1-Menstruation.csv")
id_mapping <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Original Study Mapping - Sheet3.csv", header = TRUE)
survey_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/alice_cleaned_survey_data.csv", header=TRUE)
# uminn_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Swabs with blood - Sheet1.csv", header=TRUE)
# uminn_data_qc <- read.csv("/Users/alicezhang/Desktop/microbiome_data/Swabs with blood - QC.csv", header=TRUE)
uminn_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_uminn_data.csv", header=TRUE)
samples_data <- read.csv("/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_samples.csv", header=TRUE)

### MAPPING IDS

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
# Check missing ids
missing_list <- menses_data %>%
  filter(is.na(as.numeric(biome_id)))
print(unique(missing_list$biome_id))
# convert to numeric
menses_data$biome_id <- as.numeric(menses_data$biome_id)
head(menses_data)

### DATA CLEANING

menses_data_cleaned <- menses_data %>% 
  filter(menstruation!="") %>% 
  mutate(menstruation_numeric = case_when(
    menstruation == "none" ~ -1,
    menstruation == "light" ~ 1,
    menstruation == "medium" ~ 2,
    menstruation == "heavy" ~ 3,
    menstruation == "" ~ NA_real_, # self report, but no input
    TRUE ~ NA_real_
  )) %>%
  mutate(menstruation_text=menstruation) %>% 
  mutate(menstruation=ifelse(menstruation_numeric>0, 1, menstruation_numeric)) %>% 
  # select(biome_id, logDate, menstruation) %>% 
  mutate(biome_id=as.numeric(biome_id)) %>% 
  mutate(inSelfReport=TRUE)
head(menses_data_cleaned)

### UMINN DATA INTERPOLATION

# get samples
uminn_data_subset <- uminn_data %>% 
  select(Sample.ID, Well, Special.Notes)

# Use uminn_data_subset - set IDs
# uminn_data_subset$qr <- sub("_.*", "", uminn_data_subset$Sample.ID)
uminn_data_subset$inUminn <- TRUE

# match uminn samples with blood to vaginal samples
vaginal_samples <- samples_data %>% 
  filter(sampleType=="vaginal")

# keep all the vaginal samples, mark uminn menstruation by if there is blood
uminn_data_subset <- uminn_data_subset %>% 
  left_join(vaginal_samples, by="qr") %>% 
  filter(!is.na(sampleType))

# clean UMinn samples
uminn_data_subset <- uminn_data_subset %>% 
  filter(logDate!="0000-00-00") %>% # no corresponding logDate
  mutate(uMinn_menstruation=ifelse(str_detect(Special.Notes, "Blood")==TRUE, 1, 0)) %>%  # set to menstruation true
  select(biome_id, logDate, uMinn_menstruation, inUminn, Special.Notes) %>% 
  mutate(biome_id=as.numeric(biome_id)) #%>% 
  # errors in sample processing from UMinn
  # filter(!str_detect(Special.Notes, "error")) %>% 
  # filter(!str_detect(Special.Notes, "No swab in tube"))

# error_uminn_data_subset <- uminn_data_subset %>% 
  # filter(str_detect(Special.Notes, "Processor error") | str_detect(Special.Notes, "Technical error") | str_detect(Special.Notes, "No swab in tube"))

# identical rows in UMinn
identical_uminn_data_subset <- uminn_data_subset %>% 
  group_by(across(everything())) %>%
  filter(n() > 1) %>%
  ungroup()

# filter for the identical rows
uminn_data_subset_unique <- uminn_data_subset %>%
  distinct()

# join menses data (self-report) with UMinn data (blood in sample)
full_menstruation_data <- uminn_data_subset_unique %>% 
  full_join(menses_data_cleaned) %>% 
  mutate(inUminn=ifelse(is.na(inUminn), FALSE, inUminn),
         inSelfReport=ifelse(is.na(inSelfReport), FALSE, inSelfReport)) %>% 
  select(biome_id, logDate, uMinn_menstruation, inUminn, Special.Notes, 
         menstruation, menstruation_numeric, menstruation_text, inSelfReport)

# filter for identical rows by biome_id and logDate
full_menstruation_data <- full_menstruation_data %>%
  group_by(biome_id, logDate) %>%
  filter(!(n() > 1 & (is.na(Special.Notes) | Special.Notes == ""))) %>%
  filter(!(n() > 1 & !str_detect(Special.Notes, "^Blood on swab$"))) %>%
  ungroup() %>% 
  distinct()

# check incongruity: if any days with self report no menses & sample with blood 
incong_menses <- full_menstruation_data %>% 
  filter(uMinn_menstruation==TRUE & menstruation==-1)
incong_menses

# na_menses <- full_menstruation_data %>% 
  # filter(full_menstruation_data$inSelfReport==TRUE & is.na(full_menstruation_data$menstruation))

full_menstruation_data <- full_menstruation_data %>%
  mutate(menstruation_status=case_when(
    # only UMinn sample with blood
    (uMinn_menstruation==1 & is.na(menstruation)) ~ 1,
    # only self report menstruation (no sample)
    (inUminn==FALSE & menstruation==1) ~ 2,
    # UMinn sample with blood & self report menses
    (uMinn_menstruation==1 & menstruation==1) ~ 3,
    # self report no menses (no sample)
    (inUminn==FALSE & menstruation==-1) ~ 4,
    # UMinn sample with no blood & self report no menses
    (uMinn_menstruation==0 & menstruation==-1) ~ 5,
    # sample with no blood (no self report)
    (uMinn_menstruation==0 & inSelfReport==FALSE) ~ 6,#is.na(menstruation)) ~ 6,
    # self report menstruation & no blood on sample
    (uMinn_menstruation==0 & menstruation==1) ~ 7,
    # no sample & no self-report
    (inUminn==FALSE & is.na(menstruation)) ~ 8,
    TRUE ~ NA_real_
  ))

# change NAs to 0s for menstruation (self-report) variable
full_menstruation_data <- full_menstruation_data %>% 
  mutate(menstruation=ifelse((inSelfReport==TRUE & is.na(menstruation)), 0, menstruation))

full_menstruation_data <- full_menstruation_data %>% 
  select(biome_id, logDate, menstruation_status, inSelfReport, menstruation, inUminn, uMinn_menstruation, Special.Notes)

table(full_menstruation_data$menstruation_status)
sum(is.na(full_menstruation_data$menstruation_status))

### Save final data output
write.csv(full_menstruation_data,
          file = "/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_menstruation_data.csv",
          row.names = FALSE)

# Vaginal samples from uminn (menses and non menses)
uminn_data_vaginal <- uminn_data %>% 
  mutate(qr=sub("_.*", "", Sample.ID)) %>% 
  mutate(inUminn=TRUE) %>% 
  select(Sample.ID, qr, inUminn, Special.Notes)
vaginal_samples <- samples_data %>% 
  filter(sampleType=="vaginal")
uminn_data_vaginal <- uminn_data_vaginal %>% 
  left_join(vaginal_samples, by="qr") %>% 
  filter(sampleType=="vaginal") %>% 
  filter(logDate!="0000-00-00") %>% 
  select(qr, inUminn, biome_id, logDate, sampleType, timestamp, Special.Notes) %>% 
  mutate(
    biome_id = as.numeric(biome_id),
    uMinnMenses = if_else(grepl("Blood", Special.Notes, ignore.case = TRUE), TRUE, FALSE)
  ) %>% 
  filter(!is.na(biome_id)) %>% # filter NA ids
  filter(!str_detect(Special.Notes, "error")) # filter out errors

write.csv(uminn_data_vaginal,
          file = "/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_vaginal_samples_data.csv",
          row.names = FALSE)


#### MISSINGNESS GRAPHS

############# MENSTRUATION DATA FROM PEOPLE SUBMISSION

# Get all days (including no reports)
all_days <- seq.Date(as.Date(min(menses_data_cleaned$logDate)), as.Date(max(menses_data_cleaned$logDate)), by = "day")
all_days <- as.character(all_days)

### Heatmap of missingness
menses_data_cleaned_subset <- menses_data_cleaned %>% 
  select(biome_id, logDate, menstruation_numeric)
# unique biome_id and days
complete_grid <- expand.grid(biome_id = unique(menses_data_cleaned_subset$biome_id),
                             logDate = all_days)
dim(complete_grid)
heatmap_data <- complete_grid %>%
  left_join(menses_data_cleaned_subset, by = c("biome_id", "logDate")) 

ggplot(heatmap_data, aes(x = logDate, y = factor(biome_id), fill = menstruation_numeric)) +
  geom_tile(color = "gray25") +
  scale_fill_gradientn(colors = c("lightblue", "blue", "black"), 
                       values = scales::rescale(c(0, 1, 2, 3)),  # Adjust values for your scale
                       na.value = "white") +  # Color for NA values (if any)
  labs(title = "Menstruation Heatmap", 
       x = "Date", 
       y = "Biome ID", 
       fill = "Menstruation Scale") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) 

## BW Missingness heatmap

# Reorder study id to be in order of missingness
survey_counts <- sort(table(menses_data_cleaned$biome_id), decreasing = TRUE)
ordered_study_ids <- names(survey_counts) # in descending order
num_ID <- length(unique(menses_data_cleaned$biome_id)) # one row for each ID
num_days <- length(unique(menses_data_cleaned$logDate))

grid <- matrix("", nrow = num_ID, ncol = length(all_days))
dim(grid)
rownames(grid) <- unique(menses_data_cleaned$biome_id)
colnames(grid) <- as.character(all_days)

for (i in 1:nrow(menses_data_cleaned)) {
  day <- menses_data_cleaned[i, "logDate"]
  grid[as.character(menses_data_cleaned$biome_id[i]), day] <- "X"
}

# Convert the matrix to a data frame
d_grid<- reshape2::melt(grid)
colnames(d_grid) <- c("Study_ID", "Day", "Value")
d_grid$Study_ID <- factor(d_grid$Study_ID, levels=rev(ordered_study_ids))

## All days plot
grid_plot <- ggplot(d_grid, aes(x = Day, y = reorder(Study_ID, Study_ID, decreasing=FALSE), fill = Value)) +
  geom_tile(color = "gray25") +
  labs(title = "Participant Entry vs. Day", y = "Study ID",
       x = "Date") +
  scale_fill_manual(values = c("white", "gray50"),
                    labels=c("No Response", "Responded"),
                    guide = guide_legend(reverse = TRUE)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
        text = element_text(family = "serif"),
        plot.title = element_text(hjust = 0.5))

grid_plot

### Just submission days
grid <- matrix("", nrow = num_ID, ncol = num_days)
dim(grid)
rownames(grid) <- unique(menses_data_cleaned$biome_id)
colnames(grid) <- unique(menses_data_cleaned$logDate)

for (i in 1:nrow(menses_data_cleaned)) {
  day <- menses_data_cleaned[i, "logDate"]
  grid[as.character(menses_data_cleaned$biome_id[i]), day] <- "X"
}

## Submission days plot
# Convert the matrix to a data frame
d_grid<- reshape2::melt(grid)
colnames(d_grid) <- c("Study_ID", "Day", "Value")
d_grid$Study_ID <- factor(d_grid$Study_ID, levels=rev(ordered_study_ids))

grid_plot <- ggplot(d_grid, aes(x = Day, y = reorder(Study_ID, Study_ID, decreasing=FALSE), fill = Value)) +
  geom_tile(color = "gray25") +
  labs(title = "Participant Entry vs. Day", y = "Study ID",
       x = "Date") +
  scale_fill_manual(values = c("white", "gray50"),
                    labels=c("No Response", "Responded"),
                    guide = guide_legend(reverse = TRUE)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
        text = element_text(family = "serif"),
        plot.title = element_text(hjust = 0.5))

grid_plot

###########

# ### Save final data output
# write.csv(menses_data_cleaned,
#           file = "/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_Report 1-Menstruation.csv",
#           row.names = FALSE)

### UMinn Plot
# all_days <- seq.Date(as.Date(min(uminn_data_subset_no0$logDate, na.rm=TRUE)), 
#                      as.Date(max(uminn_data_subset_no0$logDate, na.rm=TRUE)), by = "day")
# all_days <- as.character(all_days)
# uminn_data_subset <- uminn_data_subset_no0 %>% 
#   select(biome_id, logDate, menstruating)
# complete_grid <- expand.grid(biome_id = unique(uminn_data_subset$biome_id),
#                              logDate = all_days)
# heatmap_data <- complete_grid %>%
#   left_join(uminn_data_subset, by = c("biome_id", "logDate")) 

# ## All days plot
# heatmap_data$Value <- ifelse(is.na(heatmap_data$menstruating), "Not menstruating", "Menstruating")
# heatmap_data$Value <- factor(heatmap_data$Value, levels = c("Not menstruating", "Menstruating"))
# 
# ggplot(heatmap_data, aes(x = logDate, y = factor(biome_id), fill = Value)) +
#   geom_tile(color = "gray25") +
#   scale_fill_manual(values = c("white", "darkblue"), na.value = "white", 
#                     labels = c("Not menstruating", "Menstruating")) +
#   labs(title = "Menstruation Heatmap", 
#        x = "Date", 
#        y = "Biome ID", 
#        fill = "Menstruation Scale") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) 
# 
# ### UMinn & Self-submission Menses Plot
# all_days <- seq.Date(as.Date(min(full_menstruation_data$logDate, na.rm=TRUE)), 
#                      as.Date(max(full_menstruation_data$logDate, na.rm=TRUE)), by = "day")
# all_days <- as.character(all_days)
# full_menstruation_data_subset <- full_menstruation_data %>% 
#   select(biome_id, logDate, menstruating)
# complete_grid <- expand.grid(biome_id = unique(full_menstruation_data$biome_id),
#                              logDate = all_days)
# heatmap_data <- complete_grid %>%
#   left_join(full_menstruation_data_subset, by = c("biome_id", "logDate")) 
# 
# ## All days plot
# heatmap_data$Value <- ifelse(is.na(heatmap_data$menstruating), "Not menstruating", "Menstruating")
# heatmap_data$Value <- factor(heatmap_data$Value, levels = c("Not menstruating", "Menstruating"))
# 
# ggplot(heatmap_data, aes(x = logDate, y = factor(biome_id), fill = Value)) +
#   geom_tile(color = "gray25") +
#   scale_fill_manual(values = c("white", "darkblue"), na.value = "white", 
#                     labels = c("Not menstruating", "Menstruating")) +
#   labs(title = "Menstruation Heatmap", 
#        x = "Date", 
#        y = "Biome ID", 
#        fill = "Menstruation Scale") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) 
