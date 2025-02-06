library(tidyverse)
library(viridis)

uminn_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_uminn_data.csv", header=TRUE)
samples_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_samples.csv", header=TRUE)

source("~/Microbiome Thesis/functions.R")

## Sample data Cleaning
fecal_samples_data <- samples_data %>%
  filter(sampleType=="fecal")

### UMinn Data Cleaning
uminn_data_subset <- uminn_data %>% 
  select(Sample.ID, Well, Special.Notes)
# Use uminn_data_subset - set IDs
uminn_data_subset$qr <- sub("_.*", "", uminn_data_subset$Sample.ID)
uminn_data_subset$inUminn <- TRUE
# keep all the fecal samples
uminn_data_subset <- uminn_data_subset %>%
  left_join(fecal_samples_data, by="qr") %>%
  filter(!is.na(sampleType))
# clean UMinn samples
uminn_data_subset <- uminn_data_subset %>% 
  # no corresponding logDate
  filter(logDate!="0000-00-00") %>%
  select(biome_id, logDate, Special.Notes, sampleType, qr) %>% 
  mutate(biome_id=as.numeric(biome_id))

## Filtering
# identical rows in UMinn
identical_uminn_data_subset <- uminn_data_subset %>% 
  group_by(across(everything())) %>%
  filter(n() > 1) %>%
  ungroup()
# filter for the identical rows
uminn_data_subset_unique <- uminn_data_subset %>%
  distinct()

# Set submission
fecal_samples_data$submitted <- TRUE

# Save fecal data obj
write.csv(fecal_samples_data,
          file = "/Volumes/T7/microbiome_data/cleaned_data/cleaned_fecal_samples_data.csv",
          row.names = FALSE)

####################################################################
# Figures

# Get all days (including no reports)
all_days <- seq.Date(as.Date(min(fecal_samples_data$logDate)), as.Date(max(fecal_samples_data$logDate)), by = "day")
all_days <- as.character(all_days)

### Heatmap of completeness
fecal_samples_data_subset <- fecal_samples_data %>% 
  select(biome_id, logDate, submitted)

# unique biome_id and days
complete_grid <- expand.grid(biome_id = unique(fecal_samples_data_subset$biome_id),
                             logDate = all_days)
heatmap_data <- complete_grid %>%
  left_join(fecal_samples_data_subset, by = c("biome_id", "logDate")) 

# Reorder
ggplot(heatmap_data, aes(x = logDate, y = factor(biome_id), fill=factor(submitted))) +
  geom_tile(color = "gray25") +
  scale_fill_viridis_d(na.value = "white") + 
  labs(title = "Fecal Sample Submission", 
       x = "Date", 
       y = "Biome ID",
       fill="Submission") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) 

