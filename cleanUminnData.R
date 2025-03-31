library(tidyverse)

# uminn_data <- read.csv("/Volumes/T7/microbiome_data/Swabs with blood - Sheet1.csv", header=TRUE)
uminn_data <- read.csv("/Volumes/T7/microbiome_data/original_data/UMinn Samples - Sheet1.csv", header=TRUE)

dim(uminn_data)

error_data <- uminn_data %>% 
  filter(Sample.ID != "BLANK") %>%
  # errors in sample processing from UMinn
  filter(str_detect(Special.Notes, "error"))
dim(error_data)

missing_data <- uminn_data %>% 
  filter(Sample.ID != "BLANK") %>%
  filter(str_detect(Special.Notes, "No swab in tube"))
dim(missing_data)

uminn_data <- uminn_data %>% 
  filter(Sample.ID != "BLANK") %>%
  # errors in sample processing from UMinn
  filter(!str_detect(Special.Notes, "error")) %>% 
  filter(!str_detect(Special.Notes, "No swab in tube"))

# QR code
uminn_data$qr <- sub("_.*", "", uminn_data$Sample.ID)

dim(uminn_data)

write.csv(uminn_data,
          file = "/Volumes/T7/microbiome_data/cleaned_data/cleaned_uminn_data.csv",
          row.names = FALSE)

# View(uminn_data %>%  filter(Sample.ID == "BLANK"))
table_output <- table(uminn_data$Sample.ID)
subset(table_output, table_output > 1)
