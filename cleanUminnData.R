library(dplyr)

uminn_data <- read.csv("/Volumes/T7/microbiome_data/Swabs with blood - Sheet1.csv", header=TRUE)

uminn_data <- uminn_data %>% 
  filter(Sample.ID != "BLANK") %>% 
  # errors in sample processing from UMinn
  filter(!str_detect(Special.Notes, "error")) %>% 
  filter(!str_detect(Special.Notes, "No swab in tube"))

# QR code
uminn_data$qr <- sub("_.*", "", uminn_data$Sample.ID)

write.csv(uminn_data,
          file = "/Volumes/T7/microbiome_data/cleaned_data/cleaned_uminn_data.csv",
          row.names = FALSE)
