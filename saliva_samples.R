library(tidyverse)
library(viridis)

samples.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_samples.csv")

saliva_samples <- samples.data %>% 
  filter(sampleType != "vaginal" & sampleType != "fecal")

dim(saliva_samples)
participant_ids <- unique(saliva_samples$biome_id)
length(participant_ids)

# Frequency plot
saliva_samples_freq <- saliva_samples %>% 
  group_by(biome_id) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n)) %>% 
  mutate(biome_id=factor(biome_id, levels=biome_id))

ggplot(saliva_samples_freq, aes(x=biome_id, y=n)) +
  geom_bar(stat="identity", fill="blue", color="black") +
  theme_minimal() +
  ylim(0,10) +
  labs(x="Participant", y="Frequency")

all_days <- seq.Date(as.Date("2022-10-13"), as.Date("2022-12-16"), by = "day")
all_days <- data.frame(logDate=as.character(all_days))
all_days_expanded <- expand_grid(biome_id = participant_ids, logDate = all_days$logDate) %>% 
  as.data.frame()

saliva22 <- all_days_expanded %>%
  left_join(saliva_samples, by = c("logDate", "biome_id")) %>% 
  mutate(biome_id=factor(biome_id))
saliva22_joined <- saliva22 %>% 
  left_join(saliva_samples_freq, by="biome_id")

ggplot(saliva22_joined, aes(x=as.factor(logDate), y=reorder(biome_id, n), fill=as.factor(sampleType))) +
  geom_tile(color="black") +
  scale_fill_manual(
    values = c("saliva" = "blue"),
    na.value = "white",
    name = "Value"
  ) +
  labs(x = " ", y = " ", title = " ") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid = element_blank()
  )

