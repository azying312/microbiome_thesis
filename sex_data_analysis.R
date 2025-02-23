library(tidyverse)

sex_act_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 3-Sexual Activity.csv")

#### Correlate with shannon data
vaginal.microbial.menses.24 <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifetyle/vaginal.microbial.menses.24.csv")
vaginal.microbial.menses.24 <- vaginal.microbial.menses.24[,-1]

vaginal.microbial.menses.24.subset <- vaginal.microbial.menses.24 %>% 
  select(biome_id, logDate, SampleID, CST, shannon, sampleType, max_taxa, OTU, menses_day)
vaginal.microbial.menses.24.subset.summary <- vaginal.microbial.menses.24.subset %>% 
  group_by(biome_id) %>% 
  summarise(avg_shannon = sum(shannon)/n())

## analyze by sex active v. not

vaginal.microbial.menses.24.subset.summary <- vaginal.microbial.menses.24.subset.summary %>% 
  mutate(sex_active = ifelse(biome_id %in% unique(sex_act_data$biome_id), "active", "not active"))

# boxplot of sex active v. not
ggplot(vaginal.microbial.menses.24.subset.summary, aes(x = sex_active, y = avg_shannon)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color="orchid") +
  labs(
    x = " ", 
    y = "Average Shannon Diversity", 
    title = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none") 

# testing
t.test(avg_shannon ~ sex_active, data=vaginal.microbial.menses.24.subset.summary)
wilcox.test(avg_shannon ~ sex_active, data = vaginal.microbial.menses.24.subset.summary)

# analyze by those who are by contraception, num_partners, type of intercourse

vaginal.sexactive <- sex_act_data %>% 
  left_join(vaginal.microbial.menses.24.subset, by=c("biome_id", "logDate")) %>% 
  filter(!is.na(SampleID))

head(vaginal.sexactive)


##########################################################################################

# Correlate with CST
vaginal.microbial.menses.24_collapsed <- vaginal.microbial.menses.24 %>% 
  # Collapse by person (assign to most frequent CST)
  group_by(biome_id) %>% 
  count(CST, name="frequency") %>% # total 148 (participants have more than 1 CST) 
  slice_max(frequency, n=1) %>% # collapse to their most frequent CST
  rename(CST_max=CST) %>% 
  ungroup() %>%
  left_join(sex_act_data, by = "biome_id") %>%
  select(biome_id, CST_max, logDate, type_of_intercourse, gender_of_partner, num_partners, contraceptions) %>% 
  distinct(biome_id, .keep_all = TRUE) %>% 
  mutate(sex_active = ifelse(biome_id %in% unique(sex_act_data$biome_id), "active", "not active"))

# CST and sexually active
sexactive_CST <- table(vaginal.microbial.menses.24_collapsed$sex_active, vaginal.microbial.menses.24_collapsed$CST_max)
chisq.test(sexactive_CST)

vaginal.microbial.menses.24_collapsed <- vaginal.microbial.menses.24_collapsed %>% 
  mutate(sex_active_bin = ifelse(sex_active == "active", 0, 1))
glm.sexactive_CST <- glm(sex_active_bin ~ CST_max, data = vaginal.microbial.menses.24_collapsed, family = binomial)
summary(glm.sexactive_CST)


