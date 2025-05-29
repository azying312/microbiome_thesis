library(tidyverse)
library(phyloseq)

gut.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/relabeled_data/gut.lifestyle.merged.csv")
vaginal.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/vaginal.lifestyle.csv")

vaginal.data <- vaginal.data %>% 
  select(-X) %>% 
  rename(vaginal_shannon = shannon,
         # vaginal_max_taxa = max_taxa,
         vaginal_OTU = OTU,
         vaginal_sampleID = SampleID
         )

gut.data <- gut.data %>% 
  select(-c(X)) %>% 
  rename(gut_shannon=shannon,
         # gut_max_taxa = max_taxa,
         gut_OTU = OTU,
         gut_sampleID = SampleID)

cross.df <- vaginal.data %>% 
  left_join(gut.data) %>% 
  filter(!is.na(gut_shannon))

names(cross.df)
dim(cross.df)

## bacterial data - prune to one sample a day (first)

vag.bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/relabeled_data/vaginal_cleaned_max_taxa.rds")
gut.bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/relabeled_data/fecal_cleaned_max_taxa.rds")

vag_ids_to_keep <- cross.df$vaginal_sampleID
gut_ids_to_keep <- cross.df$gut_sampleID

# Subset the phyloseq dfs
vag.bacterial.subset <- prune_samples(vag_ids_to_keep, vag.bacterial.data)
gut.bacterial.subset <- prune_samples(gut_ids_to_keep, gut.bacterial.data)

# metadata df
vag.meta <- sample_data(vag.bacterial.subset)
gut.meta <- sample_data(gut.bacterial.subset)

# merge bacterial data of sites
merge.site.bacterial <- merge_phyloseq(vag.bacterial.subset, gut.bacterial.subset)

# bray curtis distance of pairs
bray_dist <- phyloseq::distance(merge.site.bacterial, method = "bray")

bray_df <- as.matrix(bray_dist) %>%
  as.data.frame() %>%
  rownames_to_column("Sample1") %>%
  pivot_longer(-Sample1, names_to = "Sample2", values_to = "bray") %>%
  filter(Sample1 < Sample2) 

# join df
meta <- as(sample_data(merge.site.bacterial), "data.frame") %>% 
  select(biome_id, SampleID, logDate, sampleType)
bray_meta <- bray_df %>%
  left_join(meta, by = c("Sample1" = "SampleID")) %>%
  rename(biome_id_1 = biome_id, site1 = sampleType, logDate_1 = logDate) %>%
  left_join(meta, by = c("Sample2" = "SampleID")) %>%
  rename(biome_id_2 = biome_id, site2 = sampleType, logDate_2 = logDate)

# Filter to within-person, cross-site comparisons
bray_cross_site <- bray_meta %>%
  filter(biome_id_1 == biome_id_2,
         site1 != site2,
         logDate_1 == logDate_2)
# Join meta data
bray_cross_site.full <- bray_cross_site %>% 
  left_join(cross.df, by = c("logDate_1" = "logDate", "biome_id_1" = "biome_id"))
names(bray_cross_site.full)

# gut vs vaginal Shannon
ggplot(bray_cross_site.full, aes(x = vaginal_shannon, y = gut_shannon, col=as.factor(biome_id_1))) +
  geom_point(alpha = 0.6) +
  # geom_smooth(method = "loess", se = FALSE) +
  geom_smooth(method = "loess", se = FALSE, color="blue") +
  labs(
    x = "Vaginal Shannon Diversity",
    y = "Gut Shannon Diversity",
    title = " "
  ) +
  theme_minimal() +
  ylim(0,5)+
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0),
        text=element_text(size=16),
        legend.position="none")

cor.test(bray_cross_site.full$vaginal_shannon, 
           bray_cross_site.full$gut_shannon, use = "complete.obs")

## Bray

# Beta Diversity Between Gut and Vaginal Microbiomes pairs
ggplot(bray_cross_site.full, aes(x = as.factor(biome_id_1), y = bray)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "steelblue") +
  labs(
    x = "Participant ID",
    y = "Bray-Curtis Dissimilarity (Gut vs Vaginal)",
    title = " "
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5),
        text = element_text(size = 14))
# Cross-site Bray-Curtis by CST
ggplot(bray_cross_site.full, aes(x = CST, y = bray, fill = CST)) +
  geom_boxplot() +
  geom_jitter(color="steelblue", width = 0.2, alpha = 0.5) +
  labs(
    x = "CST",
    y = "Bray-Curtis Dissimilarity",
    title = ""
  ) +
  theme_minimal() +
  scale_fill_viridis_d() +
  theme(text = element_text(size = 14)) +
  guides(fill = guide_legend(title = "CST"))

## Edit df
bray_cross_site.full.filtered <- bray_cross_site.full %>% 
  select(-c(Sample1, Sample2, logDate_2, biome_id_2, site1, site2, vaginal_sampleID, gut_sampleID)) %>% 
  rename(logDate = logDate_1,
         biome_id = biome_id_1)
colSums(is.na(bray_cross_site.full.filtered))

bray_cross_site.full.filtered <- bray_cross_site.full.filtered %>% 
  mutate(stress_severity = case_when(
    stress_score <= 14 ~ "Normal",
    stress_score >= 15 & stress_score <= 18 ~ "Mild",
    stress_score >= 19 & stress_score <= 25 ~ "Moderate",
    stress_score >= 26 & stress_score <= 33 ~ "Severe",
    stress_score >= 34 ~ "Extremely Severe",
    TRUE ~ NA_character_
  )) %>% 
  select(-stress_score)
colSums(is.na(bray_cross_site.full.filtered))

# filter all NA cols and all NA rows
bray_cross_site.full.filtered <- bray_cross_site.full.filtered %>%
  select(where(~ !all(is.na(.)))) %>%
  filter(if_any(everything(), ~ !is.na(.))) 
dim(bray_cross_site.full.filtered)
colSums(is.na(bray_cross_site.full.filtered))

# lmer
lmer.obj <- lmer(bray~.+(1|`biome_id`), data=bray_cross_site.full.filtered[,-3])
r2(lmer.obj)
summary(lmer.obj)

lmer.obj <- lmer(bray~gut_shannon+vaginal_shannon+as.factor(menses_day)+(1|`biome_id`), data=bray_cross_site.full.filtered)
r2(lmer.obj)
summary(lmer.obj)

ggplot(bray_cross_site.full.filtered, 
       aes(x = vaginal_shannon, y = bray, color = as.factor(menses_day))) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(color = "Menses Day",
       x = "Vaginal Shannon Diversity", 
       y = "Bray-Curtis Dissimilarity") +
  theme_minimal()

t.test(vaginal_shannon ~ menses_day, 
       data = bray_cross_site.full.filtered)

####

bray_combined <- bind_rows(
  bray_cross_site %>% mutate(comparison = "within-person"),
  bray_meta %>%
    filter(biome_id_1 != biome_id_2, site1 != site2, logDate_1 == logDate_2) %>%
    mutate(comparison = "between-person")
)

# Plot with a side-by-side boxplot
ggplot(bray_combined, aes(x = comparison, y = bray, fill = comparison)) +
  geom_boxplot(alpha = 0.7) +
  # geom_jitter(color="orchid", alpha=0.3) +
  labs(x = "", y = "Bray-Curtis Dissimilarity", title = "Cross-site Dissimilarity (Gut vs Vaginal)") +
  scale_fill_manual(values = c("within-person" = "#1f77b4", "between-person" = "#ff7f0e")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
        text=element_text(size=16))
ggplot(bray_combined, aes(x = comparison, y = bray, fill = comparison)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, alpha = 0.3) + 
  labs(x = "", y = "Bray-Curtis Dissimilarity", title = "Cross-site Dissimilarity (Gut vs Vaginal)") +
  scale_fill_manual(values = c("within-person" = "#1f77b4", "between-person" = "#ff7f0e")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        text=element_text(size=16))

## meta data

bray_cross_site.full <- bray_cross_site %>% 
  left_join(cross.df, by = c("logDate_1" = "logDate", "biome_id_1" = "biome_id"))

ggplot(bray_cross_site.full, aes(x = steps, y = bray, color = site1)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE) +
  labs(title = "Bray-Curtis Dissimilarity vs Physical Activity (Steps)",
       x = "Steps", y = "Bray-Curtis Dissimilarity") +
  theme_minimal()

lifestyle_factors <- c("stress_score", "caloriesall_avg", "cholesterol_prop", "satFat_prop", 
                       "sodium_prop", "carb_prop", "dietFib_prop", "sugar_prop", 
                       "protein_prop", "fat_prop", "fat_cal_prop", "addedSugarall_prop", 
                       "calories_burned", "steps", "distance", "minutes_sedentary", 
                       "minutes_lightly_active", "minutes_fairly_active", 
                       "minues_very_active", "activity_calories", "sport", 
                       "probiotic", "study_menstruate", "sexuallyActive")



results <- lapply(lifestyle_factors, function(factor) {
  
  # Filter rows where there are no NAs in the relevant lifestyle factor
  bray_cross_site.clean_filtered <- bray_cross_site.clean %>%
    filter(!is.na(.data[[factor]]))  # Dynamically refer to the factor column
  
  # Check if there are enough rows left after filtering
  if(nrow(bray_cross_site.clean_filtered) < 2) {
    return(c(factor, NA))  # If not enough data, return NA for p-value
  }
  
  # Compute Bray-Curtis distance on the filtered data
  bray_dist2 <- vegdist(bray_cross_site.clean_filtered$bray, method = "bray")
  
  # Build the formula dynamically for adonis2
  formula <- as.formula(paste("bray_dist2 ~", factor))
  
  # Run the adonis2 test on the filtered data
  perm_result <- adonis2(formula, data = bray_cross_site.clean_filtered, permutations = 999)
  
  # Return the factor and p-value
  return(c(factor, perm_result$`Pr(>F)`)) 
})

# Convert results into a data frame for easy viewing
results_df <- do.call(rbind, results)
colnames(results_df) <- c("Factor", "p_value")

# View the results
print(results_df)

names(bray_cross_site.clean)
perm_result <- adonis2(bray_dist2 ~ wa, data = bray_cross_site.clean, permutations = 999)
perm_result

lm_result <- lm(bray ~ stress_score, data = bray_cross_site.full)

