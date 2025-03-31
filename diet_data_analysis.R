########################
#
# Exploration of the Merged Diet Data
# Last updated: 03/18/2025
#
#########################

library(tidyverse)
library(viridis)

source("~/Microbiome Thesis/functions.R")

## Packages
packs <- c("tidyverse")
lapply(packs, require, character.only = TRUE)

## Load Data
merged_diet_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/fully_merged_diet_data.csv")

##########################################################################################
# Visualize diet completeness map

merged_diet_data <- merged_diet_data %>% 
  mutate(Date=as.Date(Date)) %>% 
  filter(Date < as.Date("2022-12-17")) 

# Get average cal intake for participants
merged_diet_data_collapsed <- merged_diet_data %>% 
  group_by(biome_id, Date) %>% 
  summarise(total_cals = sum(caloriesall),
            entries_per_day = n()) %>% 
  ungroup() %>% 
  group_by(biome_id) %>% 
  summarise(avg_study_cals = sum(total_cals / n()),
            submission_days = n(),
            entries_through_study = sum(entries_per_day),
            avg_entries = sum(entries_per_day)/n())

head(merged_diet_data_collapsed, 7)

# summary stats about diet data
summary(merged_diet_data_collapsed$submission_days)
summary(merged_diet_data_collapsed$avg_study_cals)
summary(merged_diet_data_collapsed$avg_entries)

# Figure: average caloric intake across study for participants
ggplot(merged_diet_data_collapsed, aes(x =  avg_study_cals)) + 
  geom_histogram(bins = 50, fill="orchid", color = "black") +
  # step size adj
  scale_x_continuous(
    breaks = seq(0, max(merged_diet_data_collapsed$avg_study_cals), by = 50)  
  ) +
  labs(x = "Daily Average Calories", y = "Frequency", title = "") +
  theme_minimal()

summary(merged_diet_data_collapsed$avg_cals)

all_days <- seq.Date(as.Date(min(merged_diet_data$Date, na.rm=TRUE)), 
                     as.Date(max(merged_diet_data$Date, na.rm=TRUE)), by = "day")

# Expand to include all days
expanded_data <- expand.grid(biome_id = unique(merged_diet_data$biome_id), Date = all_days)
expanded_data <- left_join(expanded_data, merged_diet_data, by = c("biome_id", "Date"))

diet_data_subset <- merged_diet_data %>% 
  dplyr::select(biome_id, Date, type, name, caloriesall, proteinall, sugarsall, fatall, food_other)

# get meals per day (only 3 meals max)
meals_per_day <- diet_data_subset %>%
  filter(type %in% c("breakfast", "lunch", "dinner")) %>% 
  group_by(biome_id, Date) %>%
  summarize(
    MealsPerDay = n_distinct(type),
    .groups = "drop"
  ) %>% 
  mutate(study_id=as.numeric(biome_id)) %>% 
  filter(!is.na(biome_id))

complete_grid <- expand.grid(biome_id = unique(meals_per_day$biome_id),
                             Date = all_days)
heatmap_data <- complete_grid %>%
  left_join(meals_per_day, by = c("biome_id", "Date"))

## All days plot
# heatmap_data$Value <- ifelse((is.na(heatmap_data$type)), "Not menstruating", "Menstruating")
# heatmap_data$Value <- factor(heatmap_data$Value, levels = c("Not menstruating", "Menstruating"))

heatmap_data_plot <- heatmap_data %>% 
  group_by(biome_id) %>%
  mutate(meals_count = sum(MealsPerDay, na.rm = TRUE)) %>%
  mutate(MealsPerDay=ifelse(is.na(MealsPerDay), 0, MealsPerDay)) %>% 
  mutate(MealsPerDay=as.factor(MealsPerDay)) %>% 
  arrange(desc(meals_count), biome_id) %>% 
  ungroup()

# all dates showing
ggplot(heatmap_data_plot, aes(x = as.factor(Date), y = reorder(factor(biome_id), meals_count), fill = MealsPerDay)) +
  geom_tile(color = "gray25") +
  scale_fill_manual(values = c("white", "lightblue", "blue", "darkblue"), 
                    na.value = "white")+
  labs(title = " ", 
       x = "Date", 
       y = "Biome ID", 
       fill = "Number of Meals") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

##########################################################################################
sum(is.na(merged_diet_data$serving.size))
dim(merged_diet_data)
sum(is.na(merged_diet_data$servings))

dim(merged_diet_data)
# glimpse(merged_diet_data)
names(merged_diet_data)

head(merged_diet_data)

merged_diet_data %>% 
  group_by(biome_id, Date) %>% 
  summarise(count=n()) %>% 
  dplyr::select(count)

#### Vegetarian
common_meats <- c("beef", "fish", "pork", "chicken", "steak", "turkey", "shrimp", "sushi", "burger king hamburger", 
                    "salmon", "tonkotsu", "mcrib", "bacon", "meatballs", "taco de birria", "meat", "sausage", "burger",
                  "tilapia", "cod", "kielbasa", "clam", "meatloaf", "chop suey", "crab", "ham", "dumplings", "lgblt",
                  "enchilada", "hawaiian pizza", "hot dog", "italian wedding soup", "jimmy johns sandwich", "lamb",
                  "leaky beaker sandwich", "london broil", "mto burrito bar", "mto baked potato bar", "mto deli bar",
                  "mto gyro bar", "pot roast", "meatball calzone", "mcdonald's hamburger", "mulligatawny", "bolognese",
                  "pepperoni", "pho", "chorizo", "ribs", "prosciutto", "scallops", "southwest breakfast wrap",
                  "salami", "taco pizza", "tuna", "wonton", "calamari", "chick fil a sandwich",
                  "duck", "el table sandwich", "wonton soup")

exclude_foods <- c("goldfish", "plant-based chicken nuggets", "goldfish, original", "plant based", "vegan", "vegetarian", 
                   "vegetarian sausage patty", "vegetarian sausage patties", "beyond burger", "no meat", "vegetable sushi",
                   "Chipotle Black Bean Burger")
# Create regex patterns for meat and exclusion lists
meat_pattern <- paste0("\\b(", paste(common_meats, collapse = "|"), ")\\b")
exclude_pattern <- paste0("\\b(", paste(exclude_foods, collapse = "|"), ")\\b")

# Check the common meats
merged_diet_data <- merged_diet_data %>%
  mutate(name_lower = str_to_lower(name),  # Convert names to lowercase
         vegetarian=TRUE,
         # Check for common meats first
         vegetarian = ifelse(str_detect(name_lower, meat_pattern), FALSE, vegetarian)) 

# Look at exclusions
merged_diet_data <- merged_diet_data %>%
  mutate(vegetarian = ifelse(str_detect(name_lower, exclude_pattern), TRUE, vegetarian)) %>% 
  dplyr::select(-name_lower)

unique_vegetarian_pairs <- merged_diet_data %>%
  dplyr::select(name, vegetarian) %>%  
  distinct() %>%                
  arrange(name) 

## Check if participant vegetarian
participant_vegetarian_status <- merged_diet_data %>%
  group_by(biome_id) %>%
  # Check if all foods for a participant are vegetarian
  summarise(is_vegetarian = all(vegetarian))  
participant_vegetarian_status$biome_id <- as.numeric(participant_vegetarian_status$biome_id)
participant_vegetarian_status <- participant_vegetarian_status %>%
  drop_na()

table(participant_vegetarian_status$is_vegetarian) # 9 vegetarian

# Add to diet data
# participant_vegetarian_status$study_id <- as.character(participant_vegetarian_status$study_id)
merged_diet_data <- merged_diet_data %>%
  left_join(participant_vegetarian_status, by = "biome_id")

##########################################################################################
merged_diet_data_vegetarian <- merged_diet_data %>% 
  dplyr::select(biome_id, is_vegetarian) %>% 
  distinct(biome_id, is_vegetarian)

## Collapse diet data to day nutrition
merged_diet_data_daily <- merged_diet_data %>% 
  group_by(biome_id, Date) %>% 
  summarise(
    caloriesall = sum(caloriesall, na.rm = TRUE),
    cholesterolall = sum(cholesterolall, na.rm = TRUE),
    saturatedFatall = sum(saturatedFatall, na.rm = TRUE),
    sodiumall = sum(sodiumall, na.rm = TRUE),
    carbohydratesall = sum(carbohydratesall, na.rm = TRUE),
    dietaryFiberall = sum(dietaryFiberall, na.rm = TRUE),
    sugarsall = sum(sugarsall, na.rm = TRUE),
    proteinall = sum(proteinall, na.rm = TRUE),
    fatall = sum(fatall, na.rm = TRUE),
    caloriesFromFat = sum(caloriesFromFat, na.rm = TRUE),
    saturatedFat = sum(saturatedFat, na.rm = TRUE),
    caloriesFromSatFat = sum(caloriesFromSatFat, na.rm = TRUE),
    transFat = sum(transFat, na.rm = TRUE),
    addedSugarall = sum(addedSugarall, na.rm = TRUE)
  )

## Change data to past 2 days
library(lubridate)
library(slider)

past_two_days_diet_data <- merged_diet_data_daily %>% 
  group_by(biome_id) %>%
  arrange(biome_id, Date) %>%
  mutate(aggregate_date = Date + 2) %>%
  mutate(
    caloriesall      = slider::slide_dbl(caloriesall, sum, .before = 2, .after = -1, .complete = TRUE),
    cholesterolall   = slider::slide_dbl(cholesterolall, sum, .before = 2, .after = -1, .complete = TRUE),
    saturatedFatall  = slider::slide_dbl(saturatedFatall, sum, .before = 2, .after = -1, .complete = TRUE),
    sodiumall        = slider::slide_dbl(sodiumall, sum, .before = 2, .after = -1, .complete = TRUE),
    carbohydratesall = slider::slide_dbl(carbohydratesall, sum, .before = 2, .after = -1, .complete = TRUE),
    dietaryFiberall  = slider::slide_dbl(dietaryFiberall, sum, .before = 2, .after = -1, .complete = TRUE),
    sugarsall        = slider::slide_dbl(sugarsall      , sum, .before = 2, .after = -1, .complete = TRUE),
    proteinall       = slider::slide_dbl(proteinall     , sum, .before = 2, .after = -1, .complete = TRUE),
    caloriesFromFat  = slider::slide_dbl(caloriesFromFat, sum, .before = 2, .after = -1, .complete = TRUE),
    saturatedFat     = slider::slide_dbl(saturatedFat   , sum, .before = 2, .after = -1, .complete = TRUE),
    addedSugarall    = slider::slide_dbl(addedSugarall  , sum, .before = 2, .after = -1, .complete = TRUE),
    fatall           = slider::slide_dbl(fatall         , sum, .before = 2, .after = -1, .complete = TRUE),
  )

head(past_two_days_diet_data)
# write.csv(past_two_days_diet_data, "/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/relabeled_data/past_two_days_diet_data_3_20.csv")

# Correlate with vegetarian status and vaginal microbiome
# vaginal.microbial.menses.24 <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifetyle/vaginal.microbial.menses.24.csv")
vaginal.microbial.menses.24 <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/relabeled_data/vaginal.microbial.menses.24.csv")
vaginal.microbial.menses.24 <- vaginal.microbial.menses.24[,-1]

vaginal.microbial.menses.24.veg <- vaginal.microbial.menses.24 %>% 
  left_join(merged_diet_data_vegetarian, by="biome_id")

vaginal.microbial.menses.24.summary <- vaginal.microbial.menses.24 %>% 
  group_by(biome_id) %>% 
  summarise(avg_shannon=sum(shannon)/n())

vaginal.microbial.menses.24.summary <- vaginal.microbial.menses.24.summary %>% 
  left_join(merged_diet_data_vegetarian, by="biome_id")

# participants with no diet data
vaginal.microbial.menses.24.summary %>% 
  filter(is.na(is_vegetarian)) # 15, 33 no diet data

vaginal.microbial.menses.24.summary <- vaginal.microbial.menses.24.summary %>% 
  filter(!is.na(is_vegetarian)) %>% 
  mutate(is_vegetarian_bin=ifelse(is_vegetarian==TRUE, "Vegetarian", "Non-vegetarian")) 

table(vaginal.microbial.menses.24.summary$is_vegetarian_bin)

ggplot(vaginal.microbial.menses.24.summary, aes(x = is_vegetarian_bin, y = avg_shannon)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.9, color="orchid") +
  labs(
    x = " ", 
    y = "Average Shannon Diversity", 
    title = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none") 

# testing
t.test(avg_shannon ~ is_vegetarian, data=vaginal.microbial.menses.24.summary)
wilcox.test(avg_shannon ~ is_vegetarian, data = vaginal.microbial.menses.24.summary)

##########################################################################################

# PERCENTAGE VEGETARIAN

veg_perc_df <- merged_diet_data %>% 
  group_by(biome_id, vegetarian) %>% 
  summarise(total_cal = sum(caloriesall)) %>% 
  ungroup()
veg_perc_df <- veg_perc_df %>% 
  group_by(biome_id) %>% 
  summarise(perc_veg = 100*sum(total_cal[vegetarian]) / sum(total_cal)) %>%
  ungroup()

vegetarian.diet.df <- merged_diet_data %>% 
  left_join(veg_perc_df, by="biome_id")

summary(vegetarian.diet.df$perc_veg)

# write.csv(
#   vegetarian.diet.df,
#   "/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/relabeled_data/nutrition_withVegPerc_data_3_20.csv"
# )

veg_perc_df <- vaginal.microbial.menses.24 %>% 
  left_join(veg_perc_df, by="biome_id")



# shannon and percent veg - cannot analyze like this
# lm.obj <- lm(shannon ~ perc_veg, data = veg_perc_df)
# summary(lm.obj)

# average shannon and percent veg
veg_perc_df_summary  <- veg_perc_df %>% 
  group_by(biome_id) %>% 
  summarise(avg_shannon = sum(shannon) / n(),
            perc_veg = first(perc_veg))
lm.obj.summary <- lm(avg_shannon ~ perc_veg, data = veg_perc_df_summary)
summary(lm.obj.summary)

#### Correlate with birth control
shannon.cst.qr.merged.24 <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/relabeled_data/shannon.cst.qr.merged.24.csv", header=TRUE)
shannon.cst.qr.merged.24 <- shannon.cst.qr.merged.24[,-1]
participant.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 9-Volunteer Medical History.csv", header = TRUE)

# select birth control
birthControl.df <- participant.data %>% 
  dplyr::select(biome_id, birthControl)
birthControl.collapsed <- birthControl.df %>% 
  count(birthControl, name="frequency")
shannon.birthControl <- shannon.cst.qr.merged.24 %>% 
  left_join(birthControl.df, by="biome_id") %>% 
  filter(!is.na(birthControl))

# add the average shannon idx for a participant to df
shannon.birthControl <- shannon.birthControl %>% 
  group_by(biome_id) %>% 
  mutate(avg_shannon=sum(shannon)/n())

# Collapse by person (assign to most frequent CST)
shannon.birthControl.collapsed <- shannon.birthControl %>%
  group_by(biome_id) %>% 
  count(CST, name="frequency") %>% # total 148 (participants have more than 1 CST) 
  slice_max(frequency, n=1) %>% # collapse to their most frequent CST
  rename(CST_max=CST) %>% 
  ungroup() %>%
  left_join(dplyr::select(shannon.birthControl, biome_id, birthControl, shannon, avg_shannon, CST), by = "biome_id") %>% 
  distinct(biome_id, .keep_all = TRUE)

shannon.birthControl.collapsed.diet <- shannon.birthControl.collapsed %>% 
  left_join(veg_perc_df_summary, by=c("biome_id", "avg_shannon"))

# Regress shannon diversity on diet + diet*birth control, if interaction terms (2) is significant
shannon.birthControl.collapsed.diet <- shannon.birthControl.collapsed.diet %>% 
  mutate(birthControl_collapsed=ifelse(birthControl=="Systemic Combined (E&P)" | (birthControl=="Systemic P only"), "Systemic",
                                       birthControl)) %>%
  mutate(birthControl=as.factor(birthControl),
         birthControl_collapsed=as.factor(birthControl_collapsed))

shannon.birthControl.collapsed.diet$birthControl <- factor(shannon.birthControl.collapsed.diet$birthControl, levels=c("None", "Local P", "Systemic P only", "Systemic Combined (E&P)"))
shannon.birthControl.collapsed.diet$birthControl_collapsed <- factor(shannon.birthControl.collapsed.diet$birthControl_collapsed, levels=c("None", "Local P", "Systemic"))

# Regressing models - birth control not collapsed
diet.birthControl <- lm(avg_shannon ~ perc_veg*birthControl, data=shannon.birthControl.collapsed.diet)
summary(diet.birthControl)

diet.birthControl2 <- lm(avg_shannon ~ perc_veg+birthControl, data=shannon.birthControl.collapsed.diet)
summary(diet.birthControl2)

anova(diet.birthControl,diet.birthControl2)

aov.obj <- aov(shannon.birthControl.collapsed.diet$avg_shannon ~ shannon.birthControl.collapsed.diet$birthControl)
summary(aov.obj)

# Regressing models - birth control collapsed
diet.birthControl_collapsed <- lm(avg_shannon ~ perc_veg*birthControl_collapsed, data=shannon.birthControl.collapsed.diet)
summary(diet.birthControl_collapsed)

diet.birthControl_collapsed2 <- lm(avg_shannon ~ perc_veg+birthControl_collapsed, data=shannon.birthControl.collapsed.diet)
summary(diet.birthControl_collapsed2)

anova(diet.birthControl_collapsed,diet.birthControl_collapsed2)

aov.diet.birthCtrl <- aov(shannon.birthControl.collapsed.diet$avg_shannon ~ shannon.birthControl.collapsed.diet$birthControl_collapsed)
summary(aov.diet.birthCtrl)

# Diet: scatter plot of percent vegetarian on shannon diversity
plot(veg_perc_df$shannon, veg_perc_df$perc_veg)
ggplot(veg_perc_df, aes(x = shannon, y = perc_veg)) +
  geom_point() +
  labs(x = "Average Shannon Diversity", y = "Percent Vegetarian", 
       title = "Scatter Plot of Percent Vegetarian vs. Shannon Diversity") +
  theme_minimal()

# Diet: scatter plot of percent vegetarian on average shannon diversity
ggplot(shannon.birthControl.collapsed.diet, aes(x = perc_veg, y = avg_shannon)) +
  geom_point() +
  labs(x = "Percent Vegetarian", y = "Average Shannon Diversity", 
       title = " ") +
  theme_minimal()

# Diet: scatter plot of percent vegetarian on average shannon diversity colored by Birth control
shannon.birthControl.collapsed.diet %>%
  filter(!is.na(birthControl)) %>% 
  ggplot(aes(x = perc_veg, y = avg_shannon, color=birthControl)) +
  geom_point() +
  labs(x = "Percent Vegetarian", y = "Average Shannon Diversity", 
       title = "Scatter Plot of Average Shannon Diversity v. Percent Vegetarian") +
  theme_minimal()

# Diet: scatter plot of percent vegetarian on average shannon diversity colored by collapsed Birth control
shannon.birthControl.collapsed.diet %>%
  filter(!is.na(birthControl_collapsed)) %>% 
  ggplot(aes(x = perc_veg, y = avg_shannon, color=birthControl_collapsed)) +
  geom_point() +
  labs(x = "Percent Vegetarian", y = "Average Shannon Diversity", 
       title = "Scatter Plot of Average Shannon Diversity v. Percent Vegetarian") +
  theme_minimal()

##########################################################################################

# vaginal microbiota and specific nutrient intake

nutrient.diet.22 <-  merged_diet_data %>% #past_two_days_diet_data %>% #
  group_by(biome_id, Date) %>% 
  summarise(cholesterol_prop = sum(cholesterolall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            satFat_prop = sum(saturatedFatall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            sodium_prop = sum(sodiumall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            carb_prop = sum(carbohydratesall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            dietFib_prop = sum(dietaryFiberall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            sugar_prop = sum(sugarsall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            protein_prop = sum(proteinall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            fat_prop = sum(fatall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            # sat fat of all fat prop
            satfat_prop_fat = sum(saturatedFatall, na.rm=TRUE) / sum(fatall, na.rm=TRUE),
            fat_cal_prop = sum(caloriesFromFat, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE)) %>% 
  ungroup() %>% 
  rename(logDate=Date) %>% 
  mutate(logDate=as.character(logDate))

shannon.cst.qr.merged.24.collapsed <- shannon.cst.qr.merged.24 %>% 
  group_by(biome_id) %>% 
  summarise(avg_shannon = sum(shannon)/n(),
            max_CST = names(sort(table(CST), decreasing = TRUE))[1])

nutrient.diet.22.microbiome <- shannon.cst.qr.merged.24 %>% 
  left_join(nutrient.diet.22, by = c("biome_id", "logDate"))

# regress the nutrients on shannon -- can't do this kind of model because participant level data
# lm.nutrients <- lm(avg_shannon ~ ., data = nutrient.diet.22.microbiome[,-c(1,3)])
# lm.nutrients <- lm(shannon ~ ., data = nutrient.diet.22.microbiome[,c(1,)])
# summary(lm.nutrients)
# visualize significant nutrients
# ggplot(nutrient.diet.22.microbiome, aes(x = cholesterol_prop, y = avg_shannon)) +
#   geom_point() +
#   geom_smooth(method = "lm", col = "blue", se=FALSE) +
#   labs(title = "", x = "Cholesterol Intake (Proportion of Total Calories)", y = "Average Shannon Diversity") +
#   theme_minimal()
# 
# ggplot(nutrient.diet.22.microbiome, aes(x = protein_prop, y = avg_shannon)) +
#   geom_point() +
#   geom_smooth(method = "lm", col = "blue", se=FALSE) +
#   labs(title = "", x = "Protein Intake (Proportion of Total Calories)", y = "Average Shannon Diversity") +
#   theme_minimal()
# 
# ggplot(nutrient.diet.22.microbiome, aes(x = satFat_prop, y = avg_shannon)) +
#   geom_point() +
#   geom_smooth(method = "lm", col = "blue", se=FALSE) +
#   labs(title = "", x = "Saturated Fat (Proportion of Total Calories)", y = "Average Shannon Diversity") +
#   theme_minimal()


##########################################################################################
# Regression analysis

# Diet nutrients: pairwise plots for correlation between nutrients
macronutrients <- nutrient.diet.22.microbiome %>% 
  dplyr::select(-c(SampleID, CST,  biome_id, shannon, qr, logDate, status, timestamp, sampleType, max_taxa, OTU, CST_max))
  # dplyr::select(!c(biome_id, avg_shannon, max_CST))

macronutrients_renamed <- macronutrients %>%
  rename(
    Cholesterol = cholesterol_prop,
    `Saturated Fat` = satFat_prop,
    Sodium = sodium_prop,
    Carbohydrates = carb_prop,
    `Dietary Fiber` = dietFib_prop,
    Sugar = sugar_prop,
    Protein = protein_prop,
    Fat = fat_prop,
  )

names(macronutrients_renamed)

# Figure: corr plots of macronutrients
pairs(macronutrients_renamed[,-c(9,10)])

# # regress the nutrients on CST
# library(nnet)
# library(glmmTMB)
# library(mme)
# multinom.nutrients.CST <- multinom(CST ~ ., data = nutrient.diet.22.microbiome[,-c(1,3,4,5,6,7,8,9,10,11)])
# summary(multinom.nutrients.CST)

# library(car)
# # Wald test on coefficients
# Anova(multinom.nutrients.CST, type = "II")
# 
# pred_probs <- predict(multinom.nutrients.CST, nutrient.diet.22.microbiome, type = "probs")
# nutrient.diet.22.microbiome_pred_probs <- cbind(nutrient.diet.22.microbiome, pred_probs)
# 
# # format into matrix
# y_matrix <- model.matrix(~ 0 + CST, data = nutrient.diet.22.microbiome)
# x_matrix <- model.matrix(~ ., data = nutrient.diet.22.microbiome[,-c(1,3,4,5,6,7,8,9,10,11,12)])
# group <- as.factor(nutrient.diet.22.microbiome$biome_id)
# 
# names(nutrient.diet.22.microbiome)
# head(nutrient.diet.22.microbiome[,-c(1,3, 4, 7:12)])
# 
# nutrient.diet.22.microbiome.subset <- 
#   nutrient.diet.22.microbiome[,-c(1,3, 4, 7:12)] %>% 
#   rename(satfatToFat_prop = satfat_prop_fat,
#          fatcal_prop = fat_cal_prop)

# Mixed effects multinomial model
# install.packages("brms")
# library(brms)
# 
# model <- brm(
#   formula = CST ~ cholesterol_prop + satFat_prop + sodium_prop + carb_prop + dietFib_prop + sugar_prop + 
#     protein_prop + fat_prop + satfatToFat_prop + fatcal_prop + (1 | biome_id),
#   data = nutrient.diet.22.microbiome.subset,
#   family = categorical(),
#   # prior = c(
#   #   set_prior("normal(0, 5)", class = "b"),    # Priors for fixed effects (predictors)
#   #   set_prior("normal(0, 5)", class = "sd"),   # Priors for random effects (biome_id)
#   #   set_prior("normal(0, 5)", class = "Intercept")  # Priors for intercepts in multinomial model
#   # ),
#   chains = 4,
#   iter = 2000
# )
# summary(model)
# 
# # attempt 2
# 
# model <- brm(
#   formula = CST ~ cholesterol_prop + satFat_prop + sodium_prop + carb_prop + dietFib_prop + sugar_prop + 
#     protein_prop + fat_prop + satfatToFat_prop + fatcal_prop + (1 | biome_id),
#   data = nutrient.diet.22.microbiome.subset,
#   family = categorical(),
#   # prior = c(
#   #   set_prior("normal(0, 5)", class = "b"),    # Priors for fixed effects (predictors)
#   #   set_prior("normal(0, 5)", class = "sd"),   # Priors for random effects (biome_id)
#   #   set_prior("normal(0, 5)", class = "Intercept")  # Priors for intercepts in multinomial model
#   # ),
#   chains = 5,
#   iter = 5000
# )
# summary(model)

head(nutrient.diet.22.microbiome.subset)

## visualize significant nutrients
# dietary fiber
ggplot(nutrient.diet.22.microbiome_pred_probs, aes(x = dietFib_prop)) +
  geom_line(aes(y = pred_probs[,1], color = "CST 1")) +
  geom_line(aes(y = pred_probs[,2], color = "CST 2")) +
  geom_line(aes(y = pred_probs[,3], color = "CST 3")) +
  geom_line(aes(y = pred_probs[,4], color = "CST 4")) +
  geom_line(aes(y = pred_probs[,5], color = "CST 5")) +
  labs(title = "", x = "Dietary Fiber Intake", y = "Predicted Probability")

# carbs
ggplot(nutrient.diet.22.microbiome_pred_probs, aes(x = carb_prop)) +
  geom_line(aes(y = pred_probs[,1], color = "CST 1")) +
  geom_line(aes(y = pred_probs[,2], color = "CST 2")) +
  geom_line(aes(y = pred_probs[,3], color = "CST 3")) +
  geom_line(aes(y = pred_probs[,4], color = "CST 4")) +
  geom_line(aes(y = pred_probs[,5], color = "CST 5")) +
  labs(title = "", x = "Carbohydrate Intake", y = "Predicted Probability")

# cholesterol_prop
ggplot(nutrient.diet.22.microbiome_pred_probs, aes(x = cholesterol_prop)) +
  geom_line(aes(y = pred_probs[,1], color = "CST 1")) +
  geom_line(aes(y = pred_probs[,2], color = "CST 2")) +
  geom_line(aes(y = pred_probs[,3], color = "CST 3")) +
  geom_line(aes(y = pred_probs[,4], color = "CST 4")) +
  geom_line(aes(y = pred_probs[,5], color = "CST 5")) +
  labs(title = "", x = "Cholesterol Intake", y = "Predicted Probability")

# sugar_prop
ggplot(nutrient.diet.22.microbiome_pred_probs, aes(x = sugar_prop)) +
  geom_line(aes(y = pred_probs[,1], color = "CST 1")) +
  geom_line(aes(y = pred_probs[,2], color = "CST 2")) +
  geom_line(aes(y = pred_probs[,3], color = "CST 3")) +
  geom_line(aes(y = pred_probs[,4], color = "CST 4")) +
  geom_line(aes(y = pred_probs[,5], color = "CST 5")) +
  labs(title = "", x = "Sugar Intake", y = "Predicted Probability")

# protein_prop
ggplot(nutrient.diet.22.microbiome_pred_probs, aes(x = protein_prop)) +
  geom_line(aes(y = pred_probs[,1], color = "CST 1")) +
  geom_line(aes(y = pred_probs[,2], color = "CST 2")) +
  geom_line(aes(y = pred_probs[,3], color = "CST 3")) +
  geom_line(aes(y = pred_probs[,4], color = "CST 4")) +
  geom_line(aes(y = pred_probs[,5], color = "CST 5")) +
  labs(title = "", x = "Protein Intake", y = "Predicted Probability")

##########################################################################################

### Mixed Effects Models
library(lme4)
library(lmerTest)
library(performance)

# Shannon diversity and nutrients
# vaginal microbiota and specific nutrient intake

nutrient.diet.22.day <- merged_diet_data %>% # past_two_days_diet_data %>% #
  group_by(biome_id, Date) %>% 
  summarise(caloriesall_avg = sum(caloriesall, na.rm=TRUE) / n(), # get avg cals per day
            cholesterol_prop = sum(cholesterolall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            satFat_prop = sum(saturatedFatall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            sodium_prop = sum(sodiumall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            carb_prop = sum(carbohydratesall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            dietFib_prop = sum(dietaryFiberall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            sugar_prop = sum(sugarsall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            protein_prop = sum(proteinall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            fat_prop = sum(fatall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            fat_cal_prop = sum(caloriesFromFat, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            addedSugarall_prop = sum(addedSugarall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE)) %>% 
  ungroup() %>% 
  rename(logDate = Date) %>% 
  mutate(logDate = as.character(logDate))

# Merge data for samples and nutrition data for the day
nutrient.diet.22.day.microbiome <- shannon.cst.qr.merged.24 %>% 
  left_join(nutrient.diet.22.day, by = c("biome_id", "logDate"))

nutrient.diet.22.day.microbiome.filtered <- nutrient.diet.22.day.microbiome %>% 
  group_by(biome_id) #%>% 
  # filter(n() > 10)
table(nutrient.diet.22.day.microbiome.filtered$biome_id)

## Different baseline, same slope, different intercept

# num participants
length(unique(nutrient.diet.22.day.microbiome.filtered$biome_id))

# Full model 
lmer.full <- lmer(shannon~caloriesall_avg+cholesterol_prop+satFat_prop+
                    sodium_prop+carb_prop+dietFib_prop+sugar_prop+
                    protein_prop+fat_prop+fat_cal_prop+addedSugarall_prop+(1|`biome_id`), 
                  data = nutrient.diet.22.day.microbiome.filtered)
summary(lmer.full)
r2(lmer.full)

# vaginal.microbial.menses.24.veg
nutrient.diet.22.day.microbiome.filtered <- nutrient.diet.22.day.microbiome.filtered %>% 
  mutate(caloriesall_avg_scale = scale(caloriesall_avg))

# Full model (no avg cal)
lmer.full <- lmer(shannon~caloriesall_avg_scale+cholesterol_prop+satFat_prop+
                    sodium_prop+carb_prop+dietFib_prop+sugar_prop+
                    protein_prop+fat_prop+fat_cal_prop+addedSugarall_prop+(1|`biome_id`), 
                  data = nutrient.diet.22.day.microbiome.filtered)
summary(lmer.full)
r2(lmer.full)

lmer.avg.cal <- lmer(shannon~caloriesall_avg+(1|`biome_id`), data=nutrient.diet.22.day.microbiome.filtered)
lmer.choles  <- lmer(shannon~cholesterol_prop+(1|`biome_id`),  data=nutrient.diet.22.day.microbiome.filtered)
lmer.satfat  <- lmer(shannon~satFat_prop+(1|`biome_id`),  data=nutrient.diet.22.day.microbiome.filtered)
lmer.sodium  <- lmer(shannon~sodium_prop+(1|`biome_id`),  data=nutrient.diet.22.day.microbiome.filtered)
lmer.carbs   <- lmer(shannon~carb_prop+(1|`biome_id`),  data=nutrient.diet.22.day.microbiome.filtered)
lmer.dietfib <- lmer(shannon~dietFib_prop+(1|`biome_id`),  data=nutrient.diet.22.day.microbiome.filtered)
lmer.sugar   <- lmer(shannon~sugar_prop+(1|`biome_id`),  data=nutrient.diet.22.day.microbiome.filtered)
lmer.prot    <- lmer(shannon~protein_prop+(1|`biome_id`),  data=nutrient.diet.22.day.microbiome.filtered)
lmer.fat     <- lmer(shannon~fat_prop+(1|`biome_id`), data=nutrient.diet.22.day.microbiome.filtered)
lmer.fat.cal <- lmer(shannon~fat_cal_prop+(1|`biome_id`), data=nutrient.diet.22.day.microbiome.filtered)
lmer.addsug  <- lmer(shannon~addedSugarall_prop+(1|`biome_id`), data=nutrient.diet.22.day.microbiome.filtered)

r2(lmer.avg.cal)
r2(lmer.choles ) # 52% of variation explained by which person is which and by the amount of cholesterol
r2(lmer.satfat )
r2(lmer.sodium )

r2(lmer.carbs  )
r2(lmer.dietfib)
r2(lmer.sugar  )
r2(lmer.prot   )
r2(lmer.fat    )
r2(lmer.fat.cal)
r2(lmer.addsug )

summary(lmer.avg.cal)
summary(lmer.choles )
summary(lmer.satfat )
summary(lmer.sodium )
summary(lmer.carbs  )
summary(lmer.dietfib)
summary(lmer.sugar  )
summary(lmer.prot   )
summary(lmer.fat    )
summary(lmer.fat.cal)
summary(lmer.addsug )

## Random slopes model, different slope, different intercept

# Full model (no avg cal)
# rs_full_model <- lmer(shannon~cholesterol_prop+(cholesterol_prop|`biome_id`) +
#                         satFat_prop+ (satFat_prop|`biome_id`) +
#                         # sodium_prop+ (sodium_prop||`biome_id`) +
#                         # carb_prop +(carb_prop|`biome_id`) +
#                         # dietFib_prop+  (dietFib_prop|`biome_id`) 
#                         #sugar_prop+ #(sugar_prop|`biome_id`)+
#                         # protein_prop + #(protein_prop|`biome_id`)+
#                         # fat_cal_prop + #(fat_cal_prop|`biome_id`)+
#                         # addedSugarall_prop+ #(addedSugarall_prop||`biome_id`)+
#                         # fat_prop+(fat_prop|`biome_id`)
#                       ,  data=nutrient.diet.22.day.microbiome.filtered)
# r2(rs_full_model)
# summary(rs_full_model)

# rs_lmer.avg.cal <- lmer(shannon~caloriesall_avg_scale+(caloriesall_avg_scale||biome_id), data=nutrient.diet.22.day.microbiome.filtered)
# rs_lmer.choles  <- lmer(shannon~cholesterol_prop+(cholesterol_prop|`biome_id`),  data=nutrient.diet.22.day.microbiome.filtered)
# rs_lmer.satfat  <- lmer(shannon~satFat_prop+(satFat_prop|`biome_id`),  data=nutrient.diet.22.day.microbiome.filtered)
# rs_lmer.sodium  <- lmer(shannon~sodium_prop+(sodium_prop||`biome_id`),  data=nutrient.diet.22.day.microbiome.filtered)
# rs_lmer.carbs   <- lmer(shannon~carb_prop+(carb_prop|`biome_id`),  data=nutrient.diet.22.day.microbiome.filtered)
# rs_lmer.dietfib <- lmer(shannon~dietFib_prop+(dietFib_prop|`biome_id`),  data=nutrient.diet.22.day.microbiome.filtered)
# rs_lmer.sugar   <- lmer(shannon~sugar_prop+(sugar_prop|`biome_id`),  data=nutrient.diet.22.day.microbiome.filtered)
# rs_lmer.prot    <- lmer(shannon~protein_prop+(protein_prop|`biome_id`),  data=nutrient.diet.22.day.microbiome.filtered)
# rs_lmer.fat     <- lmer(shannon~fat_prop+(fat_prop|`biome_id`), data=nutrient.diet.22.day.microbiome.filtered)
# rs_lmer.fat.cal <- lmer(shannon~fat_cal_prop+(fat_cal_prop|`biome_id`), data=nutrient.diet.22.day.microbiome.filtered)
# rs_lmer.addsug  <- lmer(shannon~addedSugarall_prop+(addedSugarall_prop||`biome_id`), data=nutrient.diet.22.day.microbiome.filtered)
# 
# r2(rs_lmer.avg.cal)
# r2(rs_lmer.choles )
# r2(rs_lmer.satfat )
# r2(rs_lmer.sodium )
# r2(rs_lmer.carbs  )
# r2(rs_lmer.dietfib)
# r2(rs_lmer.prot)
# r2(rs_lmer.fat    )
# r2(rs_lmer.sugar  )
# r2(rs_lmer.fat.cal)
# r2(rs_lmer.addsug )
# 
# summary(rs_lmer.avg.cal)
# summary(rs_lmer.choles )
# summary(rs_lmer.satfat )
# summary(rs_lmer.sodium )
# summary(rs_lmer.carbs  )
# summary(rs_lmer.dietfib)
# summary(rs_lmer.sugar  )
# summary(rs_lmer.prot   )
# summary(rs_lmer.fat    )
# summary(rs_lmer.fat.cal)
# summary(rs_lmer.addsug )
# 
# # anova(lmer.avg.cal, rs_lmer.avg.cal)
# anova(lmer.choles, rs_lmer.choles)
# anova(lmer.satfat, rs_lmer.satfat)
# anova(lmer.sodium, rs_lmer.sodium)
# anova(lmer.carbs, rs_lmer.carbs)
# anova(lmer.dietfib, rs_lmer.dietfib)# slightly sign 0.076
# anova(lmer.sugar, rs_lmer.sugar)
# anova(lmer.prot, rs_lmer.prot)
# anova(lmer.fat, rs_lmer.fat)
# anova(lmer.fat.cal, rs_lmer.fat.cal)
# anova(lmer.addsug, rs_lmer.addsug)

##########################################################################################

## Add time into model
nutrient.diet <- study_days(nutrient.diet.22.day.microbiome.filtered)

# Full model (no avg cal) - fixed slopes
lmer.full <- lmer(shannon~cholesterol_prop +
                    satFat_prop + sodium_prop + carb_prop +
                    dietFib_prop + sugar_prop + protein_prop +
                    fat_prop + addedSugarall_prop + caloriesall_avg +
                    (1|`biome_id`), 
                  data = nutrient.diet)
summary(lmer.full)
r2(lmer.full)

## Time
nutrients.time.model <- lmer(shannon ~ study_day + cholesterol_prop +
                               satFat_prop + sodium_prop + carb_prop +
                               dietFib_prop + sugar_prop + protein_prop +
                               fat_prop + addedSugarall_prop + caloriesall_avg +
                               (1 | biome_id),
                             data = nutrient.diet)
r2(nutrients.time.model)
summary(nutrients.time.model)
anova(lmer.full, nutrients.time.model)

## Time and Time^2
nutrients.time.model2 <- lmer(shannon ~ study_day + I(study_day^2) + cholesterol_prop +
                                satFat_prop + sodium_prop + carb_prop +
                                dietFib_prop + sugar_prop + protein_prop +
                               fat_prop + addedSugarall_prop + caloriesall_avg +
                              (1 | biome_id),
                              data = nutrient.diet)
r2(nutrients.time.model2)
summary(nutrients.time.model2)

anova(lmer.full, nutrients.time.model2)
anova(nutrients.time.model, nutrients.time.model2)

##########################################################################################
library(ggrepel)

## Genus level analyses
bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/relabeled_data/vaginal_cleaned_max_taxa.rds")


bacterial.taxa <- tax_table(bacterial.data)

sample_data_df <- as(sample_data(genus.ra), "data.frame")
length(unique(sample_data_df$biome_id))

# aggregate at genus level & relative abundances
genus.taxa <- tax_glom(bacterial.data, taxrank = "Genus")
genus.ra <- transform_sample_counts(genus.taxa, function(x) x / sum(x)) 
bactieral.meta <- as(sample_data(genus.ra), "data.frame")

# make meta data obj
veg_perc_df <- veg_perc_df %>% 
  mutate(is_vegetarian=ifelse(perc_veg ==100, TRUE, FALSE))
veg_perc_df.subset <- veg_perc_df %>%
  dplyr::select(SampleID, CST, shannon, biome_id, logDate, sampleType, perc_veg, is_vegetarian)
genus.ra.merged <- bactieral.meta %>% 
  left_join(veg_perc_df.subset, by=c("SampleID", "biome_id", "logDate", "sampleType"))
nutrient.diet.22.day.microbiome.filtered.subset <- nutrient.diet.22.day.microbiome.filtered %>% 
  dplyr::select(!c(qr, status, timestamp, sampleType))
genus.ra.nutrient.merged <- genus.ra.merged %>% 
  left_join(nutrient.diet.22.day.microbiome.filtered.subset, by=c("SampleID", "biome_id", "logDate", "CST", "shannon"))
rownames(genus.ra.nutrient.merged) <- genus.ra.nutrient.merged$SampleID
sample_data(genus.ra) <- sample_data(genus.ra.nutrient.merged)

# ordination
genus.ordination <- ordinate(genus.ra, method = "PCoA", distance = "bray")

# Diet: ordincation cluster by vegetarian
plot_ordination(genus.ra, genus.ordination, color = "is_vegetarian") + 
  geom_point(size=3) +
  labs(color="Vegetarian Status")+
  scale_color_viridis(discrete=TRUE)+
  geom_text_repel(aes(label = biome_id), size = 3) +
  theme_minimal()

# Diet: cluster by percent vegetarian
plot_ordination(genus.ra, genus.ordination, color = "perc_veg") + 
  geom_point(size=3) +
  geom_text_repel(aes(label = biome_id), size = 3) +
  labs(color="Percent Vegetarian")+
  scale_color_viridis()+
  theme_minimal()

# Diet: cluster by Cholesterol Intake
plot_ordination(genus.ra, genus.ordination, color = "cholesterol_prop") + 
  geom_point(size=3) +
  geom_text_repel(aes(label = biome_id), size = 3) +
  labs(color="Percent Vegetarian")+
  scale_color_viridis()+
  theme_minimal()

# clean df
meta_df <- as(sample_data(genus.ra), "data.frame")
meta_df_clean <- na.omit(meta_df)
genus.ra.clean <- prune_samples(rownames(meta_df_clean), genus.ra)

# analyze if microbiome composition (Bray-Curtis distances) is significantly associated with the dietary factors
vegan::adonis2(distance(genus.ra.clean, method="bray") ~ cholesterol_prop + satFat_prop, #+ 
                 # sodium_prop + carb_prop + dietFib_prop + sugar_prop +
                 # protein_prop + fat_prop + fat_cal_prop + addedSugarall_prop, 
               data = meta_df_clean)

# time and dietary category
genus.df <- psmelt(genus.ra)
genus.df.time <- study_days(genus.df)

lmer_genus <- lmer(Abundance ~ study_day + Abundance + (1 | biome_id), data = genus.df.time)
summary(lmer_genus)

##########################################################################################
merged_diet_data_daily2 <- merged_diet_data_daily %>% 
  rename(logDate=Date) %>% 
  mutate(logDate = as.character(logDate))
vaginal.shannon <- vaginal.microbial.menses.24 %>% 
  dplyr::select(SampleID, shannon, biome_id, logDate)
diet.vaginal.shannon <- vaginal.shannon %>% 
  left_join(merged_diet_data_daily2, by=c("biome_id", "logDate"))

# Figure: Log Macronutrients Plot with Shannon Diversity Index
par(mfrow = c(4,3))

plot(col = diet.vaginal.shannon$biome_id, y = diet.vaginal.shannon$shannon , x = log(diet.vaginal.shannon$caloriesall))
plot(col = diet.vaginal.shannon$biome_id, y = diet.vaginal.shannon$shannon , x = log(diet.vaginal.shannon$cholesterolall))
plot(col = diet.vaginal.shannon$biome_id, y = diet.vaginal.shannon$shannon , x = log(diet.vaginal.shannon$saturatedFatall))
plot(col = diet.vaginal.shannon$biome_id, y = diet.vaginal.shannon$shannon , x = log(diet.vaginal.shannon$sodiumall))
plot(col = diet.vaginal.shannon$biome_id, y = diet.vaginal.shannon$shannon , x = log(diet.vaginal.shannon$carbohydratesall))
plot(col = diet.vaginal.shannon$biome_id, y = diet.vaginal.shannon$shannon , x = log(diet.vaginal.shannon$dietaryFiberall))
plot(col = diet.vaginal.shannon$biome_id, y = diet.vaginal.shannon$shannon , x = log(diet.vaginal.shannon$sugarsall))
plot(col = diet.vaginal.shannon$biome_id, y = diet.vaginal.shannon$shannon , x = log(diet.vaginal.shannon$proteinall))
plot(col = diet.vaginal.shannon$biome_id, y = diet.vaginal.shannon$shannon , x = log(diet.vaginal.shannon$fatall))
plot(col = diet.vaginal.shannon$biome_id, y = diet.vaginal.shannon$shannon , x = log(diet.vaginal.shannon$saturatedFat))
plot(col = diet.vaginal.shannon$biome_id, y = diet.vaginal.shannon$shannon , x = log(diet.vaginal.shannon$transFat))
plot(col = diet.vaginal.shannon$biome_id, y = diet.vaginal.shannon$shannon , x = log(diet.vaginal.shannon$addedSugarall))

par(mfrow = c(1,1))

par(mfrow = c(2,3))
plot(y=nutrient.diet$shannon, x= log((nutrient.diet$cholesterol_prop)/(1-(nutrient.diet$cholesterol_prop)))) # nutrient.diet$cholesterol_prop)
plot(y=nutrient.diet$shannon, x= log((nutrient.diet$satFat_prop)/(1-(nutrient.diet$satFat_prop)))) #nutrient.diet$satFat_prop)
plot(y=nutrient.diet$shannon, x= log((nutrient.diet$sodium_prop)/(1-(nutrient.diet$sodium_prop)))) # nutrient.diet$sodium_prop)
plot(y=nutrient.diet$shannon, x= log((nutrient.diet$carb_prop)/(1-(nutrient.diet$carb_prop)))) # nutrient.diet$carb_prop)
plot(y=nutrient.diet$shannon, x= log((nutrient.diet$dietFib_prop)/(1-(nutrient.diet$dietFib_prop)))) # nutrient.diet$dietFib_prop)
plot(y=nutrient.diet$shannon, x= log((nutrient.diet$sugar_prop)/(1-(nutrient.diet$sugar_prop)))) # nutrient.diet$sugar_prop)
plot(y=nutrient.diet$shannon, x= log((nutrient.diet$protein_prop)/(1-(nutrient.diet$protein_prop)))) # nutrient.diet$protein_prop)
plot(y=nutrient.diet$shannon, x= log((nutrient.diet$fat_prop)/(1-(nutrient.diet$fat_prop)))) # nutrient.diet$fat_prop)
plot(y=nutrient.diet$shannon, x= log((nutrient.diet$fat_cal_prop)/(1-(nutrient.diet$fat_cal_prop)))) # nutrient.diet$fat_cal_prop)
plot(y=nutrient.diet$shannon, x=log((nutrient.diet$addedSugarall_prop)/(1-(nutrient.diet$addedSugarall_prop))))




