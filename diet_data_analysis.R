########################
#
# Exploration of the Merged Diet Data
# 11 September 2024
#
#########################

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

all_days <- seq.Date(as.Date(min(merged_diet_data$Date, na.rm=TRUE)), 
                     as.Date(max(merged_diet_data$Date, na.rm=TRUE)), by = "day")

# Expand to include all days
expanded_data <- expand.grid(biome_id = unique(merged_diet_data$biome_id), Date = all_days)
expanded_data <- left_join(expanded_data, merged_diet_data, by = c("biome_id", "Date"))

diet_data_subset <- merged_diet_data %>% 
  select(biome_id, Date, type, name, caloriesall, proteinall, sugarsall, fatall, food_other)

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
  select(-name_lower)

unique_vegetarian_pairs <- merged_diet_data %>%
  select(name, vegetarian) %>%  
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
  select(biome_id, is_vegetarian) %>% 
  distinct(biome_id, is_vegetarian)

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
  geom_jitter(width = 0.2, alpha = 0.6, color="orchid") +
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
  select(biome_id, birthControl)
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
  left_join(select(shannon.birthControl, biome_id, birthControl, shannon, avg_shannon, CST), by = "biome_id") %>% 
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

##########################################################################################


# vaginal microbiota and specific nutrient intake

nutrient.diet.22 <- merged_diet_data %>% 
  group_by(biome_id) %>% 
  summarise(cholesterol_prop = sum(cholesterolall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            satFat_prop = sum(saturatedFatall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            sodium_prop = sum(sodiumall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            carb_prop = sum(carbohydratesall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            dietFib_prop = sum(dietaryFiberall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            sugar_prop = sum(sugarsall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            protein_prop = sum(proteinall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            fat_prop = sum(fatall, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE),
            fat_cal_prop = sum(caloriesFromFat, na.rm=TRUE) / sum(caloriesall, na.rm=TRUE)) %>% 
  ungroup()

shannon.cst.qr.merged.24.collapsed <- shannon.cst.qr.merged.24 %>% 
  group_by(biome_id) %>% 
  summarise(avg_shannon = sum(shannon)/n(),
            max_CST = names(sort(table(CST), decreasing = TRUE))[1])

nutrient.diet.22.microbiome <- shannon.cst.qr.merged.24.collapsed %>% 
  left_join(nutrient.diet.22, by = "biome_id")

# regress the nutrients on avg_shannon
lm.nutrients <- lm(avg_shannon ~ ., data = nutrient.diet.22.microbiome[,-c(1,3)])
summary(lm.nutrients)

# visualize significant nutrients
ggplot(nutrient.diet.22.microbiome, aes(x = cholesterol_prop, y = avg_shannon)) +
  geom_point() +
  geom_smooth(method = "lm", col = "blue", se=FALSE) +
  labs(title = "", x = "Cholesterol Intake (Proportion of Total Calories)", y = "Average Shannon Diversity") +
  theme_minimal()

ggplot(nutrient.diet.22.microbiome, aes(x = protein_prop, y = avg_shannon)) +
  geom_point() +
  geom_smooth(method = "lm", col = "blue", se=FALSE) +
  labs(title = "", x = "Protein Intake (Proportion of Total Calories)", y = "Average Shannon Diversity") +
  theme_minimal()

ggplot(nutrient.diet.22.microbiome, aes(x = satFat_prop, y = avg_shannon)) +
  geom_point() +
  geom_smooth(method = "lm", col = "blue", se=FALSE) +
  labs(title = "", x = "Saturated Fat (Proportion of Total Calories)", y = "Average Shannon Diversity") +
  theme_minimal()

# regress the nutrients on CST
library(nnet)

multinom.nutrients.CST <- multinom(max_CST ~ ., data = nutrient.diet.22.microbiome[,-c(1,2)])
summary(multinom.nutrients.CST)

library(car)
# Wald test on coefficients
Anova(multinom.nutrients.CST, type = "II")

pred_probs <- predict(multinom.nutrients.CST, nutrient.diet.22.microbiome, type = "probs")
nutrient.diet.22.microbiome_pred_probs <- cbind(nutrient.diet.22.microbiome, pred_probs)

# visualize significant nutrients
ggplot(nutrient.diet.22.microbiome_pred_probs, aes(x = dietFib_prop)) +
  geom_line(aes(y = pred_probs[,1], color = "CST 1")) +
  geom_line(aes(y = pred_probs[,2], color = "CST 2")) +
  geom_line(aes(y = pred_probs[,3], color = "CST 3")) +
  geom_line(aes(y = pred_probs[,4], color = "CST 4")) +
  geom_line(aes(y = pred_probs[,5], color = "CST 5")) +
  labs(title = "", x = "Dietary Fiber Intake", y = "Predicted Probability")

##########################################################################################

### Mixed Effects Models
library(lme4)
library(lmerTest)

# Shannon diversity and nutrients
# vaginal microbiota and specific nutrient intake

nutrient.diet.22.day <- merged_diet_data %>% 
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
  group_by(biome_id) %>% 
  filter(n() > 10)
table(nutrient.diet.22.day.microbiome.filtered$biome_id)

## Different baseline, same slope, different intercept

# vaginal.microbial.menses.24.veg
nutrient.diet.22.day.microbiome.filtered <- nutrient.diet.22.day.microbiome.filtered %>% 
  mutate(caloriesall_avg_scale = scale(caloriesall_avg))

lmer.avg.cal <- lmer(shannon~caloriesall_avg_scale+(1|`biome_id`), data=nutrient.diet.22.day.microbiome.filtered)
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
rs_lmer.avg.cal <- lmer(shannon~caloriesall_avg_scale+(caloriesall_avg_scale||biome_id), data=nutrient.diet.22.day.microbiome.filtered)
rs_lmer.choles  <- lmer(shannon~cholesterol_prop+(cholesterol_prop|`biome_id`),  data=nutrient.diet.22.day.microbiome.filtered)
rs_lmer.satfat  <- lmer(shannon~satFat_prop+(satFat_prop|`biome_id`),  data=nutrient.diet.22.day.microbiome.filtered)
rs_lmer.sodium  <- lmer(shannon~sodium_prop+(sodium_prop||`biome_id`),  data=nutrient.diet.22.day.microbiome.filtered)
rs_lmer.carbs   <- lmer(shannon~carb_prop+(carb_prop|`biome_id`),  data=nutrient.diet.22.day.microbiome.filtered)
rs_lmer.dietfib <- lmer(shannon~dietFib_prop+(dietFib_prop|`biome_id`),  data=nutrient.diet.22.day.microbiome.filtered)
rs_lmer.sugar   <- lmer(shannon~sugar_prop+(sugar_prop|`biome_id`),  data=nutrient.diet.22.day.microbiome.filtered)
rs_lmer.prot    <- lmer(shannon~protein_prop+(protein_prop|`biome_id`),  data=nutrient.diet.22.day.microbiome.filtered)
rs_lmer.fat     <- lmer(shannon~fat_prop+(fat_prop|`biome_id`), data=nutrient.diet.22.day.microbiome.filtered)
rs_lmer.fat.cal <- lmer(shannon~fat_cal_prop+(fat_cal_prop|`biome_id`), data=nutrient.diet.22.day.microbiome.filtered)
rs_lmer.addsug  <- lmer(shannon~addedSugarall_prop+(addedSugarall_prop||`biome_id`), data=nutrient.diet.22.day.microbiome.filtered)

summary(rs_lmer.avg.cal)
summary(rs_lmer.choles )
summary(rs_lmer.satfat )
summary(rs_lmer.sodium )
summary(rs_lmer.carbs  )
summary(rs_lmer.dietfib)
summary(rs_lmer.sugar  )
summary(rs_lmer.prot   )
summary(rs_lmer.fat    )
summary(rs_lmer.fat.cal)
summary(rs_lmer.addsug )

anova(lmer.avg.cal, rs_lmer.avg.cal)
anova(lmer.choles, rs_lmer.choles)
anova(lmer.satfat, rs_lmer.satfat)
anova(lmer.sodium, rs_lmer.sodium)
anova(lmer.carbs, rs_lmer.carbs)
anova(lmer.dietfib, rs_lmer.dietfib)# slightly sign 0.076
anova(lmer.sugar, rs_lmer.sugar)
anova(lmer.prot, rs_lmer.prot)
anova(lmer.fat, rs_lmer.fat)
anova(lmer.fat.cal, rs_lmer.fat.cal)
anova(lmer.addsug, rs_lmer.addsug)



