########################
#
# ES Abstract Code
# Last updated: 03/31/2025
# data saved from microbiome_vaginal_analysis.R
#
#########################

source("~/Microbiome Thesis/functions.R")
library(tidyverse)

# RELABELED DATA
vaginal.microbial.menses.24 <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/relabeled_data/vaginal.microbial.menses.24.csv")
# vaginal.microbial.menses.24 <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifetyle/vaginal.microbial.menses.24.csv")

########################################################
## Corr with volunteer history data
participant.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 9-Volunteer Medical History.csv", header = TRUE)
shannon.cst.qr.merged.24 <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/relabeled_data/shannon.cst.qr.merged.24.csv", header=TRUE)

# filter non-hormonal & no samples
participant.data <- participant.data %>% 
  filter(birthControl!="Orilissa (Elagolix)")

# select birth control
birthControl.df <- participant.data %>% 
  dplyr::select(biome_id, birthControl)
birthControl.collapsed <- birthControl.df %>% 
  count(birthControl, name="frequency")

# Figure: barplot of frequencies of birth control methods
birthControl.collapsed %>% ggplot(aes(x=birthControl, y=frequency)) +
  geom_bar(stat = "identity", fill="orchid", color="black") +
  theme_minimal() +
  labs(
    x = "Birth Control Methods",
    y = "Frequency",
    title = "Frequency of Birth Control Methods"
  )

# merge df for shannon index with birth control
shannon.birthControl <- shannon.cst.qr.merged.24 %>% 
  left_join(birthControl.df, by="biome_id") 

table(shannon.birthControl$CST_max)

# Collapse by person 
shannon.birthControl.collapsed <- shannon.birthControl %>% 
  group_by(biome_id) %>% 
  mutate(avg_shannon=sum(shannon)/n()) %>% 
  filter(!is.na(birthControl))

table(shannon.birthControl.collapsed$CST_max)

shannon.birthControl.collapsed <- shannon.birthControl.collapsed %>% 
  group_by(biome_id) %>%
  summarise(avg_shannon = sum(shannon)/n(),
            CST_max = first(CST_max),
            birthControl = first(birthControl))

table(shannon.birthControl$CST)
table(shannon.birthControl.collapsed$CST_max)

## Heatmap of birth control frequencies to CST_max
library(pheatmap)
library(RColorBrewer)
table(shannon.birthControl.collapsed$CST_max, shannon.birthControl.collapsed$birthControl)
birthControl.CST.table <- table(shannon.birthControl.collapsed$CST_max, shannon.birthControl.collapsed$birthControl)
heatmap.mtx <- as.matrix(birthControl.CST.table)

# Birth Control: heatmap of cst assignment and birth control
pheatmap(heatmap.mtx,
         cluster_rows = FALSE,   
         cluster_cols = FALSE,   
         display_numbers = TRUE, 
         number_format = "%.0f",
         # main = "CST Assignment v. Birth Control",
         color = brewer.pal(9, "YlGnBu"),
         angle_col=0, #315,
         fontsize=18
         )

# heatmap of birth control to CST (all samples)
table(shannon.birthControl$CST, shannon.birthControl$birthControl)
birthControl.CST.table2 <- table(shannon.birthControl$CST, shannon.birthControl$birthControl)
heatmap.mtx2 <- as.matrix(birthControl.CST.table2)
pheatmap(heatmap.mtx2,
         cluster_rows = FALSE,   
         cluster_cols = FALSE,   
         display_numbers = TRUE, 
         number_format = "%.0f",
         main = "CST v. Birth Control",
         color = brewer.pal(9, "YlGnBu"))

# test normal; significant means, deviates from normal
shapiro.test(shannon.birthControl$shannon) # follows normality
# test if there are significant diffs across categories
kruskal_test <- kruskal.test(shannon ~ birthControl, data = shannon.birthControl)
# kruskal_test # ignore

# Association between shannon diversity and birth control
# library(FSA)
# dunn_result <- dunnTest(shannon ~ birthControl, data = shannon.birthControl, method = "bonferroni")
# print(dunn_result)

# Figure: boxplot of birth control and shannon diversity
ggplot(shannon.birthControl, aes(x = birthControl, y = shannon)) +
  geom_jitter(aes(color=as.factor(biome_id)), size=1, alpha=0.6) +
  # geom_point(aes(color=as.factor(biome_id)), size=1, alpha=0.7) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
  # geom_text(aes(label = as.factor(biome_id)), vjust = -0.5, size = 3) +
  scale_color_viridis_d(option="D") +
  labs(x = "Contraceptive", y = "Shannon Diversity Index", title = "",
       color = "Biome ID") +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),
        legend.position="none")
# 
# ggplot(shannon.birthControl, aes(x = birthControl, y = shannon)) +
#   geom_boxplot() +
#   labs(x = "Contraceptive", y = "Shannon Diversity Index", title = "") +
#   stat_compare_means(method = "kruskal.test", label.y = max(shannon.birthControl$shannon) * 1.5) +
#   stat_compare_means(comparisons = list(c("Local P", "None"),
#                                         c("Local P", "Systemic Combined (E&P)"),
#                                         c("None", "Systemic Combined (E&P)"),
#                                         c("None", "Systemic P only"),
#                                         c("Systemic Combined (E&P)", "Systemic P only")),
#                      method = "wilcox.test") +
#   theme_minimal()

# Average Shannon diversity per person
shannon.birthControl.avg <- shannon.birthControl %>% 
  group_by(biome_id) %>% 
  summarise(avg_shannon = sum(shannon, na.rm=TRUE)/n(),
            birthControl=first(birthControl))
  # distinct(biome_id, birthControl, shannon, .keep_all = TRUE)

shannon.birthControl.avg$birthControl <- factor(shannon.birthControl.avg$birthControl, 
                                               levels=c("None", "Local P", "Systemic P only", "Systemic Combined (E&P)"))

# Birth Control: boxplot of birth control and avg shannon diversity
shannon.birthControl.avg %>% 
  filter(!is.na(birthControl)) %>% 
  ggplot(aes(x = birthControl, y = avg_shannon)) +
  geom_jitter(aes(color=as.factor(biome_id)), size=2, alpha=1, show.legend = FALSE) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
  labs(x = "", y = "Average Shannon Diversity", title = "") +
  theme_minimal() +
  theme(legend.position = "none",
        text=element_text(size=20)) 

# Birth Control: boxplot of shannon diversity by CST cluster
ggplot(shannon.birthControl, aes(x = CST, y = shannon)) +
  geom_jitter(aes(color=as.factor(birthControl)), size=1, alpha=0.6) +
  # geom_point(aes(color=as.factor(birthControl)), size=1, alpha=0.7) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
  # geom_text(aes(label = as.factor(biome_id)), vjust = -0.5, hjust=1.5, size = 1.5) +
  # scale_color_viridis_d(option="D") +
  labs(x = "CST", y = "Shannon Index", title = "",
       color = "Birth Control") +
  theme_minimal()

########################################################
## Corr with DASS data/stress
# dass <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_dass.csv")

## WHEN RERUN - USE THIS INSTEAD 03/22
dass <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/DASS_0503_2024-final_df.csv")
dass <- dass %>%
  rename(biome_id=study_id)

## Stress: stress over time
head(dass)
dass.holder <- dass %>% 
  rename(logDate = Timestamp)
dass.holder <- study_days(dass.holder)
# join counts
dass.holder.count <- dass.holder %>% 
  group_by(biome_id) %>% 
  summarise(n = n())
dass.holder <- dass.holder %>% 
  left_join(dass.holder.count)
# plot
ggplot(dass.holder, aes(x = study_day, y = stress_score)) +
  geom_line(aes(group = biome_id), color = "orchid", alpha = 0.4) +
  geom_point(color = "orchid", alpha = 0.6, size = 2) +
  # overall
  geom_smooth(se = FALSE, color = "red", size = 1.2, method = "loess") +
  labs(
    x = "Study Day",
    y = "Stress Score",
    title = " "
  ) +
  theme_minimal() +
  theme(
    text=element_text(size=15),
    axis.text.x = element_text(angle = 0, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  scale_x_continuous(breaks = seq(0, max(dass.holder$study_day), by = 5))

# average stress score
dass.avg <- dass %>% 
  group_by(biome_id) %>% 
  summarise(
    #avg_depr=sum(depression_score)/n(),
    #avg_anx=sum(anxiety_score)/n(),
    avg_stress=sum(stress_score)/n()
  )

# merge df
# dass.participant <- merge(shannon.birthControl.collapsed, dass.avg, by=study_id)
dass.participant <- shannon.birthControl %>% 
  left_join(dass.avg, by="biome_id")
  # merge(shannon.birthControl, dass.avg, by="biome_id")
dass.participant.collapsed <- merge(shannon.birthControl.collapsed, dass.avg, by="biome_id")

length(unique(dass.participant.collapsed$biome_id))

# participants with stress score data
dim(dass)
length(unique(dass$biome_id))
dass.summary <- dass %>% 
  group_by(biome_id) %>% 
  summarise(count = n())
summary(dass.summary$count)
table(dass.summary$count)

# Stress: Boxplot of stress scores
ggplot(dass, aes(x = factor(biome_id), y = stress_score)) +
  geom_boxplot(fill="orchid") +
  # scale_fill_viridis_d(option = "plasma", guide = "none") +
  labs(x = "Participant", y = "Weekly Stress Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        text=element_text(size=16))

# Stress: stress score over time
ggplot(dass, aes(x = as.factor(Timestamp), y = stress_score, group = as.factor(biome_id), color = as.factor(biome_id))) +
  geom_line() +
  geom_point() + 
  labs(x = "Sample Submission Date", y = "Weekly Stress Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(legend.position = "none")

#### Join microbiome and stress data

shannon.birthControl.mod <- shannon.birthControl %>% 
  mutate(logDate = as.Date(logDate)) %>% 
  dplyr::select(-c(X, status))

# add stress severity
dass$stressseverity[dass$stress_score>=0 & dass$stress_score<=14] <- 0
dass$stressseverity[dass$stress_score>=15 & dass$stress_score<=18] <- 1
dass$stressseverity[dass$stress_score>=19 & dass$stress_score<=25] <- 2
dass$stressseverity[dass$stress_score>=26 & dass$stress_score<=33] <- 3
dass$stressseverity[dass$stress_score>=34] <- 4

dass2 <- dass %>% 
  mutate(Timestamp = as.Date(Timestamp)) %>% 
  rename(mood_date = Timestamp) %>% 
  dplyr::select(biome_id, week, mood_date, stress_score, stressseverity, cisWoman, sport, probiotic, sexuallyActive)

shannon.dass <- shannon.birthControl.mod %>% # issue week variable shows up twice?
  left_join(dass2, by = c("biome_id")) %>%
  mutate(
    diff_days = abs(as.numeric(difftime(logDate, mood_date, units = "days")))
  ) %>%
  group_by(biome_id, logDate) %>%
  # keep closest survey
  filter(diff_days == min(diff_days)) %>%
  ungroup() %>%
  # only keep within 7 day survey
  mutate(stress_score = ifelse(diff_days > 7, NA, stress_score))
dim(shannon.dass)

# filter for dupes
dupes <- shannon.dass %>% 
  filter(duplicated(SampleID) | duplicated(SampleID, fromLast=TRUE))
dim(dupes) # 104 dupes, surveys are same days apart

unique_pairs <- dupes %>% 
  group_by(biome_id, SampleID) %>% 
  count()
dim(unique_pairs)

dupes$sampleWeek <- rep(NA, nrow(dupes))
dupes$sampleWeek[dupes$logDate >= 19276 & dupes$logDate <= 19280] <- 1
dupes$sampleWeek[dupes$logDate >= 19281 & dupes$logDate <= 19286] <- 2
dupes$sampleWeek[dupes$logDate >= 19290 & dupes$logDate <= 19293] <- 3
dupes$sampleWeek[dupes$logDate >= 19296 & dupes$logDate <= 19300] <- 4
dupes$sampleWeek[dupes$logDate >= 19302 & dupes$logDate <= 19308] <- 5
dupes$sampleWeek[dupes$logDate >= 19312 & dupes$logDate <= 19315] <- 6
dupes$sampleWeek[dupes$logDate >= 19319 & dupes$logDate <= 19321] <- 7
dupes$sampleWeek[dupes$logDate >= 19323 & dupes$logDate <= 19329] <- 8
dupes$sampleWeek[dupes$logDate >= 19330 & dupes$logDate <= 19336] <- 9
dupes$sampleWeek[dupes$logDate >= 19339] <- 10
table(dupes$sampleWeek)

dupes.filtered <- dupes %>%
  filter(week==sampleWeek) %>%
  ungroup()
dim(dupes.filtered)

# take out dupes
shannon.dass.filtered <- shannon.dass %>%
  mutate(mood_date = as.Date(mood_date)) %>% 
  anti_join(dupes)
dim(shannon.dass.filtered)
dim(shannon.dass)
names(shannon.dass.filtered)

# add back unique
shannon.dass.filtered <- shannon.dass.filtered %>% 
  left_join(dupes.filtered)
dim(shannon.dass.filtered) # 1432   22
names(shannon.dass.filtered)

colSums(is.na(shannon.dass.filtered)) # 133 stress scores missing

# 04/17
# write.csv(shannon.dass.filtered, file="/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/shannon.dass.csv")


# Stress: shannon by stress score
ggplot(shannon.dass.filtered, aes(x = stress_score, y = shannon, col=as.factor(biome_id))) +
  geom_point(fill="orchid") +
  labs(x = "Weekly Stress Score", y = "Shannon") +
  theme_minimal() +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

length(unique(shannon.dass.filtered$biome_id))
length(!is.na(shannon.dass.filtered$stressseverity))
table(shannon.dass.filtered$stressseverity)

# Stress: shannon by stress severity category
ggplot(shannon.dass.filtered, aes(x = factor(stressseverity), y = shannon),
       col=as.factor(biome_id)) +
  geom_boxplot(fill = "skyblue", outlier.shape = NA, color = "black", alpha=0.1) + 
  geom_point(position = position_jitter(width = 0.4, height = 0), alpha = 0.7,
             aes(color=as.factor(biome_id), alpha = 0.7)) +
  # labs(x = "Stress Level", y = "Shannon Diversity") +
  labs(x = " ", y = "Shannon Diversity") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("0" = "Normal", "1" = "Mild", "2" = "Moderate", "3" = "Severe", "4" = "Extremely Severe")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        text = element_text(size=18))

# Linear regression: average stress and average 

shannon.dass.filtered.collapsed <- shannon.dass.filtered %>% 
  group_by(biome_id) %>% 
  summarise(avg_shannon=sum(shannon, na.rm=TRUE) / n(),
            avg_stress_score = sum(stress_score, na.rm=TRUE) / n()) %>% 
  mutate(stress_severity = case_when(
    avg_stress_score <= 14 ~ "Normal",
    avg_stress_score >= 15 & avg_stress_score <= 18 ~ "Mild",
    avg_stress_score >= 19 & avg_stress_score <= 25 ~ "Moderate",
    avg_stress_score >= 26 & avg_stress_score <= 33 ~ "Severe",
    avg_stress_score >= 34 ~ "Extremely Severe",
    TRUE ~ NA_character_
  ))

# differences in avg shannon for stress severity category
shannon.dass.filtered.collapsed$stress_severity <- as.factor(shannon.dass.filtered.collapsed$stress_severity)
lm.obj <- lm(avg_shannon~stress_severity, data=shannon.dass.filtered.collapsed)
summary(lm.obj)
TukeyHSD(aov(avg_shannon ~ stress_severity, data = shannon.dass.filtered.collapsed))

anova(lm(avg_shannon~stress_severity, data=shannon.dass.filtered.collapsed))
emmeans(lm.obj, pairwise ~ stress_severity, adjust = "tukey")
table(shannon.dass.filtered.collapsed$stress_severity)
shannon.dass.filtered.collapsed$stress_severity <- factor(shannon.dass.filtered.collapsed$stress_severity, levels=c("Normal", "Mild", "Moderate", "Severe", "Extremely Severe"))
# Stress: avg shannon by avg stress severity category
shannon.dass.filtered.collapsed %>% 
  filter(!is.na(stress_severity)) %>% 
  ggplot(aes(x = factor(stress_severity), y = avg_shannon),
       col=as.factor(biome_id)) +
  geom_boxplot(fill = "skyblue", outlier.shape = NA, color = "black", alpha=0.1) + 
  geom_point(position = position_jitter(width = 0.4, height = 0), alpha = 0.7,
             aes(color=as.factor(biome_id), alpha = 0.7)) +
  # labs(x = "Stress Level", y = "Shannon Diversity") +
  labs(x = " ", y = "Average Shannon Diversity") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("Normal" = "Normal", "Mild" = "Mild", "Moderate" = "Moderate", "Severe" = "Severe", "Extremely Severe" = "Extremely Severe")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        text = element_text(size=18))
table(shannon.dass.filtered$stressseverity)
## Mixed effects models - stress and shannon diversity
library(lme4)
library(lmerTest)
library(performance)

# null model
lmer.null <- lmer(shannon ~ (1 | `biome_id`), data=shannon.dass.filtered)
r2(lmer.null)
summary(lmer.null)

lmer.full <- lmer(shannon ~ as.factor(stressseverity) + (1 | `biome_id`), data=shannon.dass.filtered)
r2(lmer.full)
summary(lmer.full)

anova(lmer.null, lmer.full)

lmer.full2 <- lmer(shannon ~ as.factor(stressseverity) + (as.factor(stressseverity) | `biome_id`), data=shannon.dass.filtered)
r2(lmer.full2)
summary(lmer.full2)

# time
shannon.dass.filtered <- study_days(shannon.dass.filtered)
lmer.full.time <- lmer(shannon ~ as.factor(stressseverity) + study_day + (1 | `biome_id`), data=shannon.dass.filtered)
r2(lmer.full.time)
summary(lmer.full.time)

anova(lmer.full.time, lmer.full)

lmer.full.time2 <- lmer(shannon ~ as.factor(stressseverity) + study_day + I(study_day^2) + (1 | `biome_id`), data=shannon.dass.filtered)
r2(lmer.full.time2)
summary(lmer.full.time2)

anova(lmer.full.time, lmer.full.time2)
anova(lmer.null, lmer.full.time)
anova(lmer.full, lmer.full.time)
anova(lmer.full.time, lmer.full.time2)

## Analysis with birth control
dim(shannon.dass.filtered)
shannon.dass.filtered.na <- shannon.dass.filtered %>%
  drop_na(stress_score, shannon)
length(unique(shannon.dass.filtered.na$biome_id))
dim(shannon.dass.filtered.na)

# Stress: scatterplot stress v. shannon diversity index, birthControl
spline.obj <- shannon.dass.filtered %>%
  drop_na(stressseverity, shannon) %>%
  group_by(birthControl) %>%
  summarise(
    model = list(smooth.spline(x=stressseverity, y=shannon, spar=NULL)),
    .groups = "drop"
  ) %>%
  mutate(
    x = map(model, ~ as.numeric(.x$x)),
    y = map(model, ~ as.numeric(.x$y))
  ) %>%
  dplyr::select(birthControl, x, y) %>%
  unnest(c(x, y))

# Stress: spline of shannon diversity and stress # boxplots now
ggplot(shannon.dass.filtered, aes(x = as.factor(stressseverity), y = shannon, fill=birthControl)) +
  # geom_jitter(aes(color=as.factor(birthControl)), size=1, alpha=0.6) +
  geom_jitter(aes(colour = as.factor(biome_id)), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.4, size = 1.2) +
  geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.shape = NA) +
  # geom_boxplot(data = spline.obj, aes(x = x, y = y, color = birthControl), linewidth = 1) +
  # scale_color_viridis_d(option="D") +
  guides(colour = "none") +
  theme_minimal() +
  labs(x = " ", y = "Shannon Diversity", title = "",
       fill = "Birth Control") +
  scale_x_discrete(labels = c("0" = "Normal", "1" = "Mild", "2" = "Moderate", "3" = "Severe", "4" = "Extremely Severe")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        text = element_text(size=24),
        legend.position="right") 
  

# Stress: stress score and shannon by CST with splines
table(shannon.dass.filtered$CST)
spline.obj <- shannon.dass.filtered %>%
  drop_na(stress_score, shannon) %>%
  group_by(CST) %>%
  summarise(
    model = list(smooth.spline(x=stress_score, y=shannon, spar=NULL)),
    .groups = "drop"
  ) %>%
  mutate(
    x = map(model, ~ as.numeric(.x$x)),
    y = map(model, ~ as.numeric(.x$y))
  ) %>%
  dplyr::select(CST, x, y) %>%
  unnest(c(x, y))
ggplot(shannon.dass.filtered, aes(x = stress_score, y = shannon, col=CST)) +
  # geom_point(fill="orchid") +
  geom_jitter(aes(color=as.factor(CST)), size=1, alpha=0.5) +
  geom_line(data = spline.obj, aes(x = x, y = y, color = CST), linewidth = 1) +
  labs(x = "Weekly Stress Score", y = "Shannon",
       color="CST") +
  theme_minimal() +
  scale_color_viridis_d(option="D") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        text=element_text(size=18))

# Stress: stress severity and CST boxplots
ggplot(shannon.dass.filtered, aes(x = as.factor(stressseverity), y = shannon, fill = CST)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.shape = NA) +
  geom_jitter(aes(colour = as.factor(biome_id)), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.4, size = 1.2) +
  # scale_color_viridis_d(option = "C", guide = "none") +
  guides(colour = "none") +
  labs(
    x = " ",
    y = "Shannon Diversity",
    fill = "CST"
  ) +
  theme_minimal() +
  scale_fill_viridis_d(option="D") +
  scale_x_discrete(labels = c("0" = "Normal", "1" = "Mild", "2" = "Moderate", "3" = "Severe", "4" = "Extremely Severe")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        text = element_text(size=18),
        legend.position="right")

ggplot(shannon.dass.filtered, aes(x = CST, y = shannon)) +
  geom_boxplot() + 
  theme_minimal() +
  labs(title = "Shannon Diversity vs CST")

# CST models
lmer.null <- lmer(shannon ~ (1 | biome_id), data = shannon.dass.filtered)
summary(lmer.null)
r2(lmer.null)

lmer.obj <- lmer(shannon ~ as.factor(stressseverity) * CST + (1 | biome_id), data = shannon.dass.filtered)
summary(lmer.obj)
r2(lmer.obj)

lmer.objtime <- lmer(shannon ~ as.factor(stressseverity) * CST + study_day + (1 | biome_id), data = shannon.dass.filtered)
summary(lmer.objtime)
r2(lmer.objtime)

lmer.objtime2 <- lmer(shannon ~ as.factor(stressseverity) * CST + study_day + I(study_day)^2 + (1 | biome_id), data = shannon.dass.filtered)
summary(lmer.objtime2)
r2(lmer.objtime2)

anova(lmer.null, lmer.obj)
anova(lmer.null, lmer.objtime)
anova(lmer.null, lmer.objtime2)
anova(lmer.obj, lmer.objtime)
anova(lmer.objtime, lmer.objtime2) 

lmer.obj2 <- lmer(shannon ~ CST + (1 | biome_id), data = shannon.dass.filtered)
summary(lmer.obj2)
r2(lmer.obj2)

anova(lmer.null, lmer.obj2) # adding CST better than null
anova(lmer.obj, lmer.obj2) # adding interaction significantly improves fit than just cst

lmer.obj2time <- lmer(shannon ~ CST + study_day + (1 | biome_id), data = shannon.dass.filtered)
summary(lmer.obj2time)
r2(lmer.obj2time)

anova(lmer.obj2, lmer.obj2time) # adding linear time sig better

lmer.obj2time2 <- lmer(shannon ~ CST + study_day + I(study_day)^2 + (1 | biome_id), data = shannon.dass.filtered)
summary(lmer.obj2time2)
r2(lmer.obj2time2)

anova(lmer.obj2time, lmer.obj2time2) # adding quad time not

# Stress: stress severity and CST boxplots
ggplot(shannon.dass.filtered, aes(x = as.factor(stressseverity), y = shannon, fill = CST)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.shape = NA) +
  geom_jitter(aes(colour = as.factor(biome_id)), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.4, size = 1.2) +
  # scale_color_viridis_d(option = "C", guide = "none") +
  guides(colour = "none") +
  labs(
    x = " ",
    y = "Shannon Diversity",
    fill = "CST"
  ) +
  theme_minimal() +
  scale_fill_viridis_d(option="D") +
  scale_x_discrete(labels = c("0" = "Normal", "1" = "Mild", "2" = "Moderate", "3" = "Severe", "4" = "Extremely Severe")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        text = element_text(size=18),
        legend.position="right")

## aggregate results

spline.obj <- smooth.spline(x=dass.participant.collapsed$avg_stress, y=dass.participant.collapsed$avg_shannon)
spline.df <- data.frame(avg_stress=spline.obj$x,
                        avg_shannon=spline.obj$y)

# Figure: scatter plot of avg shannon idx by avg stress score, color by birth control
ggplot(dass.participant.collapsed, aes(x = avg_stress, y = avg_shannon)) +
  geom_point(aes(color=as.factor(birthControl)), size=1, alpha=1) +
  geom_line(data=spline.df, aes(x=avg_stress, y=avg_shannon), color="orchid", linewidth=1, alpha=0.7) +
  # geom_text(aes(label = as.factor(biome_id)), vjust = -0.5, hjust=1.5, size = 1.5) +
  scale_color_viridis_d(option="D") +
  labs(x = "Average Stress Score", y = "Average Shannon Index", title = "",
       color = "Birth Control") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

# Stress: boxplot CST (max) by avg stress score, color birth control
ggplot(dass.participant.collapsed, aes(x = CST_max, y = avg_stress)) +
  geom_jitter(aes(color=as.factor(birthControl)), size=1, alpha=0.6) +
  # geom_point(aes(color=as.factor(birthControl)), size=1, alpha=0.7) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
  labs(x = "CST", y = "Average Stress Score", title = "",
       color = "Birth Control") +
  theme_minimal()

# Figure: scatterplot average stress v. avg shannon diversity index, CST_max
spline.obj <- dass.participant %>%
  group_by(CST) %>%
  summarise(
    model = list(smooth.spline(x=avg_stress, y=avg_shannon, spar=1)),
    .groups = "drop"
  ) %>%
  mutate(
    x = map(model, ~ as.numeric(.x$x)),
    y = map(model, ~ as.numeric(.x$y))
    # x = map(model, ~ seq(min(.x$x), max(.x$x), length.out = 100)),
    # y = map2(model, x, ~ predict(.x, .y)$y)
  ) %>%
  dplyr::select(CST, x, y) %>%
  unnest(c(x, y))

# Stress: splines of avg shannon and avg stress with CST
ggplot(dass.participant, aes(x = avg_stress, y = avg_shannon)) +
  geom_jitter(aes(color=as.factor(CST)), size=1, alpha=0.6) +
  # geom_point(aes(color=as.factor(CST)), size=1, alpha=1) +
  geom_line(data = spline.obj, aes(x = x, y = y, color = CST), linewidth = 1) +
  scale_color_viridis_d(option="D") +
  labs(x = "Average Stress Score", y = "Average Shannon Diversity Index", title = "",
       color = "CST") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
  theme_minimal()

# CST max (assigned to most abundant CST from submitted samples)
spline.obj <- dass.participant.collapsed %>%
  group_by(CST_max) %>%
  filter(n_distinct(avg_stress) >= 4) %>%
  summarise(
    model = list(smooth.spline(x=avg_stress, y=avg_shannon, spar=1)),
    .groups = "drop"
  ) %>%
  mutate(
    # smoother representation
    x = map(model, ~ seq(min(.x$x), max(.x$x), length.out = 100)),
    y = map2(model, x, ~ predict(.x, .y)$y)
  ) %>%
  dplyr::select(CST_max, x, y) %>%
  unnest(c(x, y))

ggplot(dass.participant.collapsed, aes(x = avg_stress, y = avg_shannon)) +
  # geom_jitter(aes(color=as.factor(CST_max)), size=1, alpha=0.6) +
  geom_point(aes(color=as.factor(CST_max)), size=1, alpha=1) +
  geom_line(data = spline.obj, aes(x = x, y = y, color = CST_max), linewidth = 1) +
  scale_color_viridis_d(option="D") +
  labs(x = "Average Stress Score", y = "Average Shannon Diversity Index", title = "",
       color = "CST") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

# Average stress v. avg shannon diversity index, birth control
lm.results <- dass.participant %>%
  filter(!is.na(birthControl)) %>% 
  group_by(birthControl) %>%
  summarise(
    model = list(lm(avg_shannon ~ avg_stress, data = cur_data())),
    .groups = "drop"
  ) %>%
  mutate(
    slope = sapply(model, function(m) coef(m)["avg_stress"]),
    p_value = sapply(model, function(m) summary(m)$coefficients["avg_stress", "Pr(>|t|)"])
  ) %>%
  dplyr::select(birthControl, slope, p_value)

dass.participant.avg <- dass.participant %>% 
  group_by(biome_id) %>% 
  summarise(avg_shannon=sum(shannon, na.rm=TRUE)/n(),
            avg_stress=first(avg_stress),
            birthControl=first(birthControl)
            ) %>% 
  filter(!is.na(birthControl)) %>% 
  # impute mean
  mutate(avg_stress = ifelse(is.na(avg_stress), mean(avg_stress, na.rm = TRUE), avg_stress))

spline.obj <- dass.participant.avg %>%
  filter(!is.na(birthControl)) %>% 
  group_by(birthControl) %>%
  filter(n_distinct(avg_stress) >= 4) %>%
  summarise(
    model = list(smooth.spline(x=avg_stress, y=avg_shannon, spar=1)),
    .groups = "drop"
  ) %>%
  mutate(
    x = map(model, ~ seq(min(.x$x), max(.x$x), length.out = 100)),
    y = map2(model, x, ~ predict(.x, .y)$y)
  ) %>%
  dplyr::select(birthControl, x, y) %>%
  unnest(c(x, y))

ggplot(dass.participant.avg, aes(x = avg_stress, y = avg_shannon)) +
  geom_point(aes(color=as.factor(birthControl)), size=1, alpha=1) +
  geom_line(data = spline.obj, aes(x = x, y = y, color = birthControl), linewidth = 1) +
  scale_color_viridis_d(option="D") +
  labs(x = "Average Stress Score", y = "Average Shannon Diversity Index", title = "",
       color = "Birth Control") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

### Alpha diversity over time, color maybe by birth control, menstruate, or both
dass.participant$logDate <- as.Date(dass.participant$logDate)
length(unique(dass.participant$biome_id))

length(unique(shannon.cst.qr.merged.24))

# CST: shannon diversity over time for CSTs
shannon.cst.qr.merged.24 <- study_days(shannon.cst.qr.merged.24)
ggplot(shannon.cst.qr.merged.24, aes(x = as.Date(study_day), y = shannon)) +
  geom_point(aes(color=CST), alpha=0.5) +
  geom_smooth(method = "loess", se=FALSE, aes(color = CST)) +
  # geom_point(aes(color=birthControl)) +
  # geom_smooth(method = "lm", se=FALSE, aes(color = birthControl)) +
  labs(
    x = "Study Day", 
    y = "Shannon Diversity",
    title = ""
  ) +
  theme_minimal() +
  scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "1 day") +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),
        text=element_text(size=18)) +
  scale_color_viridis_d(option="D") +
  scale_x_continuous(breaks = seq(0, max(shannon.cst.qr.merged.24$study_day), by = 5))

# Menses: shannon and contraceptive by CST
ggplot(dass.participant, aes(x = birthControl, y = shannon)) +
  geom_jitter(aes(color=as.factor(CST)), size=1, alpha=0.6) +
  geom_point(aes(color=as.factor(CST)), size=1, alpha=0.7) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
  # geom_text(aes(label = as.factor(biome_id)), vjust = -0.5, hjust=1.5, size = 1.5) +
  scale_color_viridis_d(option="D") +
  labs(x = "Contraceptive", y = "Shannon Diversity Index", title = "",
       color = "CST") +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1)) + 
  theme_minimal()

ggplot(dass.participant.collapsed, aes(x = birthControl, y = avg_shannon)) +
  geom_jitter(aes(color=as.factor(CST_max)), size=1, alpha=0.6) +
  geom_point(aes(color=as.factor(CST_max)), size=1, alpha=0.7) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
  # geom_text(aes(label = as.factor(biome_id)), vjust = -0.5, hjust=1.5, size = 1.5) +
  scale_color_viridis_d(option="D") +
  labs(x = "Contraceptive", y = "Average Shannon Diversity Index", title = "",
       color = "CST") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) + 
  theme_minimal()


# Regress shannon diversity on stress + stress*birth control, if interaction terms (2) is significant

dass.participant.collapsed <- dass.participant.collapsed %>% 
  mutate(birthControl=as.character(birthControl)) %>% 
  mutate(birthControl_collapsed=ifelse(birthControl=="Systemic Combined (E&P)" | (birthControl=="Systemic P only"), "Systemic", 
                                       birthControl)) %>% 
  mutate(birthControl=as.factor(birthControl))

dass.participant.collapsed$birthControl_collapsed

dass.participant.collapsed$birthControl_collapsed <- factor(dass.participant.collapsed$birthControl_collapsed, levels=c("None", "Local P", "Systemic"))

# Regressing models
lm.obj <- lm(avg_shannon ~ avg_stress*birthControl_collapsed, data=dass.participant.collapsed)
summary(lm.obj)

lm.obj2 <- lm(avg_shannon ~ avg_stress+birthControl_collapsed, data=dass.participant.collapsed)
summary(lm.obj2)

anova(lm.obj,lm.obj2)

aov.obj <- aov(dass.participant.collapsed$avg_shannon ~ dass.participant.collapsed$birthControl_collapsed)
summary(aov.obj)

########################################################

# Get weeks for the bacterial data
shannon.birthControl$Timestamp <- as.Date(shannon.birthControl$logDate, format="%Y-%m-%d", tz="UTC")
id_values <- unique(shannon.birthControl$biome_id) 
time_values <- sort(unique(shannon.birthControl$Timestamp))
num_time_values <- as.numeric(time_values)
shannon.birthControl <- shannon.birthControl[order(shannon.birthControl$Timestamp),]

#Remove any values that occurred after 12/16 (end of semester)
dim(shannon.birthControl)
shannon.birthControl <- shannon.birthControl[shannon.birthControl$Timestamp<= 19342, ]
dim(shannon.birthControl)

shannon.birthControl$week <- rep(NA, nrow(shannon.birthControl))
shannon.birthControl$week[shannon.birthControl$Timestamp >= 19276 & shannon.birthControl$Timestamp <= 19280] <- 1
shannon.birthControl$week[shannon.birthControl$Timestamp >= 19281 & shannon.birthControl$Timestamp <= 19286] <- 2
shannon.birthControl$week[shannon.birthControl$Timestamp >= 19290 & shannon.birthControl$Timestamp <= 19293] <- 3
shannon.birthControl$week[shannon.birthControl$Timestamp >= 19296 & shannon.birthControl$Timestamp <= 19300] <- 4
shannon.birthControl$week[shannon.birthControl$Timestamp >= 19302 & shannon.birthControl$Timestamp <= 19308] <- 5
shannon.birthControl$week[shannon.birthControl$Timestamp >= 19312 & shannon.birthControl$Timestamp <= 19315] <- 6
shannon.birthControl$week[shannon.birthControl$Timestamp >= 19319 & shannon.birthControl$Timestamp <= 19321] <- 7
shannon.birthControl$week[shannon.birthControl$Timestamp >= 19323 & shannon.birthControl$Timestamp <= 19329] <- 8
shannon.birthControl$week[shannon.birthControl$Timestamp >= 19330 & shannon.birthControl$Timestamp <= 19336] <- 9
shannon.birthControl$week[shannon.birthControl$Timestamp >= 19339] <- 10
table(shannon.birthControl$week)

shannon.birthControl <- shannon.birthControl %>% 
  dplyr::select(!Timestamp)
# dass <- dass %>% 
#   dplyr::select(!biome_id) %>% 
#   rename(biome_id=study_id)

# Collapse shannon.birthControl by week
shannon.birthControl.collapsed <- shannon.birthControl %>% 
  group_by(biome_id, week) %>% 
  summarise(avg_shannon = sum(shannon)/n(),
            sampleType=first(sampleType))

dass$biome_id <- as.integer(dass$biome_id)
dass.shannon <- dass %>% 
  left_join(shannon.birthControl.collapsed, by=c("biome_id", "week"))

# Mixed effects models
library(lme4)
library(lmerTest)
library(performance)

# Full model - fixed slopes
# lmer.full <- lmer(avg_shannon~depression_score + anxiety_score + stress_score + 
#                     windDown + mouthDry + noPositiveFeeling + difficultyBreathing +
#                     initiative + overreact + trembling + nervous + panicSituation +
#                     noLookForward + agitated + difficultyRelax + downhearted +
#                     intolerant + closeToPanic + noEnthusiasm + feelWorthless +
#                     touchy + awareHeart + scared + lifeMeaningless +
#                     (1|`biome_id`), 
#                   data = dass.shannon)
# r2(lmer.full)

# Full model w time - fixed slopes
# lmer.time.full <- lmer(avg_shannon~week+depression_score + anxiety_score + stress_score + 
#                     windDown + mouthDry + noPositiveFeeling + difficultyBreathing +
#                     initiative + overreact + trembling + nervous + panicSituation +
#                     noLookForward + agitated + difficultyRelax + downhearted +
#                     intolerant + closeToPanic + noEnthusiasm + feelWorthless +
#                     touchy + awareHeart + scared + lifeMeaningless +
#                     (1|`biome_id`), 
#                   data = dass.shannon)
# r2(lmer.time.full)
# 
# anova(lmer.full, lmer.time.full)

# Full model w time - fixed slopes
# lmer.time2.full <- lmer(avg_shannon~week+I(week^2)+depression_score + anxiety_score + stress_score + 
#                          windDown + mouthDry + noPositiveFeeling + difficultyBreathing +
#                          initiative + overreact + trembling + nervous + panicSituation +
#                          noLookForward + agitated + difficultyRelax + downhearted +
#                          intolerant + closeToPanic + noEnthusiasm + feelWorthless +
#                          touchy + awareHeart + scared + lifeMeaningless +
#                          (1|`biome_id`), 
#                        data = dass.shannon)
# r2(lmer.time2.full)

# anova(lmer.time.full, lmer.time2.full)

# full model - random slopes, rnd intercepts - doesn't converge
# rnd.slope.lmer.full <- lmer(avg_shannon~depression_score + (depression_score|`biome_id`) +
#                               anxiety_score + (anxiety_score|`biome_id`) +
#                               stress_score + (stress_score|`biome_id`) +
#                               windDown + (windDown|`biome_id`) +
#                               mouthDry + (mouthDry|`biome_id`) +
#                               noPositiveFeeling + (noPositiveFeeling|`biome_id`) +
#                               difficultyBreathing + (difficultyBreathing|`biome_id`) +
#                               initiative + (initiative|`biome_id`) +
#                               overreact + (overreact|`biome_id`) +
#                               trembling + (trembling|`biome_id`) +
#                               nervous + (nervous|`biome_id`) +
#                               panicSituation + (panicSituation|`biome_id`) +
#                               noLookForward + (noLookForward|`biome_id`) +
#                               agitated + (agitated|`biome_id`) +
#                               difficultyRelax + (difficultyRelax|`biome_id`) +
#                               downhearted + (downhearted|`biome_id`) +
#                               intolerant + (intolerant|`biome_id`) +
#                               closeToPanic + (closeToPanic|`biome_id`) +
#                               noEnthusiasm + (noEnthusiasm|`biome_id`) +
#                               feelWorthless + (feelWorthless|`biome_id`) +
#                               touchy + (touchy|`biome_id`) +
#                               awareHeart + (awareHeart|`biome_id`) +
#                               scared + (scared|`biome_id`) +
#                               lifeMeaningless + (lifeMeaningless|`biome_id`) +
#                               (1|`biome_id`), 
#                             data = dass.shannon)
# r2(rnd.slope.lmer.full) 

########################################################
## Corr with menstruation
menses.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/relabeled_data/cleaned_menstruation_data.csv", header=TRUE)
menses.data$logDate <- as.Date(menses.data$logDate)

# RELABELED DATA
menses.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/relabeled_data/imputed_menstruation_data_3_11.csv")
menses.data <- menses.data %>% 
  rename_with(~gsub("X2022.", "2022.", .), starts_with("X2022.")) %>% 
  rename_with(~gsub("\\.", "-", .))

# Reshape
menses.data.long <- menses.data %>% 
  pivot_longer(cols=starts_with("2022-"), names_to="logDate", values_to="menses_status")
menses.data.long <- menses.data.long %>% 
  mutate(menstruate = ifelse(menses_status %in% c(1,2,3,7,9), "Menstruating", "Not menstruating"),
         logDate = as.Date(logDate))
menses.data.long$logDate <- as.character(menses.data.long$logDate)
shannon.dass.menses <- shannon.cst.qr.merged.24 %>%  #shannon.dass %>% 
  left_join(menses.data.long)

# shannon.cst.qr.merged.24

names(shannon.dass.menses)
names(participant.data)

dim(shannon.dass.menses)
names(shannon.dass.menses)

dim(shannon.dass.menses %>% 
      filter(!is.na(menstruate)))
shannon.dass.menses.filter <- shannon.dass.menses %>% 
  filter(!is.na(menstruate))
length(unique(shannon.dass.menses.filter$biome_id))

# Menses: shannon diversity over time on menses v not
length(unique(shannon.dass.menses.filter$biome_id))

# shannon.dass.menses.filter <- study_days(shannon.dass.menses.filter)

# Menses: shannon diversity over time on menses v not.png
ggplot(shannon.dass.menses.filter, aes(x = as.Date(logDate), y = shannon)) +
  geom_point(aes(color=as.factor(menstruate)), alpha=0.5) +
  geom_smooth(method = "lm", se=FALSE, aes(color = as.factor(menstruate))) +
  labs(
    x = "Days", 
    y = "Shannon Diversity Index",
    title = "",
    color="Menstruating v. Not menstruate"
  ) +
  theme_minimal() +
  scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "1 day") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
        text=element_text(size=18))

# join dass data
menses.data.long$logDate <- as.Date(menses.data.long$logDate)
shannon.dass.menses <- shannon.dass %>% 
  left_join(menses.data.long)

dass.participant.menses <- dass.participant %>% 
  mutate(logDate = as.Date(logDate)) %>% 
  left_join(menses.data.long)

participant.dass.collapsed <- dass.participant.menses %>% 
  distinct(biome_id, shannon, birthControl, avg_stress, menstruate)

# participant.dass.collapsed <- participant.dass.collapsed %>% 
#   group_by(biome_id) %>% 
#   summarise(avg_shannon = sum(shannon, na.rm=TRUE)/n(),
#             birthControl = first(birthControl),
#             avg_stress = first(avg_stress))

# Birth Control: Does Shannon diversity vary by menstruate/not menstruate
participant.dass.collapsed %>% 
  filter(!is.na(menstruate)) %>% 
  ggplot(aes(x = as.factor(menstruate), y = shannon)) +
  geom_jitter(aes(color=as.factor(birthControl)), size=1, alpha=0.6) +
  # geom_point(aes(color=as.factor(birthControl)), size=1, alpha=0.7) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
  theme_minimal() +
  theme(text=element_text(size=18)) +
  labs(x = " ", y = "Shannon Diversity Index", title = "",
       color = "Contraceptive")

# Menses: 8 boxplots, HBC, menses non menses
dass.participant.menses$birthControl <- factor(dass.participant.menses$birthControl, 
                                               levels=c("None", "Local P", "Systemic P only", "Systemic Combined (E&P)"))
dass.participant.menses %>% 
  filter(!is.na(menstruate), !is.na(birthControl)) %>%
  ggplot(aes(x = as.factor(birthControl), y = shannon, fill = as.factor(menstruate))) +
  geom_boxplot(position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(color = as.factor(menstruate)), 
              position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.7), 
              size = 1, alpha = 0.6) +
  theme_minimal() +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) + 
  theme(
    text = element_text(size = 16),
    legend.position = "right",
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5)
  ) +
  labs(
    x = " ",
    y = "Shannon Diversity",
    fill = "Menses Status",
    color = "Menses Status"
  )

# library(emmeans)
# model <- lm(shannon ~ birthControl * menstruate, data = dass.participant.menses)
# emmeans(model, pairwise ~ menstruate | birthControl)
# 
# anova(model)
# emmeans(model, pairwise ~ birthControl) 
# 
# # avg shannon ~ HBC
# lm.obj <- lm(avg_shannon ~ birthControl, data=participant.dass.collapsed)
# summary(lm.obj)

# Birth Control: boxplot of birth control and shannon diversity
participant.dass.collapsed$birthControl <- factor(participant.dass.collapsed$birthControl,
                                                  levels=c("None", "Local P", "Systemic P only",
                                                           "Systemic Combined (E&P)"))
participant.dass.collapsed %>% 
  filter(!is.na(birthControl)) %>% 
  ggplot(aes(x = birthControl, y = shannon)) +
  geom_jitter(aes(color=as.factor(biome_id)), size=1, alpha=0.6) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
  labs(x = " ", y = "Shannon Diversity", title = "",
       color = "Biome ID") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),
        legend.position = "None", text = element_text(size = 18)) 

# avg shannon, avg stress
participant.dass.collapsed.avg <- participant.dass.collapsed %>% 
  group_by(biome_id) %>% 
  summarise(avg_shannon = sum(shannon, na.rm=TRUE) / n(),
            birthControl = first(birthControl), 
            avg_stress_score = first(avg_stress)) %>% 
  mutate(stress_severity = case_when(
    avg_stress_score <= 14 ~ "Normal",
    avg_stress_score >= 15 & avg_stress_score <= 18 ~ "Mild",
    avg_stress_score >= 19 & avg_stress_score <= 25 ~ "Moderate",
    avg_stress_score >= 26 & avg_stress_score <= 33 ~ "Severe",
    avg_stress_score >= 34 ~ "Extremely Severe",
    TRUE ~ NA_character_
  ))

participant.dass.collapsed.avg$stress_severity <- factor(participant.dass.collapsed.avg$stress_severity, 
                                              levels=c("Normal", "Mild", "Moderate",
                                                       "Severe", "Extremely Severe"))

lm.obj <- lm(avg_shannon ~ birthControl+stress_severity, data = participant.dass.collapsed.avg)
summary(lm.obj)

lm.obj1 <- lm(avg_shannon ~ birthControl*stress_severity, data = participant.dass.collapsed.avg)
summary(lm.obj1)

anova(lm.obj, lm.obj1)

## Mixed effects models - stress and shannon diversity
library(lme4)
library(lmerTest)
library(performance)

dass.participant.menses.clean <- dass.participant.menses %>% 
  filter(!is.na(birthControl))

lmer.null <- lmer(shannon ~  (1 | `biome_id`), data=dass.participant.menses.clean)
r2(lmer.null)

lmer.full <- lmer(shannon ~ birthControl + (1 | `biome_id`), data=dass.participant.menses.clean)
r2(lmer.full)
summary(lmer.full)

anova(lmer.null, lmer.full)

# add time
dass.participant.menses.clean <- study_days(dass.participant.menses.clean)
lmer.full.time <- lmer(shannon ~ birthControl + study_day + (1 | `biome_id`), data=dass.participant.menses.clean)
r2(lmer.full.time)
summary(lmer.full.time)

anova(lmer.full, lmer.full.time)

lmer.full.time2 <- lmer(shannon ~ birthControl + study_day + I(study_day^2) + (1 | `biome_id`), data=dass.participant.menses.clean)
r2(lmer.full.time2)
summary(lmer.full.time2)

anova(lmer.full.time, lmer.full.time2)

# avg shannon ~ HBC
lm.obj <- lm(avg_shannon ~ birthControl, data=participant.dass.collapsed)
summary(lm.obj)

# interaction btw shannon stress HBC
# participant.dass.collapsed.nomenses <- 
# participant.dass.collapsed %>% 
#   group_by(biome_id) %>% 
#   summarise(avg_shannon=sum(avg_shannon, na.rm=TRUE)/n(),
#             avg_stress=sum(avg_stress, na.rm=TRUE)/n(),
#             birthControl=first(birthControl))
# shannon.dass.menses <- shannon.dass.filtered.collapsed %>% 
#   left_join(participant.dass.collapsed.nomenses) %>% 
#   mutate(stress_severity = case_when(
#     avg_stress_score <= 14 ~ "Normal",
#     avg_stress_score >= 15 & avg_stress_score <= 18 ~ "Mild",
#     avg_stress_score >= 19 & avg_stress_score <= 25 ~ "Moderate",
#     avg_stress_score >= 26 & avg_stress_score <= 33 ~ "Severe",
#     avg_stress_score >= 34 ~ "Extremely Severe",
#     TRUE ~ NA_character_
#   ))
# shannon.dass.menses$stress_severity <- factor(shannon.dass.menses$stress_severity, 
#                                               levels=c("Normal", "Mild", "Moderate",
#                                                        "Severe", "Extremely Severe"))
shannon.dass.menses <- participant.dass.collapsed.avg
lm.obj <- lm(avg_shannon~stress_severity+birthControl, data=shannon.dass.menses)
summary(lm.obj)

lm.obj2 <- lm(avg_shannon~stress_severity*birthControl, data=shannon.dass.menses)
summary(lm.obj2)

# null model
shannon.dass.filtered$birthControl <- factor(shannon.dass.filtered$birthControl, levels=c("None", "Local P", "Systemic P only", "Systemic Combined (E&P)"))

lmer.null <- lmer(shannon ~ (1|`biome_id`), data=shannon.dass.filtered)
r2(lmer.null)
summary(lmer.null)

# mixed effects
# shannon.dass.filtered$birthControl <- factor(shannon.dass.filtered$birthControl, levels=c("None", "Local P", "Systemic P only", "Systemic Combined (E&P)"))
lmer.obj <- lmer(shannon~as.factor(stressseverity) + birthControl + (1|`biome_id`), data=shannon.dass.filtered)
r2(lmer.obj)
summary(lmer.obj)

anova(lmer.null, lmer.obj)

# interaction with the stress severity and HBC
lmer.obj2 <- lmer(shannon~as.factor(stressseverity) * birthControl + (1|`biome_id`), data=shannon.dass.filtered)
r2(lmer.obj2)
summary(lmer.obj2)

anova(lmer.obj, lmer.obj2)

lmer.obj3 <- lmer(shannon~as.factor(stressseverity) + (1|`biome_id`), data=shannon.dass.filtered)
r2(lmer.obj3)
summary(lmer.obj3)

lmer.obj.null <- lmer(shannon~(1|`biome_id`), data=shannon.dass.filtered)
r2(lmer.obj.null)
summary(lmer.obj.null)

anova(lmer.obj.null, lmer.obj3)

lmer.obj4 <- lmer(shannon~as.factor(birthControl) + (1|`biome_id`), data=shannon.dass.filtered)
r2(lmer.obj4)
summary(lmer.obj4)

# 04/02 uncomment to save
# write.csv(shannon.dass.filtered, "/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/dass.participant.csv")

##########################################################################################



