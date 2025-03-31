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
  left_join(birthControl.df, by="biome_id") %>% 
  filter(!is.na(birthControl))
# add the average shannon idx for a participant to df
shannon.birthControl <- shannon.birthControl %>% 
  group_by(biome_id) %>% 
  mutate(avg_shannon=sum(shannon)/n())

# Collapse by person 
shannon.birthControl.collapsed <- shannon.birthControl %>% 
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

# heatmap of cst assignment and birth control
pheatmap(heatmap.mtx,
         cluster_rows = FALSE,   
         cluster_cols = FALSE,   
         display_numbers = TRUE, 
         number_format = "%.0f",
         main = "CST Assignment v. Birth Control",
         color = brewer.pal(9, "YlGnBu"),
         angle_col=315)

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
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
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
  distinct(biome_id, birthControl, avg_shannon, .keep_all = TRUE)

# Figure: boxplot of birth control and avg shannon diversity
ggplot(shannon.birthControl.avg, aes(x = birthControl, y = avg_shannon)) +
  geom_jitter(aes(color=as.factor(biome_id)), size=1, alpha=0.6) +
  geom_point(aes(color=as.factor(biome_id)), size=1, alpha=0.7) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
  # geom_text(aes(label = as.factor(biome_id)), vjust = -0.5, hjust=1.5, size = 1.5) +
  scale_color_viridis_d(option="D") +
  labs(x = "Contraceptive", y = "Average Shannon Diversity Index", title = "",
       color = "Biome ID") +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1))

# Figure: boxplot of shannon diversity by CST cluster
ggplot(shannon.birthControl, aes(x = CST, y = shannon)) +
  geom_jitter(aes(color=as.factor(birthControl)), size=1, alpha=0.6) +
  geom_point(aes(color=as.factor(birthControl)), size=1, alpha=0.7) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
  # geom_text(aes(label = as.factor(biome_id)), vjust = -0.5, hjust=1.5, size = 1.5) +
  scale_color_viridis_d(option="D") +
  labs(x = "CST", y = "Shannon Index", title = "",
       color = "Birth Control") #+

########################################################
## Corr with DASS data/stress
# dass <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_dass.csv")

## WHEN RERUN - USE THIS INSTEAD 03/22
dass <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/DASS_0503_2024-final_df.csv")
dass <- dass %>%
  rename(biome_id=study_id)

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
dass.participant <- merge(shannon.birthControl, dass.avg, by="biome_id")
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

# Stress: stress score over time
ggplot(dass, aes(x = factor(biome_id), y = stress_score)) +
  geom_boxplot(fill="orchid") +
  # scale_fill_viridis_d(option = "plasma", guide = "none") +
  labs(x = "Participant", y = "Weekly Stress Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

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
  select(-c(X, status))

# add stress severity
dass$stressseverity[dass$stress_score>=0 & dass$stress_score<=14] <- 0
dass$stressseverity[dass$stress_score>=15 & dass$stress_score<=18] <- 1
dass$stressseverity[dass$stress_score>=19 & dass$stress_score<=25] <- 2
dass$stressseverity[dass$stress_score>=26 & dass$stress_score<=33] <- 3
dass$stressseverity[dass$stress_score>=34] <- 4

dass2 <- dass %>% 
  mutate(Timestamp = as.Date(Timestamp)) %>% 
  rename(mood_date = Timestamp) %>% 
  select(biome_id, week, mood_date, stress_score, stressseverity, cisWoman, sport, probiotic, sexuallyActive)

shannon.dass <- shannon.birthControl.mod %>%
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

# add week - join sulogDate# add week - join survey for that given week
# shannon.dass$sampleWeek <- rep(NA, nrow(shannon.dass))
# shannon.dass$sampleWeek[shannon.dass$logDate >= 19276 & shannon.dass$logDate <= 19280] <- 1
# shannon.dass$sampleWeek[shannon.dass$logDate >= 19281 & shannon.dass$logDate <= 19286] <- 2
# shannon.dass$sampleWeek[shannon.dass$logDate >= 19290 & shannon.dass$logDate <= 19293] <- 3
# shannon.dass$sampleWeek[shannon.dass$logDate >= 19296 & shannon.dass$logDate <= 19300] <- 4
# shannon.dass$sampleWeek[shannon.dass$logDate >= 19302 & shannon.dass$logDate <= 19308] <- 5
# shannon.dass$sampleWeek[shannon.dass$logDate >= 19312 & shannon.dass$logDate <= 19315] <- 6
# shannon.dass$sampleWeek[shannon.dass$logDate >= 19319 & shannon.dass$logDate <= 19321] <- 7
# shannon.dass$sampleWeek[shannon.dass$logDate >= 19323 & shannon.dass$logDate <= 19329] <- 8
# shannon.dass$sampleWeek[shannon.dass$logDate >= 19330 & shannon.dass$logDate <= 19336] <- 9
# shannon.dass$sampleWeek[shannon.dass$logDate >= 19339] <- 10
# table(shannon.dass$sampleWeek)

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
  anti_join(dupes)
dim(shannon.dass.filtered)
dim(shannon.dass)

# add back unique
shannon.dass.filtered <- shannon.dass.filtered %>% 
  left_join(dupes.filtered)
dim(shannon.dass.filtered) # 1432   22

# Stress: shannon by stress score
ggplot(shannon.dass.filtered, aes(x = stress_score, y = shannon, col=as.factor(biome_id))) +
  geom_point(fill="orchid") +
  labs(x = "Weekly Stress Score", y = "Shannon") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

length(unique(shannon.dass.filtered$biome_id))
length(!is.na(shannon.dass.filtered$stressseverity))
table(shannon.dass.filtered$stressseverity)

# Stress: shannon by stress severity category
ggplot(shannon.dass.filtered, aes(x = factor(stressseverity), y = shannon, col=as.factor(biome_id))) +
  geom_boxplot(fill = "white", color = "black") + 
  geom_point(position = position_jitter(width = 0.4, height = 0), alpha = 0.7) +
  labs(x = "Stress Level", y = "Shannon Diversity") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("0" = "Normal", "1" = "Mild", "2" = "Moderate", "3" = "Severe", "4" = "Extremely Severe")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

## Mixed effects models
library(lme4)
library(lmerTest)
library(performance)

lmer.full <- lmer(shannon ~ stress_score + (1 | `biome_id`), data=shannon.dass.filtered)
r2(lmer.full)
summary(lmer.full)

lmer.full <- lmer(shannon ~ stress_score + (stress_score || `biome_id`), data=shannon.dass.filtered)
r2(lmer.full)
summary(lmer.full)

## Analysis with birth control
dim(shannon.dass.filtered)
shannon.dass.filtered.na <- shannon.dass.filtered %>%
  drop_na(stress_score, shannon)
length(unique(shannon.dass.filtered.na$biome_id))
dim(shannon.dass.filtered.na)

# Stress: scatterplot stress v. shannon diversity index, birthControl
spline.obj <- shannon.dass.filtered %>%
  drop_na(stress_score, shannon) %>%
  group_by(birthControl) %>%
  summarise(
    model = list(smooth.spline(x=stress_score, y=shannon, spar=NULL)),
    .groups = "drop"
  ) %>%
  mutate(
    x = map(model, ~ as.numeric(.x$x)),
    y = map(model, ~ as.numeric(.x$y))
  ) %>%
  dplyr::select(birthControl, x, y) %>%
  unnest(c(x, y))

ggplot(shannon.dass.filtered, aes(x = stress_score, y = shannon)) +
  geom_jitter(aes(color=as.factor(birthControl)), size=1, alpha=0.6) +
  geom_line(data = spline.obj, aes(x = x, y = y, color = birthControl), linewidth = 1) +
  scale_color_viridis_d(option="D") +
  labs(x = "Weekly Stress Score", y = "Shannon Diversity Index", title = "",
       color = "Birth Control") +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1))

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
  geom_jitter(aes(color=as.factor(CST)), size=1, alpha=0.6) +
  geom_line(data = spline.obj, aes(x = x, y = y, color = CST), linewidth = 1) +
  labs(x = "Weekly Stress Score", y = "Shannon",
       color="CST") +
  theme_minimal() +
  scale_color_viridis_d(option="D") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

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

# Figure: boxplot CST (max) by avg stress score, color birth control
ggplot(dass.participant.collapsed, aes(x = CST_max, y = avg_stress)) +
  geom_jitter(aes(color=as.factor(birthControl)), size=1, alpha=0.6) +
  geom_point(aes(color=as.factor(birthControl)), size=1, alpha=0.7) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
  # geom_text(aes(label = as.factor(biome_id)), vjust = -0.5, hjust=1.5, size = 1.5) +
  scale_color_viridis_d(option="D") +
  labs(x = "CST", y = "Average Stress Score", title = "",
       color = "Birth Control")

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

ggplot(dass.participant, aes(x = avg_stress, y = avg_shannon)) +
  geom_jitter(aes(color=as.factor(CST)), size=1, alpha=0.6) +
  # geom_point(aes(color=as.factor(CST)), size=1, alpha=1) +
  geom_line(data = spline.obj, aes(x = x, y = y, color = CST), linewidth = 1) +
  scale_color_viridis_d(option="D") +
  labs(x = "Average Stress Score", y = "Average Shannon Diversity Index", title = "",
       color = "CST") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

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

spline.obj <- dass.participant %>%
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

ggplot(dass.participant, aes(x = avg_stress, y = avg_shannon)) +
  geom_point(aes(color=as.factor(birthControl)), size=1, alpha=1) +
  geom_line(data = spline.obj, aes(x = x, y = y, color = birthControl), linewidth = 1) +
  scale_color_viridis_d(option="D") +
  labs(x = "Average Stress Score", y = "Average Shannon Diversity Index", title = "",
       color = "Birth Control") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

### Alpha diversity over time, color maybe by birth control, menstruate, or both
dass.participant$logDate <- as.Date(dass.participant$logDate)

# Figure: shannon diversity over time for CSTs
ggplot(dass.participant, aes(x = logDate, y = shannon)) +
  geom_point(aes(color=CST)) +
  geom_smooth(method = "lm", se=FALSE, aes(color = CST)) +
  # geom_point(aes(color=birthControl)) +
  # geom_smooth(method = "lm", se=FALSE, aes(color = birthControl)) +
  labs(
    x = "Days", 
    y = "Shannon Diversity Index",
    title = ""
  ) +
  theme_minimal() +
  scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "1 day") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
  scale_color_viridis_d(option="D")

# Menses: shannon and contraceptive by CST
ggplot(dass.participant, aes(x = birthControl, y = shannon)) +
  geom_jitter(aes(color=as.factor(CST)), size=1, alpha=0.6) +
  geom_point(aes(color=as.factor(CST)), size=1, alpha=0.7) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
  # geom_text(aes(label = as.factor(biome_id)), vjust = -0.5, hjust=1.5, size = 1.5) +
  scale_color_viridis_d(option="D") +
  labs(x = "Contraceptive", y = "Shannon Diversity Index", title = "",
       color = "CST") +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1))

ggplot(dass.participant.collapsed, aes(x = birthControl, y = avg_shannon)) +
  geom_jitter(aes(color=as.factor(CST_max)), size=1, alpha=0.6) +
  geom_point(aes(color=as.factor(CST_max)), size=1, alpha=0.7) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
  # geom_text(aes(label = as.factor(biome_id)), vjust = -0.5, hjust=1.5, size = 1.5) +
  scale_color_viridis_d(option="D") +
  labs(x = "Contraceptive", y = "Average Shannon Diversity Index", title = "",
       color = "CST") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

# Regress shannon diversity on stress + stress*birth control, if interaction terms (2) is significant

dass.participant.collapsed <- dass.participant.collapsed %>% 
  mutate(birthControl_collapsed=ifelse(birthControl=="Systemic Combined (E&P)" | (birthControl=="Systemic P only"), "Systemic", 
                                       birthControl)) %>% 
  mutate(birthControl=as.factor(birthControl))

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
dass <- dass %>% 
  dplyr::select(!biome_id) %>% 
  rename(biome_id=study_id)

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
lmer.full <- lmer(avg_shannon~depression_score + anxiety_score + stress_score + 
                    windDown + mouthDry + noPositiveFeeling + difficultyBreathing +
                    initiative + overreact + trembling + nervous + panicSituation +
                    noLookForward + agitated + difficultyRelax + downhearted +
                    intolerant + closeToPanic + noEnthusiasm + feelWorthless +
                    touchy + awareHeart + scared + lifeMeaningless +
                    (1|`biome_id`), 
                  data = dass.shannon)
r2(lmer.full)

# Full model w time - fixed slopes
lmer.time.full <- lmer(avg_shannon~week+depression_score + anxiety_score + stress_score + 
                    windDown + mouthDry + noPositiveFeeling + difficultyBreathing +
                    initiative + overreact + trembling + nervous + panicSituation +
                    noLookForward + agitated + difficultyRelax + downhearted +
                    intolerant + closeToPanic + noEnthusiasm + feelWorthless +
                    touchy + awareHeart + scared + lifeMeaningless +
                    (1|`biome_id`), 
                  data = dass.shannon)
r2(lmer.time.full)

anova(lmer.full, lmer.time.full)

# Full model w time - fixed slopes
lmer.time2.full <- lmer(avg_shannon~week+I(week^2)+depression_score + anxiety_score + stress_score + 
                         windDown + mouthDry + noPositiveFeeling + difficultyBreathing +
                         initiative + overreact + trembling + nervous + panicSituation +
                         noLookForward + agitated + difficultyRelax + downhearted +
                         intolerant + closeToPanic + noEnthusiasm + feelWorthless +
                         touchy + awareHeart + scared + lifeMeaningless +
                         (1|`biome_id`), 
                       data = dass.shannon)
r2(lmer.time2.full)

anova(lmer.time.full, lmer.time2.full)

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
# participant.dass <- merge(participant.data, dass.avg, by = "biome_id")
# participant.dass <- participant.dass %>% 
#   dplyr::select(!logDate)
# shannon.birthControl.subset<- shannon.birthControl %>% 
#   dplyr::select(shannon, CST, biome_id, avg_shannon, logDate)
# participant.dass <- shannon.birthControl.subset %>% 
#   left_join(participant.dass, by = "biome_id")
# 
# participant.dass$logDate <- as.Date(participant.dass$logDate)
# # person 31 menstruates, must've not submitted any samples
# participant.dass[which(is.na(participant.dass$survey_menstruate)),]$survey_menstruate <- 1
# participant.dass <- participant.dass %>% 
#   filter(biome_id!="31")
# 
# participant.dass$survey_menstruate_text <- as.factor(ifelse(participant.dass$survey_menstruate == 1, "Menstruate", "Do not menstruate"))
# participant.dass$survey_menstruate_text <- factor(participant.dass$survey_menstruate_text, levels=c("Menstruate", "Do not menstruate"))
# 
shannon.dass.menses <- shannon.dass %>% 
  left_join(menses.data.long, by=c("biome_id", "logDate"))

dim(shannon.dass.menses)
names(shannon.dass.menses)

dim(shannon.dass.menses %>% 
      filter(!is.na(menstruate)))
shannon.dass.menses.filter <- shannon.dass.menses %>% 
  filter(!is.na(menstruate))
length(unique(shannon.dass.menses.filter$biome_id))

# Figure: shannon diversity over time on menses v not
ggplot(shannon.dass.menses.filter, aes(x = logDate, y = shannon)) +
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
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

participant.dass.collapsed <- participant.dass %>% 
  distinct(biome_id, avg_shannon, birthControl, avg_stress, survey_menstruate)

# Figure: Does Shannon diversity vary by menstruate/not menstruate?
ggplot(participant.dass.collapsed, aes(x = as.factor(survey_menstruate), y = avg_shannon)) +
  geom_jitter(aes(color=as.factor(birthControl)), size=1, alpha=0.6) +
  geom_point(aes(color=as.factor(birthControl)), size=1, alpha=0.7) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
  # geom_text(aes(label = as.factor(biome_id)), vjust = -0.5, hjust=1.5, size = 1.5) +
  # scale_color_viridis_d(option="D") +
  labs(x = "Menstruate v. No menstruate", y = "Average Shannon Diversity Index", title = "",
       color = "Contraceptive") #+
# theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

t.test(participant.dass.collapsed$avg_shannon ~ participant.dass.collapsed$survey_menstruate)


