########################
#
# ES Abstract Code
# Last updated: 03/18/2025
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
  select(biome_id, birthControl)
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
  summarise()

table(shannon.birthControl$CST)
table(shannon.birthControl.collapsed$CST_max)

## Heatmap of birth control frequencies to CST_max
library(pheatmap)
library(RColorBrewer)
table(shannon.birthControl.collapsed$CST_max, shannon.birthControl.collapsed$birthControl)
birthControl.CST.table <- table(shannon.birthControl.collapsed$CST_max, shannon.birthControl.collapsed$birthControl)
heatmap.mtx <- as.matrix(birthControl.CST.table)

pheatmap(heatmap.mtx,
         cluster_rows = FALSE,   
         cluster_cols = FALSE,   
         display_numbers = TRUE, 
         number_format = "%.0f",
         main = "CST Assignment v. Birth Control",
         color = brewer.pal(9, "YlGnBu"))
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
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

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
dass <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/DASS-21-Cleaning - CLEAN.csv")

# ------------------------------------------------------------------------------
# calculate depression, anxiety, and stress scores 
# ------------------------------------------------------------------------------
dass$Timestamp <- as.Date(dass$Timestamp, format="%m/%d/%y", tz="UTC")
id_values <- unique(dass$study_id) 
time_values <- sort(unique(dass$Timestamp))
num_time_values <- as.numeric(time_values)
dass <- dass[order(dass$Timestamp),]

#Remove any values that occurred after 12/16 (end of semester)
dass <- dass[dass$Timestamp<= 19342, ]

dass$week <- rep(NA, nrow(dass))
dass$week[dass$Timestamp >= 19276 & dass$Timestamp <= 19280] <- 1
dass$week[dass$Timestamp >= 19281 & dass$Timestamp <= 19286] <- 2
dass$week[dass$Timestamp >= 19290 & dass$Timestamp <= 19293] <- 3
dass$week[dass$Timestamp >= 19296 & dass$Timestamp <= 19300] <- 4
dass$week[dass$Timestamp >= 19302 & dass$Timestamp <= 19308] <- 5
dass$week[dass$Timestamp >= 19312 & dass$Timestamp <= 19315] <- 6
dass$week[dass$Timestamp >= 19319 & dass$Timestamp <= 19321] <- 7
dass$week[dass$Timestamp >= 19323 & dass$Timestamp <= 19329] <- 8
dass$week[dass$Timestamp >= 19330 & dass$Timestamp <= 19336] <- 9
dass$week[dass$Timestamp >= 19339] <- 10

dass <- dass %>% 
  filter(!is.na(week))

dass$depression_score <- rep(NA, nrow(dass))
dass$anxiety_score <- rep(NA, nrow(dass))
dass$stress_score <- rep(NA, nrow(dass))
id_values <- unique(dass$study_id) 
for(id in id_values){
  for(week in 1:10){
    dass$depression_score[dass$study_id==id & dass$week==week] <- 2*sum(dass$noPositiveFeeling[dass$study_id==id & dass$week==week] + dass$initiative[dass$study_id==id & dass$week==week] + dass$noLookForward[dass$study_id==id & dass$week==week] + 
                                                                          dass$downhearted[dass$study_id==id & dass$week==week] + dass$noEnthusiasm[dass$study_id==id & dass$week==week] + dass$feelWorthless[dass$study_id==id & dass$week==week] + 
                                                                          dass$lifeMeaningless[dass$study_id==id & dass$week==week], na.rm=TRUE)
    
    dass$anxiety_score[dass$study_id==id & dass$week==week] <- 2*sum(dass$mouthDry[dass$study_id==id & dass$week==week] + dass$difficultyBreathing[dass$study_id==id & dass$week==week] + dass$trembling[dass$study_id==id & dass$week==week] + 
                                                                       dass$panicSituation[dass$study_id==id & dass$week==week] + dass$closeToPanic[dass$study_id==id & dass$week==week] + 
                                                                       dass$awareHeart[dass$study_id==id & dass$week==week] + dass$scared[dass$study_id==id & dass$week==week], na.rm=TRUE)
    
    dass$stress_score[dass$study_id==id & dass$week==week] <- 2*sum(dass$windDown[dass$study_id==id  & dass$week==week] + dass$overreact[dass$study_id==id & dass$week==week] + dass$nervous[dass$study_id==id & dass$week==week] + 
                                                                      dass$agitated[dass$study_id==id & dass$week==week] + dass$difficultyRelax[dass$study_id==id & dass$week==week] + 
                                                                      dass$intolerant[dass$study_id==id & dass$week==week] + dass$touchy[dass$study_id==id & dass$week==week], na.rm=TRUE)
  }
}

#creation of variables to categorize scores into level of severity 
dass$depressionseverity <- rep(NA, nrow(dass))
dass$anxietyseverity <- rep(NA, nrow(dass))
dass$stressseverity<- rep(NA, nrow(dass))

dass$depressionseverity[dass$depression_score>=0 & dass$depression_score<=9] <- 0 
dass$depressionseverity[dass$depression_score>=10 & dass$depression_score<=13] <- 1 
dass$depressionseverity[dass$depression_score>=14 & dass$depression_score<=20] <- 2
dass$depressionseverity[dass$depression_score>=21 & dass$depression_score<=27] <- 3
dass$depressionseverity[dass$depression_score>=28] <- 4

dass$anxietyseverity[dass$anxiety_score>=0 & dass$anxiety_score<=7] <- 0
dass$anxietyseverity[dass$anxiety_score>=8 & dass$anxiety_score<=9] <- 1
dass$anxietyseverity[dass$anxiety_score>=10 & dass$anxiety_score<=14] <- 2
dass$anxietyseverity[dass$anxiety_score>=15 & dass$anxiety_score<=19] <- 3
dass$anxietyseverity[dass$anxiety_score>=20] <- 4

dass$stressseverity[dass$stress_score>=0 & dass$stress_score<=14] <- 0
dass$stressseverity[dass$stress_score>=15 & dass$stress_score<=18] <- 1
dass$stressseverity[dass$stress_score>=19 & dass$stress_score<=25] <- 2
dass$stressseverity[dass$stress_score>=26 & dass$stress_score<=33] <- 3
dass$stressseverity[dass$stress_score>=34] <- 4

# average stress score
dass.avg <- dass %>% 
  group_by(study_id) %>% 
  summarise(
    avg_depr=sum(depression_score)/n(),
    avg_anx=sum(anxiety_score)/n(),
    avg_stress=sum(stress_score)/n()
  ) %>% 
  rename(biome_id=study_id)

# merge df
# dass.participant <- merge(shannon.birthControl.collapsed, dass.avg, by=study_id)
dass.participant <- merge(shannon.birthControl, dass.avg, by="biome_id")
dass.participant.collapsed <- merge(shannon.birthControl.collapsed, dass.avg, by="biome_id")

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
  select(CST, x, y) %>%
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
    # x = map(model, ~ as.numeric(.x$x)),
    # y = map(model, ~ as.numeric(.x$y))
    # smoother representation
    x = map(model, ~ seq(min(.x$x), max(.x$x), length.out = 100)),
    y = map2(model, x, ~ predict(.x, .y)$y)
  ) %>%
  select(CST_max, x, y) %>%
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
  select(birthControl, slope, p_value)

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
  select(birthControl, x, y) %>%
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

ggplot(dass.participant, aes(x = birthControl, y = shannon)) +
  geom_jitter(aes(color=as.factor(CST)), size=1, alpha=0.6) +
  geom_point(aes(color=as.factor(CST)), size=1, alpha=0.7) +
  geom_boxplot(fill="skyblue", outlier.shape = NA, alpha = 0.1) +
  # geom_text(aes(label = as.factor(biome_id)), vjust = -0.5, hjust=1.5, size = 1.5) +
  scale_color_viridis_d(option="D") +
  labs(x = "Contraceptive", y = "Shannon Diversity Index", title = "",
       color = "CST") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

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
  select(!Timestamp)
dass <- dass %>% 
  select(!biome_id) %>% 
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
participant.dass <- merge(participant.data, dass.avg, by = "biome_id")
participant.dass <- participant.dass %>% 
  select(!logDate)
shannon.birthControl.subset<- shannon.birthControl %>% 
  select(shannon, CST, biome_id, avg_shannon, logDate)
participant.dass <- shannon.birthControl.subset %>% 
  left_join(participant.dass, by = "biome_id")

participant.dass$logDate <- as.Date(participant.dass$logDate)
# person 31 menstruates, must've not submitted any samples
participant.dass[which(is.na(participant.dass$survey_menstruate)),]$survey_menstruate <- 1
participant.dass <- participant.dass %>% 
  filter(biome_id!="31")

ggplot(participant.dass, aes(x = logDate, y = shannon)) +
  geom_point(aes(color=as.factor(survey_menstruate))) +
  geom_smooth(method = "lm", se=FALSE, aes(color = as.factor(survey_menstruate))) +
  # geom_point(aes(color=birthControl)) +
  # geom_smooth(method = "lm", se=FALSE, aes(color = birthControl)) +
  labs(
    x = "Days", 
    y = "Shannon Diversity Index",
    title = "",
    color="Menstruate v. No menstruate"
  ) +
  theme_minimal() +
  scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "1 day") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

participant.dass.collapsed <- participant.dass %>% 
  distinct(biome_id, avg_shannon, birthControl, avg_stress, survey_menstruate)

# Does Shannon diversity vary by menstruate/not menstruate?
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


