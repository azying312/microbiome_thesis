########################
#
# Clean DASS Data
# Last updated: 03/21/2025
# Cleaning code from: https://github.com/jyontika/Tetel-Lab/blob/main/DASS/DASS_Cleaning.R
#
#########################

source("~/Microbiome Thesis/functions.R")

library(tidyverse)

dass <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/DASS.csv")
id_mapping <- read.csv("/Volumes/T7/microbiome_data/original_data/Original Study Mapping - Sheet3.csv", header = TRUE)

dass$Timestamp <- as.Date(dass$Timestamp, format="%m/%d/%y", tz="UTC")
id_values <- unique(dass$study_id) 
time_values <- sort(unique(dass$Timestamp))
num_time_values <- as.numeric(time_values)
dass <- dass[order(dass$Timestamp),]

# Remove any values that occurred after 12/16 (end of semester)
dass <- dass[dass$Timestamp<= 19342, ]

# Remove people with unidentifiable study ID
dass <- dass[dass$study_id!="#N/A", ] 
dass$study_id <- as.factor(dass$study_id)

# Weeks variable
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

table(dass$week)

### Map uid to study id
# dass <- dass %>% 
#   rename(biome_id=Your.Biome.Health.App.ID)

# Duplicated surveys
total_num_dups <- sum(duplicated(dass[c("study_id", "week")]))
print(total_num_dups)
list_dups <- dass[duplicated(dass[c("study_id", "week")]), c("study_id", "week")]
print(list_dups)
df_dups <- as.data.frame(list_dups)

# Randomly choose from one of the duplicated rows
set.seed(489)
for (i in 1:nrow(df_dups)) {
  study_id <- df_dups$study_id[i]
  week <- df_dups$week[i]
  
  duplicate_rows <- dass[dass$study_id == study_id & dass$week == week, ]
  random_row_index <- sample(1:2, 1) # Choose random row
  
  random_row <- data.frame(study_id = study_id, week = week) 
  
  duplicate_rows <- dass[dass$study_id == study_id & dass$week == week, ]
  random_row_index <- sample(1:2, 1) # Choose random row
  
  random_row <- duplicate_rows[random_row_index, ] # Select random row
  
  # Convert "biome_id" to character to avoid data type issues
  random_row$biome_id <- as.character(random_row$biome_id)
  
  # Redo total column to be sum of the random values (excluding biome_id and study_id)
  
  random_row$TOTAL <- sum(random_row[, 4:24], na.rm = TRUE)
  
  # Keep only the random rows, remove the others
  dass <- bind_rows(dass[!(dass$study_id == study_id & dass$week == week), ], random_row)
}

# ------------------------------------------------------------------------------
# address internal missing values
# ------------------------------------------------------------------------------
# Locate missing values across columns
sapply(dass, function(x) sum(is.na(x))) # 8 missing values total

# Creating a vector of strings for question column names
question_names <- c("windDown", "mouthDry", "noPositiveFeeling",
                    "difficultyBreathing", "initiative",
                    "overreact", "trembling", "nervous", 
                    "panicSituation", "noLookForward", 
                    "agitated", "difficultyRelax",
                    "downhearted", "intolerant", "closeToPanic", 
                    "noEnthusiasm","feelWorthless", "touchy", 
                    "awareHeart", "scared", "lifeMeaningless")

#Create new column in df that sums num. of NAs across each row
dass$num_NA_pre <- rep(NA, nrow(dass))
sum_NA <- function(x) sum(is.na(x))
dass$num_NA_pre <- apply(dass[, question_names], 1, sum_NA)
dass[dass$num_NA_pre > 0, c("week","study_id",question_names)]

# NA imputation -----------------------------------------------------------

## random imputation because we assume MCAR

stress_vars <- c("windDown", "overreact", "nervous", "agitated", "difficultyRelax", "intolerant", "touchy")
anxiety_vars <- c("mouthDry", "difficultyBreathing", "trembling", "panicSituation", "closeToPanic", "awareHeart", "scared")
depression_vars <- c("noPositiveFeeling", "initiative", "noLookForward", "downhearted", "noEnthusiasm", "feelWorthless", "lifeMeaningless")

set.seed(7)

#missing values
missing_vars <- c("overreact", "difficultyBreathing", "trembling", "noPositiveFeeling", "agitated", 
                  "difficultyRelax", "scared")

#impute NA value for participant 46 -- week 4, overreact, stress
ID_46_impute <- dass[dass$week == 4 & dass$study_id == 46, stress_vars[c(1, 3:7)]]
dass[dass$week == 4 & dass$study_id == 46, stress_vars[2]] <- sample(ID_46_impute, 1)

#impute NA values for participant 57, two values for this participant
##week 6, difficulty breathing, anxiety
ID_57_impute <- dass[dass$week == 6 & dass$study_id == 57, anxiety_vars[c(1, 3:7)]]
dass[dass$week == 6 & dass$study_id == 57, anxiety_vars[2]] <- sample(ID_57_impute, 1)

##week 6, trembling, anxiety
ID_57_impute2 <- dass[dass$week == 6 & dass$study_id == 57, anxiety_vars[c(1:2, 4:7)]]
dass[dass$week == 6 & dass$study_id == 57, anxiety_vars[3]] <- sample(ID_57_impute2, 1)

#impute NA value for participant 44 -- week 7, no positive feeling, depression
ID_44_impute <- dass[dass$week == 7 & dass$study_id == 44, depression_vars[c(2:7)]]
dass[dass$week == 7 & dass$study_id == 44, depression_vars[1]] <- sample(ID_44_impute, 1)

#impute NA value for participant 32 -- week 5, agitated, stress
ID_32_impute <- dass[dass$week == 5 & dass$study_id == 32, stress_vars[c(1:3, 5:7)]]
dass[dass$week == 5 & dass$study_id == 32, stress_vars[4]] <- sample(ID_32_impute, 1)

#impute NA value for participant 61 -- week 9, difficulty relax, stress
ID_61_impute <- dass[dass$week == 9 & dass$study_id == 61, stress_vars[c(1:4, 6:7)]]
dass[dass$week == 9 & dass$study_id == 61, stress_vars[5]] <- sample(ID_61_impute, 1)

#impute NA value for participant 18 -- week 5, scared, anxiety
ID_18_impute <- dass[dass$week == 5 & dass$study_id == 18, anxiety_vars[c(1:6)]]
dass[dass$week == 5 & dass$study_id == 18, anxiety_vars[7]] <- sample(ID_18_impute, 1)

#impute NA value for participant 49 -- week 10, scared, anxiety
ID_49_impute <- dass[dass$week == 10 & dass$study_id == 49, anxiety_vars[c(1:6)]]
dass[dass$week == 10 & dass$study_id == 49, anxiety_vars[7]] <- sample(ID_49_impute, 1)

# no more NA values
dass$num_NA_post <- apply(dass[, question_names], 1, sum_NA)
dass[dass$num_NA_post > 0, c("week","study_id",question_names)] 

# ------------------------------------------------------------------------------
# calculate depression, anxiety, and stress scores 
# ------------------------------------------------------------------------------

dass$depression_score <- rep(NA, nrow(dass))
dass$anxiety_score <- rep(NA, nrow(dass))
dass$stress_score <- rep(NA, nrow(dass))

for(id in id_values){
  for(week in 1:10){
    dass$depression_score[dass$study_id==id & dass$week==week] <- 2*sum(dass$noPositiveFeeling[dass$study_id==id & dass$week==week] + dass$initiative[dass$study_id==id & dass$week==week] + dass$noLookForward[dass$study_id==id & dass$week==week] + 
                                                                          dass$downhearted[dass$study_id==id & dass$week==week] + dass$noEnthusiasm[dass$study_id==id & dass$week==week] + dass$feelWorthless[dass$study_id==id & dass$week==week] + 
                                                                          dass$lifeMeaningless[dass$study_id==id & dass$week==week])
    
    dass$anxiety_score[dass$study_id==id & dass$week==week] <- 2*sum(dass$mouthDry[dass$study_id==id & dass$week==week] + dass$difficultyBreathing[dass$study_id==id & dass$week==week] + dass$trembling[dass$study_id==id & dass$week==week] + 
                                                                       dass$panicSituation[dass$study_id==id & dass$week==week] + dass$closeToPanic[dass$study_id==id & dass$week==week] + 
                                                                       dass$awareHeart[dass$study_id==id & dass$week==week] + dass$scared[dass$study_id==id & dass$week==week])
    
    dass$stress_score[dass$study_id==id & dass$week==week] <- 2*sum(dass$windDown[dass$study_id==id  & dass$week==week] + dass$overreact[dass$study_id==id & dass$week==week] + dass$nervous[dass$study_id==id & dass$week==week] + 
                                                                      dass$agitated[dass$study_id==id & dass$week==week] + dass$difficultyRelax[dass$study_id==id & dass$week==week] + 
                                                                      dass$intolerant[dass$study_id==id & dass$week==week] + dass$touchy[dass$study_id==id & dass$week==week])
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

dass <- dass %>% 
  select(!biome_id) %>% 
  rename(biome_id=study_id)

# ------------------------------------------------------------------------------
# calculate depression, anxiety, and stress scores 
# ------------------------------------------------------------------------------

dass$depression_score <- rep(NA, nrow(dass))
dass$anxiety_score <- rep(NA, nrow(dass))
dass$stress_score <- rep(NA, nrow(dass))

for(id in id_values){
  for(week in 1:10){
    dass$depression_score[dass$study_id==id & dass$week==week] <- 2*sum(dass$noPositiveFeeling[dass$study_id==id & dass$week==week] + dass$initiative[dass$study_id==id & dass$week==week] + dass$noLookForward[dass$study_id==id & dass$week==week] + 
                                                                          dass$downhearted[dass$study_id==id & dass$week==week] + dass$noEnthusiasm[dass$study_id==id & dass$week==week] + dass$feelWorthless[dass$study_id==id & dass$week==week] + 
                                                                          dass$lifeMeaningless[dass$study_id==id & dass$week==week])
    
    dass$anxiety_score[dass$study_id==id & dass$week==week] <- 2*sum(dass$mouthDry[dass$study_id==id & dass$week==week] + dass$difficultyBreathing[dass$study_id==id & dass$week==week] + dass$trembling[dass$study_id==id & dass$week==week] + 
                                                                       dass$panicSituation[dass$study_id==id & dass$week==week] + dass$closeToPanic[dass$study_id==id & dass$week==week] + 
                                                                       dass$awareHeart[dass$study_id==id & dass$week==week] + dass$scared[dass$study_id==id & dass$week==week])
    
    dass$stress_score[dass$study_id==id & dass$week==week] <- 2*sum(dass$windDown[dass$study_id==id  & dass$week==week] + dass$overreact[dass$study_id==id & dass$week==week] + dass$nervous[dass$study_id==id & dass$week==week] + 
                                                                      dass$agitated[dass$study_id==id & dass$week==week] + dass$difficultyRelax[dass$study_id==id & dass$week==week] + 
                                                                      dass$intolerant[dass$study_id==id & dass$week==week] + dass$touchy[dass$study_id==id & dass$week==week])
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

dim(dass)
head(dass)

### Save final data output
write.csv(dass,
          file = "/Volumes/T7/microbiome_data/cleaned_data/cleaned_dass.csv",
          row.names = FALSE)


