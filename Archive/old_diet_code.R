
# Note: cut from diet_data_analysis.R

# For a participant on a given day, get cals that are vegetarian v. not vegetarian
# collapse data into foods for that participant by day and the percent of food that is vegetarian

# vegetarian_percentage <- merged_diet_data %>% 
#   mutate(study_id=as.numeric(study_id)) %>% 
#   filter(!is.na(study_id)) %>% 
#   group_by(study_id, Date) %>% 
#   summarise(full_day_cal=sum(caloriesall, na.rm = TRUE),
#             vegetarian_perc=sum(caloriesall[is_vegetarian==TRUE], na.rm=TRUE)) %>% 
#   ungroup()
# 
# View(vegetarian_percentage)

### Save final data output
# write.csv(merged_diet_data,
#           file = "/Users/alicezhang/Desktop/microbiome_data/cleaned_data/cleaned_diet.csv",
#           row.names = FALSE)

###############

## Creating dataset with participant, date, item, and nutritional info
big <- merged_diet_data

### Collapse by Day
small <- big %>%
  group_by(study_id, Date) %>%
  summarise(across(
    c(caloriesall,cholesterolall,carbohydratesall,sodiumall,
      proteinall,transFat,caloriesFromSatFat,caloriesFromFat,
      sugarsall,fatall,addedSugarall,dietaryFiberall), sum
  ))



# Week Variable
small$Timestamp <- as.Date(small$Date, format="%m/%d/%y", tz="UTC")
time_values <- sort(unique(small$Timestamp))
num_time_values <- as.numeric(time_values)

## Figuring out numbers and dates
as.numeric(as.Date("2022-10-12", format="%Y-%m-%d"))


small$week <- rep(NA, nrow(small))
# Week 1: October 12-October 18
small$week[small$Timestamp >= 19272 & small$Timestamp <= 19280] <- 1
small$week[small$Timestamp >= 19281 & small$Timestamp <= 19286] <- 2
small$week[small$Timestamp >= 19287 & small$Timestamp <= 19293] <- 3
small$week[small$Timestamp >= 19294 & small$Timestamp <= 19300] <- 4
small$week[small$Timestamp >= 19301 & small$Timestamp <= 19308] <- 5
small$week[small$Timestamp >= 19309 & small$Timestamp <= 19315] <- 6
small$week[small$Timestamp >= 19316 & small$Timestamp <= 19321] <- 7
small$week[small$Timestamp >= 19322 & small$Timestamp <= 19329] <- 8
small$week[small$Timestamp >= 19330 & small$Timestamp <= 19336] <- 9
small$week[small$Timestamp >= 19337] <- 10

table(small$week)
sum(is.na(small$week)) # 10

small_naweek10 <- small %>% 
  filter(is.na(week))
dim(small_naweek10)
#  1   2   3   4   5   6   7   8   9  10 
# 150 239 216 184 177 120  49 111  89  57 

# Creating a variable with weeks
small$week <- rep(NA, nrow(small))
small$week[small$Date >= "2022-10-08" & small$Date <= "2022-10-14"] <- 1
small$week[small$Date >= "2022-10-15" & small$Date <= "2022-10-21"] <- 2
small$week[small$Date >= "2022-10-22" & small$Date <= "2022-10-27"] <- 3
small$week[small$Date >= "2022-10-28" & small$Date <= "2022-11-03"] <- 4
small$week[small$Date >= "2022-11-04" & small$Date <= "2022-11-10"] <- 5
small$week[small$Date >= "2022-11-11" & small$Date <= "2022-11-17"] <- 6
small$week[small$Date >= "2022-11-18" & small$Date <= "2022-11-24"] <- 7
small$week[small$Date >= "2022-11-25" & small$Date <= "2022-12-01"] <- 8
small$week[small$Date >= "2022-12-02" & small$Date <= "2022-12-08"] <- 9
small$week[small$Date >= "2022-12-09" & small$Date <= "2022-12-15"] <- 10
small$week[small$Date >= "2022-12-15" & small$Date <= "2022-12-21"] <- 11

# 1   2   3   4   5   6   7   8   9  10  11 
# 109 280 190 184 173 120  76  89  91  73   4 
# 13 NA
sum(is.na(small$week)) # 13
small_naweek <- small %>% 
  filter(is.na(week))
small_week11<- small %>% 
  filter(week==11)
dim(small_week11)

#########################

# Collapse by Meal
small2<-big%>%
  group_by(study_id, Date,type)%>%
  summarise(across(c(caloriesall, cholesterolall, carbohydratesall,
                     sodiumall, proteinall, transFat, caloriesFromSatFat,
                     caloriesFromFat, sugarsall, fatall, addedSugarall, dietaryFiberall), sum))


# Contigency table for Study Participant and Meal Type frequency
table(small2$study_id, small2$type)

#### TO DO: beverage and dessert need to be added to the collapse by meal

# Look at beverages
table(merged_diet_data$type)

beverage_data <- merged_diet_data %>% 
  filter(type == "beverage")
dim(beverage_data)
head(beverage_data)

table(beverage_data$name)

# Look at other
other_data <- merged_diet_data %>% 
  filter(type == "other")
dim(other_data)
head(other_data)

table(other_data$name)
# These are drinks? Tequila, 100% Honeycrisp Apple Cider, Apple Cider, Beer, Boba, Brewed Tea, Hot chocolate, latte, Nesquik Chocolate Powder, Pina Colada, Rose wine, sangria, shirley temple, tea, vodka
