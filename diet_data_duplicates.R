## Packages
packs <- c("tidyverse")
lapply(packs, require, character.only = TRUE)

## Load Data
merged_diet_data <- read_csv("/Users/alicezhang/Desktop/microbiome_data/alice_fully_merged_diet_data.csv")
names(merged_diet_data)

# food_other is true
food_other <- merged_diet_data %>%
  filter(food_other == TRUE) %>%
  select(study_id, Date, type)

# extract duplicate check data
check_dupes <- merged_diet_data %>%
  semi_join(food_other, by = c("study_id", "Date", "type"))

dim(check_dupes)

# randomize groups
set.seed(360)
grouped_check_dupes <- check_dupes %>%
  distinct(study_id, Date, type)
shuffled_groups <- grouped_check_dupes %>%
  slice_sample(n = nrow(grouped_check_dupes))
randomized_check_dupes <- check_dupes %>%
  semi_join(shuffled_groups, by = c("study_id", "Date", "type")) %>%
  arrange(match(paste(study_id, Date, type), paste(shuffled_groups$study_id, shuffled_groups$Date, shuffled_groups$type)))

### Save final data output
write.csv(randomized_check_dupes,
          file = "/Users/alicezhang/Desktop/microbiome_data/check_dupes_diet_data.csv",
          row.names = FALSE)
