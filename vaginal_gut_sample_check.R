library(tidyverse)
library(phyloseq)

source("~/Microbiome Thesis/functions.R")

bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/02_11/bacteria_intermediary2.rds")
metadata.22 <- sample_data(bacterial.data)

vaginal.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/vaginal_cleaned_max_taxa.rds")
vaginal.microbial.menses.24 <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/vaginal.microbial.menses.24.csv")

gut.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/fecal_bacteria_filteredv2.rds")
gut.microbial.menses.24 <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/gut.microbial.menses.24.csv")

#### Participant 1
# 2022-12-01
bacterial_plotting(vaginal.data, 1, "2022-12-01")

#### Participant 3
# 2022-12-05
bacterial_plotting(vaginal.data, 3, "2022-12-05") # S9620_V3V5_S105 - likely fecal
bacterial_plotting(gut.data, 3, "2022-12-10") # F1849_V3V5_S1633 - likely vaginal

## Plot gut
gut.data.3 <- prune_samples(sample_data(gut.data)$biome_id == 3, gut.data)
bacterial_plotting2(gut.data.3)

## Plot vaginal
vaginal.data.3 <- prune_samples(sample_data(vaginal.data)$biome_id == 3, vaginal.data)
bacterial_plotting2(vaginal.data.3)

# Swap swabs
sample_data(bacterial.data)["S9620_V3V5_S105", "sampleType"] <- "fecal"
sample_data(bacterial.data)["F1849_V3V5_S1633", "sampleType"] <- "vaginal"

#### Participant 4
# 2022-10-25
bacterial_plotting(vaginal.data, 4, c("2022-10-26", "2022-11-03")) # S8993_V3V5_S1237, S951_V3V5_S1184
bacterial_plotting(gut.data, 4, "2022-12-10") # 

## Plot gut
gut.data.4 <- prune_samples(sample_data(gut.data)$biome_id == 4, gut.data)
bacterial_plotting2(gut.data.4)

## Plot vaginal
vaginal.data.4 <- prune_samples(sample_data(vaginal.data)$biome_id == 4, vaginal.data)
bacterial_plotting2(vaginal.data.4)

# Swap swabs
sample_data(bacterial.data)["S8993_V3V5_S1237", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S951_V3V5_S1184", "sampleType"] <- "fecal"

#### Participant 5

## Plot gut
gut.data.5 <- prune_samples(sample_data(gut.data)$biome_id == 5, gut.data)
bacterial_plotting2(gut.data.5)

## Plot vaginal
vaginal.data.5 <- prune_samples(sample_data(vaginal.data)$biome_id == 5, vaginal.data)
bacterial_plotting2(vaginal.data.5)


#### Participant 6

## Plot gut
gut.data.6 <- prune_samples(sample_data(gut.data)$biome_id == 6, gut.data)
bacterial_plotting2(gut.data.6)

## Plot vaginal
vaginal.data.6 <- prune_samples(sample_data(vaginal.data)$biome_id == 6, vaginal.data)
bacterial_plotting2(vaginal.data.6)

#### Participant 7

## Plot gut
gut.data.7 <- prune_samples(sample_data(gut.data)$biome_id == 7, gut.data) # F1594_V3V5_S1096, S1587_2_V3V5_S1535 - CANT CONFIRM
bacterial_plotting2(gut.data.7)

## Plot vaginal
vaginal.data.7 <- prune_samples(sample_data(vaginal.data)$biome_id == 7, vaginal.data) # S1853_V3V5_S384 - CANT CONFIRM
bacterial_plotting2(vaginal.data.7)

#### Participant 10

## Plot gut
gut.data.10 <- prune_samples(sample_data(gut.data)$biome_id == 10, gut.data) # S3977_2_V3V5_S1523, S8810_thawed_V3V5_S1751, S8812_V3V5_S142
bacterial_plotting2(gut.data.10)

## Plot vaginal
vaginal.data.10 <- prune_samples(sample_data(vaginal.data)$biome_id == 10, vaginal.data) # S3988_V3V5_S1650, S8839_V3V5_S1497
bacterial_plotting2(vaginal.data.10)

# Swap swabs
sample_data(bacterial.data)["S3977_2_V3V5_S1523", "sampleType"] <- "vaginal"
sample_data(bacterial.data)["S3977_V3V5_S1235", "sampleType"] <- "vaginal"
sample_data(bacterial.data)["S8810_thawed_V3V5_S1751", "sampleType"] <- "vaginal"
sample_data(bacterial.data)["S8812_V3V5_S142", "sampleType"] <- "vaginal"
sample_data(bacterial.data)["S3988_V3V5_S1650", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S8839_V3V5_S1497", "sampleType"] <- "fecal"

#### Participant 11
bacterial_plotting(vaginal.data, 11, "2022-11-14") # F1442_V3V5_S1520

bacterial_plotting(gut.data, 11, c("2022-10-21", "2022-10-27"))

## Plot gut
gut.data.11 <- prune_samples(sample_data(gut.data)$biome_id == 11, gut.data)
bacterial_plotting2(gut.data.11)

## Plot vaginal
vaginal.data.11 <- prune_samples(sample_data(vaginal.data)$biome_id == 11, vaginal.data)
bacterial_plotting2(vaginal.data.11)

# Swap swabs
sample_data(bacterial.data)["F1442_V3V5_S1520", "sampleType"] <- "fecal"

#### Participant 12

bacterial_plotting(vaginal.data, 12, "2022-10-29") # S8403_V3V5_S1122
bacterial_plotting(vaginal.data, 12, "2022-11-02") # S8371_V3V5_S1697

## Plot gut
gut.data.12 <- prune_samples(sample_data(gut.data)$biome_id == 12, gut.data)
bacterial_plotting2(gut.data.12)

## Plot vaginal
vaginal.data.12 <- prune_samples(sample_data(vaginal.data)$biome_id == 12, vaginal.data)
bacterial_plotting2(vaginal.data.12)

# Swap swabs
sample_data(bacterial.data)["S8403_V3V5_S1122", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S8371_V3V5_S1697", "sampleType"] <- "fecal"

#### Participant 13

## Plot gut
gut.data.13 <- prune_samples(sample_data(gut.data)$biome_id == 13, gut.data)
bacterial_plotting2(gut.data.13)

## Plot vaginal
vaginal.data.13 <- prune_samples(sample_data(vaginal.data)$biome_id == 13, vaginal.data)
bacterial_plotting2(vaginal.data.13)

#### Participant 14

bacterial_plotting(vaginal.data, 14, "2022-10-15") # S18_V3V5_S1536
bacterial_plotting(vaginal.data, 14, "2022-10-30") # S7953_V3V5_S1713

## Plot gut
gut.data.14 <- prune_samples(sample_data(gut.data)$biome_id == 14, gut.data)
bacterial_plotting2(gut.data.14)

## Plot vaginal
vaginal.data.14 <- prune_samples(sample_data(vaginal.data)$biome_id == 14, vaginal.data)
bacterial_plotting2(vaginal.data.14)

# Swap swabs
sample_data(bacterial.data)["S18_V3V5_S1536", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S7953_V3V5_S1713", "sampleType"] <- "fecal"

#### Participant 15
#### Participant 16

#### Participant 17

bacterial_plotting(vaginal.data, 17, "2022-10-31") # S9274_V3V5_S1542
bacterial_plotting(vaginal.data, 17, "2022-11-22") # F264_V3V5_S1347

## Plot gut
gut.data.17 <- prune_samples(sample_data(gut.data)$biome_id == 17, gut.data)
bacterial_plotting2(gut.data.17)

## Plot vaginal
vaginal.data.17 <- prune_samples(sample_data(vaginal.data)$biome_id == 17, vaginal.data)
bacterial_plotting2(vaginal.data.17)

# Swap swabs
sample_data(bacterial.data)["S9274_V3V5_S1542", "sampleType"] <- "fecal"
sample_data(bacterial.data)["F264_V3V5_S1347", "sampleType"] <- "fecal"


#### Participant 18

bacterial_plotting(vaginal.data, 18, "2022-11-07") # S7293_V3V5_S153
bacterial_plotting(gut.data, 18, "2022-10-18") # S1814_V3V5_S388

## Plot gut
gut.data.18 <- prune_samples(sample_data(gut.data)$biome_id == 18, gut.data)
bacterial_plotting2(gut.data.18)

## Plot vaginal
vaginal.data.18 <- prune_samples(sample_data(vaginal.data)$biome_id == 18, vaginal.data)
bacterial_plotting2(vaginal.data.18)

# Swap swabs
sample_data(bacterial.data)["S7293_V3V5_S153", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S1814_V3V5_S388", "sampleType"] <- "vaginal"

#### Participant 20
#### Participant 21

#### Participant 22

bacterial_plotting(gut.data, 22, "2022-10-26") # S2169_V3V5_S1589

## Plot gut
gut.data.22 <- prune_samples(sample_data(gut.data)$biome_id == 22, gut.data)
bacterial_plotting2(gut.data.22)

## Plot vaginal
vaginal.data.22 <- prune_samples(sample_data(vaginal.data)$biome_id == 22, vaginal.data)
bacterial_plotting2(vaginal.data.22)

# Swap swabs
sample_data(bacterial.data)["S2169_V3V5_S1589", "sampleType"] <- "vaginal"

#### Participant 23

bacterial_plotting(vaginal.data, 23, "2022-10-15") # S2341_V3V5_S1164
bacterial_plotting(vaginal.data, 23, "2022-12-03") # F1781_V3V5_S1084

## Plot gut
gut.data.23 <- prune_samples(sample_data(gut.data)$biome_id == 23, gut.data)
bacterial_plotting2(gut.data.23)

## Plot vaginal
vaginal.data.23 <- prune_samples(sample_data(vaginal.data)$biome_id == 23, vaginal.data)
bacterial_plotting2(vaginal.data.23)

# Swap swabs
sample_data(bacterial.data)["S2341_V3V5_S1164", "sampleType"] <- "fecal"
sample_data(bacterial.data)["F1781_V3V5_S1084", "sampleType"] <- "fecal"

#### Participant 24

#### Participant 25

bacterial_plotting(vaginal.data, 25, "2022-11-03") # S7323_V3V5_S172

## Plot gut
gut.data.25 <- prune_samples(sample_data(gut.data)$biome_id == 25, gut.data)
bacterial_plotting2(gut.data.25)

## Plot vaginal
vaginal.data.25 <- prune_samples(sample_data(vaginal.data)$biome_id == 25, vaginal.data)
bacterial_plotting2(vaginal.data.25)

# Swap swabs
sample_data(bacterial.data)["S7323_V3V5_S172", "sampleType"] <- "fecal"

#### Participant 26
#### Participant 27
#### Participant 29
#### Participant 30
#### Participant 31

#### Participant 32

bacterial_plotting(gut.data, 32, "2022-10-19") # S1007_V3V5_S499

## Plot gut
gut.data.32 <- prune_samples(sample_data(gut.data)$biome_id == 32, gut.data)
bacterial_plotting2(gut.data.32)

## Plot vaginal
vaginal.data.32 <- prune_samples(sample_data(vaginal.data)$biome_id == 32, vaginal.data)
bacterial_plotting2(vaginal.data.32)

# Swap swabs
sample_data(bacterial.data)["S1007_V3V5_S499", "sampleType"] <- "vaginal"

#### Participant 33

#### Participant 34

bacterial_plotting(vaginal.data, 34, "2022-10-14") # S4220_V3V5_S837
bacterial_plotting(vaginal.data, 34, "2022-10-19") # S4348_V3V5_S1149

## Plot gut
gut.data.34 <- prune_samples(sample_data(gut.data)$biome_id == 34, gut.data)
bacterial_plotting2(gut.data.34)

## Plot vaginal
vaginal.data.34 <- prune_samples(sample_data(vaginal.data)$biome_id == 34, vaginal.data)
bacterial_plotting2(vaginal.data.34)

# Swap swabs
sample_data(bacterial.data)["S4220_V3V5_S837", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S4348_V3V5_S1149", "sampleType"] <- "fecal"

#### Participant 35

#### Participant 36

bacterial_plotting(vaginal.data, 36, "2022-11-27") # S7498_V3V5_S1678

## Plot gut
gut.data.36 <- prune_samples(sample_data(gut.data)$biome_id == 36, gut.data)
bacterial_plotting2(gut.data.36)

## Plot vaginal
vaginal.data.36 <- prune_samples(sample_data(vaginal.data)$biome_id == 36, vaginal.data)
bacterial_plotting2(vaginal.data.36)

# Swap swabs
sample_data(bacterial.data)["S7498_V3V5_S1678", "sampleType"] <- "fecal"

#### Participant 37

## Plot gut
gut.data.37 <- prune_samples(sample_data(gut.data)$biome_id == 37, gut.data)
bacterial_plotting2(gut.data.37)

## Plot vaginal
vaginal.data.37 <- prune_samples(sample_data(vaginal.data)$biome_id == 37, vaginal.data)
bacterial_plotting2(vaginal.data.37)

#### Participant 38

bacterial_plotting(gut.data, 38, "2022-12-05") # S8947_V3V5_S104

## Plot gut
gut.data.38 <- prune_samples(sample_data(gut.data)$biome_id == 38, gut.data)
bacterial_plotting2(gut.data.38)

## Plot vaginal
vaginal.data.38 <- prune_samples(sample_data(vaginal.data)$biome_id == 38, vaginal.data)
bacterial_plotting2(vaginal.data.38)

# Swap swabs
sample_data(bacterial.data)["S8947_V3V5_S104", "sampleType"] <- "vaginal"

#### Participant 39

bacterial_plotting(vaginal.data, 39, "2022-10-15") # S3244_V3V5_S1187
bacterial_plotting(vaginal.data, 39, "2022-11-07") # S2866_thawed_V3V5_S672

## Plot gut
gut.data.39 <- prune_samples(sample_data(gut.data)$biome_id == 39, gut.data)
bacterial_plotting2(gut.data.39)

## Plot vaginal
vaginal.data.39 <- prune_samples(sample_data(vaginal.data)$biome_id == 39, vaginal.data)
bacterial_plotting2(vaginal.data.39)

# Swap swabs
sample_data(bacterial.data)["S3244_V3V5_S1187", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S3244_2_V3V5_S1475", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S2866_thawed_V3V5_S672", "sampleType"] <- "fecal"

#### Participant 40

bacterial_plotting(vaginal.data, 40, "2022-10-14") # S1747_V3V5_S860
bacterial_plotting(gut.data, 40, "2022-11-10") # S8467_thawed_V3V5_S244

## Plot gut
gut.data.40 <- prune_samples(sample_data(gut.data)$biome_id == 40, gut.data)
bacterial_plotting2(gut.data.40)

## Plot vaginal
vaginal.data.40 <- prune_samples(sample_data(vaginal.data)$biome_id == 40, vaginal.data)
bacterial_plotting2(vaginal.data.40)

# Swap swabs
sample_data(bacterial.data)["S1747_V3V5_S860", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S8467_thawed_V3V5_S244", "sampleType"] <- "vaginal"

#### Participant 41
#### Participant 42

#### Participant 43

bacterial_plotting(gut.data, 43, "2022-11-03") # S3854_V3V5_S185

## Plot gut
gut.data.43 <- prune_samples(sample_data(gut.data)$biome_id == 43, gut.data)
bacterial_plotting2(gut.data.43)

## Plot vaginal
vaginal.data.43 <- prune_samples(sample_data(vaginal.data)$biome_id == 43, vaginal.data)
bacterial_plotting2(vaginal.data.43)

# Swap swabs
sample_data(bacterial.data)["S3854_V3V5_S185", "sampleType"] <- "vaginal"

#### Participant 44

bacterial_plotting(gut.data, 44, "2022-10-16") # S3195_V3V5_S191 - CAN'T CONFIRM

#### Participant 45

## Plot gut
gut.data.45 <- prune_samples(sample_data(gut.data)$biome_id == 45, gut.data)
bacterial_plotting2(gut.data.45)

## Plot vaginal
vaginal.data.45 <- prune_samples(sample_data(vaginal.data)$biome_id == 45, vaginal.data)
bacterial_plotting2(vaginal.data.45)

#### Participant 46

bacterial_plotting(vaginal.data, 46, "2022-10-15") # S1352_V3V5_S1453

## Plot gut
gut.data.46 <- prune_samples(sample_data(gut.data)$biome_id == 46, gut.data)
bacterial_plotting2(gut.data.46)

## Plot vaginal
vaginal.data.46 <- prune_samples(sample_data(vaginal.data)$biome_id == 46, vaginal.data)
bacterial_plotting2(vaginal.data.46)

# Swap swabs
sample_data(bacterial.data)["S1352_V3V5_S1453", "sampleType"] <- "fecal"

#### Participant 47

#### Participant 48

bacterial_plotting(vaginal.data, 48, "2022-10-16") # S1831_V3V5_S1221
bacterial_plotting(vaginal.data, 48, "2022-11-05") # S9507_thawed_V3V5_S653
bacterial_plotting(vaginal.data, 48, "2022-12-01") # F475_V3V5_S1700
bacterial_plotting3(gut.data, 48, "2022-10-16") # S1825_V3V5_S335
bacterial_plotting(gut.data, 48, "2022-12-01") # F473_blood_V3V5_S1898

## Plot gut
gut.data.48 <- prune_samples(sample_data(gut.data)$biome_id == 48, gut.data)
bacterial_plotting2(gut.data.48)

## Plot vaginal
vaginal.data.48 <- prune_samples(sample_data(vaginal.data)$biome_id == 48, vaginal.data)
bacterial_plotting2(vaginal.data.48)

# Swap swabs
sample_data(bacterial.data)["S1831_V3V5_S1221", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S9507_thawed_V3V5_S653", "sampleType"] <- "fecal"
sample_data(bacterial.data)["F475_V3V5_S1700", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S1825_V3V5_S335", "sampleType"] <- "vaginal"
sample_data(bacterial.data)["F473_blood_V3V5_S1898", "sampleType"] <- "vaginal"

#### Participant 49

bacterial_plotting(gut.data, 49, "2022-10-16") # S3791_2_V3V5_S1511

## Plot gut
gut.data.49 <- prune_samples(sample_data(gut.data)$biome_id == 49, gut.data)
bacterial_plotting2(gut.data.49)

## Plot vaginal
vaginal.data.49 <- prune_samples(sample_data(vaginal.data)$biome_id == 49, vaginal.data)
bacterial_plotting2(vaginal.data.49)

sample_data(bacterial.data)["S3791_2_V3V5_S1511", "sampleType"] <- "vaginal"
sample_data(bacterial.data)["S3791_V3V5_S1223", "sampleType"] <- "vaginal"

#### Participant 50

bacterial_plotting(vaginal.data, 50, "2022-10-26") # S2112_V3V5_S1205
bacterial_plotting(vaginal.data, 50, "2022-11-08") # S9212_V3V5_S136

## Plot gut
gut.data.50 <- prune_samples(sample_data(gut.data)$biome_id == 50, gut.data)
bacterial_plotting2(gut.data.50)

## Plot vaginal
vaginal.data.50 <- prune_samples(sample_data(vaginal.data)$biome_id == 50, vaginal.data)
bacterial_plotting2(vaginal.data.50)

# Swap swabs
sample_data(bacterial.data)["S2112_V3V5_S1205", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S9212_V3V5_S136", "sampleType"] <- "fecal"

#### Participant 51

bacterial_plotting3(vaginal.data, 51, "2022-11-30") # F2432_V3V5_S1909 - CAN'T TELL
bacterial_plotting3(vaginal.data, 51, "2022-10-14")

## Plot gut
gut.data.51 <- prune_samples(sample_data(gut.data)$biome_id == 51, gut.data)
bacterial_plotting2(gut.data.51)

## Plot vaginal
vaginal.data.51 <- prune_samples(sample_data(vaginal.data)$biome_id == 51, vaginal.data)
bacterial_plotting2(vaginal.data.51)

#### Participant 52

bacterial_plotting3(vaginal.data, 52, "2022-10-16") # S4103_V3V5_S1175

## Plot gut
gut.data.52 <- prune_samples(sample_data(gut.data)$biome_id == 52, gut.data)
bacterial_plotting2(gut.data.52)

## Plot vaginal
vaginal.data.52 <- prune_samples(sample_data(vaginal.data)$biome_id == 52, vaginal.data)
bacterial_plotting2(vaginal.data.52)

# Swap swabs
sample_data(bacterial.data)["S4103_V3V5_S1175", "sampleType"] <- "fecal"

#### Participant 53

bacterial_plotting(vaginal.data, 53, "2022-10-20") # S1719_V3V5_S755
bacterial_plotting3(gut.data, 53, "2022-10-20") # S1721_V3V5_S1432
bacterial_plotting3(gut.data, 53, "2022-10-17") # "S2558_2_V3V5_S1499

## Plot gut
gut.data.53 <- prune_samples(sample_data(gut.data)$biome_id == 53, gut.data)
bacterial_plotting2(gut.data.53)

## Plot vaginal
vaginal.data.53 <- prune_samples(sample_data(vaginal.data)$biome_id == 53, vaginal.data)
bacterial_plotting2(vaginal.data.53)

# Swap swabs
sample_data(bacterial.data)["S1719_V3V5_S755", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S1721_V3V5_S1432", "sampleType"] <- "vaginal"
sample_data(bacterial.data)["S2558_2_V3V5_S1499", "sampleType"] <- "vaginal"
sample_data(bacterial.data)["S2558_V3V5_S1211", "sampleType"] <- "vaginal"

#### Participant 54

bacterial_plotting(vaginal.data, 54, "2022-12-04") # S8491_V3V5_S418

## Plot gut
gut.data.54 <- prune_samples(sample_data(gut.data)$biome_id == 54, gut.data)
bacterial_plotting2(gut.data.54)

## Plot vaginal
vaginal.data.54 <- prune_samples(sample_data(vaginal.data)$biome_id == 54, vaginal.data)
bacterial_plotting2(vaginal.data.54)

# Swap swabs
sample_data(bacterial.data)["S8491_V3V5_S418", "sampleType"] <- "fecal"

#### Participant 55

bacterial_plotting(gut.data, 55, "2022-10-15") # S2792_V3V5_S761
bacterial_plotting(vaginal.data, 55, "2022-11-03") # S2448_thawed_V3V5_S225

## Plot gut
gut.data.55 <- prune_samples(sample_data(gut.data)$biome_id == 55, gut.data) 
bacterial_plotting2(gut.data.55)

## Plot vaginal
vaginal.data.55 <- prune_samples(sample_data(vaginal.data)$biome_id == 55, vaginal.data)
bacterial_plotting2(vaginal.data.55)

# Swap swabs
sample_data(bacterial.data)["S2448_thawed_V3V5_S225", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S2792_V3V5_S761", "sampleType"] <- "vaginal"

#### Participant 56

bacterial_plotting(vaginal.data, 56, "2022-10-16") # S4693_V3V5_S1185
bacterial_plotting(vaginal.data, 56, "2022-10-23") # S4659_V3V5_S1596
bacterial_plotting(vaginal.data, 56, "2022-10-25") # S4691_V3V5_S1256
bacterial_plotting(vaginal.data, 56, "2022-10-30") # S6378_V3V5_S1722
bacterial_plotting(vaginal.data, 56, "2022-11-02") # S6409_V3V5_S1224
bacterial_plotting(vaginal.data, 56, "2022-11-03") # S6415_V3V5_S1662
bacterial_plotting(vaginal.data, 56, "2022-11-05") # S6397_thawed_V3V5_S580
bacterial_plotting(vaginal.data, 56, "2022-11-06") # S6468_thawed_V3V5_S671
bacterial_plotting(vaginal.data, 56, "2022-11-17") # S6723_V3V5_S1738
bacterial_plotting(vaginal.data, 56, "2022-11-22") # S6712_V3V5_S1135
bacterial_plotting(vaginal.data, 56, "2022-11-27") # S6716_V3V5_S1666
bacterial_plotting3(gut.data, 56, "2022-12-05") # S6403_V3V5_S1726

## Plot gut
gut.data.56 <- prune_samples(sample_data(gut.data)$biome_id == 56, gut.data) 
bacterial_plotting2(gut.data.56)

## Plot vaginal
vaginal.data.56 <- prune_samples(sample_data(vaginal.data)$biome_id == 56, vaginal.data)
bacterial_plotting2(vaginal.data.56)

# Swap swabs
sample_data(bacterial.data)["S4693_V3V5_S1185", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S4659_V3V5_S1596", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S4691_V3V5_S1256", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S6378_V3V5_S1722", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S6409_V3V5_S1224", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S6415_V3V5_S1662", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S6397_thawed_V3V5_S580", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S6468_thawed_V3V5_S671", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S6723_V3V5_S1738", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S6712_V3V5_S1135", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S6716_V3V5_S1666", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S6403_V3V5_S1726", "sampleType"] <- "vaginal"

#### Participant 58

#### Participant 59

bacterial_plotting(vaginal.data, 59, "2022-10-14") # S2364_V3V5_S850
bacterial_plotting(vaginal.data, 59, "2022-11-05") # S2420_V3V5_S1735
bacterial_plotting3(gut.data, 59, "2022-11-05") # S6306_V3V5_S604

## Plot gut
gut.data.59 <- prune_samples(sample_data(gut.data)$biome_id == 59, gut.data) 
bacterial_plotting2(gut.data.59)

## Plot vaginal
vaginal.data.59 <- prune_samples(sample_data(vaginal.data)$biome_id == 59, vaginal.data)
bacterial_plotting2(vaginal.data.59)

# Swap swabs
sample_data(bacterial.data)["S2364_V3V5_S850", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S2420_V3V5_S1735", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S6306_V3V5_S604", "sampleType"] <- "vaginal"

#### Participant 60

bacterial_plotting(vaginal.data, 60, "2022-11-30") # F240_V3V5_S1698

## Plot gut
gut.data.60 <- prune_samples(sample_data(gut.data)$biome_id == 60, gut.data) 
bacterial_plotting2(gut.data.60)

## Plot vaginal
vaginal.data.60 <- prune_samples(sample_data(vaginal.data)$biome_id == 60, vaginal.data)
bacterial_plotting2(vaginal.data.60)

# Swap swabs
sample_data(bacterial.data)["F240_V3V5_S1698", "sampleType"] <- "fecal"

#### Participant 61

## Plot gut
gut.data.61 <- prune_samples(sample_data(gut.data)$biome_id == 61, gut.data) 
bacterial_plotting2(gut.data.61)

## Plot vaginal
vaginal.data.61 <- prune_samples(sample_data(vaginal.data)$biome_id == 61, vaginal.data)
bacterial_plotting2(vaginal.data.61)

#### Participant 62

bacterial_plotting(vaginal.data, 62, "2022-10-27") # S3319_V3V5_S1621

## Plot gut
gut.data.62 <- prune_samples(sample_data(gut.data)$biome_id == 62, gut.data) 
bacterial_plotting2(gut.data.62)

## Plot vaginal
vaginal.data.62 <- prune_samples(sample_data(vaginal.data)$biome_id == 62, vaginal.data)
bacterial_plotting2(vaginal.data.62)

# Swap swabs
sample_data(bacterial.data)["S3319_V3V5_S1621", "sampleType"] <- "fecal"

#### Participant 63

bacterial_plotting(vaginal.data, 63, "2022-11-07") # S8297_thawed_V3V5_S668

## Plot gut
gut.data.63 <- prune_samples(sample_data(gut.data)$biome_id == 63, gut.data) 
bacterial_plotting2(gut.data.63)

## Plot vaginal
vaginal.data.63 <- prune_samples(sample_data(vaginal.data)$biome_id == 63, vaginal.data)
bacterial_plotting2(vaginal.data.63)

# Swap swabs
sample_data(bacterial.data)["S8297_thawed_V3V5_S668", "sampleType"] <- "fecal"

#### Participant 64

bacterial_plotting(vaginal.data, 64, "2022-10-23") # S3326_V3V5_S1182
bacterial_plotting3(gut.data, 64, "2022-11-26") # F74_capoff_V3V5_S1673

## Plot gut
gut.data.64 <- prune_samples(sample_data(gut.data)$biome_id == 64, gut.data) 
bacterial_plotting2(gut.data.64)

## Plot vaginal
vaginal.data.64 <- prune_samples(sample_data(vaginal.data)$biome_id == 64, vaginal.data)
bacterial_plotting2(vaginal.data.64)

# Swap swabs
sample_data(bacterial.data)["S3326_V3V5_S1182", "sampleType"] <- "fecal"
sample_data(bacterial.data)["F74_capoff_V3V5_S1673", "sampleType"] <- "vaginal"

#### Participant 64
#### Participant 65

#### Participant 66

bacterial_plotting(vaginal.data, 66, "2022-11-10") # S8913_thawed_V3V5_S267
bacterial_plotting(vaginal.data, 66, "2022-11-21") # S5726_V3V5_S1138
bacterial_plotting(vaginal.data, 66, "2022-12-12") # F1379_thawed_V3V5_S631
bacterial_plotting3(gut.data, 66, "2022-11-08") # S8887_V3V5_S141
bacterial_plotting3(gut.data, 66, "2022-12-04") # F1612_V3V5_S99

## Plot gut
gut.data.66 <- prune_samples(sample_data(gut.data)$biome_id == 66, gut.data) 
bacterial_plotting2(gut.data.66) # 2

## Plot vaginal
vaginal.data.66 <- prune_samples(sample_data(vaginal.data)$biome_id == 66, vaginal.data)
bacterial_plotting2(vaginal.data.66) # 3

# Swap swabs
sample_data(bacterial.data)["S8913_thawed_V3V5_S267", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S5726_V3V5_S1138", "sampleType"] <- "fecal"
sample_data(bacterial.data)["F1379_thawed_V3V5_S631", "sampleType"] <- "fecal"
sample_data(bacterial.data)["S8887_V3V5_S141", "sampleType"] <- "vaginal"
sample_data(bacterial.data)["F1612_V3V5_S99", "sampleType"] <- "vaginal"

#### Participant 67

bacterial_plotting(gut.data, 67, "2022-10-14") # S2964_V3V5_S752
bacterial_plotting3(vaginal.data, 67, "2022-12-14") # S9450_thawed_V3V5_S990

## Plot gut
gut.data.67 <- prune_samples(sample_data(gut.data)$biome_id == 67, gut.data) 
bacterial_plotting2(gut.data.67) # 1

## Plot vaginal
vaginal.data.67 <- prune_samples(sample_data(vaginal.data)$biome_id == 67, vaginal.data)
bacterial_plotting2(vaginal.data.67) # 1

# Swap swabs
sample_data(bacterial.data)["S2964_V3V5_S752", "sampleType"] <- "vaginal"
sample_data(bacterial.data)["S9450_thawed_V3V5_S990", "sampleType"] <- "fecal"

#### Participant 68

#### Participant 69

bacterial_plotting(gut.data, 69, "2022-11-10") # S6023_V3V5_S1645
bacterial_plotting3(vaginal.data, 69, "2022-11-05") # S6015_thawed_V3V5_S629

## Plot gut
gut.data.69 <- prune_samples(sample_data(gut.data)$biome_id == 69, gut.data) 
bacterial_plotting2(gut.data.69) # 1

## Plot vaginal
vaginal.data.69 <- prune_samples(sample_data(vaginal.data)$biome_id == 69, vaginal.data)
bacterial_plotting2(vaginal.data.69) # 1

# Swap swabs
sample_data(bacterial.data)["S6023_V3V5_S1645", "sampleType"] <- "vaginal"
sample_data(bacterial.data)["S6015_thawed_V3V5_S629", "sampleType"] <- "fecal"

#### Participant 70

#### Participant 71

bacterial_plotting3(vaginal.data, 71, "2022-10-30") # S353_V3V5_S1124

## Plot gut
gut.data.71 <- prune_samples(sample_data(gut.data)$biome_id == 71, gut.data) 
bacterial_plotting2(gut.data.71) # 

## Plot vaginal
vaginal.data.71 <- prune_samples(sample_data(vaginal.data)$biome_id == 71, vaginal.data)
bacterial_plotting2(vaginal.data.71) # 1

# Swap swabs
sample_data(bacterial.data)["S353_V3V5_S1124", "sampleType"] <- "fecal"

#### Participant 73

## Plot gut
gut.data.73 <- prune_samples(sample_data(gut.data)$biome_id == 73, gut.data) 
bacterial_plotting2(gut.data.73) 

## Plot vaginal
vaginal.data.73 <- prune_samples(sample_data(vaginal.data)$biome_id == 73, vaginal.data)
bacterial_plotting2(vaginal.data.73)

# Save new obj
save.image("/Volumes/T7/microbiome_data/R_environments/microbiome_vaginal_gut_sample_check.RData")

otu_table(bacterial.data) <- t(otu_table(bacterial.data))

# Save new obj
saveRDS(bacterial.data, file = "/Volumes/T7/microbiome_data/sequenced_data/02_11/bacteria_cleanedv2.rds")

####################################################################################
# fix data

samples.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_samples.csv")

metadata.22 <- as(sample_data(bacterial.data), "data.frame")
metadata.22.subset <- metadata.22 %>% 
  select(qr, sampleType) %>% 
  rename(sampleType_fixed = sampleType)

samples.data <- samples.data %>% 
  left_join(metadata.22.subset, by="qr")

samples.data <- samples.data %>% 
  mutate(sampleType=ifelse(is.na(sampleType), sampleType, sampleType_fixed)) %>% 
  select(!sampleType_fixed) %>% 
  distinct()

dupes <- names(table(samples.data$qr))[table(samples.data$qr) > 1]
dupes

write.csv(samples.data,
          file = "/Volumes/T7/microbiome_data/cleaned_data/cleaned_samplesv2.csv",
          row.names = FALSE)
