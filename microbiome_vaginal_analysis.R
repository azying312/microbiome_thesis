library(vegan)
library(phyloseq)
library(tidyverse)

source("~/Microbiome Thesis/functions.R")

bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/relabeled_data/vaginal_cleaned_max_taxa.rds")
# bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/vaginal_cleaned_max_taxa.rds")

###############################################################################################

bacteria_taxa_df <- tax_table(bacterial.data)
bacteria_metadata_df <- sample_data(bacterial.data)

length(unique(bacteria_metadata_df$biome_id))

#### Exploratory Data Analysis

# Relative abundances
vaginal_relative_abundances <- transform_sample_counts(bacterial.data, function(x) x/sum(x))
relative_abundance_otu <- as.data.frame(otu_table(vaginal_relative_abundances))
relative_abundance_otu_t <- t(relative_abundance_otu) %>% as.data.frame()

# Add sample ID
relative_abundance_otu$SampleID <- rownames(relative_abundance_otu)
bacteria_metadata_df <- sample_data(bacterial.data)
bacteria_metadata_df <- bacteria_metadata_df[,-2]
otu_with_participant <-  relative_abundance_otu %>%
  left_join(bacteria_metadata_df, by="SampleID")
rownames(otu_with_participant) <- otu_with_participant$SampleID
relative_abundance_otu <- relative_abundance_otu %>%
  dplyr::select(!SampleID)

dim(bacteria_metadata_df)

###################################################################################################

###### Five Community State Type (CSTs) Clustering

# species to CST mapping
species_to_cst <- data.frame(
  Species = c("crispatus", "gasseri", "iners", "jensenii"),
  CST = c("I", "II", "III", "V") # IV is anaerobic/diverse cluster
)

# Get most abundant OTU per sample
# max_taxa <- apply(relative_abundance_otu_t, 2, function(sample) { 
  max_taxa <- apply(relative_abundance_otu_t, 1, function(sample) { # on relabeled data
  taxa_idx <- which.max(sample)
  taxa_names(vaginal_relative_abundances)[taxa_idx]
})

# Map most abundant OTU to the sample data
bacteria_metadata_df$max_taxa <- max_taxa

# Get max taxa names
# bacteria_metadata_df$OTU <- as.character(bacteria_taxa_df[bacteria_metadata_df$max_taxa, "Species_exact"]) # 2039 samples
bacteria_metadata_df$OTU <- as.character(bacteria_taxa_df[bacteria_metadata_df$max_taxa, "BLAST_species"])

# Set sample data in phyloseq obj
sample_data(bacterial.data) <- bacteria_metadata_df
# bacteria_taxa_df <- as.data.frame(bacteria_taxa_df)

# Assign CST col
bacteria_metadata_df <- as(bacteria_metadata_df, "data.frame")
class(bacteria_metadata_df)

bacteria_metadata_df <- bacteria_metadata_df %>%
  mutate(
    CST = case_when(
      str_detect(OTU, "gasseri") ~ "II",
      str_detect(OTU, "crispatus") ~ "I",
      str_detect(OTU, "iners") ~ "III",
      str_detect(OTU, "jensenii") ~ "V",
      TRUE ~ "IV" # Assign "IV" for diverse/anaerobic or unclassified species
    )
  )

table(bacteria_metadata_df$CST)
cst_summary <- bacteria_metadata_df %>% 
  count(CST) %>% 
  mutate(Percentage=100*(n/sum(n)))
print(cst_summary)

sample.cst <- bacteria_metadata_df$CST %>% 
  as.data.frame()
colnames(sample.cst) <- "CST"
rownames(sample.cst) <- rownames(bacteria_metadata_df)
sample.cst$SampleID <- rownames(sample.cst)

dim(bacteria_metadata_df)

################################################################################

head(bacteria_metadata_df)
vaginal_relative_abundances <- transform_sample_counts(bacterial.data, function(x) x/sum(x))

################################################################################

bacteria_taxa_table <- tax_table(bacterial.data)
bacteria_taxa_df <- as.data.frame(bacteria_taxa_table)

########################################################

bacteria_metadata_df <- sample_data(bacterial.data)
# bacteria_metadata_df <- as.data.frame(as(otu_table(bacterial.data_subset), "matrix"))
otu_table_df <- as(otu_table(bacterial.data), "matrix")

# Aggregate data by participants by mean relative abundance for a given OTU
# participant_otu <- tapply(sample_names(bacteria_metadata_df),
#                           sample_data(bacterial.data)$biome_id,
#                           function(samples) rowMeans(t(otu_table(vaginal_relative_abundances))[, samples, drop = FALSE]))

participant_otu <- tapply(sample_names(bacteria_metadata_df),
                          sample_data(bacterial.data)$biome_id,
                          function(samples) rowMeans(otu_table(vaginal_relative_abundances)[, samples, drop = FALSE]))

# 2015 1571
participant_otu <- do.call(cbind, participant_otu)
rownames(participant_otu) <- taxa_names(vaginal_relative_abundances)

## Alpha Div - Shannon Index
shannon.24 <- vegan::diversity(t(otu_table_df), "shannon")

# Add participant IDs from sample data | Merge the calculated Shannon diversity values with metadata
bacteria_metadata_df <- as(bacteria_metadata_df, "data.frame")
shannon.df.24 <- data.frame("SampleID"=names(shannon.24), "shannon"=shannon.24)
shannon.qr.merged.24 <- merge(shannon.df.24, bacteria_metadata_df, by="SampleID")
shannon.cst.qr.merged.24 <- merge(sample.cst, shannon.qr.merged.24, by="SampleID") %>%
  mutate(biome_id=as.integer(biome_id)) %>% 
  filter(!is.na(biome_id)) %>% 
  # Filter the data to be within study days: 10-13 to 12-14
  filter(logDate > "2022-10-11" & logDate < "2022-12-16") # 2038 to 1971 --> 1571 in relabeled

head(shannon.cst.qr.merged.24)

CST_max <- shannon.cst.qr.merged.24 %>% 
  group_by(biome_id) %>% 
  summarise(CST_max = names(sort(table(CST), decreasing = TRUE)[1]))

shannon.cst.qr.merged.24 <- shannon.cst.qr.merged.24 %>% 
  left_join(CST_max, by = "biome_id")
# Save R environment - did not resave after relabeling data
# save.image("/Volumes/T7/microbiome_data/R_environments/vaginal_microbiome_relAbundance.RData")

###########################################################################
library(tidyverse)
library(viridis)

### Check UMinn Spreadsheet v. Sequenced Data
uminn_data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_uminn_data.csv")
samples.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/relabeled_data/cleaned_samplesv2.csv")
# samples.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_samples.csv")

uminn_data <- uminn_data %>% 
  dplyr::select(Sample.ID, Special.Notes) %>% 
  filter(Sample.ID != "BLANK") %>% 
  # errors in sample processing from UMinn
  filter(!str_detect(Special.Notes, "error")) %>% 
  filter(!str_detect(Special.Notes, "No swab in tube")) #2880 -> 2791
uminn_data$qr <- sub("_.*", "", uminn_data$Sample.ID)

vaginal_data <- samples.data %>% 
  filter(sampleType=="vaginal")

uminn_data_vaginal <- uminn_data %>% 
  filter(qr %in% vaginal_data$qr)
uminn_data_vaginal <- uminn_data_vaginal %>% 
  left_join(vaginal_data, by="qr")

length(setdiff(shannon.cst.qr.merged.24$qr, uminn_data_vaginal$qr)) # 0 -- 386; turned into 522? -> relabeled: 5
length(setdiff(uminn_data_vaginal$qr, shannon.cst.qr.merged.24$qr)) # 57 -- 19; turned into 18 - less concerned -> relabeled: 57

diff.qr <- setdiff(shannon.cst.qr.merged.24$qr, uminn_data_vaginal$qr)

setdiff(unique(shannon.cst.qr.merged.24$biome_id), unique(uminn_data_vaginal$biome_id)) # 68 ---> none
setdiff(unique(uminn_data_vaginal$biome_id), unique(shannon.cst.qr.merged.24$biome_id)) # 20 ---> none

# write.csv(shannon.cst.qr.merged.24, "/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/shannon.cst.qr.merged.24.csv")
write.csv(shannon.cst.qr.merged.24, "/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/relabeled_data/shannon.cst.qr.merged.24.csv")

################################################################################

### Compare results with 2017-18 findings

# Section: Community state types of vaginal microbiota
shannon.cst.summary <- shannon.cst.qr.merged.24 %>% 
  group_by(biome_id, CST) %>% 
  mutate(count=n(), .groups="drop") %>% 
  group_by(biome_id) %>% 
  mutate(prop=count/sum(count)) %>% 
  dplyr::select(!.groups) %>% 
  ungroup() %>% 
  filter(!is.na(biome_id)) %>% 
  mutate(biome_id=as.numeric(biome_id))

dominant_cst <- shannon.cst.summary %>%
  group_by(biome_id) %>%
  filter(prop == max(prop)) %>%
  slice(1) %>% # person 52 has 1 sample in 2 CSTs
  dplyr::select(biome_id, CST) %>%
  distinct()

id_order <- dominant_cst %>% 
  arrange(CST, biome_id)

# number of indvs with CST info
length(unique(shannon.cst.summary$biome_id))

# FIGURE
cst_plt <- ggplot(shannon.cst.summary, aes(x=factor(biome_id, levels=id_order$biome_id), y=prop, fill=CST)) +
  geom_bar(stat="identity") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal() +
  labs(x = "Participant", y = "Proportion of Samples", fill = "CST") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  theme(legend.position = "bottom")
# CST label overlay
# geom_tile(data = dominant_cst, aes(x = factor(biome_id, levels = id_order$biome_id), y = -0.1, fill = CST), height = 0.05) +
# theme(legend.position = "right")
cst_plt

## Heatmap of CSTs for each participant
CST_tbl <- as.data.frame(table(shannon.cst.qr.merged.24$biome_id, shannon.cst.qr.merged.24$CST))
colnames(CST_tbl) <- c("biome_id", "CST", "Frequency")
ggplot(CST_tbl, aes(y=as.factor(biome_id), x=as.factor(CST), fill=Frequency)) +
  geom_tile(color="black") +
  scale_fill_viridis(option="C", direction=1) +
  theme_minimal() +
  labs(y="Participant", x="CST", title="", fill="Frequency")

## Species abundance visualization -- come back to this
# find which ppl have most samples submitted - try person 17

# Section: Menstrual fluctuations of vaginal microbiota.
# OLD DATA: menses.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/imputed_menstruation_data.csv")
# menses.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/imputed_menstruation_data_2_12.csv")
# RELABELED DATA
menses.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/relabeled_data/imputed_menstruation_data_3_11.csv")
menses.data <- menses.data %>% 
  rename_with(~gsub("X2022.", "2022.", .), starts_with("X2022.")) %>% 
  rename_with(~gsub("\\.", "-", .))

# Reshape
menses.data.long <- menses.data %>% 
  pivot_longer(cols=starts_with("2022-"), names_to="logDate", values_to="menses_status")

vaginal.microbial.menses.24 <- shannon.cst.qr.merged.24 %>% 
  left_join(menses.data.long, by=c("biome_id", "logDate"))

# make df of the Wilcox test comparing participants shannon index when they are on menses v. not
# vaginal.microbial.menses.24 <-  vaginal.microbial.menses.24 %>% 
#   mutate(menses_day = ifelse(menses_status %in% c(1,2,3,7,9,78), "menses", "not_menses"))

vaginal.microbial.menses.24 <- vaginal.microbial.menses.24 %>% 
  mutate(menses_day = ifelse(menses_status %in% c(1,2,3,7,9,78), "menses", 
                             ifelse(menses_status %in% c(4,5,6,10), "not_menses", NA)))

length(unique(menses.data$biome_id))
table(vaginal.microbial.menses.24$menses_day)

# write.csv(vaginal.microbial.menses.24, file="/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifetyle/vaginal.microbial.menses.24.csv")
# uncomment to save file
# write.csv(vaginal.microbial.menses.24, "/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/relabeled_data/vaginal.microbial.menses.24.csv")

############################################# MISSING DATA ISSUES [RESOLVED - 2/12]
## Check missing data
volunteer.df <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 9-Volunteer Medical History.csv")
# na.vaginal.microbial.menses.24 <- vaginal.microbial.menses.24 %>% 
#   filter(is.na(menses_day)) %>% 
#   select(c(SampleID, CST, shannon, biome_id, logDate, survey_menstruate, menses_status, menses_day, sampleType))
# 
# # 68 not in the uminn returned data
# no.survey.ids <- setdiff(na.vaginal.microbial.menses.24$biome_id, uminn_data_vaginal$biome_id)
# dim(na.vaginal.microbial.menses.24 %>% filter(biome_id %in% no.survey.ids))
# 
# # 70 71 69 73 68 don't have volunteer surveys
# no.survey.ids <- setdiff(na.vaginal.microbial.menses.24$biome_id, volunteer.df$biome_id)
# dim(na.vaginal.microbial.menses.24 %>% filter(biome_id %in% no.survey.ids))
# 
# impute.data <- na.vaginal.microbial.menses.24 %>% 
#   filter(!(biome_id %in% no.survey.ids))
# impute.data$qr <- sub("_.*", "", impute.data$SampleID)
# 
# # check sample ids with uminn_data_vaginal
# notUminn.qr <- setdiff(impute.data$qr, uminn_data$qr)
# length(setdiff(impute.data$qr, uminn_data$qr)) # 309 qr codes that were sequenced but not in Uminn dataset?
# # length(setdiff(uminn_data$qr, impute.data$qr))
# 
# impute.data <- impute.data %>% 
#   filter(!(qr %in% notUminn.qr))
# dim(impute.data)

na.menses.day <- vaginal.microbial.menses.24 %>% 
  filter(is.na(menses_day)) %>% 
  filter(biome_id %in% volunteer.df$biome_id)
# View(na.menses.day)
##########################################################################################

menses.table.df <- vaginal.microbial.menses.24 %>% 
  dplyr::select(biome_id, logDate, CST, shannon, qr, max_taxa, OTU, menses_status, menses_day)

# Wilcox Sign Rank test on menses v. not menses day
table(menses.table.df$biome_id, menses.table.df$menses_day)

wilcox.test(filter_id_data(menses.table.df, 17)$shannon ~ filter_id_data(menses.table.df, 17)$menses_day)
wilcox.test(filter_id_data(menses.table.df, 17)$shannon ~ filter_id_data(menses.table.df, 17)$menses_day)

###

# Plot Shannon diversity over time for each participant
participant_ids <- unique(vaginal.microbial.menses.24$biome_id)
all_days <- seq.Date(as.Date("2022-10-13"), as.Date("2022-12-16"), by = "day")
all_days <- data.frame(logDate=as.character(all_days))
file_path <- "/Volumes/T7/microbiome_data/graphics/Results from 2022 (compare to 2017-18)/relabeled_data/vaginal_shannon_diversity_logDates/"
for(id in participant_ids) {
  # print(id)
  file_name_id <- paste0(file_path, id, "_id.png")
  vaginal.microbial.menses.24.participant <- vaginal.microbial.menses.24 %>% 
    filter(biome_id==id)
  vaginal.microbial.menses.24.participant <- all_days %>%
    left_join(vaginal.microbial.menses.24.participant, by = "logDate")
  vaginal.microbial.menses.24.participant <- vaginal.microbial.menses.24.participant %>% 
    mutate(menses_day=ifelse(!is.na(menses_day), menses_day, 
                             ifelse(!is.na(shannon), "not_recorded", NA)))
  dupe.day <- vaginal.microbial.menses.24.participant %>% 
    group_by(logDate) %>% 
    filter(n() > 1) %>% 
    distinct(logDate)
  vaginal.microbial.menses.24.participant <- study_days(vaginal.microbial.menses.24.participant)
  
  shannon_plt <- ggplot(vaginal.microbial.menses.24.participant, 
                        aes(x=(study_day), y=(shannon), color=as.factor(menses_day), 
                            shape=as.factor(menses_day))) +
    geom_point() +
    scale_color_manual(name = "Menses Status",
                       values = c("not_menses" = "black", "menses" = "red", "not_recorded" = "purple"),
                       labels = c("not_menses" = "Not Menstruating", 
                                  "menses" = "Menstruating", 
                                  "not_recorded" = "Not Recorded")) +  
    scale_shape_manual(name = "Menses Status",
                       values = c("not_menses" = 16, "menses" = 17, "not_recorded"=18),
                       labels = c("not_menses" = "Not Menstruating", 
                                  "menses" = "Menstruating", 
                                  "not_recorded" = "Not Recorded")) +
    ylim(0,3)+
    # Add lines for the days with multiple samples
    geom_segment(data = dupe.day, 
                 aes(x = logDate, xend = logDate, y = -Inf, yend = Inf), 
                 linetype = "dashed", color = "gray50", inherit.aes = FALSE) +  
    theme_minimal() +
    labs(x = "Study Day", y = "Shannon Diversity", 
         # title=paste(id, "Shannon diversity days"),
         ) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1),
          text=element_text(size=16)) +
    scale_x_continuous(breaks = seq(0, max(vaginal.microbial.menses.24.participant$study_day), by = 5))
  
  ggsave(file_name_id, shannon_plt, width = 8, height = 6, dpi = 300)
  
  print(shannon_plt)
}

##########################################################################################

# Figure: plot of average shannon for each participant on menses vs. not on menses
vaginal.microbial.menses.24.summary <- vaginal.microbial.menses.24 %>% 
  group_by(biome_id, menses_day) %>% 
  summarise(avg_shannon=mean(shannon), .groups="drop") %>% 
  pivot_wider(names_from = menses_day, values_from = avg_shannon, names_prefix = "menses_day_") %>% 
  filter(!is.na(menses_day_not_menses) & !is.na(menses_day_menses)) %>% 
  dplyr::select(!menses_day_NA)
### Wilcox test
wilcox.test(vaginal.microbial.menses.24.summary$menses_day_not_menses,
            vaginal.microbial.menses.24.summary$menses_day_menses, paired=TRUE)
t.test(vaginal.microbial.menses.24.summary$menses_day_not_menses,
       vaginal.microbial.menses.24.summary$menses_day_menses, paired=TRUE)

# Menses: Average Shannon Diversity for Participants on menses v. not
length(unique(vaginal.microbial.menses.24.summary$biome_id))
ggplot(
  vaginal.microbial.menses.24.summary,
  aes(x = menses_day_not_menses, y = menses_day_menses)) + #,
     # color=as.factor(biome_id))) +
  geom_point(color="orchid") +
  xlim(0,4)+
  ylim(0,4)+
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  theme_minimal() +
  theme(
    text=element_text(size=15),
  ) +
  labs(x="Not Menstruating (Shannon)", y="Menstruating (Shannon)",
       color="Participant ID")

##########################################################################################

length(unique(vaginal.microbial.menses.24$biome_id))

# Long term CST transition - participant shifts from one CST to another and stays in the new CST for 10 days or more
longterm_cst_vaginal.microbial.menses.24 <- vaginal.microbial.menses.24 %>% 
  dplyr::select(SampleID, CST, shannon, qr, biome_id, logDate, sampleType, max_taxa, OTU, menses_day) %>% 
  # filter(!is.na(menses_day)) %>% # only if you need menses data
  arrange(biome_id, logDate) %>% 
  group_by(biome_id) %>% 
  mutate(
    logDate = as.Date(logDate),
    CST_previous = lag(CST),
    transition = CST != CST_previous,
    transition_start = ifelse(transition, logDate, NA)
  ) %>% 
  fill(transition_start, .direction = "down") %>%
  group_by(biome_id, transition_start) %>%
  mutate(
    transition_duration = as.numeric(difftime(max(logDate), min(logDate), units = "days"))
  ) %>%
  ungroup()

study_end <- max(longterm_cst_vaginal.microbial.menses.24$logDate)

long_term_transitions <- longterm_cst_vaginal.microbial.menses.24 %>%
  filter(transition == TRUE & (transition_duration >= 10 | (logDate >= study_end - 10)))

table(long_term_transitions$transition_duration)
length(unique(long_term_transitions$biome_id))

# Count transitions during and outside of menses
transition_counts <- longterm_cst_vaginal.microbial.menses.24 %>%
  group_by(menses_day) %>%
  summarise(count = n())

transition_counts

# test if transition more likely in menses
binom.test(x = transition_counts$count[transition_counts$menses_day == "menses"][1], 
           n = sum(transition_counts$count), 
           p = 0.5)

length(unique(longterm_cst_vaginal.microbial.menses.24$biome_id))

##########################################################################################

## Analyze Lactobacillus abundance throughout study

bacterial_ra <- transform_sample_counts(bacterial.data, function(x) x / sum(x))  
lactobacillus_phyloseq <- subset_taxa(bacterial_ra, Genus %in% c("Lactobacillus"))
lactobacillus_df_ra <- psmelt(lactobacillus_phyloseq)

lactobacillus_df_ra_subset <- lactobacillus_df_ra %>% 
  rename(LactobacillusRelativeAbundance=Abundance) %>% 
  select(c(SampleID, LactobacillusRelativeAbundance, biome_id, Kingdom, Phylum, Class, Order, Family, Genus, Species, Species_exact))

# Merge with sample metadata
meta_df <- as.data.frame(sample_data(bacterial.data))
lactobacillus_df <- left_join(lactobacillus_df_ra_subset, meta_df, by = c("SampleID", "biome_id"))
lactobacillus_df <- study_days(lactobacillus_df)

length(unique(lactobacillus_df$biome_id))

# Lactobacillus: Logged Relative Abundance of Lactobacillus Abundance Over Time for participants
ggplot(lactobacillus_df, aes(x = as.numeric(study_day), y = LactobacillusRelativeAbundance, group=biome_id, color=as.factor(biome_id))) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5)+ #, span=0.5) +  
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, max(lactobacillus_df$study_day), by = 1)) +
  # scale_y_continuous(breaks = seq(0, 1, by = 0.1))+
  theme_minimal() +
  # facet_wrap(~ Species_exact) + 
  theme(legend.position="none") +
  labs(x = "Study Day", y = "Log Relative Abundance of Lactobacillus",
       color="Participant ID",
       ylim=c(0,1))

# Join menses and vaginal microbiome data
vaginal.microbial.menses.24.subset <- vaginal.microbial.menses.24 %>% 
  select(biome_id, logDate, CST, SampleID, shannon, timestamp, menses_day)
lactobacillus.menses <- lactobacillus_df %>% 
  left_join(vaginal.microbial.menses.24.subset, by=c("SampleID", "biome_id", "logDate", "timestamp"))

library(lme4)
library(lmerTest)
library(performance)

# Add a small constant to avoid zero values
lactobacillus.menses$Lactobacillus_log <- log(lactobacillus.menses$LactobacillusRelativeAbundance + 1e-6)

lm.model <- lm(Lactobacillus_log~menses_day*study_day, data=lactobacillus.menses)
summary(lm.model)

lacto.menses.model <- lmer(Lactobacillus_log ~ menses_day + (1|biome_id), data = lactobacillus.menses)
r2(lacto.menses.model)
summary(lacto.menses.model)

lacto.menses.model <- lmer(Lactobacillus_log ~ menses_day + study_day + (1|biome_id), data = lactobacillus.menses)
r2(lacto.menses.model)

# Menses: Shannon Diversity Over Time by menses day
longterm_cst_vaginal.microbial.menses.24 <- study_days(longterm_cst_vaginal.microbial.menses.24)
longterm_cst_vaginal.microbial.menses.24 %>% 
  filter(!is.na(menses_day)) %>% 
  ggplot(aes(x = as.numeric(study_day), 
                                          y = shannon, 
                                          color=as.factor(menses_day))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) + 
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, max(lactobacillus_df$study_day), by = 5)) +
  theme_minimal() +
  scale_color_discrete(labels = c("Menstruating", "Not Menstruating")) + 
  theme(text=element_text(size=16))+
  labs(x = "Study Day", y = "Shannon Diversity", color="Menstruation Status") 

# Lactobacillus: Relative Abundance of Lactobacillus Abundance Over Time by menses day
lactobacillus.menses.filtered <- lactobacillus.menses %>% 
  filter(!is.na(menses_day))
table(lactobacillus.menses.filtered$menses_day)
ggplot(lactobacillus.menses.filtered, aes(x = as.numeric(study_day), 
                                          y = LactobacillusRelativeAbundance, 
             color=as.factor(menses_day))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) + 
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, max(lactobacillus_df$study_day), by = 5)) +
  theme_minimal() +
  scale_color_discrete(labels = c("Menstruating", "Not Menstruating")) + 
  theme(text=element_text(size=16))+
  labs(x = "Study Day", y = "Log Relative Abundance of Lactobacillus", color="Menstruation Status") 

# Figure: plot of average shannon for each participant on menses vs. not on menses
lactobacillus.menses.summary <- lactobacillus.menses %>% 
  group_by(biome_id, menses_day) %>% 
  summarise(avg_lacto = mean(LactobacillusRelativeAbundance, na.rm = TRUE)) %>%
  pivot_wider(names_from = menses_day, values_from = avg_lacto, names_prefix = "menses_day_") %>%
  filter(!is.na(menses_day_not_menses) & !is.na(menses_day_menses)) %>%
  select(!menses_day_NA)

### Wilcox test
wilcox.test(lactobacillus.menses.summary$menses_day_not_menses,
            lactobacillus.menses.summary$menses_day_menses, paired=TRUE)
t.test(lactobacillus.menses.summary$menses_day_not_menses,
       lactobacillus.menses.summary$menses_day_menses, paired=TRUE)

## Individual Analysis for genus level (case study: participant 61 and 65)

# Menses: participant 61 menses and not menses graph
lactobacillus.menses.filtered.61 <- lactobacillus.menses.filtered %>% 
  filter(biome_id == 61)
ggplot(lactobacillus.menses.filtered.61, aes(x = as.numeric(study_day), 
                                          y = LactobacillusRelativeAbundance, 
                                          color=as.factor(menses_day))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) + 
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, max(lactobacillus_df$study_day), by = 5)) +
  theme_minimal() +
  scale_color_discrete(labels = c("Menstruating", "Not Menstruating")) + 
  theme(text=element_text(size=16))+
  labs(x = "Study Day", y = "Log Relative Abundance of Lactobacillus", color="Menstruation Status") 
# participant 61: top lactobacillus species over time
ggplot(lactobacillus.menses.filtered.61, 
       aes(x = as.numeric(study_day), 
           y = LactobacillusRelativeAbundance, 
           color = OTU, 
           group = OTU)) +
  geom_line() +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, max(lactobacillus_df$study_day), by = 1)) +
  labs(title = "Relative Abundance of Lactobacillus Species Over Time (Participant 61)",
       x = "Study Day",
       y = "Relative Abundance",
       color = "Lactobacillus Species") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values=c("Lactobacillus crispatus" = "darkblue",
                     "Lactobacillus iners"="darkgreen")) +
  theme(legend.position = "right")

# Menses: participant 65 menses and not menses graph
lactobacillus.menses.filtered.65 <- lactobacillus.menses.filtered %>% 
  filter(biome_id == 65)
ggplot(lactobacillus.menses.filtered.65, aes(x = as.numeric(study_day), 
                                             y = LactobacillusRelativeAbundance, 
                                             color=as.factor(menses_day))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) + 
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, max(lactobacillus_df$study_day), by = 5)) +
  theme_minimal() +
  theme(text=element_text(size=16))+
  scale_color_discrete(labels = c("Menstruating", "Not Menstruating")) + 
  labs(x = "Study Day", y = "Log Relative Abundance of Lactobacillus", color="Menstruation Status") 

# participant 65: top lactobacillus species over time
ggplot(lactobacillus.menses.filtered.65, 
       aes(x = as.numeric(study_day), 
           y = LactobacillusRelativeAbundance, 
           color = OTU, 
           group = OTU)) +
  geom_line() +
  geom_point(alpha = 0.7) +
  scale_x_continuous(breaks = seq(0, max(lactobacillus_df$study_day), by = 1)) +
  theme_minimal() +
  labs(title = "Relative Abundance of Lactobacillus Species Over Time (Participant 61)",
       x = "Study Day",
       y = "Relative Abundance",
       color = "Lactobacillus Species") +
  scale_color_manual(values=c("Lactobacillus crispatus" = "darkblue",
                              "Lactobacillus iners"="darkgreen")) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = "right")

# Analyze Lactobacillus species
lactobacillus_species_df <- lactobacillus_df_ra_subset %>%
  group_by(Species_exact, biome_id) %>%
  summarise(mean_abundance = mean(LactobacillusRelativeAbundance, na.rm = TRUE),
            median_abundance = median(LactobacillusRelativeAbundance, na.rm = TRUE),
            n = n()) %>%
  ungroup()

# add menses
lactobacillus_df_ra_subset_menses <- lactobacillus_df_ra_subset %>% 
  left_join(longterm_cst_vaginal.microbial.menses.24) %>% 
  filter(!is.na(menses_day)) %>% 
  filter(Species_exact!="")
dim(lactobacillus_df_ra_subset_menses)[1]/dim(lactobacillus_df_ra_subset)[1]

##########################################################################################

## Changes in Genus level analysis in Menses

# relative abundance
bacterial.ra <- transform_sample_counts(bacterial.data, function(x) x / sum(x))  

# aggregate to genus level
genus.ra <- tax_glom(bacterial.ra, taxrank = "Genus")
genus.df <- psmelt(genus.ra) %>%
  rename(RelativeAbundance = Abundance)

dim(genus.df) # 611119 samples and otu combined

# merge with menses data
vaginal.microbial.menses.24.subset <- vaginal.microbial.menses.24 %>% 
  select(SampleID, CST, shannon, biome_id, logDate, sampleType, regular_periods, menses_status, menses_day)
genus.df.menses <- genus.df %>% 
  left_join(vaginal.microbial.menses.24.subset)

# collapse data to be paired
genus.df.menses.paired <- genus.df.menses %>% 
  filter(!is.na(menses_day)) %>% 
  group_by(biome_id, menses_day, Genus) %>% 
  summarise(avg_shannon = sum(shannon) / n(),
            avg_ra = sum(RelativeAbundance) / n(),
            .groups="drop") %>% 
  group_by(biome_id, Genus) %>% 
  filter(n_distinct(menses_day) == 2) %>% 
  ungroup() %>%
  # reshape df
  pivot_wider(names_from = menses_day, values_from = c(avg_ra, avg_shannon),
              names_sep="_") %>%
  drop_na()

# Menses: Log fold change - volcano plot
genus_stats <- genus.df.menses.paired %>%
  group_by(Genus) %>%
  summarise(
    # Log fold change
    fold_change = log2(mean(avg_ra_menses) / mean(avg_ra_not_menses)),  
    p_value = t.test(avg_ra_menses, avg_ra_not_menses, paired = TRUE)$p.value
  ) %>%
  # bonferonni correction
  mutate(
    p_adj = p.adjust(p_value, method = "bonferroni"),
    significance = ifelse(p_adj < 0.05, "Significant", "Not Significant")
  ) %>% 
  arrange(p_adj)

genus_stats.filtered <- genus_stats %>% 
  filter(!is.na(p_adj)) %>% 
  filter(fold_change != -Inf & fold_change != Inf) %>% 
  mutate(
    significance = ifelse(p_adj < 1, "Significant", "Not Significant"),
    color = case_when(
      significance == "Significant" & fold_change > 0 ~ "darkgray",    
      significance == "Significant" & fold_change < 0 ~ "blue",  
      TRUE ~ "gray" 
    )
  )

genus_stats.filtered <- genus_stats.filtered %>%
  mutate(hjust_text = ifelse(fold_change < 0, 1.1, -0.1))

# Menses: log fold change genus plot
ggplot(genus_stats.filtered, aes(x = fold_change, y = -log10(p_adj), color = color)) +
  geom_point(alpha = 0.7) +
  # geom_segment(data = subset(genus_stats.filtered, significance == "Significant"), 
               # aes(xend = fold_change, yend = -log10(p_value)), 
               # linetype = "dashed", color = "black", size = 0.5) + 
  geom_text(data = subset(genus_stats.filtered, significance == "Significant"), 
            aes(label = Genus,hjust = hjust_text), 
            vjust = -0.5, 
            # hjust = 1.2, 
            size = 5) +  
  geom_hline(yintercept = -log10(0.05), linetype = "solid", color = "black") + 
  geom_vline(xintercept = 0, linetype = "solid", color = "black") + 
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (Menses vs. Not Menses)",
    y = "Adjusted p-value (-log10)",
    title = " ",
    color = "Significance"
  ) +
  xlim(-10, 10) +
  scale_color_identity() +
  theme(text=element_text(size=18))

# top genus change
top_genera <- genus_stats %>%
  filter(p_value < 0.05) %>%
  arrange(p_value)

# Menses: boxplot of change in genus levels on and not on menses
genus.df.menses.paired %>%
  filter(Genus %in% top_genera$Genus) %>%
  pivot_longer(cols = starts_with("avg_ra_"), names_to = "Condition", values_to = "Abundance") %>%
  mutate(Condition = ifelse(Condition == "avg_ra_menses", "Menses", "Not Menses")) %>%
  ggplot(aes(x = Condition, y = Abundance, fill = Condition)) +
  geom_boxplot() +
  facet_wrap(~Genus, scales = "free_y") +
  theme_minimal() +
  labs(
    x = "Condition",
    y = "Relative Abundance",
    title = " "
  ) +
  scale_fill_manual(values = c("Menses" = "red", "Not Menses" = "blue"))

# heatmap of all genus up and down in menses or not
library(pheatmap)
genus.df.menses.paired.collapsed <- genus.df.menses.paired %>% 
  group_by(Genus) %>% 
  summarise(avg_ra_menses = sum(avg_ra_menses) / n(),
            avg_ra_not_menses = sum(avg_ra_not_menses) / n(),
            avg_shannon_menses = sum(avg_shannon_menses) / n(),
            avg_shannon_not_menses = sum(avg_shannon_not_menses) / n())
genus_matrix <- genus.df.menses.paired.collapsed %>%
  # filter(Genus %in% top_genera$Genus) %>%
  filter(avg_ra_menses > 0.001 & avg_ra_not_menses > 0.001) %>%
  select(Genus, avg_ra_menses, avg_ra_not_menses) %>%
  column_to_rownames("Genus") %>%
  as.matrix()

pheatmap(log10(genus_matrix + 1),  # Log10 transform for better visualization
         cluster_rows = TRUE, 
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Genus Relative Abundance During Menses")

##########################################################################################
library(ggrepel)

## Changes in Genus level analysis 

# read in stress and participant data: ES_2025_analysis.R line 718
shannon.dass.filtered <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/dass.participant.csv")
shannon.dass.filtered <- shannon.dass.filtered[,-1]
shannon.dass.filtered <- shannon.dass.filtered %>% 
  select(SampleID, CST, shannon, qr, biome_id, logDate, CST_max, stress_score, stressseverity, sport)

# Corr with volunteer history data
participant.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 9-Volunteer Medical History.csv", header = TRUE)
# filter non-hormonal & no samples
participant.data <- participant.data %>% 
  filter(birthControl!="Orilissa (Elagolix)")
# select birth control
birthControl.df <- participant.data %>% 
  dplyr::select(biome_id, birthControl)

# aggregate at genus level & relative abundances
genus.taxa <- tax_glom(bacterial.data, taxrank = "Genus")
genus.ra <- transform_sample_counts(genus.taxa, function(x) x / sum(x)) 
bactieral.meta <- as(sample_data(genus.ra), "data.frame")

# add dass data
bactieral.meta.dass <- bactieral.meta %>% 
  left_join(shannon.dass.filtered)
# add menses data
vaginal.microbial.menses.24.subset2 <- vaginal.microbial.menses.24.subset %>% 
  select(SampleID, biome_id, menses_status, menses_day)
bactieral.meta.dass.menses <- bactieral.meta.dass %>% 
  left_join(vaginal.microbial.menses.24.subset2, 
            by=c("SampleID", "biome_id"))
# add birth control data
bactieral.meta.dass.menses.bc <- bactieral.meta.dass.menses %>% 
  left_join(birthControl.df, by="biome_id")

# set sample data
rownames(bactieral.meta.dass.menses.bc) <- bactieral.meta.dass.menses.bc$SampleID
sample_data(genus.ra) <- bactieral.meta.dass.menses.bc
names(bactieral.meta.dass.menses.bc)

# ordination
genus.ordination <- ordinate(genus.ra, method = "PCoA", distance = "bray")

# vaginal: ordination cluster by CST
plot_ordination(genus.ra, genus.ordination, color = "CST") + 
  geom_point(size=1, alpha=0.6) +
  # geom_jitter(aes(color=as.factor(CST)), size=1, alpha=0.6) +
  labs(color="CST")+
  scale_color_viridis(discrete=TRUE)+
  geom_text_repel(aes(label = biome_id), size = 3,
                  show.legend=FALSE) +
  theme_minimal()

# vaginal: ordination cluster by birth control
plot_ordination(genus.ra, genus.ordination, color = "birthControl") + 
  geom_point(size=1, alpha=0.6) +
  labs(color="birthControl")+
  scale_color_viridis(discrete=TRUE)+
  geom_text_repel(aes(label = biome_id), size = 3,
                  show.legend=FALSE) +
  theme_minimal()

# Figure: plot of average shannon for each participant on birth controls
lactobacillus.menses.bc <- lactobacillus.menses %>% 
  left_join(birthControl.df, by="biome_id") %>% 
  filter(!is.na(birthControl))

# relevel
lactobacillus.menses.bc$birthControl <- 
  factor(lactobacillus.menses.bc$birthControl, 
         levels=c("None", "Local P", "Systemic P only", "Systemic Combined (E&P)"))

# Birth Control: boxplot of Lactobacillus log relative abundance
ggplot(lactobacillus.menses.bc, aes(x = birthControl, y = LactobacillusRelativeAbundance, fill = birthControl)) +
  geom_jitter(aes(color=as.factor(biome_id)), size=1, alpha=0.6,
              show.legend = FALSE) +
  geom_boxplot(alpha=0.5) +
  scale_y_log10() +
  # facet_wrap(~Genus, scales = "free_y") +
  theme_minimal() +
  labs(
    x = "Birth Control",
    y = "Lactobacillus Log Relative Abundance",
    title = " ",
    fill = "Birth Control"
  ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        text=element_text(size=18))
################################################################################ 

######### Merge dataframes
dim(longterm_cst_vaginal.microbial.menses.24)
longterm_cst_vaginal.microbial.menses.24.subset <- longterm_cst_vaginal.microbial.menses.24 %>% 
  select(CST, shannon, biome_id, logDate, OTU, menses_day, SampleID) %>% 
  group_by(biome_id, logDate) %>% 
  summarise(CST=first(CST),
            shannon=first(shannon),
            OTU=first(OTU),
            menses_day=first(menses_day),
            SampleID = first(SampleID))
dim(longterm_cst_vaginal.microbial.menses.24.subset) # filter out repeat days

shannon.dass <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/shannon.dass.csv")
shannon.dass.subset <- shannon.dass %>% 
  select(biome_id, logDate, birthControl, stress_score) %>% 
  mutate(logDate=as.Date(logDate)) %>% 
  group_by(biome_id, logDate) %>% 
  summarise(birthControl = first(birthControl),
            stress_score = first(stress_score))
dim(shannon.dass.subset)
colSums(is.na(shannon.dass.subset)) # 115 stress scores missing

vaginal.merged.df <- longterm_cst_vaginal.microbial.menses.24.subset %>% 
  left_join(shannon.dass.subset, by=c("biome_id", "logDate")) 
colSums(is.na(vaginal.merged.df)) # 115 stress scores missing

nutrient.diet <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/nutrient.diet.csv")
nutrient.diet.subset <- nutrient.diet %>% 
  select(c(biome_id, logDate, caloriesall_avg, cholesterol_prop, satFat_prop, sodium_prop,
           carb_prop, dietFib_prop, sugar_prop, protein_prop, fat_prop, fat_cal_prop,
           addedSugarall_prop)) %>% 
  mutate(logDate=as.Date(logDate)) %>% 
  group_by(biome_id, logDate) %>% 
  summarise(caloriesall_avg = first(caloriesall_avg),
            cholesterol_prop = first(cholesterol_prop),
            satFat_prop = first(satFat_prop),
            sodium_prop = first(sodium_prop),
            carb_prop = first(carb_prop),
            dietFib_prop = first(dietFib_prop),
            sugar_prop = first(sugar_prop),
            protein_prop = first(protein_prop),
            fat_prop = first(fat_prop),
            fat_cal_prop = first(fat_cal_prop),
            addedSugarall_prop = first(addedSugarall_prop)
  )
vaginal.merged.df <- vaginal.merged.df %>% 
  left_join(nutrient.diet.subset)

colSums(is.na(vaginal.merged.df))

activity.shannon <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/activity.data.shannon.csv")
activity.shannon.subset <- activity.shannon %>% 
  select(c(biome_id, logDate, calories_burned, steps, distance, minutes_sedentary,
           minutes_lightly_active, minutes_fairly_active, minues_very_active, activity_calories,
           survey_menstruate)) %>% 
  mutate(logDate=as.Date(logDate)) %>% 
  group_by(biome_id, logDate) %>% 
  summarise(calories_burned = first(calories_burned),
            steps = first(steps),
            distance = first(distance),
            minutes_sedentary = first(minutes_sedentary),
            minutes_lightly_active = first(minutes_lightly_active),
            minutes_fairly_active = first(minutes_fairly_active),
            minues_very_active = first(minues_very_active),
            activity_calories = first(activity_calories))

vaginal.merged.df <- vaginal.merged.df %>% 
  left_join(activity.shannon.subset, by=c("biome_id", "logDate"))

participant.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 9-Volunteer Medical History.csv", header = TRUE)
participant.data.subset <- participant.data %>% 
  select(biome_id, cisWoman, sport, probiotic, birthControl, study_menstruate,
         sexuallyActive)
vaginal.merged.df <- vaginal.merged.df %>% 
  left_join(participant.data.subset)

# 04/13 uncomment to save - full vaginal df
# write.csv(vaginal.merged.df, "/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/vaginal.lifestyle.csv")

# Aggregate model
vaginal.merged.df.summary <- vaginal.merged.df %>% 
  group_by(biome_id) %>%
  summarise(avg_shannon=sum(shannon, na.rm=TRUE) / n(),
            avg_stress_score=sum(stress_score, na.rm=TRUE) / n(),
            cisWoman = first(cisWoman),
            sport = first(sport),
            probiotic = first(probiotic),
            birthControl = first(birthControl),
            sexuallyActive = first(sexuallyActive),
            caloriesall_avg = sum(caloriesall_avg, na.rm=TRUE) / n(),
            cholesterol_prop = sum(cholesterol_prop, na.rm=TRUE) / n(),
            satFat_prop = sum(satFat_prop, na.rm=TRUE) / n(),
            sodium_prop = sum(sodium_prop, na.rm=TRUE) / n(),
            carb_prop = sum(carb_prop, na.rm=TRUE) / n(),
            dietFib_prop = sum(dietFib_prop, na.rm=TRUE) / n(),
            sugar_prop = sum(sugar_prop, na.rm=TRUE) / n(),
            protein_prop = sum(protein_prop, na.rm=TRUE) / n(),
            fat_prop = sum(fat_prop, na.rm=TRUE) / n(),
            fat_cal_prop = sum(fat_cal_prop, na.rm=TRUE) / n(),
            addedSugarall_prop = sum(addedSugarall_prop, na.rm=TRUE) / n(),
            study_menstruate = first(study_menstruate),
            calories_burned = sum(calories_burned, na.rm=TRUE) / n(),
            steps = sum(steps, na.rm=TRUE) / n(),
            distance = sum(distance, na.rm=TRUE) / n(),
            minutes_sedentary = sum(minutes_sedentary, na.rm=TRUE) / n(),
            minutes_lightly_active = sum(minutes_lightly_active, na.rm=TRUE) / n(),
            minutes_fairly_active = sum(minutes_fairly_active, na.rm=TRUE) / n(),
            minues_very_active = sum(minues_very_active, na.rm=TRUE) / n(),
            activity_calories = sum(activity_calories, na.rm=TRUE) / n()) %>% 
  mutate(stress_severity = case_when(
    avg_stress_score <= 14 ~ "Normal",
    avg_stress_score >= 15 & avg_stress_score <= 18 ~ "Mild",
    avg_stress_score >= 19 & avg_stress_score <= 25 ~ "Moderate",
    avg_stress_score >= 26 & avg_stress_score <= 33 ~ "Severe",
    avg_stress_score >= 34 ~ "Extremely Severe",
    TRUE ~ NA_character_
  )) %>% 
  select(-avg_stress_score) 

colSums(is.na(vaginal.merged.df.summary))

vaginal.merged.df.summary.filtered <- 
  vaginal.merged.df.summary[complete.cases(vaginal.merged.df.summary),]

lm.obj.orig <- lm(avg_shannon~., data=vaginal.merged.df.summary.filtered[,-1])
summary(lm.obj.orig)

# menses
lm.obj <- lm(avg_shannon~., data=vaginal.merged.df.summary.filtered[,-c(1,19)])
summary(lm.obj)
# 0.2864
anova(lm.obj.orig, lm.obj)

# stress
lm.obj <- lm(avg_shannon~., data=vaginal.merged.df.summary.filtered[,-c(1, 28)])
summary(lm.obj)

anova(lm.obj.orig, lm.obj)

# HBC
lm.obj <- lm(avg_shannon~., data=vaginal.merged.df.summary.filtered[,-c(1, 6)])
summary(lm.obj)

anova(lm.obj.orig, lm.obj)

# nutrition
lm.obj <- lm(avg_shannon~., data=vaginal.merged.df.summary.filtered[,-c(1, 8:18)])
summary(lm.obj)
# 0.2864

anova(lm.obj.orig, lm.obj)

# activity
lm.obj <- lm(avg_shannon ~ ., 
             data=vaginal.merged.df.summary.filtered[,-c(1, 20:27)])
summary(lm.obj)
# 0.2864

anova(lm.obj.orig, lm.obj)

# sport
lm.obj <- lm(avg_shannon~., data=vaginal.merged.df.summary.filtered[,-c(1, 4)])
summary(lm.obj)
# 0.2864
anova(lm.obj.orig, lm.obj)

# probiotic
lm.obj <- lm(avg_shannon~., data=vaginal.merged.df.summary.filtered[,-c(1, 5)])
summary(lm.obj)
# 0.2864
anova(lm.obj.orig, lm.obj)

# sex act
lm.obj <- lm(avg_shannon~., data=vaginal.merged.df.summary.filtered[,-c(1, 7)])
summary(lm.obj)
# 0.2864
anova(lm.obj.orig, lm.obj)

# gender - cisWoman
lm.obj <- lm(avg_shannon~., data=vaginal.merged.df.summary.filtered[,-c(1, 3)])
summary(lm.obj)
# 0.2864
anova(lm.obj.orig, lm.obj)

## full model with time

vaginal.merged.df2 <- vaginal.merged.df %>% 
  mutate(stress_severity = case_when(
    stress_score <= 14 ~ "Normal",
    stress_score >= 15 & stress_score <= 18 ~ "Mild",
    stress_score >= 19 & stress_score <= 25 ~ "Moderate",
    stress_score >= 26 & stress_score <= 33 ~ "Severe",
    stress_score >= 34 ~ "Extremely Severe",
    TRUE ~ NA_character_
  )) %>% 
  select(-stress_score) 

colSums(is.na(vaginal.merged.df2))
lmer.obj <- lmer(shannon~. - biome_id +(1|`biome_id`), data=vaginal.merged.df2[,-c(2,7)])
r2(lmer.obj)
summary(lmer.obj)





