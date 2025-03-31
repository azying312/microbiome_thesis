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
write.csv(vaginal.microbial.menses.24, "/Volumes/T7/microbiome_data/cleaned_data/microbiome_lifestyle/relabeled_data/vaginal.microbial.menses.24.csv")

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
  
  shannon_plt <- ggplot(vaginal.microbial.menses.24.participant, 
                        aes(x=as.factor(logDate), y=(shannon), color=as.factor(menses_day), 
                            shape=as.factor(menses_day))) +
    geom_point() +
    scale_color_manual(name = "Menses Status",
                       values = c("not_menses" = "black", "menses" = "red", "not_recorded" = "purple")) +  
    scale_shape_manual(name = "Menses Status",
                       values = c("not_menses" = 16, "menses" = 17, "not_recorded"=18)) +
    ylim(0,3)+
    # Add lines for the days with multiple samples
    geom_segment(data = dupe.day, 
                 aes(x = logDate, xend = logDate, y = -Inf, yend = Inf), 
                 linetype = "dashed", color = "gray50", inherit.aes = FALSE) +  
    theme_minimal() +
    labs(x = "Date", y = "Shannon Index", title=paste(id, "Shannon diversity days"),
         ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
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

# Figure: Average Shannon Diversity for Participants on menses v. not
ggplot(
  vaginal.microbial.menses.24.summary,
  aes(x = menses_day_not_menses, y = menses_day_menses,
      color=as.factor(biome_id))) +
  geom_point() +
  xlim(0,3)+
  ylim(0,3)+
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  theme_minimal() +
  labs(x="Not Menses (Shannon)", y="Menses (Shannon)",
       color="Participant ID")

##########################################################################################

length(unique(vaginal.microbial.menses.24$biome_id))

# Long term CST transition - participant shifts from one CST to another and stays in the new CST for 10 days or more
longterm_cst_vaginal.microbial.menses.24 <- vaginal.microbial.menses.24 %>% 
  dplyr::select(SampleID, CST, shannon, qr, biome_id, logDate, sampleType, max_taxa, OTU, menses_day) %>% 
  # filter(!is.na(menses_day)) %>% # only if you need menses data
  arrange(biome_id, logDate) %>% 
  group_by(biome_id) %>% 
  mutate(CST_previous = lag(CST),   
         transition = CST != CST_previous)

length(unique(longterm_cst_vaginal.microbial.menses.24$biome_id))

longterm_cst_vaginal.microbial.menses.24 <- longterm_cst_vaginal.microbial.menses.24 %>%
  mutate(logDate = as.Date(logDate)) %>% 
  group_by(biome_id, CST) %>%
  mutate(duration = as.numeric(difftime(max(logDate), min(logDate), units = "days"))) %>%
  ungroup()

table(longterm_cst_vaginal.microbial.menses.24$duration)

## Max gap in samples
max_gaps <- longterm_cst_vaginal.microbial.menses.24 %>%
  arrange(biome_id, logDate) %>%
  group_by(biome_id) %>%
  mutate(day_gap = as.numeric(difftime(logDate, lag(logDate), units = "days"))) %>%
  summarise(max_gap = max(day_gap, na.rm = TRUE))
table(max_gaps$max_gap)

# 10 day length
long_term_transitions <- longterm_cst_vaginal.microbial.menses.24 %>%
  filter(transition == TRUE, duration >= 10)

length(unique(longterm_cst_vaginal.microbial.menses.24$biome_id))
length(unique(long_term_transitions$biome_id))

# Count transitions during and outside of menses
transition_counts <- long_term_transitions %>%
  group_by(menses_day) %>%
  summarise(count = n())

# test if transition more likely in menses
binom.test(x = transition_counts$count[transition_counts$menses_day == "menses"][1], 
           n = sum(transition_counts$count), 
           p = 0.5)

length(unique(long_term_transitions$biome_id))

##########################################################################################

## Analyze Lactobacillus abundance throughout study

# bacteria_otu_table <- as.data.frame(otu_table(bacterial.data))
# bacteria_taxa_table <- as.data.frame(tax_table(bacterial.data))
# otu_taxa_table <- cbind(bacteria_otu_table, bacteria_taxa_table)
# # get Lactobacillus abundances
# lactobacillus_df <- otu_taxa_table %>% filter(Genus == "Lactobacillus")
# relative abundance
# lactobacillus_df_ra <- sweep(lactobacillus_df[, 1:ncol(bacteria_otu_table)], 2, colSums(bacteria_otu_table), FUN = "/")

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

# Figure: Logged Relative Abundance of Lactobacillus Abundance Over Time
ggplot(lactobacillus_df, aes(x = as.numeric(study_day), y = LactobacillusRelativeAbundance, group=biome_id, color=as.factor(biome_id))) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) +  
  # geom_smooth(method = "loess", span = 0.2, se = FALSE, color = "orchid") +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, max(lactobacillus_df$study_day), by = 1)) +
  theme_minimal() +
  labs(x = "Study Day", y = "Logged Relative Abundance of Lactobacillus",
       color="Participant ID")

# Join menses and vaginal microbiome data
vaginal.microbial.menses.24.subset <- vaginal.microbial.menses.24 %>% 
  select(biome_id, logDate, CST, SampleID, shannon, timestamp, menses_day)
lactobacillus.menses <- lactobacillus_df %>% 
  left_join(vaginal.microbial.menses.24.subset, by=c("SampleID", "biome_id", "logDate", "timestamp"))

library(lme4)
library(lmerTest)
library(performance)

lm.model <- lm(Lactobacillus_log~menses_day*study_day, data=lactobacillus.menses)
summary(lm.model)

# Add a small constant to avoid zero values
lactobacillus.menses$Lactobacillus_log <- log(lactobacillus.menses$LactobacillusRelativeAbundance + 1e-6)

lacto.menses.model <- lmer(Lactobacillus_log ~ menses_day + (1|biome_id), data = lactobacillus.menses)
r2(lacto.menses.model)
summary(lacto.menses.model)

lacto.menses.model <- lmer(Lactobacillus_log ~ menses_day + study_day + (1|biome_id), data = lactobacillus.menses)
r2(lacto.menses.model)

# Figure: Relative Abundance of Lactobacillus Abundance Over Time by menses day
lactobacillus.menses %>% 
  filter(!is.na(menses_day)) %>% 
  ggplot(aes(x = as.numeric(study_day), y = LactobacillusRelativeAbundance, group=menses_day, color=as.factor(menses_day))) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "loess", se = FALSE, size = 0.5) + 
    scale_x_continuous(breaks = seq(0, max(lactobacillus_df$study_day), by = 1)) +
    theme_minimal() +
    scale_y_log10() +
    labs(x = "Study Day", y = "Logged Relative Abundance of Lactobacillus",
         color="Participant ID")


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

##########################################################################################

