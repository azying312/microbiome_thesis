library(dplyr)
library(phyloseq)
library(vegan)
library(pheatmap)
library(tidyverse)
library(Matrix)

library(ggplot2)

library(cluster)
library("igraph")
library("markovchain")
library("RColorBrewer")
library("gridExtra")

# Use corrected data after
# bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/vaginal_bacteria_cleaned.rds")
# Old data - ampliseq
# bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/old_data/vaginal_bacteria_cleaned.rds")

###############################################################################################

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
  select(!SampleID)
# bray_curtis_dist <- vegdist(t(relative_abundance_otu), method = "bray")
# specnumber.ra <- specnumber(relative_abundance_otu_t)

###### Five Community State Type (CSTs) Clustering

# species to CST mapping
species_to_cst <- data.frame(
  Species = c("crispatus", "gasseri", "iners", "jensenii"),
  CST = c("I", "II", "III", "V") # IV is anaerobic/diverse cluster
)

# Get most abundant OTU per sample
max_taxa <- apply(relative_abundance_otu_t, 2, function(sample) {
  taxa_idx <- which.max(sample)
  taxa_names(vaginal_relative_abundances)[taxa_idx]
})

# Map most abundant OTU to the sample data
bacteria_metadata_df$max_taxa <- max_taxa

###################################################################################################

### Get Species variable cleaned in new data

bacteria_taxa_table <- tax_table(bacterial.data)
bacteria_taxa_df <- as.data.frame(bacteria_taxa_table)

dim(bacteria_taxa_df) # 92499    10

length(bacteria_taxa_df$Species) # 92499
sum(bacteria_taxa_df$Species=="") # 40397
sum(bacteria_taxa_df$Species_exact=="") # 90484

bacteria_taxa_df <- bacteria_taxa_df %>% 
  mutate(Species_exact=ifelse(Species_exact=="", Species, Species_exact))
sum(bacteria_taxa_df$Species_exact=="") # 39845
sum(bacteria_taxa_df$Species_exact!="") # 52654

# Subset taxa; Filter ASVs w/o Species_exact assignment
bacterial.data_subset <- subset_taxa(
  bacterial.data,
  Species_exact != ""
)

bacteria_taxa_table <- tax_table(bacterial.data_subset)
bacteria_taxa_df <- as.data.frame(bacteria_taxa_table)
dim(bacteria_taxa_df) # 2015 OTU <-- 52654 OTUs  10 cols

# Get max taxa names
bacteria_metadata_df$OTU <- as.character(bacteria_taxa_df[bacteria_metadata_df$max_taxa, "Species_exact"]) # 2039 samples

# Set sample data in phyloseq obj
sample_data(bacterial.data) <- bacteria_metadata_df
bacteria_taxa_table <- as.data.frame(bacteria_taxa_table)

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
## Other extreme
# bacteria_metadata_df <- bacteria_metadata_df %>%
#   mutate(
#     CST = case_when(
#       str_detect(OTU, "^crispatus$") ~ "I",
#       str_detect(OTU, "^gasseri$") ~ "II",
#       str_detect(OTU, "^iners$") ~ "III",
#       str_detect(OTU, "^jensenii$") ~ "V",
#       TRUE ~ "IV" # Assign "IV" for diverse/anaerobic or unclassified species
#     )
#   )

########################################################
# PNAS Clustering

# bacterial.data_subset
ps <- bacterial.data_subset
tt <- data.frame(tax_table(bacterial.data_subset))
# ps <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))

# Relative abundances
ps <- transform_sample_counts(ps, function(x) x/sum(x))

# MDS (or PCoA) ordination
braydist <- phyloseq::distance(ps, method="bray")
ord = ordinate(ps, method = "MDS", distance = braydist)
plot_scree(ord) + xlim(as.character(seq(1,12))) + ggtitle("MDS-bray ordination eigenvalues")

evs <- ord$value$Eigenvalues
print(evs[1:20])
print(tail(evs))

# Denoise distance matrix
h_sub5 <- hist(evs[6:length(evs)], 100)
plot(h_sub5$mids, h_sub5$count, log="y", type='h', lwd=10, lend=2)

# Determine number of clusters
NDIM <- 10
x <- ord$vectors[,1:NDIM]  # rows=sample, cols=MDS axes, entries = value
pamPCoA = function(x, k) {
  list(cluster = pam(x[,1:2], k, cluster.only = TRUE))
}
gs = clusGap(x, FUN = pamPCoA, K.max = 12, B = 50)
plot_clusgap(gs) + scale_x_continuous(breaks=c(seq(0, 12, 2)))

# Cluster into CSTs
K <- 5
x <- ord$vectors[,1:NDIM]
clust <- as.factor(pam(x, k=K, cluster.only=T))

# SWAPPING THE ASSIGNMENT OF 2 AND 3 TO MATCH RAVEL CST ENUMERATION
clust[clust==2] <- NA
clust[clust==3] <- 2
clust[is.na(clust)] <- 3
sample_data(ps)$CST <- clust
CSTs <- as.character(seq(K))

# Evaluate clustering
CSTColors <- brewer.pal(6,"Paired")[c(1,3,2,5,4,6)] # Length 6 for consistency with pre-revision CST+ coloration
names(CSTColors) <- CSTs
CSTColorScale <- scale_colour_manual(name = "CST", values = CSTColors[1:5])
CSTFillScale <- scale_fill_manual(name = "CST", values = CSTColors[1:5])
p1 <- plot_ordination(ps, ord, color = "CST") + CSTColorScale
p2 <- plot_ordination(ps, ord, axes = c(3, 4), color = "CST") + CSTColorScale

grid.arrange(p1, p2, top ="Ordination by Cluster")

# Graph to show clusters
plot_ordination(ps, ordinate(ps, method="NMDS", distance=braydist), color="CST") + CSTColorScale + ggtitle("NMDS -- bray -- By Cluster")

# Heatmaps of Clusters - Don't seem right...
taxa.order <- names(sort(taxa_sums(ps)))
for(CST in CSTs) {
  pshm <- prune_taxa(names(sort(taxa_sums(ps), T))[1:25], ps)
  pshm <- prune_samples(sample_data(pshm)$CST == CST, pshm)
  print(plot_heatmap(pshm, taxa.label="Species", taxa.order=taxa.order) + ggtitle(paste("CST:", CST)))
}

########################################################

bacteria_metadata_df <- sample_data(bacterial.data_subset)
# bacteria_metadata_df <- as.data.frame(as(otu_table(bacterial.data_subset), "matrix"))
otu_table_df <- as(otu_table(bacterial.data_subset), "matrix")

# Aggregate data by participants by mean relative abundance for a given OTU
participant_otu <- tapply(sample_names(bacteria_metadata_df), 
                          sample_data(bacterial.data_subset)$biome_id, 
                          function(samples) rowMeans(t(otu_table(vaginal_relative_abundances))[, samples, drop = FALSE]))
participant_otu <- do.call(cbind, participant_otu)
rownames(participant_otu) <- taxa_names(vaginal_relative_abundances)

## Alpha Div - Shannon Index
shannon.24 <- vegan::diversity(otu_table_df, "shannon")

# Add participant IDs from sample data | Merge the calculated Shannon diversity values with metadata
bacteria_metadata_df <- as(bacteria_metadata_df, "data.frame")
shannon.df.24 <- data.frame("SampleID"=names(shannon.24), "shannon"=shannon.24)
shannon.qr.merged.24 <- merge(shannon.df.24, bacteria_metadata_df, by="SampleID")
shannon.cst.qr.merged.24 <- merge(sample.cst, shannon.qr.merged.24, by="SampleID") %>% 
  mutate(biome_id=as.integer(biome_id))

# Save R environment
save.image("/Volumes/T7/microbiome_data/R_environments/vaginal_microbiome_relAbundance.RData")

# Define unique user IDs and labels
# unique.uid.24 <- sort(unique(shannon.qr.merged.24$biome_id))
# alph <- LETTERS[1:length(unique.uid.24)] 

# par(mfrow = c(2, 2), mai = c(.6, .6, .3, .3))

# Loop through each user and plot diversity trends
# i <- 1
# for (i.uid in unique.uid.24[1:4]) {
#   # Extract data for the current user
#   current.div <- shannon.qr.merged.24$shannon[shannon.qr.merged.24$biome_id == i.uid]
#   current.dates <- shannon.qr.merged.24$logDate[shannon.qr.merged.24$biome_id == i.uid]
#   
#   # Remove NA values from both dates and diversity scores
#   valid_indices <- !is.na(current.dates) & !is.na(current.div)
#   current.dates <- current.dates[valid_indices]
#   current.div <- current.div[valid_indices]
#   
#   # Plot diversity trend
#   plot(current.dates, current.div,
#        xlim = range(as.Date(shannon.qr.merged.24$logDate)),
#        ylim = c(0, max(shannon.qr.merged.24$shannon, na.rm = TRUE)),
#        xlab = "Date",
#        ylab = "Shannon Index",
#        main = paste("User", i.uid),
#        pch = 16, col = "blue")
#   
#   i <- i + 1
# }

load("/Volumes/T7/microbiome_data/R_environments/vaginal_microbiome_relAbundance.RData")

########################################################
## Corr with diet
diet.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/fully_merged_diet_data.csv", header = TRUE)

diet.data.collapsed <- diet.data %>% 
  group_by(study_id, Date) %>% 
  summarise(total_cal=sum(caloriesall)) %>% 
  rename(biome_id=study_id,
         logDate=Date)

ggplot(diet.data.collapsed, aes(x=total_cal, color=factor(biome_id))) +
  geom_histogram(binwidth=100)

# remove outliers with IQR
diet.data.collapsed.filtered <- diet.data.collapsed %>% 
  mutate(Q1=quantile(total_cal, 0.25),
         Q3=quantile(total_cal, 0.75),
         IQR=Q3-Q1) %>%
  filter(total_cal >= (Q1 - 1.5 * IQR) & total_cal <= (Q3 + 1.5 * IQR)) %>%
  select(-Q1, -Q3, -IQR)

ggplot(diet.data.collapsed.filtered, aes(x=total_cal, color=factor(biome_id))) +
  geom_histogram(binwidth=100)

shannon.diet <- shannon.cst.qr.merged.24 %>% 
  left_join(diet.data.collapsed.filtered, by=c("biome_id", "logDate")) %>%
  filter(!is.na(total_cal)) # 2039 -> 1196 = 843 no diet data | before IQR remove

ggplot(shannon.diet, aes(x = total_cal, y = shannon, color=factor(biome_id))) +
  geom_point() +
  geom_smooth(method = "lm", col = "blue") +
  labs(
    title = "",
    x = "Total Calories",
    y = "Shannon Diversity Index"
  ) +
  theme_minimal()

# shannon.diet %>% ggplot(aes(x=CST, y=total_cal) )+
#   geom_bar(stat = "identity", fill="orchid", color="black") +
#   theme_minimal() +
#   labs(
#     x = "Total Cals",
#     y = "Frequency",
#     title = "Frequency of Birth Control Methods"
#   )

cor.test(shannon.diet$shannon, shannon.diet$total_cal, method = "spearman")

# group_by biome_id
participant.spearman <- shannon.diet %>% 
  group_by(biome_id) %>% 
  summarise(spearman=cor(shannon, total_cal, method="spearman"))

hist(x=participant.spearman$spearman)

########################################################
## Corr with volunteer history data
participant.data <- read.csv("/Volumes/T7/microbiome_data/cleaned_data/cleaned_Report 9-Volunteer Medical History.csv", header = TRUE)

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

# Collapse by person (assign to most frequent CST)
shannon.birthControl.collapsed <- shannon.birthControl %>%
  group_by(biome_id) %>% 
  count(CST, name="frequency") %>% # total 148 (participants have more than 1 CST) 
  slice_max(frequency, n=1) %>% # collapse to their most frequent CST
  rename(CST_max=CST) %>% 
  ungroup() %>%
  left_join(select(shannon.birthControl, biome_id, birthControl, shannon, avg_shannon, CST), by = "biome_id") %>% 
  distinct(biome_id, .keep_all = TRUE)

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
shapiro.test(shannon.birthControl$shannon)
# test if there are significant diffs across categories
kruskal_test <- kruskal.test(shannon ~ birthControl, data = shannon.birthControl)
kruskal_test

# Association between shannon diversity and birth control
library(FSA)
dunn_result <- dunnTest(shannon ~ birthControl, data = shannon.birthControl, method = "bonferroni")
print(dunn_result)

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

# Figure: boxplot of shannon diversity by cluster
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
# dass.participant.binary <- dass.participant %>% 
#   mutate(birthControl.binary)
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



# lm.results <- participant.dass %>%
#   group_by(survey_menstruate) %>%
#   summarise(
#     model = list(lm(shannon ~ logDate, data = cur_data())),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     slope = sapply(model, function(m) coef(m)["logDate"]),
#     p_value = sapply(model, function(m) summary(m)$coefficients["logDate", "Pr(>|t|)"])
#   ) %>%
#   select(survey_menstruate, slope, p_value)

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

