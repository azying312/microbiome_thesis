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
bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/vaginal_bacteria_cleaned.rds")
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
# 
# table(bacteria_metadata_df$CST)
# cst_summary <- bacteria_metadata_df %>% 
#   count(CST) %>% 
#   mutate(Percentage=100*(n/sum(n)))
# print(cst_summary)

# Expand meta data for CST assignment
# match_species <- species_to_cst$Species

# bacteria_metadata_df_long <- bacteria_metadata_df %>% 
#   mutate(separate_flag = grepl(paste(match_species, collapse = "|"), OTU))  %>%
#   # Separate rows only for flagged rows
#   separate_rows(OTU, sep = "/") %>%
#   filter(separate_flag | (!separate_flag & !grepl("/", OTU))) %>%
#   mutate(CST_species=OTU %in% match_species) %>% 
#   group_by(across(-OTU)) %>%
#   # Separate CST I, II, III, V species
#   reframe(
#     OTU = c(
#       OTU[CST_species],                     
#       paste(OTU[!CST_species], collapse = "/") # IV Species
#     )) %>% 
#   filter(!(OTU=="")) %>% 
#   select(-c(separate_flag, CST_species))

# bacteria_metadata_df_long$CST <- ifelse(
#   bacteria_metadata_df_long$OTU %in% species_to_cst$Species,
#   species_to_cst$CST[match(bacteria_metadata_df_long$OTU, species_to_cst$Species)],
#   # Assign to CST IV if not in I, II, III, V
#   "IV"  
# )
# 
# table(bacteria_metadata_df_long$CST)
# table(bacteria_metadata_df_long$OTU)

## Check repeats
# freq_table <- table(bacteria_metadata_df_long$SampleID)
# freq_df <- as.data.frame(freq_table)
# colnames(freq_df) <- c("SampleID", "Freq")
# summary(freq_df$Freq)
# 
# freq_table <- table(bacteria_metadata_df_long$max_taxa)
# freq_df <- as.data.frame(freq_table)
# colnames(freq_df) <- c("OTU", "Freq")
# summary(freq_df$Freq)
# 
# freq_table <- table(bacteria_metadata_df$max_taxa)
# freq_df <- as.data.frame(freq_table)
# colnames(freq_df) <- c("OTU", "Freq")
# summary(freq_df$Freq)

# bacteria_metadata_df_long_filter <- bacteria_metadata_df_long %>% 
#   filter(max_taxa=="4554fad0a0d39e7fe9a9d3ed5686ec7e")

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
# 
# # Heatmap helper functions
# make_hcb <- function(data, var, name = NULL, fillScale = NULL, ...) {
#   hcb <- ggplot(data=data, aes_string(x="index", y=1, fill=var)) + 
#     geom_raster() +
#     scale_y_continuous(expand=c(0,0), breaks=1, labels=name) + 
#     scale_x_continuous(expand=c(0,0)) +
#     xlab(NULL) + ylab(NULL) +
#     theme(axis.title=element_blank(), axis.ticks=element_blank()) +
#     theme(axis.text.x=element_blank()) +
#     theme(axis.text.y=element_text(size=8, face="bold")) +
#     theme(plot.margin=unit(c(0,0,0,0),"lines"), 
#           axis.ticks.margin = unit(0,"null"), ...) +
#     guides(fill=F)
#   if(!is.null(fillScale)) hcb <- hcb + fillScale
#   return(hcb)
# }
# 
# plot_heatmap.2 <- function(ps, sample.label=NULL, taxa.label=NULL, ...) {
#   hm <- plot_heatmap(ps, taxa.label="Species", sample.order=sample.order, taxa.order = taxa.order)
#   hm <- hm + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
#   low = "#000033"; high = "#66CCFF"; trans = scales::log_trans(4); na.value = "black" # From plot_heatmap defaults
#   new_gradient <- scale_fill_gradient(low = low, high = high, trans = trans, na.value = na.value, breaks = c(0.001, 0.01, 0.1, 1), name="Relative\nabundance")
#   hm <- hm + theme(plot.margin=unit(c(0,0.5,0.5,0.5),"lines"))
#   hm <- hm + new_gradient
#   hm <- hm + geom_raster() #
#   hm <- hm + ylab("Taxa")
#   hm$layers <- hm$layers[2] #
#   return(hm)
# }
# 
# mush <- function(hmap, hcbs) {
#   cbgs <- lapply(hcbs, ggplotGrob)
#   hmg <- ggplotGrob(hmap)
#   # Make sure both plots have the same width in our final output
#   cbWidths <- lapply(cbgs, function(x) x$widths[1:4])
#   maxWidth <- do.call(unit.pmax, cbWidths)
#   maxWidth <- unit.pmax(hmg$widths[1:4], maxWidth)
#   
#   # For visibility, set to the maximum width
#   hmg$widths[1:4] <- as.list(maxWidth)
#   for(i in seq_along(cbgs)) {
#     cbgs[[i]]$widths[1:5] <- as.list(unit.c(maxWidth, hmg$widths[5]+hmg$widths[6]))
#   }
#   heights <- unit.c(unit(rep(1,length(cbgs)), "lines"), unit(1, "null"))
#   rval <- do.call(arrangeGrob, args = c(cbgs, list(hmg), ncol=1, heights=list(heights)))
#   return(rval)
# }
# 
# # Taxa sample image
# top25 <- names(sort(taxa_sums(ps), decreasing=T))[1:25]
# pshm <- prune_taxa(top25,ps)
# taxa.order <- names(sort(taxa_sums(pshm)))
# 
# sample.order <- rownames(sample_data(pshm)[order(get_variable(pshm, "CST"))])
# hm <- plot_heatmap.2(pshm, taxa.label="Species", sample.order=sample.order, taxa.order=taxa.order)
# hm <- hm + theme(axis.title.x = element_text(size=10),
#                  axis.title.y = element_text(size=10),
#                  axis.text.x = element_text(size=7),
#                  axis.text.y = element_text(size=7),
#                  plot.title = element_text(size=8),
#                  legend.text = element_text(size=7),
#                  legend.title = element_text(size=8),
#                  #                legend.margin = unit(c(0.1,0.1,0.1,0.1),"mm"),
#                  #                legend.key.height = unit(1, "in"),
#                  legend.key.width = unit(0.15, "in"),
#                  plot.margin=unit(c(0,0,0,0),"mm"))
# 
# ### CHANGING SPECIES TO TAXA ON YLABEL
# labvec <- as(tax_table(pshm)[, "Species"], "character")
# names(labvec) <- taxa_names(pshm)
# labvec <- labvec[taxa.order]
# labvec[is.na(labvec)] <- ""
# labvec[which(labvec == "Lactobacillus reuteri-vaginalis")] <- "L. reuteri-vaginalis"
# hm <- hm + scale_y_discrete("Taxa", labels = labvec)
# hm <- hm + theme(axis.title = element_text(size=10))
# 
# hcbdf <- data.frame(sample_data(pshm))[sample.order,]
# hcbdf$index <- seq(1,nsamples(pshm))
# hcb <- make_hcb(hcbdf, "CST", name="CST", fillScale = CSTFillScale)
# hcb <- hcb + annotate("text", x=tapply(hcbdf$index, hcbdf[,"CST",drop=T], mean), y=1, label=levels(hcbdf[,"CST",drop=T]), size=2)
# 
# hcbPreterm <- make_hcb(hcbdf, "Outcome", name="Very Pre Term",
#                        fillScale = scale_fill_manual(values=c("Term"="grey60", "Preterm"="maroon", "VeryPreterm"="magenta2", "Marginal"="white")))
# #hcbPreterm <- hcbPreterm + theme(axis.text.y = element_text(size=8, face="bold", color="magenta2"))
# #hcbPreterm <- hcbPreterm + theme(axis.text.y = element_text(size=8, face="bold", color="maroon"))
# hcbPreterm <- hcbPreterm + theme(axis.text.y = element_text(size=8, face="bold", color="grey60"))
# 
# Fig2 <- mush(hm, list(hcbPreterm, hcb))
# print(Fig2)

########################################################

# Aggregate data by participants by mean relative abundance for a given OTU
participant_otu <- tapply(sample_names(bacteria_metadata_df), 
                          sample_data(bacteria_metadata_df)$biome_id, 
                          function(samples) rowMeans(t(otu_table(vaginal_relative_abundances))[, samples, drop = FALSE]))
participant_otu <- do.call(cbind, participant_otu)
rownames(participant_otu) <- taxa_names(vaginal_relative_abundances)

## Alpha Div - Shannon Index
shannon.24 <- diversity((otu_table(bacterial.data)), index="shannon")

# Add participant IDs from sample data | Merge the calculated Shannon diversity values with metadata
bacteria_metadata_df <- as(bacteria_metadata_df, "data.frame")
shannon.df.24 <- data.frame("SampleID"=names(shannon.24), "shannon"=shannon.24)
shannon.qr.merged.24 <- merge(shannon.df.24, bacteria_metadata_df[,1:4], by="SampleID") %>% 
  mutate(biome_id=as.integer(biome_id))
# sample_data_df <- as.data.frame(sample_data(fungal_physeq_filtered))
# alpha_diversity_df$biome_id <- sample_data_df$biome_id

# Define unique user IDs and labels
unique.uid.24 <- sort(unique(shannon.qr.merged.24$biome_id))
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

shannon.diet <- shannon.qr.merged.24 %>% 
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

cor.test(shannon.diet$shannon, shannon.diet$total_cal, method = "spearman")

# group_by biome_id
participant.spearman <- shannon.diet %>% 
  group_by(biome_id) %>% 
  summarise(spearman=cor(shannon, total_cal, method="spearman"))

hist(x=participant.spearman$spearman)

