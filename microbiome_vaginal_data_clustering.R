# Cite: https://susan.su.domains/papers/Pregnancy/PNAS_Vaginal_Analysis.html

library(phyloseq)
library(tidyverse)
library(vegan)
library(cluster)
library(igraph)
library(markovchain)
library(RColorBrewer)
library(gridExtra)

# Data created in microbiome_vaginal_data_cleaning.R
# bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/vaginal_bacteria_cleanedv3.rds")
# bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/vaginal_cleaned_max_taxa.rds")
# RELABELED DATA
bacterial.data <- readRDS("/Volumes/T7/microbiome_data/sequenced_data/relabeled_data/vaginal_cleaned_max_taxa.rds")

otu.22 <- otu_table(bacterial.data)
tt.22 <- tax_table(bacterial.data)
tt.df <- as.data.frame(tt.22)

# transform to proportions
ps <- transform_sample_counts(bacterial.data, function(OTU) OTU/sum(OTU))
## Cluster into CSTs

# Genetically similar != functionally similar - use a non-phylogenetically aware distance measure
braydist <- distance(ps, method="bray") # BC dissimilarity matrix
ord = ordinate(ps, method = "MDS", distance = braydist) # reduce the dim

# Figure: plot first 12 eigenvalues
plot_scree(ord) + xlim(as.character(seq(1,12))) +
  ggtitle("MDS-bray ordination eigenvalues")
evs <- ord$value$Eigenvalues
print(evs[1:20])

# Denoise distance matrix - get dims that actually mean something
h_sub5 <- hist(evs[6:length(evs)], 100)

# don't really get what this is
plot(h_sub5$mids, h_sub5$count, log="y", type='h', lwd=10, lend=2, xlab ="Midpoint of eigenvalue bins (by 100)",
       ylab ="Count (log scale)")

## Determine number of clusters
NDIM <- 7
x <- ord$vectors[,1:NDIM]  # rows=sample, cols=MDS axes, entries = value
pamPCoA = function(x, k) {
  list(cluster = pam(x[,1:2], k, cluster.only = TRUE))
}
gs = clusGap(x, FUN = pamPCoA, K.max = 12, B = 50)
# use gap statistic to choose the number of clusters

# Cluster: Gap Statistic for Determining Optimal Clusters
plot_clusgap(gs) + 
  scale_x_continuous(breaks=c(seq(0, 12, 2))) +
  theme_minimal() +
  labs(title = "", 
       x = "Number of Clusters (K)", 
       y = "Gap Statistic")

## Cluster into CSTs: PAM-5 clustering
K <- 6# 5
x <- ord$vectors[,1:NDIM]
clust <- as.factor(pam(x, k=K, cluster.only=T))

# SWAPPING THE ASSIGNMENT OF 2 AND 3 TO MATCH RAVEL CST ENUMERATION
# clust[clust==2] <- NA
# clust[clust==4] <- 2
# clust[is.na(clust)] <- 4

# clust[clust==3] <- NA
# clust[clust==4] <- 3
# clust[is.na(clust)] <- 4
sample_data(ps)$CST <- clust
CSTs <- as.character(seq(K))

# Evaluate Clustering
CSTColors <- brewer.pal(6,"Paired")[c(1,3,2,5,4,6)] # Length 6 for consistency with pre-revision CST+ coloration
names(CSTColors) <- CSTs
CSTColorScale <- scale_colour_manual(name = "CST", values = CSTColors[1:5])
CSTFillScale <- scale_fill_manual(name = "CST", values = CSTColors[1:5])

plt1 <- plot_ordination(ps, ord, color="CST")
plt2 <- plot_ordination(ps, ord, axes=c(3,4), color="CST")
grid.arrange(plt1 + CSTColorScale, plt2 + CSTColorScale, top="Ordination by Cluster")

plot_ordination(ps, ordinate(ps, method = "NMDS", distance = braydist), color =
                  "CST") + CSTColorScale + ggtitle("NMDS -- bray -- By Cluster")

# Heatmaps of clustering
taxa.order <- names(sort(taxa_sums(ps)))

i <- 1
plt <- list()
for(CST in CSTs) {
  pshm <- prune_taxa(names(sort(taxa_sums(ps), T))[1:25], ps)
  pshm <- prune_samples(sample_data(pshm)$CST == CST, pshm)
  colnames(tax_table(pshm))[colnames(tax_table(pshm)) == "BLAST_species"] <- "Species (BLAST)"
  
  # print(plot_heatmap(pshm, taxa.label="BLAST_species", taxa.order=taxa.order))
  plt[[i]] <- plot_heatmap(pshm,
               taxa.label = "Species (BLAST)")  +
    theme(axis.text.x = element_blank(),
          # plot.margin = margin(1, 1, 1, 1, "cm"),
          axis.text.y = element_text(size = 6),
          legend.key.size = unit(0.5, "cm")) +
    xlab("Samples")
  i <- i+1
  # print(plt)
}

grid.arrange(grobs = plt, ncol=2)

# grid.arrange(grobs = plt, layout_matrix = rbind(c(1,2), c(3,4), c(5, NA)))


####
# Save R Environment - on not relabeled data (did not resave after relabeling)
# save.image("/Volumes/T7/microbiome_data/R_environments/vaginal_microbiome_CST_clustering.RData")
# load("/Volumes/T7/microbiome_data/R_environments/vaginal_microbiome_CST_clustering.RData")


