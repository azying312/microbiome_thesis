library(tidyverse)

bacterial.data <- readRDS("/Users/alicezhang/Desktop/microbiome_data/sequenced_data/bacteria_cleaned.rds")

# Get vaginal swabs
vaginal_metadata <- bacteria_metadata_df %>% 
  filter(sampleType=="vaginal")

# Prune taxa and OTU table by metadata
bacteria_subset <- prune_samples(samples=vaginal_metadata$SampleID, x=bacterial.data)

#### Taxa table
bacteria_taxa_table <- as.data.frame(tax_table(bacteria_subset))
bacteria_otu_table <- otu_table(bacteria_subset_Species) # OTU row, sample col

bacteria_taxa_table_df <- as.data.frame(bacteria_taxa_table)
# bacteria_taxa_table_df <- bacteria_taxa_table_df %>% 
#   filter(Species=="Lactobacillus amylovorus")
# 
# write.csv(bacteria_taxa_table_df, "/Users/alicezhang/Desktop/lamylovorous.csv", row.names = TRUE)
dim(bacteria_taxa_table_df)

fasta_df <- data.frame(seq_id=rownames(bacteria_taxa_table),
                       seq=bacteria_taxa_table$sequence,
                       stringsAsFactors = FALSE)
write_fasta <- function(df, file_name) {
  con <- file(file_name, "w")
  # Loop through each row and write to the file
  for (i in 1:nrow(df)) {
    # Write the FASTA header
    writeLines(paste(">", df$seq_id[i]), con)
    # Write the sequence
    writeLines(df$seq[i], con)
  }
  close(con)
}

write_fasta(fasta_df, "/Users/alicezhang/Desktop/microbiome_data/sequenced_data/vaginal_bacteria_sequences.fasta")