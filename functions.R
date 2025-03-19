### Helper Functions

#' map biome health app IDs to study ID
#'
#' @param data input - make health app ID "biome_id" first
#' @param mapping input
#' 
#' @export
study_mapping <- function(data, mapping){
  # Get unique study and biome health pairings
  study_and_u_id <- unique(mapping %>% select(STUDY.ID, Biome.Health.App.ID))
  # Match and join columns
  study_and_u_id <- study_and_u_id %>% 
    rename("study_id" = "STUDY.ID") %>% 
    rename("biome_id" = "Biome.Health.App.ID")
  study_and_u_id$study_id <- as.character(study_and_u_id$study_id)
  # Map ids
  data <- data %>%
    left_join(study_and_u_id, by = "biome_id") %>%
    mutate(biome_id = coalesce(study_id, biome_id)) %>%
    select(-study_id)
  # Check missing ids
  missing_list <- data %>%
    filter(is.na(as.numeric(biome_id)))
  print(unique(missing_list$biome_id))
  
  # convert to numeric
  data$biome_id <- as.numeric(data$biome_id)
  # return data
  data
}

#' map biome health app IDs to study ID
#'
#' @param data input - make health app ID "biome_id" first
#' @param mapping input
#' 
#' @export
filter_id_data <- function(data, id){
  data <- data %>% 
    filter(biome_id == id)
  # return data
  data
}

#' filter to days of study
#'
#' @param data input
#' 
#' @export
filter_days <- function(data){
  
  data <- data %>% 
    filter(logDate < as.Date("2022-12-17") & logDate > as.Date("2022-10-13"))
  
  # return data
  data
}

#' add study days
#'
#' @param data input # needs logDate variable
#' @param mapping input
#' 
#' @export
study_days <- function(data){
  data <- data %>% 
    mutate(Date=as.Date(logDate, format = "%Y-%m-%d"),
           study_day = as.numeric(Date - min(Date, na.rm = TRUE)) + 1)
  # return data
  data
}

#' menstruation heat map plot
#'
#' @param data input
#' 
#' @export
heatmap_plot <- function(data){
  
  all_days <- seq.Date(
    as.Date(min(names(data)[grepl("^2022-", names(data))])), 
    as.Date(max(names(data)[grepl("^2022-", names(data))])), 
    by = "day"
  )
  all_days <- as.character(all_days)
  
  data <- data %>%
    mutate(
      submission_days_count = rowSums(
        across(starts_with("2022-"), ~ !is.na(.), .names = "temp")
        , na.rm = TRUE)
    )
  
  heatmap_data <- data %>%
    select(biome_id, submission_days_count, starts_with("2022-")) %>%
    pivot_longer(
      cols = starts_with("2022-"),
      names_to = "logDate",
      values_to = "value"
    ) %>%
    mutate(
      logDate = as.character(logDate)
    )
  
  # return plot
  ggplot(heatmap_data, aes(x = logDate, y = reorder(factor(biome_id), submission_days_count), fill = factor(value))) +
    geom_tile(color = "black") +
    scale_fill_manual(
      values = c("1" = "red", "0" = "black", "99"="blue", "NA" = "white"),
      na.value = "white",
      name = "Value"
    ) +
    labs(
      x = " ",
      y = " ",
      title = " "
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.grid = element_blank()
    )
  
}

#' menstruation heat map plot
#'
#' @param data input
#' 
#' @export
heatmap_plot1 <- function(data){
  
  all_days <- seq.Date(
    as.Date(min(names(data)[grepl("^2022-", names(data))])), 
    as.Date(max(names(data)[grepl("^2022-", names(data))])), 
    by = "day"
  )
  all_days <- as.character(all_days)
  
  data <- data %>%
    mutate(
      menstruating_days_count = rowSums(
        across(starts_with("2022-"), ~ . == 1, .names = "temp")
        , na.rm = TRUE)
    )

  heatmap_data <- data %>%
    select(biome_id, menstruating_days_count, starts_with("2022-")) %>%
    pivot_longer(
      cols = starts_with("2022-"),
      names_to = "logDate",
      values_to = "value"
    ) %>%
    mutate(
      logDate = as.character(logDate)
    )
  
  # return plot
  ggplot(heatmap_data, aes(x = logDate, y = reorder(factor(biome_id), menstruating_days_count), fill = factor(value))) +
    geom_tile(color = "black") +
    scale_fill_manual(
      values = c("1" = "red", "0" = "black", "99"="blue", "NA" = "white"),
      na.value = "white",
      name = "Value"
    ) +
    labs(
      x = " ",
      y = " ",
      title = " "
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.grid = element_blank()
    )

}

#' menstruation heat map plot - 10 level version
#'
#' @param data input
#' 
#' @export
heatmap_plot2 <- function(data){
  
  all_days <- seq.Date(
    as.Date(min(names(data)[grepl("^2022-", names(data))])), 
    as.Date(max(names(data)[grepl("^2022-", names(data))])), 
    by = "day"
  )
  all_days <- as.character(all_days)
  
  data <- data %>%
    mutate(
      menstruating_days_count = rowSums(
        across(starts_with("2022-"), ~ . == 1, .names = "temp")
        , na.rm = TRUE)
    )
  
  heatmap_data <- data %>%
    select(biome_id, menstruating_days_count, starts_with("2022-")) %>%
    pivot_longer(
      cols = starts_with("2022-"),
      names_to = "logDate",
      values_to = "value"
    ) %>%
    mutate(
      logDate = as.character(logDate)
    )
  
  # return plot
  ggplot(heatmap_data, aes(x = logDate, y = reorder(factor(biome_id), menstruating_days_count), fill = factor(value))) +
    geom_tile(color = "black") +
    scale_fill_manual(
      values = c("1" = "red3", "2" = "pink", "3"="darkred", "4"="darkblue", 
                 "5"="black","6"="blue",
                 "7"="red2", "8" = "white", # 8 is turned into NA
                 "9"="orchid", "10"="slateblue1"),
      na.value = "white",
      name = "Value"
    ) +
    labs(
      x = " ",
      y = " ",
      title = " "
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.grid = element_blank()
    )
  
}

#' menstruation heat map plot - 10 level version & fecal samples
#'
#' @param data input
#' 
#' @export
heatmap_plot_fecal <- function(data){
  
  all_days <- seq.Date(
    as.Date(min(names(data)[grepl("^2022-", names(data))])), 
    as.Date(max(names(data)[grepl("^2022-", names(data))])), 
    by = "day"
  )
  all_days <- as.character(all_days)
  
  data <- data %>%
    mutate(
      menstruating_days_count = rowSums(
        across(starts_with("2022-"), ~ . == 1, .names = "temp")
        , na.rm = TRUE)
    )
  
  heatmap_data <- data %>%
    select(biome_id, menstruating_days_count, starts_with("2022-")) %>%
    pivot_longer(
      cols = starts_with("2022-"),
      names_to = "logDate",
      values_to = "value"
    ) %>%
    mutate(
      logDate = as.character(logDate)
    )
  
  # return plot
  ggplot(heatmap_data, aes(x = logDate, y = reorder(factor(biome_id), menstruating_days_count), fill = factor(value))) +
    geom_tile(color = "black") +
    scale_fill_manual(
      values = c("1" = "red3", "2" = "pink", "3"="darkred", "4"="darkblue", 
                 "5"="black","6"="blue",
                 "7"="red2", "8" = "white", # 8 is turned into NA
                 "9"="orchid", "10"="slateblue1", "50"="yellow"),
      na.value = "white",
      name = "Value"
    ) +
    labs(
      x = " ",
      y = " ",
      title = " "
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.grid = element_blank()
    )
  
}

#' menstruation heat map plot - Imputation
#'
#' @param data input
#' 
#' @export
heatmap_plot_imputation <- function(data){
  
  all_days <- seq.Date(
    as.Date(min(names(data)[grepl("^2022-", names(data))])), 
    as.Date(max(names(data)[grepl("^2022-", names(data))])), 
    by = "day"
  )
  all_days <- as.character(all_days)
  
  data <- data %>%
    mutate(
      menstruating_days_count = rowSums(
        across(starts_with("2022-"), ~ . == 1, .names = "temp")
        , na.rm = TRUE)
    )
  
  heatmap_data <- data %>%
    select(biome_id, menstruating_days_count, starts_with("2022-")) %>%
    pivot_longer(
      cols = starts_with("2022-"),
      names_to = "logDate",
      values_to = "value"
    ) %>%
    mutate(
      logDate = as.character(logDate)
    )
  
  # return plot
  ggplot(heatmap_data, aes(x = logDate, y = reorder(factor(biome_id), menstruating_days_count), fill = factor(value))) +
    geom_tile(color = "black") +
    scale_fill_manual(
      values = c("1" = "red3", "2" = "pink", "3"="darkred", "4"="darkblue", 
                 "5"="black","6"="blue",
                 "7"="red2", "8" = "white", # 8 is turned into NA
                 "9"="orchid", "10"="slateblue1", "50"="yellow",
                 "78"="purple"),
      na.value = "white",
      name = "Value"
    ) +
    labs(
      x = " ",
      y = " ",
      title = " "
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.grid = element_blank()
    )
  
}

#' specific bacterial data filtering
#'
#' @param data, participantID, logDate_list input
#' 
#' @export
#' # vaginal.data, 4, "2022-10-25"
bacterial_plotting <- function(data, participantID, logDate_list){
  
  metadata.22 <- as(sample_data(data), "data.frame")
  
  participant <- metadata.22 %>% 
    filter(biome_id %in% participantID & logDate %in% logDate_list)
  
  participant_sample <- prune_samples(sample_names(data) %in% participant$SampleID, data)
  print(participant$SampleID)
  
  sample_aggregate <- tax_glom(participant_sample, "Genus")
  sample_aggregate_rel <- transform_sample_counts(sample_aggregate, function(x) x / sum(x))
  sample_aggregate_rel_df <- psmelt(sample_aggregate_rel)
  
  top20 <- sample_aggregate_rel_df %>%
    group_by(Genus) %>%
    summarize(TotalAbundance = sum(Abundance)) %>%
    top_n(20, TotalAbundance) %>%
    pull(Genus)
  sample_aggregate_rel_df_top20 <- sample_aggregate_rel_df %>% filter(Genus %in% top20)
  
  # ggplot(sample_aggregate_rel_df_top20, aes(x = Sample, y = Abundance, fill = Genus)) +
  #   geom_bar(stat = "identity", position = position_dodge()) +
  #   theme_minimal() +
  #   labs(title = "Microbial Composition",
  #        x = "Sample ID",
  #        y = "Relative Abundance") +
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # 
  ggplot(sample_aggregate_rel_df_top20, aes(x = Sample, y = Abundance, fill = Genus)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme_minimal() +
    scale_fill_manual(values = c(
      "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
      "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
      "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
      "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3"
    )) +
    labs(title = "Microbial Composition",
         x = "Sample ID",
         y = "Relative Abundance") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

}

#' microbial sample plots
#'
#' @param bacterial.data input
#' 
#' @export
bacterial_plotting2 <- function(bacterial.data){
  
  participant_sample <- prune_samples(sample_names(bacterial.data), bacterial.data)
  
  # aggregate
  sample_aggregate <- tax_glom(participant_sample, "Genus")
  sample_aggregate_rel <- transform_sample_counts(sample_aggregate, function(x) x / sum(x))
  sample_aggregate_rel_df <- psmelt(sample_aggregate_rel)
  
  top20 <- sample_aggregate_rel_df %>%
    group_by(Genus) %>%
    summarize(TotalAbundance = sum(Abundance)) %>%
    top_n(20, TotalAbundance) %>%
    pull(Genus)
  sample_aggregate_rel_df_top20 <- sample_aggregate_rel_df %>% filter(Genus %in% top20)
  
  ggplot(sample_aggregate_rel_df_top20, aes(x = Sample, y = Abundance, fill = Genus)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_minimal() +
    scale_fill_manual(values = c(
      "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
      "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
      "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
      "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3"
    )) +
    labs(title = "Microbial Composition",
         x = "Sample ID",
         y = "Relative Abundance") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
}

#' specific bacterial data filtering
#'
#' @param data, participantID, logDate_list input
#' 
#' @export
#' # vaginal.data, 4, "2022-10-25"
bacterial_plotting3 <- function(data, participantID, logDate_list){
  
  metadata.22 <- as(sample_data(data), "data.frame")
  
  participant <- metadata.22 %>% 
    filter(biome_id %in% participantID & logDate %in% logDate_list)
  
  participant_sample <- prune_samples(sample_names(data) %in% participant$SampleID, data)
  print(participant$SampleID)
  
  sample_aggregate <- tax_glom(participant_sample, "Genus")
  sample_aggregate_rel <- transform_sample_counts(sample_aggregate, function(x) x / sum(x))
  sample_aggregate_rel_df <- psmelt(sample_aggregate_rel)
  
  top20 <- sample_aggregate_rel_df %>%
    group_by(Genus) %>%
    summarize(TotalAbundance = sum(Abundance)) %>%
    top_n(20, TotalAbundance) %>%
    pull(Genus)
  sample_aggregate_rel_df_top20 <- sample_aggregate_rel_df %>% filter(Genus %in% top20)

  ggplot(sample_aggregate_rel_df_top20, aes(x = Sample, y = Abundance, fill = Genus)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme_minimal() +
    scale_fill_manual(values = c(
      "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
      "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
      "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
      "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3",
      "#7B4173", "#D6616B"
    )) +
    labs(title = "Microbial Composition",
         x = "Sample ID",
         y = "Relative Abundance") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

}

#' microbial sample plots
#'
#' @param bacterial.data input
#' 
#' @export
bacterial_plotting4 <- function(bacterial.data){
  
  participant_sample <- prune_samples(sample_names(bacterial.data), bacterial.data)
  
  # aggregate
  sample_aggregate <- tax_glom(participant_sample, "Genus")
  sample_aggregate_rel <- transform_sample_counts(sample_aggregate, function(x) x / sum(x))
  sample_aggregate_rel_df <- psmelt(sample_aggregate_rel)
  
  top20 <- sample_aggregate_rel_df %>%
    group_by(Genus) %>%
    summarize(TotalAbundance = sum(Abundance)) %>%
    top_n(20, TotalAbundance) %>%
    pull(Genus)
  sample_aggregate_rel_df_top20 <- sample_aggregate_rel_df %>% filter(Genus %in% top20)
  
  ggplot(sample_aggregate_rel_df_top20, aes(x = Sample, y = Abundance, fill = Genus)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_minimal() +
    scale_fill_manual(values = c(
      "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
      "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
      "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
      "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3",
      "#7B4173", "#D6616B"
    )) +
    labs(title = "Microbial Composition",
         x = "Sample ID",
         y = "Relative Abundance") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
}


#' mixed effects model, fixed slope
#'
#' @param data input
#' @param thres input
#' @param response input
#' @param fixed input
#' @export lmer_models
mixed_effects_fixed_slope <- function(data, thres, response, fixed_list){
  
  data.filtered <- data %>% 
    group_by(biome_id) %>% 
    filter(n() > thres) %>% 
    ungroup()
  # table(data$biome_id)
  
  lmer_models <- list()
  
  for(fixed in fixed_list) {
    lmer.formula <- as.formula(paste(response, "~", fixed, "+ (1|biome_id)"))
    lmer.obj <- lmer(lmer.formula, data=data.filtered)
    
    lmer_models[[fixed]] <- lmer.obj
  }
  
  return(lmer_models)
  
}

#' mixed effects model, random slope
#'
#' @param data input
#' @param thres input
#' @param response input
#' @param fixed input
#' @export lmer_models
mixed_effects_rnd_slope <- function(data, thres, response, fixed_list){
  
  data.filtered <- data %>% 
    group_by(biome_id) %>% 
    filter(n() > thres) %>% 
    ungroup()
  # table(data$biome_id)
  
  lmer_models <- list()
  
  for(fixed in fixed_list) {
    lmer.formula <- as.formula(paste(response, "~", fixed, "+ (", fixed, "|biome_id)"))
    lmer.obj <- lmer(lmer.formula, data=data.filtered)
    
    lmer_models[[fixed]] <- lmer.obj
  }
  
  return(lmer_models)
  
}
