clean_taxa_names <- function(taxa) {
  taxa <- gsub(".", " ", taxa, fixed = TRUE)
  taxa <- gsub("  ", " ", taxa, fixed = TRUE)
  taxa <- gsub("U29 B03", "U29-B03", taxa, fixed = TRUE)
  taxa <- gsub("wb1 A12", "wb1-A12", taxa, fixed = TRUE)
  taxa <- gsub("CSP1 2", "CSP1-2", taxa, fixed = TRUE)
  taxa <- gsub("OM60 NOR5", "OM60(NOR5)", taxa, fixed = TRUE)
  return(taxa)
}

# Wrapper to process and plot top 20 taxa
plot_top_taxa <- function(rf_model, comm_matrix, stations, title, palette) {
  importance_values <- as.data.frame(rf_model$importance$Genus)
  importance_values <- importance_values[order(importance_values$MeanDecreaseAccuracy, decreasing = TRUE), ]

  top20 <- importance_values[1:20, , drop = FALSE]
  rownames(top20) <- clean_taxa_names(rownames(top20))

  # Match cleaned names with community matrix
  matching_taxa <- rownames(top20)[rownames(top20) %in% colnames(comm_matrix)]
  top_comm <- comm_matrix[, matching_taxa, drop = FALSE]
  top_comm$habitat <- as.factor(stations$Type.1)

  # Convert to long format
  long_data <- top_comm %>%
    pivot_longer(cols = -habitat, names_to = "bacterium", values_to = "abundance")

  # Plot
  ggplot(long_data, aes(x = habitat, y = abundance, fill = habitat)) +
    geom_boxplot() +
    facet_wrap(~ bacterium, scales = "free_y") +
    scale_fill_manual(values = palette) +
    theme_minimal() +
    labs(title = title, y = "Abundance", x = "Habitat")
}
