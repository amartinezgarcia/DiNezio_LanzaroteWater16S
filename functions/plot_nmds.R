plot_nmds <- function(comm_matrix, stations, palette, plot_title,
                      distance = "jaccard",
                      point_size = 1.2,
                      label_size = 2.3,
                      label_vjust = -0.5) {
  # Reorder rows to match station metadata
  comm_matrix <- comm_matrix[match(rownames(stations), rownames(comm_matrix)), , drop = FALSE]
  stopifnot(all(rownames(comm_matrix) == rownames(stations)))

  # Colors by Type
  stations$Type <- factor(stations$Type, levels = names(palette))
  col_vec <- palette[as.character(stations$Type)]

  # Run NMDS
  nmds <- vegan::metaMDS(comm_matrix, distance = distance, trymax = 100, autotransform = FALSE)
  scores <- as.data.frame(vegan::scores(nmds, display = "sites"))
  scores$Sample <- rownames(scores)
  scores$Type <- stations$Type

  # Plot
  ggplot2::ggplot(scores, ggplot2::aes(x = NMDS1, y = NMDS2, color = Type)) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::geom_text(ggplot2::aes(label = Sample), vjust = label_vjust, size = label_size) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::labs(title = plot_title) +
    ggplot2::theme_minimal()
}