#' Plot β-diversity densities within vs across habitat types
#'
#' @param comm   community matrix/data.frame (samples in rows, taxa in cols)
#' @param habitat factor/character vector of habitat labels (length = nrow(comm)).
#'                If named, names must match rownames(comm); otherwise assumed in row order.
#' @param beta_obj list of distance objects from BAT::beta() (Btotal, Brepl, Brich)
#'                 or betapart::beta.pair() (beta.sor, beta.sim, beta.sne).
#' @param show_medians logical; add median vlines per group.
#' @param adjust numeric; density bandwidth adjust passed to geom_density().
#' @return ggplot object (invisible list with data attached as attr "data")
beta_density_plot <- function(comm, habitat, beta_obj,
                              show_medians = TRUE, 
                              adjust = 1.2,
                              fill_values = c("Within habitat" = "#1b9e77",
                                              "Across habitats" = "#d95f02"),
                              line_values = fill_values,
                              alpha = 0.35) {
  stopifnot(nrow(comm) > 1)
  
  # Ensure sample names
  if (is.null(rownames(comm))) {
    rownames(comm) <- sprintf("S%03d", seq_len(nrow(comm)))
  }
  sample_ids <- rownames(comm)
  
  # Normalize habitat vector to named factor aligned to comm
  if (is.null(names(habitat))) {
    if (length(habitat) != nrow(comm)) stop("Un-named 'habitat' must have length nrow(comm).")
    habitat <- factor(habitat)
    names(habitat) <- sample_ids
  } else {
    habitat <- factor(habitat[sample_ids])  # reorder by comm
  }
  
  # Helper to melt a dist/matrix to long upper-triangle
  melt_dist <- function(d, label) {
    # Accept 'dist' or matrix/data.frame
    if (inherits(d, "dist")) {
      # Try to recover labels; fall back to comm rownames
      labs <- attr(d, "Labels")
      if (is.null(labs)) labs <- sample_ids
      m <- as.matrix(d)
      rownames(m) <- colnames(m) <- labs
    } else {
      m <- as.matrix(d)
      if (is.null(rownames(m)) || is.null(colnames(m))) {
        rownames(m) <- colnames(m) <- sample_ids
      }
    }
    
    # Align & subset to the samples in 'comm'
    # Check coverage
    if (!all(sample_ids %in% rownames(m)) || !all(sample_ids %in% colnames(m))) {
      stop("Distance component '", label, "' does not contain all samples from 'comm'.")
    }
    m <- m[sample_ids, sample_ids, drop = FALSE]
    
    idx <- which(upper.tri(m), arr.ind = TRUE)
    tibble::tibble(
      sample1 = rownames(m)[idx[, 1]],
      sample2 = colnames(m)[idx[, 2]],
      value   = m[idx],
      metric  = label
    )
  }
  
  # Detect available components (BAT / betapart)
  have <- names(beta_obj)
  pieces <- list()
  
  if (any(c("Btotal", "Brepl", "Brich") %in% have)) {
    if ("Btotal" %in% have) pieces[["Total β"]] <- beta_obj$Btotal
    if ("Brepl"  %in% have) pieces[["Turnover"]] <- beta_obj$Brepl/beta_obj$Btotal
    if ("Brich"  %in% have) pieces[["Nestedness"]] <- beta_obj$Brich/beta_obj$Btotal
  } else if (any(c("beta.sor", "beta.sim", "beta.sne") %in% have)) {
    # Baselga family (betapart)
    if ("beta.sor" %in% have) pieces[["Total β (Sørensen)"]] <- beta_obj$beta.sor
    if ("beta.sim" %in% have) pieces[["Turnover (Simpson)"]] <- beta_obj$beta.sim/beta_obj$beta.sor
    if ("beta.sne" %in% have) pieces[["Nestedness"]] <- beta_obj$beta.sne/beta_obj$beta.sor
  } else {
    stop("Could not recognize components in 'beta_obj'. Pass output from BAT::beta() or betapart::beta.pair().")
  }
  
  # Build long table
  beta_long <- dplyr::bind_rows(
    lapply(names(pieces), function(lbl) melt_dist(pieces[[lbl]], lbl))
  ) |>
    dplyr::mutate(
      h1 = habitat[sample1],
      h2 = habitat[sample2],
      pair_type = dplyr::if_else(h1 == h2, "Within habitat", "Across habitats"),
      metric = factor(metric, levels = names(pieces))
    )
  
  # Medians for lines
  meds <- beta_long |>
    dplyr::group_by(metric, pair_type) |>
    dplyr::summarise(median = stats::median(value, na.rm = TRUE), .groups = "drop")
  
  # Plot
  p <- ggplot2::ggplot(beta_long, ggplot2::aes(x = value, fill = pair_type)) +
    ggplot2::geom_density(alpha = 0.35, adjust = adjust, na.rm = TRUE) +
    { if (show_medians)
      ggplot2::geom_vline(data = meds,
                          ggplot2::aes(xintercept = median, linetype = pair_type),
                          linewidth = 0.4, show.legend = TRUE)
    } +
    ggplot2::facet_wrap(~ metric, ncol = 1, scales = "free_y") +
    ggplot2::scale_fill_manual(values = fill_values, drop = FALSE) +
    ggplot2::labs(
      x = "Pairwise β-diversity", y = "Density",
      fill = "", linetype = "",
      title = "β-diversity within vs across habitat types"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "top")
  
  attr(p, "data") <- beta_long
  invisible(p)
}