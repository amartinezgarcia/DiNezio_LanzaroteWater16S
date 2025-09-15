#' Plot β-diversity densities within vs across habitat types
#'
#' @param comm   community matrix/data.frame (samples in rows, taxa in cols)
#' @param habitat factor/character vector of habitat labels (length = nrow(comm)).
#' @param beta_obj list of distance objects from BAT::beta() (Btotal, Brepl, Brich)
#'                 or betapart::beta.pair() (beta.sor, beta.sim, beta.sne).
#' @param show_medians logical; add median vlines per group.
#' @param adjust numeric; density bandwidth adjust passed to geom_density().
#' @param facet_cols integer; number of columns for facet layout (1 => 1×k, k => k×1).
#' @param include_total logical; include Total β panel.
#' @param include_nestedness logical; include Nestedness panel.
#' @param fill_values named colors for "Within habitat" and "Across habitats".
#' @param alpha numeric; fill alpha for densities.
#' @param title character; plot title.
#' @param subtitle character; optional subtitle.
#' @return ggplot object (with raw long data attached as attr "data")
beta_density_plot <- function(
  comm, habitat, beta_obj,
  show_medians = TRUE,
  adjust = 1.2,
  facet_cols = 1,
  include_total = TRUE,
  include_nestedness = TRUE,
  fill_values = c("Within habitat" = "#1b9e77",
                  "Across habitats" = "#d95f02"),
  alpha = 0.35,
  title = "β-diversity within vs across habitat",
  subtitle = NULL
) {
  stopifnot(nrow(comm) > 1)

  # Ensure sample names
  if (is.null(rownames(comm))) rownames(comm) <- sprintf("S%03d", seq_len(nrow(comm)))
  sample_ids <- rownames(comm)

  # Align habitat to comm
  if (is.null(names(habitat))) {
    if (length(habitat) != nrow(comm)) stop("Un-named 'habitat' must have length nrow(comm).")
    habitat <- factor(habitat); names(habitat) <- sample_ids
  } else {
    habitat <- factor(habitat[sample_ids])
  }

  melt_dist <- function(d, label) {
    if (inherits(d, "dist")) {
      labs <- attr(d, "Labels"); if (is.null(labs)) labs <- sample_ids
      m <- as.matrix(d); rownames(m) <- colnames(m) <- labs
    } else {
      m <- as.matrix(d)
      if (is.null(rownames(m)) || is.null(colnames(m))) {
        rownames(m) <- colnames(m) <- sample_ids
      }
    }
    if (!all(sample_ids %in% rownames(m)) || !all(sample_ids %in% colnames(m))) {
      stop("Distance component '", label, "' does not contain all samples from 'comm'.")
    }
    m <- m[sample_ids, sample_ids, drop = FALSE]
    idx <- which(upper.tri(m), arr.ind = TRUE)
    tibble::tibble(
      sample1 = rownames(m)[idx[,1]],
      sample2 = colnames(m)[idx[,2]],
      value   = m[idx],
      metric  = label
    )
  }

  have <- names(beta_obj)
  pieces_all <- list()

  if (any(c("Btotal","Brepl","Brich") %in% have)) {
    if (include_total && "Btotal" %in% have) pieces_all[["Total β"]] <- beta_obj$Btotal
    if ("Brepl" %in% have && "Btotal" %in% have) pieces_all[["Turnover"]] <- beta_obj$Brepl / beta_obj$Btotal
    if (include_nestedness && "Brich" %in% have && "Btotal" %in% have) pieces_all[["Nestedness"]] <- beta_obj$Brich / beta_obj$Btotal
  } else if (any(c("beta.sor","beta.sim","beta.sne") %in% have)) {
    if (include_total && "beta.sor" %in% have) pieces_all[["Total β (Sørensen)"]] <- beta_obj$beta.sor
    if ("beta.sim" %in% have && "beta.sor" %in% have) pieces_all[["Turnover (Simpson)"]] <- beta_obj$beta.sim / beta_obj$beta.sor
    if (include_nestedness && "beta.sne" %in% have && "beta.sor" %in% have) pieces_all[["Nestedness"]] <- beta_obj$beta.sne / beta_obj$beta.sor
  } else {
    stop("Unrecognized 'beta_obj'. Pass output from BAT::beta() or betapart::beta.pair().")
  }

  if (length(pieces_all) == 0) stop("No metrics to plot. Check include_* flags and beta_obj contents.")

  beta_long <- dplyr::bind_rows(lapply(names(pieces_all), function(lbl) melt_dist(pieces_all[[lbl]], lbl))) |>
    dplyr::mutate(
      h1 = habitat[sample1],
      h2 = habitat[sample2],
      pair_type = dplyr::if_else(h1 == h2, "Within habitat", "Across habitats"),
      metric = factor(metric, levels = c(
        grep("^Total", names(pieces_all), value = TRUE),
        grep("Turnover", names(pieces_all), value = TRUE),
        grep("Nestedness", names(pieces_all), value = TRUE)
      ))
    )

  meds <- beta_long |>
    dplyr::group_by(metric, pair_type) |>
    dplyr::summarise(median = stats::median(value, na.rm = TRUE), .groups = "drop")

  p <- ggplot2::ggplot(beta_long, ggplot2::aes(x = value, fill = pair_type)) +
    ggplot2::geom_density(alpha = alpha, adjust = adjust, na.rm = TRUE, color = NA) +
    { if (show_medians)
      ggplot2::geom_vline(data = meds,
                          ggplot2::aes(xintercept = median, linetype = pair_type),
                          linewidth = 0.4, show.legend = TRUE)
    } +
    ggplot2::facet_wrap(~ metric, ncol = facet_cols, scales = "free_y") +
    ggplot2::scale_fill_manual(values = fill_values, drop = FALSE) +
    ggplot2::scale_linetype_manual(values = c("Within habitat" = "dashed",
                                              "Across habitats" = "solid")) +
    ggplot2::labs(
      x = "Pairwise β-diversity", y = "Density",
      fill = "", linetype = "",
      title = title, subtitle = subtitle
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "top")

  attr(p, "data") <- beta_long
  invisible(p)
}