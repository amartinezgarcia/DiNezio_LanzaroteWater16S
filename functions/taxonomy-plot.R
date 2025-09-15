# ==========================
# Libraries
# ==========================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(forcats)
  library(ggplot2); library(scales)
  library(pals)      # install.packages("pals")
  library(grid)      # for unit()
  library(tibble)
})

# ==========================
# Helper 1: tidy + taxonomy join + rel. abundance
# ==========================
prep_comm_long <- function(comm, taxonomy, ranks = c("Family","Genus")) {
  ranks <- base::intersect(ranks, colnames(taxonomy))
  if (!all(c("ID", ranks) %in% colnames(taxonomy))) {
    stop("`taxonomy` must contain at least columns: ID and one of: ", paste(ranks, collapse=", "))
  }

  tax_use <- taxonomy |>
    dplyr::select(ID, dplyr::all_of(ranks)) |>
    distinct()

  # warn if duplicated IDs in taxonomy (would double-count if not unique)
  dup_ids <- tax_use |> count(ID) |> filter(n > 1)
  if (nrow(dup_ids) > 0) {
    warning("Taxonomy has duplicated ASV IDs (n = ", nrow(dup_ids),
            "). Keeping all rows; consider deduplicating to avoid double counting.")
  }

  comm_long <- comm |>
    tibble::rownames_to_column("Sample") |>
    pivot_longer(-Sample, names_to = "ID", values_to = "abund") |>
    mutate(abund = suppressWarnings(as.numeric(abund))) |>
    replace_na(list(abund = 0)) |>
    left_join(tax_use, by = "ID") |>
    group_by(Sample) |>
    mutate(rel_abund = ifelse(sum(abund, na.rm = TRUE) > 0,
                              abund / sum(abund, na.rm = TRUE), 0)) |>
    ungroup()

  comm_long
}

# ==========================
# Helper 2: faceted composition plot (Family/Genus)
# ==========================
plot_comp_faceted <- function(
  dat_long,
  tax_rank      = c("Family","Genus"),
  top_n         = 20,
  sample_order,                 # character vector of sample names (desired order)
  groups,                       # character vector, same length as sample_order
  group_levels  = c("cave","pool","pond","sea","salt","well"),
  use_rel_abund = TRUE          # TRUE = relative (recommended), FALSE = raw counts
) {
  tax_rank <- match.arg(tax_rank)
  if (length(sample_order) != length(groups)) {
    stop("`sample_order` and `groups` must have the same length (1 group per sample).")
  }

  # build mapping table and warn if any plotting samples are unmapped
  map_df <- tibble(Sample = sample_order, Group = groups)

  yvar <- "abund"
  unl  <- if (tax_rank == "Family") "Unclassified_Family" else "Unclassified_Genus"
  other_lab <- if (tax_rank == "Family") "Other families" else "Other genera"

  # ensure label columns exist and are not NA/empty
  dat_rank <- dat_long |>
    mutate(
      Sample = fct_relevel(Sample, sample_order),
      .rank  = !! rlang::sym(tax_rank),
      .rank  = if_else(is.na(.rank) | .rank == "", unl, .rank)
    )

  # warn & label unmapped samples
  unmapped <- setdiff(unique(dat_rank$Sample), map_df$Sample)
  if (length(unmapped) > 0) {
    warning("Some samples have no group mapping: ", paste(unmapped, collapse=", "),
            ". They will be assigned to 'other'.")
  }

  # aggregate; lump to top-N across all samples
  top_ids <- dat_rank |>
    group_by(.rank) |>
    summarise(tot = sum(.data[[yvar]], na.rm = TRUE), .groups = "drop") |>
    slice_max(order_by = tot, n = top_n, with_ties = TRUE) |>
    pull(.rank)

  dat_plot <- dat_rank |>
    mutate(rank_plot = if_else(.rank %in% top_ids, .rank, other_lab)) |>
    group_by(Sample, rank_plot) |>
    summarise(val = sum(.data[[yvar]], na.rm = TRUE), .groups = "drop") |>
    left_join(map_df, by = "Sample") |>
    mutate(
      Group = dplyr::coalesce(Group, "other"),
      Group = factor(Group, levels = c(group_levels, "other")),
      rank_plot = fct_relevel(rank_plot, other_lab, after = Inf)
    )

  # palette for exactly the levels present
  levs <- levels(dat_plot$rank_plot)
  core <- setdiff(levs, other_lab)
  set.seed(1)
  pal  <- setNames(glasbey(length(core)), core)
  pal  <- c(pal, setNames("#BFBFBF", other_lab))  # grey for "Other"

  p <- ggplot(dat_plot, aes(Sample, val, fill = rank_plot)) +
    geom_col(position = if (use_rel_abund) "fill" else "stack", alpha = 0.95) +
    { if (use_rel_abund) scale_y_continuous(labels = percent) else scale_y_continuous() } +
    scale_fill_manual(values = pal, drop = FALSE, name = tax_rank) +
    facet_grid(~ Group, scales = "free_x", space = "free_x") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.spacing.x = unit(8, "mm"),
          strip.background = element_rect(fill = "grey95", color = NA),
          strip.text = element_text(face = "bold")) +
    labs(x = "Sample",
         y = if (use_rel_abund) "Relative abundance" else "Abundance",
         fill = tax_rank,
         title = paste("16S community composition by", tax_rank))

  list(plot = p, data = dat_plot, palette = pal)
}



plot_comp_faceted_compact <- function(
  dat_long,
  tax_rank      = c("Family","Genus"),
  top_n         = 20,
  sample_order,
  groups,
  set_title = NULL,
  group_levels  = c("cave","pool","pond","sea","salt","well"),
  use_rel_abund = TRUE,
  # --- compact knobs ---
  base_size     = 8,
  x_text_size   = 5,
  strip_size    = 8,
  bar_alpha     = 0.95,
  bar_width     = 0.85,
  panel_spacing_mm = 2,
  # --- legend controls ---
  legend_rows   = 4,          # << set number of legend ROWS (e.g., 4)
  legend_cols   = NULL,       #    or set columns instead (leave NULL to use rows)
  legend_key_mm = 3,
  legend_text   = 6,
  legend_title  = 7,
  shorten_labels = FALSE
){
  tax_rank <- match.arg(tax_rank)
  stopifnot(length(sample_order) == length(groups))

  map_df <- tibble::tibble(Sample = sample_order, Group = groups)

  yvar <- "abund"

  unl  <- if (tax_rank == "Family") "Unclassified_Family" else "Unclassified_Genus"
  other_lab <- if (tax_rank == "Family") "Other families" else "Other genera"

  dat_rank <- dat_long |>
    dplyr::mutate(
      Sample = forcats::fct_relevel(Sample, sample_order),
      .rank  = !! rlang::sym(tax_rank),
      .rank  = dplyr::if_else(is.na(.rank) | .rank == "", unl, .rank)
    )

  # top-N by abundance used for plotting
  top_ids <- dat_rank |>
    dplyr::group_by(.rank) |>
    dplyr::summarise(tot = sum(.data[[yvar]], na.rm = TRUE), .groups = "drop") |>
    dplyr::slice_max(order_by = tot, n = top_n, with_ties = TRUE) |>
    dplyr::pull(.rank)

  dat_plot <- dat_rank |>
    dplyr::mutate(rank_plot = dplyr::if_else(.rank %in% top_ids, .rank, other_lab)) |>
    dplyr::group_by(Sample, rank_plot) |>
    dplyr::summarise(val = sum(.data[[yvar]], na.rm = TRUE), .groups = "drop") |>
    dplyr::left_join(map_df, by = "Sample") |>
    dplyr::mutate(
      Group      = dplyr::coalesce(Group, "other"),
      Group      = factor(Group, levels = c(group_levels, "other")),
      Sample_lab = if (shorten_labels) base::abbreviate(Sample, minlength = 3) else Sample
    )

  # ---- enforce order so "Unclassified" and "Other" are the last two ----
  current_levels <- levels(dat_plot$rank_plot)
  if (is.null(current_levels)) current_levels <- sort(unique(dat_plot$rank_plot))
  core_levels <- setdiff(current_levels, c(other_lab, unl))
  target_levels <- c(core_levels, other_lab, unl)

  dat_plot <- dat_plot |>
    dplyr::mutate(rank_plot = forcats::fct_relevel(rank_plot, target_levels))

  # palette: core colours + black for Unclassified + grey for Other
  set.seed(1)
  pal_core <- setNames(pals::glasbey(length(core_levels)), core_levels)
  pal <- c(
    pal_core,
    setNames("#BFBFBF", other_lab),       # Other = grey
    setNames("#000000", unl)  # Unclassified = black
  )

  # build plot
  p <- ggplot2::ggplot(dat_plot, ggplot2::aes(Sample_lab, val, fill = rank_plot)) +
    ggplot2::geom_col(position = if (use_rel_abund) "fill" else "stack",
                      alpha = bar_alpha, width = bar_width) +
    { if (use_rel_abund)
        ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                    breaks = c(0, .25, .5, .75, 1), expand = c(0,0))
      else
        ggplot2::scale_y_continuous(expand = c(0,0))
    } +
    ggplot2::scale_fill_manual(values = pal, drop = FALSE, name = tax_rank) +
    ggplot2::facet_grid(~ Group, scales = "free_x", space = "free_x") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(size = base_size + 1, face = "bold"),
      axis.title.x    = ggplot2::element_text(size = base_size),
      axis.title.y    = ggplot2::element_text(size = base_size),
      axis.text.x     = ggplot2::element_text(size = x_text_size, angle = 90, hjust = 1, vjust = 0.5),
      strip.text      = ggplot2::element_text(size = strip_size, face = "bold"),
      panel.spacing.x = grid::unit(panel_spacing_mm, "mm"),
      legend.position = "bottom",
      legend.text     = ggplot2::element_text(size = legend_text),
      legend.title    = ggplot2::element_text(size = legend_title),
      legend.key.size = grid::unit(legend_key_mm, "mm")
    ) +
    ggplot2::labs(
      x = "Sample",
      y = if (use_rel_abund) "Relative abundance" else "Abundance",
      title = set_title
    )

  # legend rows/cols control
  p <- p + guides(fill = ggplot2::guide_legend(
    nrow = if (!is.null(legend_rows)) legend_rows else NULL,
    ncol = if (!is.null(legend_cols)) legend_cols else NULL
  ))

  list(plot = p, data = dat_plot, palette = pal, levels = target_levels)
}

# saver unchanged
save_compact <- function(p, file, width = 7.5, height = 10, dpi = 300){
  ggplot2::ggsave(file, p, width = width, height = height, dpi = dpi, limitsize = FALSE)
}