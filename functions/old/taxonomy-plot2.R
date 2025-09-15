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
  base_size     = 8,           # overall font sizing
  x_text_size   = 5,           # size of x tick labels
  strip_size    = 8,           # facet strip text
  bar_alpha     = 0.95,
  bar_width     = 0.85,        # bar width (1 = touching)
  panel_spacing_mm = 2,        # gap between facet panels (mm)
  legend_cols   = 8,           # legend columns (increase to compress)
  legend_key_mm = 3,           # key size
  legend_text   = 6,           # legend text size
  legend_title  = 7,           # legend title size
  shorten_labels = FALSE       # TRUE: abbreviate sample names
){
  tax_rank <- match.arg(tax_rank)
  stopifnot(length(sample_order) == length(groups))

  map_df <- tibble::tibble(Sample = sample_order, Group = groups)

  yvar <- "abund"
  unl  <- if (tax_rank == "Family") "Unclassified_Family" else "Unclassified_Genus"
  other_lab <- if (tax_rank == "Family") "Other families" else "Other genera"

  dat_rank <- dat_long |>
    mutate(
      Sample = forcats::fct_relevel(Sample, sample_order),
      .rank  = !! rlang::sym(tax_rank),
      .rank  = if_else(is.na(.rank) | .rank == "", unl, .rank)
    )

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
      Group    = dplyr::coalesce(Group, "other"),
      Group    = factor(Group, levels = c(group_levels, "other")),
      rank_plot = forcats::fct_relevel(rank_plot, other_lab, after = Inf),
      Sample_lab = if (shorten_labels) abbreviate(Sample, minlength = 3) else Sample
    )

  levs <- levels(dat_plot$rank_plot)
  core <- setdiff(levs, other_lab)

  # palette sized exactly to levels present
  set.seed(1)
  pal  <- setNames(pals::glasbey(length(core)), core)
  pal  <- c(pal, setNames("#BFBFBF", other_lab))

  p <- ggplot(dat_plot, aes(Sample_lab, val, fill = rank_plot)) +
    geom_col(position = if (use_rel_abund) "fill" else "stack",
             alpha = bar_alpha, width = bar_width) +
    { if (use_rel_abund) scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                            breaks = c(0, .25, .5, .75, 1), expand = c(0,0))
      else scale_y_continuous(expand = c(0,0)) } +
    scale_fill_manual(values = pal, drop = FALSE, name = tax_rank) +
    facet_grid(~ Group, scales = "free_x", space = "free_x") +
    theme_minimal(base_size = base_size) +
    theme(
      plot.title      = element_text(size = base_size + 1, face = "bold"),
      axis.title.x    = element_text(size = base_size),
      axis.title.y    = element_text(size = base_size),
      axis.text.x     = element_text(size = x_text_size, angle = 90, hjust = 1, vjust = 0.5),
      strip.text      = element_text(size = strip_size, face = "bold"),
      panel.spacing.x = grid::unit(panel_spacing_mm, "mm"),
      legend.position = "bottom",
      legend.text     = element_text(size = legend_text),
      legend.title    = element_text(size = legend_title),
      legend.key.size = grid::unit(legend_key_mm, "mm")
    ) +
    guides(fill = guide_legend(ncol = legend_cols)) +
    labs(x = "Sample",
         y = if (use_rel_abund) "Relative abundance" else "Abundance",
         title = set_title)

  list(plot = p, data = dat_plot, palette = pal)
}

# ---------- convenience saver for tight layouts ----------
save_compact <- function(p, file, width = 7.5, height = 10, dpi = 300){
  ggsave(file, p, width = width, height = height, dpi = dpi, limitsize = FALSE)
}