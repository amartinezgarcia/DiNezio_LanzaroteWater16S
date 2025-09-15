
####### Venn diagram part ---------------

asv_flags <- cat_map %>%
  group_by(ID) %>%
  summarise(
    Environmental = any(origin == "Environmental", na.rm = TRUE),
    Anthropogenic = any(origin == "Anthropogenic", na.rm = TRUE),
    Pathogenic    = any(is_patho),
    .groups = "drop"
  )

# helpers to build safe inputs for eulerr + “unassigned” totals
make_region_label <- function(df_logical) {
  rn <- c("Environmental","Anthropogenic","Pathogenic")
  apply(as.matrix(df_logical[, rn]), 1, function(v){
    m <- rn[as.logical(v)]
    paste(m, collapse = if (length(m) > 1) "&" else "")
  })
}

counts_for_eulerr <- function(flags_df, weights = NULL) {
  reg <- make_region_label(flags_df)
  unassigned <- if (is.null(weights)) sum(reg == "") else sum(weights[reg == ""], na.rm = TRUE)
  if (is.null(weights)) {
    cnt_df <- as.data.frame(table(reg), stringsAsFactors = FALSE)
    names(cnt_df) <- c("region","n")
  } else {
    cnt_df <- aggregate(weights ~ reg, FUN = sum, na.rm = TRUE)
    names(cnt_df) <- c("region","n")
  }
  cnt_df <- subset(cnt_df, nzchar(region))
  v <- cnt_df$n; names(v) <- cnt_df$region
  list(counts = v, unassigned = unassigned)
}


####### Histogram part ---------------


# 2a) PREP: collapse to per-sample × combo with a 'val' column
make_combo_table <- function(comm_long, richness = FALSE) {
  if (richness) {
    comm_long %>%
      mutate(present = as.integer(reads > 0)) %>%
      group_by(Sample, combo) %>%
      summarise(val = sum(present), .groups = "drop")
  } else {
    comm_long %>%
      group_by(Sample, combo) %>%
      summarise(val = sum(reads, na.rm = TRUE), .groups = "drop")
  }
}

# 2b) PLOT: compact faceted composition, mirroring plot_comp_faceted_compact()
plot_combo_faceted_compact <- function(
  dat_combo,
  sample_order, groups,
  set_title = NULL,
  group_levels  = c("cave","pool","pond","sea","salt","well"),
  use_rel_abund = TRUE,
  base_size     = 8,  x_text_size = 5,  strip_size = 8,
  bar_alpha     = 0.95, bar_width = 0.85,
  panel_spacing_mm = 2,
  legend_rows   = 4, legend_cols = NULL,
  legend_key_mm = 3, legend_text = 6, legend_title = 7,
  shorten_labels = FALSE
){
  stopifnot(length(sample_order) == length(groups))
  map_df <- tibble::tibble(Sample = sample_order, Group = groups)

  dat_plot <- dat_combo %>%
    mutate(Sample = forcats::fct_relevel(Sample, sample_order),
           Sample_lab = if (shorten_labels) abbreviate(Sample, minlength = 3) else Sample) %>%
    left_join(map_df, by = "Sample") %>%
    mutate(Group = coalesce(Group, "other"),
           Group = factor(Group, levels = c(group_levels, "other")),
           combo = fct_relevel(combo, names(cols_combo)))

  p <- ggplot(dat_plot, aes(Sample_lab, val, fill = combo)) +
    geom_col(position = if (use_rel_abund) "fill" else "stack",
             alpha = bar_alpha, width = bar_width, color = NA) +
    { if (use_rel_abund)
        scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                           breaks = c(0,.25,.5,.75,1), expand = c(0,0))
      else
        scale_y_continuous(expand = c(0,0))
    } +
    scale_fill_manual(values = cols_combo, drop = FALSE, name = "Origin × Pathogenicity") +
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
    labs(x = "Sample",
         y = if (use_rel_abund) "Relative abundance" else "Abundance",
         title = set_title)

  p + guides(fill = guide_legend(
    nrow = if (!is.null(legend_rows)) legend_rows else NULL,
    ncol = if (!is.null(legend_cols)) legend_cols else NULL
  ))
}
