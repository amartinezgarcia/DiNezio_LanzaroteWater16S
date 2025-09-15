mk_heat <- function(mod, title = "",
                    ord = c("cave","pond","pool","salt","sea","well"),
                    digits = 2, alpha_sig = 1,
                    text_white_threshold = 0.06,
                    label_size = 3,
                    limit = 0.3,
                    tile_ratio = NULL) {   # NEW: NULL = free aspect; number = y/x

  emm <- emmeans::emmeans(mod, ~ Type)
  pr  <- summary(emm, type = "response") |> as.data.frame()
  prob_col <- base::intersect(c("prob","p","response","rate","emmean"), names(pr))
  stopifnot(length(prob_col) > 0)
  probs <- pr |> dplyr::select(Type, prob = dplyr::all_of(prob_col[1]))

  ct <- summary(pairs(emm, adjust = "tukey")) |>
    tidyr::separate(contrast, c("i","j"), sep = " - ") |>
    dplyr::left_join(probs, by = c("i" = "Type")) |> dplyr::rename(p_i = prob) |>
    dplyr::left_join(probs, by = c("j" = "Type")) |> dplyr::rename(p_j = prob) |>
    dplyr::mutate(diff_prob = p_i - p_j, sig = p.value < 0.05)

  ct$i <- factor(ct$i, levels = ord)
  ct$j <- factor(ct$j, levels = ord)
  upper <- dplyr::filter(ct, as.integer(i) < as.integer(j))

  fmt <- paste0("%.", digits, "f")
  upper <- upper |>
    dplyr::mutate(
      label = ifelse(sig, sprintf(fmt, diff_prob), "-"),
      text_col_key = dplyr::case_when(
        !sig ~ "nsg",
        abs(diff_prob) > text_white_threshold ~ "white",
        TRUE ~ "black"
      )
    )

  p <- ggplot2::ggplot(upper, ggplot2::aes(i, j, fill = diff_prob)) +
    ggplot2::geom_tile(color = "grey90", alpha = alpha_sig) +
    ggplot2::geom_text(ggplot2::aes(label = label, color = text_col_key),
                       fontface = "bold", size = label_size) +
    ggplot2::scale_fill_gradient2(
      name   = "\u0394 prob.", limits = c(-limit, limit), midpoint = 0,
      breaks = c(-limit, 0, limit), labels = scales::number_format(accuracy = 0.01),
      oob = scales::squish
    ) +
    ggplot2::scale_color_manual(values = c(white="white", black="black", nsg="grey30"),
                                guide = "none") +
    ggplot2::labs(x = "", y = "", title = title) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "right")

  # Aspect control
  if (!is.null(tile_ratio)) p <- p + ggplot2::coord_fixed(ratio = tile_ratio)
  p
}