# Map boolean membership per ASV to a region label like "Environmental&Pathogenic"
make_region_label <- function(df_logical) {
  stopifnot(all(region_names %in% names(df_logical)))
  mat <- as.matrix(df_logical[, region_names])
  apply(mat, 1, function(v) {
    members <- region_names[as.logical(v)]
    paste(members, collapse = if (length(members) > 1) "&" else "")
  })
}

# Build a named counts vector for eulerr from boolean flags (+ optional weights)
counts_for_eulerr <- function(flags_df, weights = NULL) {
  # region label per ASV
  reg <- make_region_label(flags_df)
  # "outside" = no membership in any set (empty string)
  unassigned_n <- if (is.null(weights)) sum(reg == "") else sum(weights[reg == ""], na.rm = TRUE)

  # counts for non-empty regions
  if (is.null(weights)) {
    cnt_df <- as.data.frame(table(reg), stringsAsFactors = FALSE)
    names(cnt_df) <- c("region","n")
  } else {
    cnt_df <- aggregate(weights ~ reg, FUN = sum, na.rm = TRUE)
    names(cnt_df) <- c("region","n")
  }
  cnt_df <- subset(cnt_df, nzchar(region))              # drop empty region (unassigned)
  cnt_vec <- cnt_df$n
  names(cnt_vec) <- cnt_df$region

  list(counts = cnt_vec, unassigned = unassigned_n)
}