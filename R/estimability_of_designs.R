assess_estimability <- function(full_design, spec, budgets = NULL) {
  #browser()
  if (is.null(budgets)) {
    budgets <- ceiling(seq(from = 10, to = nrow(full_design), by = 5))
  }

  estimability_table <- do.call(
    rbind,
    lapply(budgets, function(b) {
      est <- check_estimability(
        full_design[1:b, , drop = FALSE],
        spec$factor_names
      )
      est$budget <- b
      est[, c(
        "budget",
        "model",
        "n_runs",
        "n_params",
        "qr_rank",
        "estimable",
        "residual_df"
      )]
    })
  )

  return(data.table(estimability_table))
}
