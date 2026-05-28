# design_diagnostics.R
# -----------------------------------------------------------------
# Summaries for a tiered augmented fractional-factorial design.
# -----------------------------------------------------------------

#' Summarize an augmented design
#'
#' @param design A data frame returned by [build_augmented_design()].
#' @param factor_names Optional character vector of factor columns to
#'   tabulate; defaults to everything except source, tier, run_order.
#' @return A named list with run counts, tier/source breakdown,
#'   cumulative budget table, and level-frequency tables.
summarize_design <- function(design, factor_names = NULL) {
  meta_cols <- c("source", "tier", "run_order")
  if (is.null(factor_names)) {
    factor_names <- setdiff(names(design), meta_cols)
  }

  # -- tier breakdown -----------------------------------------------
  tier_labels <- c(
    "1" = "base OA (res III)",
    "2" = "foldover (~res IV)",
    "3" = "replicates (pure error)",
    "4" = "corner augment (2FI)",
    "5" = "full-factorial filler"
  )

  tier_tab <- NULL
  budget_tab <- NULL

  if ("tier" %in% names(design)) {
    tier_tab <- table(design$tier)

    # cumulative budget table: for each tier cutoff, how many runs?
    tiers_present <- sort(unique(design$tier))
    cum_n <- vapply(
      tiers_present,
      function(t) {
        sum(design$tier <= t)
      },
      integer(1)
    )
    budget_tab <- data.frame(
      tier = tiers_present,
      description = unname(tier_labels[as.character(tiers_present)]),
      runs_added = as.integer(table(design$tier)[as.character(tiers_present)]),
      cumulative = cum_n,
      stringsAsFactors = FALSE
    )
  }

  list(
    n_total = nrow(design),
    n_by_source = if ("source" %in% names(design)) {
      table(design$source)
    } else {
      NULL
    },
    n_by_tier = tier_tab,
    budget_table = budget_tab,
    level_counts = lapply(design[, factor_names, drop = FALSE], table),
    head_of_plan = utils::head(design, 10)
  )
}
