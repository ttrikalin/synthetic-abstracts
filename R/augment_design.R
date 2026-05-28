# augment_design.R
# -----------------------------------------------------------------
# Tiered augmentation of a base fractional-factorial design.
#
# Tier 1 -- base OA fraction          (resolution III main effects)
# Tier 2 -- foldover of the base      (promotes to ~res IV, de-aliases 2FIs)
# Tier 3 -- replicate runs            (pure error / lack-of-fit d.f.)
# Tier 4 -- corner/edge augmentation  (targeted 2FI estimation)
# Tier 5 -- full-factorial filler      (every remaining unique cell)
#
# With randomize = FALSE, rows are returned in tier order so that
# head(plan, budget) always gives the best design for that budget.
# With randomize = TRUE, the selected rows are shuffled and a
# run_order column records the execution sequence.
# -----------------------------------------------------------------

# -- Tier 2 helpers ------------------------------------------------

#' Compute foldover (level complement) for a single factor column
#' @keywords internal
.complement_level <- function(x) {
  lev <- levels(x)
  k <- length(lev)
  idx <- match(as.character(x), lev)
  factor(lev[k + 1L - idx], levels = lev)
}

#' Build the foldover of a base design
#'
#' Each run is replaced by its level-complement twin: every factor
#' level i is mapped to level (k + 1 - i).
#'
#' @param base_design A design from [build_base_design()].
#' @param factor_names Character vector of factor columns.
#' @return A data frame with tier = 2 and source = "foldover".
make_foldover <- function(base_design, factor_names) {
  fold <- as.data.frame(base_design[, factor_names, drop = FALSE])
  for (nm in factor_names) {
    fold[[nm]] <- .complement_level(fold[[nm]])
  }
  fold$source <- "foldover"
  fold$tier <- 2L
  fold
}

# -- Tier 3 helpers ------------------------------------------------

#' Sample replicate runs from a base design
#'
#' @param base_design A data frame / design from [build_base_design()].
#' @param n_replicates Number of replicate rows to draw.
#' @return A data frame with tier = 3 and source = "replicate".
make_replicates <- function(base_design, n_replicates = 6) {
  stopifnot(n_replicates >= 0)
  if (n_replicates == 0) {
    out <- base_design[0, , drop = FALSE]
  } else {
    replace <- n_replicates > nrow(base_design)
    idx <- sample(seq_len(nrow(base_design)), n_replicates, replace = replace)
    out <- as.data.frame(base_design[idx, , drop = FALSE])
  }
  out$source <- "replicate"
  out$tier <- 3L
  out
}

# -- Tier 4 helpers ------------------------------------------------

#' Build a set of corner-augmenting runs
#'
#' For every factor the extreme levels (first and last) are taken;
#' for the sweep_factor *all* levels are taken.
#'
#' @param nlevels Integer vector of levels per factor.
#' @param factor_names Character vector of factor names.
#' @param level_labels Named list of level labels.
#' @param n_extra Number of augmenting runs to keep.
#' @param sweep_factor Name (or index) of the factor whose levels are
#'   fully swept.  Defaults to the factor with the most levels.
#' @return A data frame with tier = 4 and source = "augment".
make_corner_augment <- function(
  nlevels,
  factor_names,
  level_labels,
  n_extra = 12,
  sweep_factor = NULL
) {
  stopifnot(length(nlevels) == length(factor_names))
  if (n_extra == 0) {
    empty <- as.data.frame(matrix(ncol = length(factor_names), nrow = 0))
    names(empty) <- factor_names
    empty$source <- character(0)
    empty$tier <- integer(0)
    return(empty)
  }

  if (is.null(sweep_factor)) {
    sweep_idx <- which.max(nlevels)
  } else if (is.character(sweep_factor)) {
    sweep_idx <- match(sweep_factor, factor_names)
  } else {
    sweep_idx <- as.integer(sweep_factor)
  }
  stopifnot(!is.na(sweep_idx), sweep_idx >= 1, sweep_idx <= length(nlevels))

  grid_levels <- vector("list", length(nlevels))
  names(grid_levels) <- factor_names
  for (i in seq_along(nlevels)) {
    labs <- level_labels[[i]]
    if (i == sweep_idx) {
      grid_levels[[i]] <- factor(labs, levels = labs)
    } else {
      extremes <- unique(c(labs[1], labs[length(labs)]))
      grid_levels[[i]] <- factor(extremes, levels = labs)
    }
  }

  full <- do.call(expand.grid, c(grid_levels, list(KEEP.OUT.ATTRS = FALSE)))
  n_keep <- min(n_extra, nrow(full))
  out <- full[sample(nrow(full), n_keep), , drop = FALSE]
  out$source <- "augment"
  out$tier <- 4L
  out
}

# -- Tier 5 helpers ------------------------------------------------

#' Build full-factorial filler (remaining unique cells)
#'
#' Generates every cell of the full factorial that is not already
#' present in existing_runs.
#'
#' @param factor_names Character vector of factor names.
#' @param level_labels Named list of level-label vectors.
#' @param existing_runs Data frame containing already-planned runs.
#' @return A data frame with tier = 5 and source = "filler".
make_filler <- function(factor_names, level_labels, existing_runs) {
  grid_levels <- lapply(factor_names, function(nm) {
    factor(level_labels[[nm]], levels = level_labels[[nm]])
  })
  names(grid_levels) <- factor_names
  full <- do.call(expand.grid, c(grid_levels, list(KEEP.OUT.ATTRS = FALSE)))

  # Unique combination key for anti-join
  key_full <- do.call(paste, c(full[factor_names], list(sep = "|")))
  key_exist <- do.call(
    paste,
    c(
      lapply(existing_runs[factor_names], as.character),
      list(sep = "|")
    )
  )
  remaining <- full[!key_full %in% key_exist, , drop = FALSE]

  # Shuffle so that order within tier 5 is arbitrary
  if (nrow(remaining) > 0) {
    remaining <- remaining[sample(nrow(remaining)), , drop = FALSE]
  }
  remaining$source <- "filler"
  remaining$tier <- 5L
  remaining
}

# -- Main entry point ----------------------------------------------

#' Build a tiered augmented design with optional budget and randomization
#'
#' Rows are assembled in tiers of decreasing information-per-run.
#' When randomize = FALSE (the default) rows stay in tier order,
#' so head(plan, budget) always gives the best design for that
#' budget.  When randomize = TRUE, rows are shuffled but the tier
#' column is preserved.
#'
#' @param base_design A design from [build_base_design()].
#' @param n_replicates Number of tier-3 replicate runs (default 6).
#' @param n_extra Number of tier-4 corner runs (default 12).
#' @param sweep_factor Factor whose levels are fully swept in tier 4.
#' @param include_foldover Logical; include the tier-2 foldover block.
#' @param include_filler Logical; include tier-5 full-factorial filler.
#' @param budget Optional integer. If supplied, only the first budget
#'   rows (in tier order) are kept.
#' @param randomize Logical; shuffle the (budget-selected) rows.
#'   Default FALSE so sequential execution follows tier order.
#' @param seed Optional integer seed for reproducibility.
#' @return A data frame with columns for each factor plus source,
#'   tier, and (if randomize = TRUE) run_order.
build_augmented_design <- function(
  base_design,
  n_replicates = 6,
  n_extra = 12,
  sweep_factor = NULL,
  include_foldover = TRUE,
  include_filler = TRUE,
  budget = NULL,
  randomize = FALSE,
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # -- discover factor metadata from the base design ----------------
  factor_names <- names(base_design)[seq_len(ncol(base_design))]
  level_labels <- lapply(base_design[factor_names], function(col) {
    if (is.factor(col)) levels(col) else sort(unique(col))
  })
  names(level_labels) <- factor_names
  nlevels_vec <- vapply(level_labels, length, integer(1))

  # -- Tier 1: base OA ---------------------------------------------
  tier1 <- as.data.frame(base_design[, factor_names, drop = FALSE])
  tier1$source <- "base"
  tier1$tier <- 1L

  # -- Tier 2: foldover --------------------------------------------
  if (include_foldover) {
    tier2 <- make_foldover(base_design, factor_names)
    # Remove foldover rows that duplicate a base row
    key_t1 <- do.call(
      paste,
      c(lapply(tier1[factor_names], as.character), list(sep = "|"))
    )
    key_t2 <- do.call(
      paste,
      c(lapply(tier2[factor_names], as.character), list(sep = "|"))
    )
    tier2 <- tier2[!key_t2 %in% key_t1, , drop = FALSE]
  } else {
    tier2 <- tier1[0, , drop = FALSE]
  }

  # -- Tier 3: replicates ------------------------------------------
  tier3 <- make_replicates(base_design, n_replicates = n_replicates)
  tier3 <- tier3[, names(tier1), drop = FALSE]

  # -- Tier 4: corner augmentation ----------------------------------
  tier4 <- make_corner_augment(
    nlevels = nlevels_vec,
    factor_names = factor_names,
    level_labels = level_labels,
    n_extra = n_extra,
    sweep_factor = sweep_factor
  )

  # -- Tier 5: full-factorial filler --------------------------------
  if (include_filler) {
    existing <- rbind(tier1, tier2, tier4)
    tier5 <- make_filler(factor_names, level_labels, existing)
  } else {
    tier5 <- tier1[0, , drop = FALSE]
  }

  # -- Align factor columns and stack in tier order -----------------
  all_tiers <- list(tier1, tier2, tier3, tier4, tier5)
  for (i in seq_along(all_tiers)) {
    df <- all_tiers[[i]]
    for (nm in factor_names) {
      lev <- as.character(level_labels[[nm]])
      df[[nm]] <- factor(as.character(df[[nm]]), levels = lev)
    }
    all_tiers[[i]] <- df
  }
  plan <- do.call(rbind, all_tiers)
  rownames(plan) <- NULL

  # -- Apply budget -------------------------------------------------
  if (!is.null(budget)) {
    stopifnot(is.numeric(budget), budget >= 1)
    budget <- min(budget, nrow(plan))
    plan <- plan[seq_len(budget), , drop = FALSE]
  }

  # -- Randomize or keep deterministic tier order -------------------
  if (randomize) {
    plan$run_order <- sample(nrow(plan))
    plan <- plan[order(plan$run_order), , drop = FALSE]
    rownames(plan) <- NULL
  }

  plan
}
