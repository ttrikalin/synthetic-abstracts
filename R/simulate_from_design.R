# simulate_from_design.R
# -----------------------------------------------------------------
# Wrapper that turns one row of an augmented design plan into a
# simulated random-effects meta-analysis dataset (OR metric).
#
# Maps the categorical design factors to the numeric arguments of
# simulate_meta_or():
#
#   n_studies     -> n_studies   (parsed integer)
#   sample_sizes  -> sample_size (per-study draws from a range)
#   effect_sizes  -> mu_or       (the OR value itself)
#   heterogeneity -> tau         (log-OR scale, via lookup)
#   stat_sig_re   -> ignored (it is an OUTCOME, not a generating input)
# -----------------------------------------------------------------

# Requires simulate_meta_or() from R/simulate_meta.R
if (!exists("simulate_meta_or")) {
  source("R/simulate_meta.R")
}

# -- Level-to-parameter lookups ------------------------------------

#' Default mapping from heterogeneity label to tau (log-OR SD)
#' @keywords internal
.default_tau_map <- c(
  none = 0.0,
  small = 0.2,
  medium = 0.4,
  large = 0.7
)

#' Default FIXED sample size (total n per study) for each label.
#'
#' We do NOT guess a size from the range in the label; each category
#' is assigned one fixed total n. "mix" is handled separately in
#' [.parse_sample_sizes()].
#' @keywords internal
.default_ss_map <- c(
  "20-50" = 30L,
  "51-100" = 80L,
  "100-1000" = 750L,
  small = 30L,
  medium = 80L,
  large = 750L
)

#' Build a per-study sample-size vector from a label
#'
#' For a non-"mix" label every study gets the same fixed total n from
#' `ss_map`. For "mix" the composition is: 1 study at the "large"
#' size, and the remaining studies split roughly half "medium" and
#' half "small" (the medium half takes the extra study when the
#' remainder is odd).
#'
#' @param label Sample-size factor level (character).
#' @param n_studies Number of studies to generate sizes for.
#' @param ss_map Named integer vector of fixed sizes per label.
#' @return Integer vector of length n_studies (even totals).
#' @keywords internal
.parse_sample_sizes <- function(label, n_studies, ss_map = .default_ss_map) {
  label <- trimws(as.character(label))

  if (tolower(label) == "mix") {
    sizes <- integer(0)
    if (n_studies >= 1L) {
      sizes <- ss_map[["large"]] # 1 large study
      remaining <- n_studies - 1L
      n_med <- ceiling(remaining / 2) # medium takes the extra
      n_sml <- remaining - n_med
      sizes <- c(
        sizes,
        rep(ss_map[["medium"]], n_med),
        rep(ss_map[["small"]], n_sml)
      )
    }
  } else {
    stopifnot(label %in% names(ss_map))
    sizes <- rep(ss_map[[label]], n_studies)
  }

  sizes <- as.integer(sizes)
  # Force even totals so the two arms split cleanly
  sizes + (sizes %% 2L)
}


# -- Main wrapper --------------------------------------------------

#' Simulate one meta-analysis dataset from a design-plan row
#'
#' @param design_row A one-row data frame / data.table from a design
#'   plan, containing the factor columns listed below.
#' @param tau_map Named numeric vector mapping heterogeneity labels to
#'   tau values on the log-OR scale. Defaults to
#'   c(none = 0, small = 0.1, medium = 0.3, large = 0.6).
#' @param baseline_risk Control-arm event probability passed through.
#' @param mode "binary" (simulate 2x2 counts) or "normal".
#' @param seed Optional integer seed for this dataset.
#' @param ... Additional arguments forwarded to simulate_meta_or().
#' @return A data frame of simulated studies (one row per study) with
#'   a design attribute recording the originating design-row factors.
simulate_from_design_row <- function(
  design_row,
  tau_map = .default_tau_map,
  ss_map = .default_ss_map,
  baseline_risk = 0.2,
  mode = "binary",
  seed = NULL,
  ...
) {
  #browser()
  if (!is.null(seed)) {
    set.seed(seed)
  }
  row <- as.data.frame(design_row, stringsAsFactors = FALSE)
  stopifnot(nrow(row) == 1)

  get_lvl <- function(nm) trimws(as.character(row[[nm]]))

  # n_studies: parse integer
  n_studies <- as.integer(get_lvl("n_studies"))

  # effect_sizes: the OR value itself
  mu_or <- as.numeric(get_lvl("effect_sizes"))

  # heterogeneity: lookup tau
  het_lvl <- tolower(get_lvl("heterogeneity"))
  stopifnot(het_lvl %in% names(tau_map))
  tau <- unname(tau_map[het_lvl])

  # sample_sizes: per-study vector
  ss_lvl <- tolower(get_lvl("sample_sizes"))
  #stopifnot(ss_lvl %in% names(ss_map))
  sample_sizes <- .parse_sample_sizes(
    label = ss_lvl,
    n_studies = n_studies,
    ss_map = .default_ss_map
  )

  out <- simulate_meta_or(
    n_studies = n_studies,
    sample_sizes = sample_sizes,
    mu_or = mu_or,
    tau = tau,
    baseline_risk = baseline_risk,
    mode = mode,
    ...
  )

  # Record the originating design factors for traceability
  attr(out, "design") <- row[,
    intersect(
      c(
        "id",
        "tier",
        "source",
        "n_studies",
        "sample_sizes",
        "effect_sizes",
        "heterogeneity",
        "stat_sig_re"
      ),
      names(row)
    ),
    drop = FALSE
  ]

  out
}

#' Simulate datasets for many design-plan rows
#'
#' Applies [simulate_from_design_row()] to each row of a plan,
#' returning a named list of simulated datasets. Seeds are derived
#' deterministically from base_seed + the row id (or row index) so the
#' whole batch is reproducible.
#'
#' @param plan A design plan (data frame / data.table) with the factor
#'   columns and, optionally, an id column.
#' @param rows Optional integer vector selecting which rows to simulate
#'   (defaults to all rows).
#' @param base_seed Integer base seed; each row gets base_seed + id.
#' @param ... Additional arguments forwarded to
#'   [simulate_from_design_row()].
#' @return A named list of simulated datasets, one per selected row.
simulate_from_design <- function(plan, rows = NULL, base_seed = 1000L, ...) {
  #browser()
  plan <- as.data.frame(plan, stringsAsFactors = FALSE)
  if (is.null(rows)) {
    rows <- seq_len(nrow(plan))
  }

  ids <- if ("id" %in% names(plan)) plan$id[rows] else rows

  datasets <- vector("list", length(rows))
  for (k in seq_along(rows)) {
    cat(paste0(k, "."))
    i <- rows[k]
    tmp <- NULL
    j <- 1
    while (TRUE) {
      tmp <- simulate_from_design_row(
        plan[i, , drop = FALSE],
        seed = j * base_seed + ids[k],
        ...
      )
      tmp <- get_meta_analysis(ma_data = tmp)
      tmp[, id := ids[k]]
      if (is_acceptable(ma_results = tmp, plan = plan[i, , drop = FALSE])) {
        tmp[, successful_simulation := TRUE]
        break
      } else {
        tmp[, successful_simulation := FALSE]
      }
      j <- j + 1
      if (j > 20000) {
        break
      }
    }
    datasets[[k]] <- tmp
  }
  names(datasets) <- paste0("design_", ids)
  datasets
}
