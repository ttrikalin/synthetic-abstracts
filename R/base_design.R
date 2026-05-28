# base_design.R
# -----------------------------------------------------------------
# Build a base fractional-factorial design from an orthogonal array
# for a user-supplied set of (possibly mixed) factor levels.
# -----------------------------------------------------------------

#' List candidate orthogonal arrays for a mixed-level design
#'
#' @param nlevels Integer vector of the number of levels per factor,
#'   e.g. `c(4, 4, 5, 4)`.
#' @param showmetrics Logical; passed to [DoE.base::show.oas()].
#' @return A data frame of candidate OAs (invisibly also printed).
list_candidate_oas <- function(nlevels, showmetrics = TRUE) {
  stopifnot(is.numeric(nlevels), length(nlevels) >= 1)
  tab <- table(nlevels)
  factors_spec <- list(
    nlevels = as.integer(names(tab)),
    number = as.integer(unname(tab))
  )
  DoE.base::show.oas(factors = factors_spec, showmetrics = showmetrics)
}

#' Build a fractional-factorial base design from an OA
#'
#' @param nlevels Integer vector of levels per factor.
#' @param factor_names Optional character vector of factor names.
#'   Defaults to `"F1", "F2", ...`.
#' @param level_labels Optional named list of level labels, one entry
#'   per factor. Defaults to `1:nlevels[i]`.
#' @param randomize Logical; randomize run order within the OA.
#' @param ordered_factors List of ordered factor names.
#' @param seed Integer seed for reproducibility.
#' @return A `design` object from [DoE.base::oa.design()].
build_base_design <- function(
  nlevels,
  factor_names = NULL,
  level_labels = NULL,
  randomize = FALSE,
  ordered_factors = NULL,
  seed = NULL
) {
  stopifnot(is.numeric(nlevels), length(nlevels) >= 1)
  k <- length(nlevels)

  if (is.null(factor_names)) {
    stopifnot(is.null(ordered_factors))
    factor_names <- paste0("F", seq_len(k))
  }
  stopifnot(length(factor_names) == k)

  if (is.null(level_labels)) {
    level_labels <- lapply(nlevels, seq_len)
  }
  names(level_labels) <- factor_names

  base_design <- DoE.base::oa.design(
    nlevels = nlevels,
    factor.names = level_labels,
    randomize = randomize,
    seed = seed
  )

  if (!is.null(ordered_factors)) {
    stopifnot(all(ordered_factors %in% factor_names))
    for (factor in ordered_factors) {
      base_design[[factor]] <- ordered(
        base_design[[factor()]],
        levels = levels(base_design[[factor]])
      )
    }
  }

  base_design
}
