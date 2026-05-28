# simulate_meta.R
# -----------------------------------------------------------------
# Simulate studies under a normal-normal random-effects
# meta-analysis model with the odds ratio (OR) as the effect metric.
#
# Model (log-OR scale):
#   theta_i ~ Normal(mu, tau^2)          # true study log-OR
#   y_i     ~ Normal(theta_i, v_i)       # observed study log-OR
#
# mu  = overall mean log-OR  (log of the pooled OR)
# tau = between-study SD     (heterogeneity, on the log-OR scale)
#
# Two data-generating modes:
#   "binary"  - simulate 2x2 event counts per study from binomial
#               draws given a baseline (control) risk, then compute
#               the empirical log-OR and its variance (with optional
#               continuity correction).  Most faithful to real data.
#   "normal"  - draw y_i directly from Normal(theta_i, v_i) using an
#               analytic large-sample variance.  Fast, no zero cells.
# -----------------------------------------------------------------

#' Simulate a random-effects meta-analysis dataset (OR metric)
#'
#' @param n_studies Integer. Number of studies to generate.
#' @param sample_sizes Per-study total sample size. Either a single
#'   value (recycled to all studies) or a length-`n_studies` vector.
#'   Split evenly between treatment and control arms unless
#'   `allocation` is given.
#' @param mu_or Overall pooled odds ratio (on the natural scale).
#'   Internally converted to `mu = log(mu_or)`. Default 1 (no effect).
#' @param tau Between-study heterogeneity SD on the log-OR scale
#'   (>= 0). `tau = 0` gives a fixed-effect model. Default 0.
#' @param baseline_risk Control-arm event probability used to
#'   generate 2x2 tables (mode = "binary") or to derive the analytic
#'   variance (mode = "normal"). Default 0.2.
#' @param allocation Fraction of each study allocated to the
#'   treatment arm. Default 0.5 (balanced).
#' @param mode Data-generating mechanism: "binary" (default) or
#'   "normal".
#' @param cc Continuity correction added to all cells of a 2x2 table
#'   when any cell is zero (mode = "binary"). Default 0.5.
#' @param seed Optional integer seed for reproducibility.
#' @return A data frame with one row per study containing the true
#'   and observed log-OR, the sampling variance/SE, the OR, and (for
#'   mode = "binary") the underlying 2x2 cell counts.
simulate_meta_or <- function(
  n_studies,
  sample_sizes,
  mu_or = 1,
  tau = 0,
  baseline_risk = 0.2,
  allocation = 0.5,
  mode = c("binary", "normal"),
  cc = 0.5,
  seed = NULL
) {
  mode <- match.arg(mode)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  #browser()

  stopifnot(
    n_studies >= 1,
    tau >= 0,
    baseline_risk > 0,
    baseline_risk < 1,
    allocation > 0,
    allocation < 1
  )

  # -- recycle / validate sample sizes ------------------------------
  if (length(sample_sizes) == 1L) {
    n_total <- rep(sample_sizes, n_studies)
  } else if (length(sample_sizes) == n_studies) {
    n_total <- sample_sizes
  } else {
    stop("`sample_sizes` must have length 1 or `n_studies`.")
  }
  n_trt <- round(n_total * allocation)
  n_ctl <- n_total - n_trt

  # -- true study effects (log-OR scale) ----------------------------
  mu <- log(mu_or)
  theta_i <- rnorm(n_studies, mean = mu, sd = tau) # true log-OR

  # control-arm odds and risk
  p_ctl <- baseline_risk
  odds_c <- p_ctl / (1 - p_ctl)

  if (mode == "binary") {
    # treatment-arm risk implied by this study's true OR
    odds_t <- odds_c * exp(theta_i)
    p_trt <- odds_t / (1 + odds_t)

    # draw event counts
    a <- rbinom(n_studies, size = n_trt, prob = p_trt) # trt events
    c <- rbinom(n_studies, size = n_ctl, prob = p_ctl) # ctl events
    b <- n_trt - a # trt non-events
    d <- n_ctl - c # ctl non-events

    # continuity correction where any cell is zero
    zero_cell <- (a == 0 | b == 0 | c == 0 | d == 0)
    a_cc <- a + ifelse(zero_cell, cc, 0)
    b_cc <- b + ifelse(zero_cell, cc, 0)
    c_cc <- c + ifelse(zero_cell, cc, 0)
    d_cc <- d + ifelse(zero_cell, cc, 0)

    y_i <- log((a_cc * d_cc) / (b_cc * c_cc)) # observed log-OR
    v_i <- 1 / a_cc + 1 / b_cc + 1 / c_cc + 1 / d_cc # sampling variance

    out <- data.frame(
      study = seq_len(n_studies),
      n_total = n_total,
      n_trt = n_trt,
      n_ctl = n_ctl,
      events_trt = a,
      events_ctl = c,
      theta_true = theta_i,
      yi = y_i,
      vi = v_i,
      sei = sqrt(v_i),
      or = exp(y_i),
      cc_applied = zero_cell
    )
  } else {
    # analytic large-sample variance of the log-OR at the implied
    # treatment risk, evaluated at the expected cell counts
    odds_t <- odds_c * exp(theta_i)
    p_trt <- odds_t / (1 + odds_t)

    a_exp <- n_trt * p_trt
    b_exp <- n_trt * (1 - p_trt)
    c_exp <- n_ctl * p_ctl
    d_exp <- n_ctl * (1 - p_ctl)
    v_i <- 1 / a_exp + 1 / b_exp + 1 / c_exp + 1 / d_exp

    y_i <- rnorm(n_studies, mean = theta_i, sd = sqrt(v_i))

    out <- data.frame(
      study = seq_len(n_studies),
      n_total = n_total,
      n_trt = n_trt,
      n_ctl = n_ctl,
      theta_true = theta_i,
      yi = y_i,
      vi = v_i,
      sei = sqrt(v_i),
      or = exp(y_i)
    )
  }

  attr(out, "params") <- list(
    mu_or = mu_or,
    mu = mu,
    tau = tau,
    baseline_risk = baseline_risk,
    allocation = allocation,
    mode = mode
  )
  out
}


get_meta_analysis <- function(ma_data) {
  res <- metafor::rma.uni(
    data = ma_data,
    yi = yi,
    sei = sei,
    method = "REML",
    measure = "GEN"
  )
  meta <- data.table::data.table(ma_data)
  meta[, re_mean := exp(res$b)]
  meta[, re_mean_lb := exp(res$ci.lb)]
  meta[, re_mean_ub := exp(res$ci.ub)]
  meta[, re_mean_pval := res$pval]
  meta[, het_I2 := res$I2]
  meta[, het_pval := res$QEp]

  return(meta)
}


is_acceptable <- function(ma_results, plan) {
  browser()

  if (plan$stat_sig_re == "yes" & ma_results$re_mean_pval[1] >= 0.05) {
    return(FALSE)
  }
  if (plan$stat_sig_re == "no" & ma_results$re_mean_pval[1] < 0.05) {
    return(FALSE)
  }

  if (plan$heterogeneity == "large" & ma_results$het_I2[1] < 2 / 3) {
    return(FALSE)
  }
  if (plan$heterogeneity == "none" & ma_results$het_I2[1] >= 0.1) {
    return(FALSE)
  }
  if (
    plan$heterogeneity == "small" &
      (ma_results$het_I2[1] < 0.1 | ma_results$het_I2[1] >= 1 / 3)
  ) {
    return(FALSE)
  }
  if (
    plan$heterogeneity == "medium" &
      (ma_results$het_I2[1] < 1 / 3 | ma_results$het_I2[1] >= 2 / 3)
  ) {
    return(FALSE)
  }

  return(TRUE)
}
