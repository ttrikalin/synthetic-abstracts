if (!requireNamespace("DoE.base", quietly = TRUE)) {
  install.packages("DoE.base")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
library(data.table)

source("R/base_design.R")
source("R/augment_design.R")
source("R/design_diagnostics.R")

spec <- list(
  nlevels = c(4, 4, 5, 4, 2),
  factor_names = c(
    "n_studies",
    "sample_sizes",
    "effect_sizes",
    "heterogeneity",
    "stat_sig_re"
  ),
  level_labels = list(
    n_studies = c("5", "10", "15", "20"),
    sample_sizes = c("20-50", "51-100", "100-1000", "mix"),
    effect_sizes = c("0.5", "0.67", "1", "1.5", "2"),
    heterogeneity = c("none", "small", "medium", "large"),
    stat_sig_re = c("no", "yes")
  ),
  seed = 12345
)

# Build the base fractional design ------------------------------
base_design <- build_base_design(
  nlevels = spec$nlevels,
  factor_names = spec$factor_names,
  level_labels = spec$level_labels,
  ordered_factors = NULL,
  randomize = FALSE,
  seed = spec$seed
)

# Full tiered plan (all 5 tiers, no budget, deterministic order)
full_plan <- build_augmented_design(
  base_design = base_design,
  n_replicates = 0,
  n_extra = NULL,
  sweep_factor = "effect_sizes",
  include_foldover = TRUE,
  include_filler = TRUE,
  budget = NULL,
  randomize = FALSE,
  seed = NULL
)

# See the tier breakdown and cumulative budget table
summarize_design(full_plan)


# Budget levels to compare
budgets <- ceiling(seq(from = 80, to = prod(spec$nlevels), by = 5))


estimability_table <- do.call(
  rbind,
  lapply(budgets, function(b) {
    plan <- build_augmented_design(
      base_design,
      n_replicates = 0,
      n_extra = NULL,
      sweep_factor = "effect_sizes",
      include_foldover = TRUE,
      include_filler = TRUE,
      budget = b,
      randomize = FALSE,
      seed = NULL
    )

    tiers <- paste(sort(unique(plan$tier)), collapse = ",")

    est <- check_estimability(plan, spec$factor_names)
    est$budget <- b
    est$tiers <- tiers
    est[, c(
      "budget",
      "tiers",
      "model",
      "n_runs",
      "n_params",
      "qr_rank",
      "estimable",
      "residual_df"
    )]
  })
)
estimability_table <- data.table(estimability_table)

# minimum rows to estimate 2FI
estimability_table[model == "two_fi" & (estimable), first(budget)]

# minimum rows to estimate 3FI
estimability_table[model == "three_fi" & (estimable), first(budget)]


# save
full_plan <- as.data.table(full_plan)
full_plan <- full_plan[, id := 1:.N][,
  c("id", "tier", "source", spec$factor_names),
  with = FALSE
]

fwrite(estimability_table, "designs/estimability_design1.csv")
fwrite(as.data.table(full_plan), "designs/design1.csv")
saveRDS(
  list(
    "full_plan" = full_plan,
    "estimability_table" = estimability_table,
    "spec" = spec
  ),
  "designs/design1.rds"
)
