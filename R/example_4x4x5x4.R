# example_4x4x5x4.R
# -----------------------------------------------------------------
# Example: tiered augmented fractional factorial for a 4x4x5x4
# experiment, demonstrating budget-based selection and the
# randomize flag.
# -----------------------------------------------------------------

if (!requireNamespace("DoE.base", quietly = TRUE)) {
  install.packages("DoE.base")
}

source("R/base_design.R")
source("R/augment_design.R")
source("R/design_diagnostics.R")

# 1. Inspect candidate orthogonal arrays (informational) -----------
candidates <- list_candidate_oas(nlevels = c(4, 4, 5, 4))

# 2. Build the base fractional design ------------------------------
base_design <- build_base_design(
  nlevels = c(4, 4, 5, 4),
  factor_names = c("A", "B", "C", "D"),
  level_labels = list(A = 1:4, B = 1:4, C = 1:5, D = 1:4),
  randomize = TRUE,
  seed = 42
)

# 3. Full tiered plan (all 5 tiers, no budget, deterministic order)
full_plan <- build_augmented_design(
  base_design = base_design,
  n_replicates = 6,
  n_extra = 12,
  sweep_factor = "C",
  include_foldover = TRUE,
  include_filler = TRUE,
  budget = NULL,
  randomize = FALSE,
  seed = 42
)

# See the tier breakdown and cumulative budget table
summarize_design(full_plan)

# 4. Budget scenarios: pick the first N rows in tier order ----------
#    Because randomize = FALSE, head(full_plan, N) is optimal for
#    any budget N.

# Budget = 80  -> tier 1 only (base OA)
plan_80 <- build_augmented_design(
  base_design,
  sweep_factor = "C",
  budget = 80,
  randomize = FALSE,
  seed = 42
)

# Budget = 150 -> tiers 1 + 2 (base + foldover)
plan_150 <- build_augmented_design(
  base_design,
  sweep_factor = "C",
  budget = 150,
  randomize = FALSE,
  seed = 42
)

# Budget = 170 -> tiers 1 + 2 + 3 + 4 (add replicates + corners)
plan_170 <- build_augmented_design(
  base_design,
  sweep_factor = "C",
  budget = 170,
  randomize = FALSE,
  seed = 42
)

# Compare
budget_comparison <- data.frame(
  budget = c(80, 150, 170, nrow(full_plan)),
  tiers_used = c(
    paste(sort(unique(plan_80$tier)), collapse = ","),
    paste(sort(unique(plan_150$tier)), collapse = ","),
    paste(sort(unique(plan_170$tier)), collapse = ","),
    paste(sort(unique(full_plan$tier)), collapse = ",")
  ),
  n_runs = c(nrow(plan_80), nrow(plan_150), nrow(plan_170), nrow(full_plan))
)
budget_comparison

# 5. Randomized plan for lab execution -----------------------------
#    Same design as full_plan but shuffled with run_order column
randomized_plan <- build_augmented_design(
  base_design = base_design,
  n_replicates = 6,
  n_extra = 12,
  sweep_factor = "C",
  include_foldover = TRUE,
  include_filler = TRUE,
  budget = 170,
  randomize = TRUE,
  seed = 42
)

# tier column preserved; run_order column added
utils::head(randomized_plan, 10)
