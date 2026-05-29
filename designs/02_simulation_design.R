library(data.table)
library(metafor)

source("R/simulate_meta.R")
source("R/simulate_from_design.R")
source("R/estimability_of_designs.R")

base_results <- readRDS("designs/design1.rds")
full_plan <- base_results$full_plan

#### Simulate from design and assemble the results
#results_list <- simulate_from_design(plan = full_plan, rows = NULL)
#saveRDS(results_list, "designs/study_results_design1.rds")
results_list <- readRDS("designs/study_results_design1.rds")
results <- rbindlist(results_list, idcol = "file")
# reorder vars
results <- results[, .(
  file,
  id,
  successful_simulation,
  study,
  n_total,
  n_trt,
  n_ctl,
  events_trt,
  events_ctl,
  or,
  or_lb,
  or_ub,
  or_pval,
  yi,
  vi,
  sei,
  theta_true,
  re_mean,
  re_mean_lb,
  re_mean_ub,
  re_mean_pval,
  het_I2,
  het_pval,
  cc_applied
)]


### add the full plan info, keep the successful simulation only
final_results <- results[full_plan, on = "id"][(successful_simulation)]
setkey(final_results, id, study)

### what is the estimability of the final_design?
final_full_plan <- final_results[,
  .(
    tier = tier[1],
    source = source[1],
    n_studies = n_studies[1],
    sample_sizes = sample_sizes[1],
    effect_sizes = effect_sizes[1],
    heterogeneity = heterogeneity[1],
    stat_sig_re = stat_sig_re[1]
  ),
  by = id
]


# Estimability
final_estimability_table <-
  assess_estimability(full_design = final_full_plan, spec = base_results$spec)

# minimum rows to estimate main effects: 40
row_main_effects <- final_estimability_table[
  model == "main_effects" & (estimable),
  first(budget)
]

# minimum rows to estimate 2FI: 140
row_2FI <- final_estimability_table[
  model == "two_fi" & (estimable),
  first(budget)
]

# minimum rows to estimate 3FI: Never, as expected
row_3FI <- final_estimability_table[
  model == "three_fi" & (estimable),
  first(budget)
]

### Save for the generation of abstracts:
saveRDS(
  list(
    feasible_plan = final_full_plan,
    spec = spec,
    data = final_results,
    min_row_main_effects = row_main_effects,
    min_row_2FI = row_2FI,
    min_row_3FI = row_3FI
  ),
  "designs/feasible_design1.rds"
)

fwrite(final_full_plan, "designs/feasible_plan.csv")
fwrite(final_results, "designs/feasible_study_data.csv")
