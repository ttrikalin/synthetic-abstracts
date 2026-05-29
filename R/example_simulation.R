library(data.table)
library(metafor)

source("R/simulate_meta.R")
source("R/simulate_from_design.R")

full_plan <- readRDS("designs/design1.rds")[[1]]

#results_list <- simulate_from_design(plan = full_plan, rows = 2:3)
results_list <- simulate_from_design(plan = full_plan, rows = NULL)

saveRDS(results_list, "designs/study_results_design1.rds")
