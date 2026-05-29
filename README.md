# synthetic-abstracts
Background data for the construction of synthetic abstract -- a design of experiments approach


# Feasible Design Generation

The file `02_simulation_design.R` generates feasible designs for meta-analyses by systematically exploring design options within budget constraints.

## Design Specifications

We generate designs with the following specifications:
- **Design Matrix**: 4×4×5×4x2 factorial design
- **Factors**: Four factors with 4, 4, 5, 4, and 2 levels respectively
- **Budget Constraints**: Various budget levels ranging from 80 to total possible runs (640)
- **Hierarchical Structure**: Tiered design with base orthogonal array (tier 1) and additional augmentations including foldover (tier 2), replicates (tier 3), and corner points (tier 4); then fillers through the full factorial design 

## Design Feasibility

Not all design combinations are feasible. We attempt to generate valid designs by:
1. Running up to 20,000 random tries to generate meta-analyses/trials for each design row; if the RE D-L meta-analysis does not fit the design row description it is rejected. If no example is found the design row is deemed infeasible.  
2. Excluding design rows that are infeasible
3. Assess the estimability of the final feasible plan for budget sizes in batches of 5.  

The resulting feasible design is prioiritized
- **Main Effects**: 40 first rows (best -- all tier one)
- **Two-Factor Interactions**: 140 first rows  
- **Three-Factor Interactions**: Never estimable because the feasible design is not full/orthogonal

## Output Files

### `feasible_plan.csv`
Contains design specifications by run ID
- `id`: Unique identifier for row/run ID
- `tier`: Tier of the design (1-4)
- `source`: Source of the design
- `n_studies`: Number of studies
- `sample_sizes`: Sample size specifications
- `effect_sizes`: Effect size categories
- `heterogeneity`: Heterogeneity levels
- `stat_sig_re`: Statistical significance of random effects

### `feasible_study_data.csv` 
Contains meta-analysis studies and results by design row/run ID:
- `file`: Source file identifier (1 to 1 with `id`)
- `id`: design row/ run ID identifier
- `successful_simulation`: Boolean indicating simulation success (all `TRUE` in this dataset)
- `study`: Study identifier
- `n_total`, `n_trt`, `n_ctl`: Total, treatment, and control group sizes
- `events_trt`, `events_ctl`: Event counts in treatment and control
- `or`, `or_lb`, `or_ub`: Odds ratio with confidence intervals
- `or_pval`: Odds ratio p-value
- `yi`, `vi`, `sei`: Effect size, variance, and standard error
- `theta_true`: True effect size
- `re_mean`, `re_mean_lb`, `re_mean_ub`, `re_mean_pval`: Random effects mean with confidence intervals and p-value
- `het_I2`, `het_pval`: Heterogeneity measures
- `cc_applied`: Whether correction was applied
