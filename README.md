# Diagnostic Accuracy Analysis for a Cluster-Randomized Trial
R-based analysis of diagnostic accuracy outcomes in a cluster-randomized clinical trial, using GEE models with Kauermann-Carroll small-sample corrections.

## Overview
This repository contains an R code template used to analyze diagnostic accuracy outcomes in a cluster-randomized trial. The analysis evaluates whether an intervention program improves healthcare workers' diagnostic accuracy for specific clinical conditions compared to a standard control group.

## Key Features
- **Cluster-level metrics:** Calculates sensitivity, specificity, PPV, NPV, and accuracy by cluster and trial arm.
- **GEE with small-sample corrections:** Implements Kauermann-Carroll (KC) corrected standard errors for clustered data.
- **Joint hypothesis tests:** Multivariate Wald tests for sensitivity/specificity and PPV/NPV pairs.
- **Link function comparison:** Supports both logit (odds ratios) and identity (risk differences) link functions for overall accuracy.
- **Covariance-based contrasts:** Retains the KC covariance matrix to calculate standard errors and confidence intervals for derived sensitivity, specificity, PPV, and NPV comparisons.
- **Comprehensive output:** Creates model results, diagnostic metric summaries, and cluster-level summary objects corresponding to the report's tables. Plot generation is optional.

## Requirements
- R (version 4.2.1 or higher)

**Required packages:**

```r
install.packages(c("foreign", "tidyverse", "geepack", "geesmv", "broom", "stringr", "ggplot2"))
```

## File Structure
```text
Trial-Diagnostic-Analysis/
├── README.md
└── Primary Analytic Code.Rmd
```

## Code Structure
The script is organized into 12 sections corresponding to the statistical analysis report:

| Section | Description | Outputs |
| :--- | :--- | :--- |
| 1–2 | Data preparation and descriptive statistics | `cond1_data`, `cond2_data`, `table4` |
| 3 | Cluster-level metrics | `cond1_cluster`, `cond2_cluster` |
| 4 | GEE function with KC correction and contrast calculations | `run_gee_kc()`, `calculate_kc_contrast()` |
| 5 | Overall accuracy models (logit) | `all_accuracy`, `table7` |
| 6 | Sensitivity/specificity models | `all_sens_spec` |
| 7 | PPV/NPV models | `all_ppv_npv` |
| 8 | Joint hypothesis tests | `all_joint` |
| 9 | Identity-link models for overall accuracy | `all_identity_acc`, Table 11 |
| 10 | Cluster line plots | Plotting function defined; generation loop commented out |
| 11 | Arm-level summaries | `cond1_summary`, `cond2_summary` |
| 12 | Cluster-level summaries by condition and arm | Appendix Tables 16–19 |

## Key Functions

### `run_gee_kc()`
Core GEE function with Kauermann–Carroll (KC) small-sample correction.

- Fits GEE models using `geepack`.
- Applies KC correction via `GEE.var.kc()` from `geesmv`.
- Returns parameter estimates, KC-corrected standard errors, test statistics, p-values, confidence intervals, and degrees of freedom.
- Retains the KC covariance matrix for calculating standard errors and confidence intervals for derived linear contrasts.

Uses adjusted degrees of freedom:

`df = (Number of Clusters) - (Number of Parameters) - 1`

### `calculate_kc_contrast()`
Calculates estimates, standard errors, confidence intervals, and p-values for linear combinations of model coefficients using the retained KC covariance matrix. This is used for the sensitivity, specificity, PPV, and NPV comparisons.

### `calculate_cluster_metrics()`
Computes diagnostic accuracy metrics for each cluster:

- Sensitivity = TP / (TP + FN)
- Specificity = TN / (TN + FP)
- PPV = TP / (TP + FP)
- NPV = TN / (TN + FN)
- Accuracy = (TP + TN) / N

## Output Objects

| Object | Description | Corresponding Table |
| :--- | :--- | :--- |
| `table4` | Descriptive statistics by condition and trial arm | Table 4 |
| `all_accuracy` | Overall accuracy GEE results | Table 7 |
| `table7` | Formatted overall accuracy model results | Table 7 |
| `all_sens_spec` | Sensitivity/specificity model results | Table 8 |
| `all_ppv_npv` | PPV/NPV model results | Table 9 |
| `all_joint` | Joint hypothesis test results | Table 10 |
| `all_identity_acc` | Identity-link accuracy results | Table 11 |
| `cond1_summary`, `cond2_summary` | Arm-level diagnostic accuracy summaries | Tables 5–6 |
| `cond1_control_summary` | Condition 1 control-arm cluster summary | Appendix Table 16 |
| `cond1_intervention_summary` | Condition 1 intervention-arm cluster summary | Appendix Table 17 |
| `cond2_control_summary` | Condition 2 control-arm cluster summary | Appendix Table 18 |
| `cond2_intervention_summary` | Condition 2 intervention-arm cluster summary | Appendix Table 19 |

## Methodological Notes
- **Clustering:** GEE models account for correlated outcomes within clusters, the randomization units.
- **Small-sample correction:** Kauermann–Carroll (KC) correction is used for individual coefficient tests and derived linear contrasts to address potential downward bias in variance estimates.
- **Joint tests:** Joint hypotheses are evaluated using Wald tests with the robust sandwich covariance matrices from the fitted GEE models. These use a different covariance estimator from the KC-corrected individual and contrast tests.
- **Link functions:**
  - Logit (primary analysis): estimates odds ratios for overall accuracy and for the modeled diagnostic accuracy comparisons.
  - Identity (secondary analysis): estimates risk differences for overall accuracy.

**Reproducibility note:** The analysis requires `analysis_data.dta` in the working directory. Verify that the KC covariance matrix returned by `geesmv` aligns with the fitted model's coefficient order, and validate the identity-link KC correction with the installed package version and analysis data before interpreting results.
