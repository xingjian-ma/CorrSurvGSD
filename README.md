
# CorrSurvGSD

## Overview

- `CorrSurvGSD` evaluates correlated PFS and OS group sequential designs under the Fleischer model.
- The package calculates calendar cutoffs, per-subject moments, the joint correlation matrix, group-sequential boundaries, theoretical power, and optional Monte Carlo summaries.
- The package supports endpoint-specific efficacy alpha spending, optional futility boundaries, and a configurable PFS/OS gatekeeping order.

## Installation

- Install the development version from GitHub with `remotes::install_github("xingjian-ma/CorrSurvGSD")`.
- From a local package checkout, run `devtools::install(".")` in the package root.
- During package development, load the current source with `devtools::load_all(".")` in the package root.

## Main analysis interface

- Run an analysis with `CorrSurvGSD::run_pipeline()`.
- Supply treatment and control sample sizes, planned PFS event counts, accrual rates, control PFS/OS medians, and either treatment medians or endpoint-specific hazard ratios.
- The returned `trial_state` contains validated inputs in `state$design`, numerical settings in `state$options`, and closed-form results in `state$theoretical_results`.
- Set `simulation = TRUE` to append Monte Carlo summaries to `state$empirical_results`.
- Use `integration_seed` for reproducible multivariate-normal integration and `simulation_seed` for reproducible Monte Carlo results.

```r
state <- CorrSurvGSD::run_pipeline(
  n_T = 100,
  n_C = 100,
  d_PFS_vec = 50,
  v_vec = 200 / 12,
  median_PFS_C = log(2) / 0.15,
  median_OS_C = log(2) / 0.05,
  median_PFS_T = log(2) / 0.15,
  median_OS_T = log(2) / 0.05,
  alpha_spending_PFS = "OF",
  alpha_spending_OS = "Pocock",
  efficacy_looks = list(PFS = 1, OS = 1),
  futility_looks = list(PFS = integer(0), OS = integer(0))
)

state$theoretical_results$joint_power
```

## Shiny application

- After installation, launch the bundled application with `CorrSurvGSD::run_app()`.
- The application calls the installed package API and does not require source files from the development checkout.

## Documentation

- The [methodology vignette](vignettes/closed-form-methodology.Rmd) describes the model, closed-form calculations, boundaries, and simulation workflow.
- The methodology vignette is available after installation with `vignette("closed-form-methodology", package = "CorrSurvGSD")`.

## Development and release

- Run tests from the package root with `R -q -e 'testthat::test_local()'`.
- Run a package check from the package root with `R CMD check .`.
- This package is synchronized from the parent workspace with Git subtree using the prefix `code/CorrSurvGSD`.
