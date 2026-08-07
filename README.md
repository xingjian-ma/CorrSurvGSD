
# CorrSurvGSD

## Overview

- `CorrSurvGSD` provides closed-form calculations for correlated PFS and OS
group sequential designs under the Fleischer model.
- The package computes calendar cutoffs, per-subject moments, joint
correlations, theoretical power, and optional Monte Carlo summaries.
- The statistical definitions and derivations are documented in
`../Blueprint.md` and the package vignette.

## Installation

- From the repository root, load the package during development with
`devtools::load_all("code/CorrSurvGSD")`.
- Install the local package with `devtools::install("code/CorrSurvGSD")`.

## Main interface

- The main public entry point is `CorrSurvGSD::run_pipeline()`.
- It returns a `trial_state` object containing the validated design and
theoretical results.
- Set `simulation = TRUE` to add Monte Carlo results to the returned state.

## Development

- Run tests with `devtools::test("code/CorrSurvGSD")`.
- Run package checks with `devtools::check("code/CorrSurvGSD")`.
- The package is synchronized to its GitHub repository with Git subtree using
the prefix `code/CorrSurvGSD/`.
