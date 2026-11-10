# dwlsmv

`dwlsmv` is a Stata command that fits confirmatory factor and structural equation models to
ordinal and binary indicators. It estimates by diagonally weighted least squares on thresholds and
polychoric correlations, with robust standard errors and the mean- and variance-adjusted
scaled-and-shifted (T3) test statistic, the approach known as WLSMV. Stata's `sem` and `gsem` do
not offer this estimator.

Current version: **0.3.0**. Requires Stata 16 or later. Tested on Stata 16.1 and 18 on Windows;
other platforms have not been checked.

## Install

```stata
net install dwlsmv, from("https://raw.githubusercontent.com/bbachilles/dwlsmv/v0.3.0/dist/") replace
```

To get the example datasets, change to the folder where you want them and run:

```stata
net get dwlsmv, from("https://raw.githubusercontent.com/bbachilles/dwlsmv/v0.3.0/dist/")
```

Then see `help dwlsmv` and `help difftest`.

## Quick start

Models use `sem`'s path syntax.

```stata
use dwlsmv_example.dta, clear
dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)
dwlsmv, standardized
```

On the example data (N = 5,000), the first model gives T3 χ²(8) = 5.154, p = 0.741.

Measurement invariance across groups, testing one level against another:

```stata
use dwlsmv_example_mg.dta, clear
dwlsmv (F -> y1-y5), group(grp) invariance(configural) metric(unitvariance)
estimates store configural
dwlsmv (F -> y1-y5), group(grp) invariance(scalar) metric(unitvariance)
estimates store scalar
difftest configural scalar
```

Exploratory factor analysis with oblique geomin rotation:

```stata
use dwlsmv_example.dta, clear
dwlsmv (F1 F2 -> y1-y6), efa(F1 F2)
```

## What it does

- **CFA and SEM with latent variables:** factors measured by ordinal or binary indicators, with
  regressions and covariances among the factors. Factors are scaled by a marker indicator or by unit
  variance.
- **Fit statistics:** T3 with its degrees of freedom and p-value, plus scaled RMSEA, CFI and TLI.
- **Standardized solutions:** delta-method standard errors, and residual variances for the
  indicators.
- **Multiple groups:** configural, thresholds, metric and scalar invariance levels, for models with
  no structural paths between factors.
- **Nested-model tests:** `difftest` runs the scaled-and-shifted difference test for two nested
  models or a chain of them, single-group or multiple-group. It first refits each model and refuses
  the comparison if a refit doesn't reproduce the stored results.
- **Exploratory factor analysis:** oblique geomin rotation with point estimates and fit statistics.
- **Sparse data:** a polychoric pair that can't be estimated stops the fit with a message naming the
  pair. A binary pair with exactly one empty cell is repaired, as lavaan and Mplus do, and listed in
  the output.

## What it does not do

- Continuous or mixed continuous and ordinal indicators.
- Observed variables as predictors, including observed covariates of a factor.
- Covariances among indicator residuals, and second-order factors.
- Sampling weights, clustering, or full-information handling of missing data. Estimation uses
  listwise deletion.
- Partial invariance, or `invariance(scalar)` when any factor has only binary indicators.
- Standard errors for rotated EFA estimates, or `difftest` on EFA fits.

Unsupported models are refused with a message rather than fitted differently. `help dwlsmv` has the
full list, including `sem` options that are not accepted.

For `difftest`, fit the models and run the test in the same Stata version and edition. For
multiple-group models, keep the number of processors the same too: that's the setting under which
its check of the refit has been validated.

## How it has been checked

Estimates, standard errors and test statistics are compared with lavaan 0.7.2 on committed reference
data, and with Mplus on one reference model. These comparisons cover specific models and data, not
every model, and `help dwlsmv` describes known differences from both programs. `dwlsmv` is pre-1.0
software, so please check results that matter against another program and
[report problems](https://github.com/bbachilles/dwlsmv/issues).

## Citation and license

`CITATION.cff` has the citation details, and GitHub's "Cite this repository" button uses it.

The software is licensed under GPL-3.0-or-later (see `LICENSE`). The specification, the prose and
mathematics of the derivation, and the documentation are licensed under CC BY 4.0; the derivation's
code chunks are GPL-3.0-or-later. `CITATION.cff` is CC0 1.0. Source files carry SPDX headers where they can, and `REUSE.toml`
covers every other file. The license texts are in `LICENSES/`.

## Development

The specification (`spec/stata-dwls-wlsmv-spec.md`), the derivations with their executable checks
(`spec/derivations/`), the reference fixtures and their generators (`fixtures/`) and the test suite
(`tests/`) are in this repository. The specification's *Evidence status* paragraph and its
Section 7.6 say what the tests establish and what supporting material is not included.
