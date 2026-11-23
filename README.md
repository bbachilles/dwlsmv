# dwlsmv

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24118857.svg)](https://doi.org/10.5281/zenodo.24118857)

`dwlsmv` is a Stata command that fits confirmatory factor and structural equation models to
ordinal and binary indicators. It estimates by diagonally weighted least squares on thresholds and
polychoric correlations, with robust standard errors and the mean- and variance-adjusted
scaled-and-shifted (T3) test statistic, the approach known as WLSMV. Stata's `sem` and `gsem` do
not offer this estimator.

Current version: **0.7.0**. Requires Stata 16 or later. Tested on Stata 16.1 and 18 on Windows;
other platforms have not been checked.

## Install

```stata
net install dwlsmv, from("https://raw.githubusercontent.com/bbachilles/dwlsmv/v0.7.0/dist/") replace
```

To get the example datasets, change to the folder where you want them and run:

```stata
net get dwlsmv, from("https://raw.githubusercontent.com/bbachilles/dwlsmv/v0.7.0/dist/")
```

Then see `help dwlsmv`, `help difftest` and `help efacompare`.

### If an older version keeps running

`dwlsmv`'s Mata code ships as a compiled library, `ldwlsmv.mlib`. Stata indexes libraries when it
starts, so after installing or upgrading, **restart Stata**. Otherwise run
`mata: mata mlib index`.

If two copies of `ldwlsmv.mlib` are on the ado-path, Stata uses the first one it finds. Stray
copies end up in `PERSONAL` after a manual copy, or in the current working directory if you start
Stata in a folder that contains the package files. The usual sign is that the command errors with
a Mata "function not found" message, or is missing features of the version you installed. To check:

```stata
which dwlsmv
findfile ldwlsmv.mlib, all
```

`which` shows the version line of the `dwlsmv.ado` that Stata runs. `findfile ..., all` lists
every `ldwlsmv.mlib` on the ado-path, in search order (see `adopath`). After `net install`, there
should be exactly one copy, in your `PLUS` folder (`.../plus/l/ldwlsmv.mlib`).

If there is more than one, delete every copy except the one in `PLUS`, and restart Stata. If the
copy in `PLUS` is itself old, run `ado uninstall dwlsmv`, install again with the command above, and
restart Stata.

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

A binary covariate: the factor regressed on an indicator coded 0 and 1, then a direct effect of the
covariate on one item:

```stata
use dwlsmv_example_mg.dta, clear
generate byte grp2 = (grp == 2)
dwlsmv (F -> y1-y5) (F <- grp2)
dwlsmv (F -> y1-y5) (F <- grp2) (y5 <- grp2)
```

On the example data (N = 4,000), the first model gives T3 χ²(9) = 27.238, p = 0.001, and a slope of
0.192 for `F <- grp2`.

Exploratory factor analysis with oblique geomin rotation:

```stata
use dwlsmv_example.dta, clear
dwlsmv (F1 F2 -> y1-y6), efa(F1 F2)
```

## What it does

- **CFA and SEM with latent variables:** factors measured by ordinal or binary indicators, with
  regressions and covariances among the factors. Factors are scaled by a marker indicator or by unit
  variance.
- **Second-order factors:** a factor measured by other factors, scaled by its variance, with a check
  of the first-order factors' disturbance variances. Single-group models only.
- **Fit statistics:** T3 with its degrees of freedom and p-value, plus scaled RMSEA, CFI and TLI.
- **Residual diagnostics:** SRMR in the header, and the correlation residuals and normalized
  residuals stored as matrices, for single-group fits without a covariate. `efacompare` adds an SRMR
  column and lists the largest normalized residuals.
- **Standardized solutions:** delta-method standard errors, and residual variances for the
  indicators.
- **Multiple groups:** configural, thresholds, metric and scalar invariance levels, for models with
  no structural paths between factors.
- **Nested-model tests:** `difftest` runs the scaled-and-shifted difference test for two nested
  models or a chain of them, single-group or multiple-group. It first refits each model and refuses
  the comparison if a refit doesn't reproduce the stored results.
- **Exploratory factor analysis:** oblique geomin rotation with point estimates, fit statistics and
  the eigenvalues for a scree plot. Standard errors and Wald intervals for the rotated loadings and
  factor correlations are shown on request (`se`, `ci`). They are delta-method approximations,
  conditional on the reported rotation, and were tested in a limited set of simulations; `help dwlsmv`
  describes what has and has not been checked.
- **Fit by number of factors:** `efacompare` fits the exploratory model with one factor, two, and so on,
  and shows one row per number of factors: the fit and some over-extraction diagnostics of the
  rotated solution. It names no best number of factors and applies no cutoff.
- **One binary covariate:** a 0/1 observed variable predicting factors and items, with the model
  estimated conditional on it; thresholds are those at covariate = 0. Single-group fits only.
- **Sparse data:** a polychoric pair that can't be estimated stops the fit with a message naming the
  pair. A binary pair with exactly one empty cell is repaired, as lavaan and Mplus do, and listed in
  the output.

## What it does not do

- Continuous or mixed continuous and ordinal indicators.
- Observed predictors other than one binary covariate (a continuous covariate, or several), and a
  covariate together with `group()`. With a covariate: `standardized`, `difftest`, `efa()`,
  `keepbound`, and CFI, TLI and RMSEA.
- Covariances among indicator residuals.
- Second-order factors together with `group()`, `efa()` or a covariate, and `difftest` between a
  second-order model and its correlated-factor parent. Both are refused with a message.
- Sampling weights, clustering, or full-information handling of missing data. Estimation uses
  listwise deletion.
- Partial invariance, or `invariance(scalar)` when any factor has only binary indicators.
- `difftest` or `standardized` on EFA fits. The rotated loadings are already standardized.
- Standard errors of rotated loadings known to be reliable near zero. Those for loadings close to zero can
  undercover, and the package does not flag which; `help dwlsmv` gives the cases tested.

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
Releases are archived on Zenodo. The DOI
[10.5281/zenodo.24118857](https://doi.org/10.5281/zenodo.24118857) always resolves to the latest
version. Each version also has its own DOI, for citing that version exactly. These are listed in
`CITATION.cff` (under `identifiers`) and in the Versions list of the Zenodo record.

The software is licensed under GPL-3.0-or-later (see `LICENSE`). The specification, the prose and
mathematics of the derivation, and the documentation are licensed under CC BY 4.0; the derivation's
code chunks are GPL-3.0-or-later. `CITATION.cff` is CC0 1.0. Source files carry SPDX headers where they can, and `REUSE.toml`
covers every other file. The license texts are in `LICENSES/`.

## Development

The specification (`spec/stata-dwls-wlsmv-spec.md`), the derivations with their executable checks
(`spec/derivations/`), the reference fixtures and their generators (`fixtures/`) and the test suite
(`tests/`) are in this repository. The specification's *Evidence status* paragraph and its
Section 7.6 say what the tests establish and what supporting material is not included.
