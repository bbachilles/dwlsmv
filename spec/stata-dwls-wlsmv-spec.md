# Stata DWLS/WLSMV Estimator (ordinal/binary v1, with future continuous and ESEM extensions) — Build Spec v9

<!-- SPDX-License-Identifier: CC-BY-4.0 -->
<!-- SPDX-FileCopyrightText: 2026 Erik Ruzek -->

**Purpose of this document:** the formal specification of a Stata command that fits SEM/CFA
models to ordinal and binary data using diagonally weighted least squares (DWLS) with the
Muthén-style mean- and variance-adjusted ("scaled and shifted") WLSMV correction — the
estimator `sem`/`gsem` currently lack. The architecture is deliberately laid out so
continuous indicators (mixed continuous/ordinal models) and exploratory structural equation
modeling (ESEM) can both be added later as extensions of the same backend.

**Version 9.** Gate closures are labelled by the spec version that closed them, so those
labels remain meaningful; the running change history that used to sit here has moved to
the private development record, which is where corrections and superseded decisions belong. This
document states what the estimator must do, not how it came to say so.

**Evidence status.** Development validation, preregistered numerical studies, and internal
decision records supporting this specification are maintained in the private development
repository and are not distributed with this release. The public tests provide executable checks
of the released implementation but do not constitute a reproduction of the complete
development-validation record.

Citations of the form E-001 name one item of that record: a decision, a preregistered study, or
its result. An ID is stable, so a citation can be matched to the item if it is released later;
the item itself is not in this release, and the claim it supports stands on the text here. Words
such as *preregistered*, *chosen after measuring* or *decided (user)* describe how the development
record reached a rule; the rule is the one stated. A passage that calls a script "a scratch run,
discarded" had nothing retained to cite. Section 7.6 says what the tests in this release do
establish.

---

## 1. Motivation / gap being filled

Stata's `sem` and `gsem` support ML and (for `gsem`) generalized-linear-response models for
ordinal indicators (`family(ordinal) link(probit)`), but nothing corresponding to
Mplus's/lavaan's WLSMV: limited-information estimation off polychoric correlations, with a diagonal
weight matrix for point estimation and the full asymptotic covariance matrix used only for robust
SEs and a corrected test statistic. WLSMV is the standard estimator for ordinal CFA/SEM in
psychometrics because it scales to many items/factors where full-information ML with numerical
integration becomes impractical. There is no general-purpose Stata implementation of this currently
on SSC or in official Stata.

## 2. Prior art / references to build against

- Muthén, B. (1984). *A general structural equation model with dichotomous, ordered categorical,
  and continuous latent variable indicators.* Psychometrika, 49(1), 115–132. — the underlying
  latent-response-variable framework, the two-stage limited-information estimator, and the origin
  of the threshold-propagation issue in Section 4.5.
- Muthén, B., & Satorra, A. (1995). *Technical aspects of Muthén's LISCOMP approach to estimation
  of latent variable relations with a comprehensive measurement model.* Psychometrika. — extends
  the asymptotic covariance derivation to the full mixed-type case (relevant once continuous
  indicators are added post-v1).
- Muthén, B., du Toit, S. H. C., & Spisic, D. (1997). *Robust inference using weighted least
  squares and quadratic estimating equations in latent variable modeling with categorical and
  continuous outcomes.* Unpublished manuscript. — the direct source for the WLSMV mean- and
  variance-adjusted (scaled-and-shifted) test statistic and the robust-SE sandwich formula used by
  Mplus. **Note: this paper does NOT contain the scaled-and-shifted statistic.**
  It gives the mean-adjusted $T_1$ (its eq. 12–13), the residual-weight matrix $U$ (eq. 14), and
  a Satterthwaite-style mean- and variance-adjusted $T_2$ that estimates the degrees of freedom
  (eq. 15–16). $U$ and the trace machinery come from here; the shift does not.
- Asparouhov, T., & Muthén, B. (2010). *Simple Second Order Chi-Square Correction.* Mplus
  technical note, May 3 2010
  ([PDF](https://www.statmodel.com/download/WLSMV_new_chi21.pdf)). **The source for the v1
  target statistic.** Defines $T_3 = aT + b$ (its eq. 4) keeping the usual degrees of freedom,
  and its Table 1 confirms the naming: Mplus WLSMV with `Satterthwaite=OFF` — the V6+ default —
  is $T_3$, which is lavaan's `scaled.shifted`.
- Asparouhov, T., & Muthén, B. (2006). *Robust Chi Square Difference Testing with Mean and
  Variance Adjusted Test Statistics.* Mplus Web Notes: No. 10
  ([PDF](http://statmodel.com/download/webnotes/webnote10.pdf)). The source for the DIFFTEST
  construction (Module 12). **Note its eq. (9) as printed does not conform dimensionally** — see
  Module 12.
- Asparouhov, T., & Muthén, B. (2010). *Weighted least squares estimation with missing data.*
  Mplus technical appendix (see [statmodel.com/wlscv.shtml](https://www.statmodel.com/wlscv.shtml)).
- Olsson, U. (1979). *Maximum likelihood estimation of the polychoric correlation coefficient.*
  Psychometrika, 44(4), 443–460. — bivariate-probit ML for a single polychoric correlation; the
  starting point for Module 2.
- Olsson, U. (1979). *On the estimation of polychoric correlations and their asymptotic covariance
  matrix.* Psychometrika. — the closed-form single-pair variance used to validate Module 4's
  diagonal; Module 4 as a whole generalizes this to the joint covariance across pairs, including
  the two-stage threshold-propagation term (Section 4.5).
- Christoffersson, A. (1975). *Factor analysis of dichotomized variables.* Psychometrika. — early
  limited-information approach for the binary case, useful for sanity-checking the tetrachoric
  special case of Modules 2–4.
- lavaan source and documentation — the working reference implementation to numerically match
  against. Specific objects to pull, confirmed directly against lavaan's documentation and source
  (not just inferred from a paraphrase — see Section 4.9's history for why that distinction
  matters here):
  - [`lavOptions`](https://rdrr.io/cran/lavaan/man/lavOptions.html) — valid `estimator`/`test`
    values.
  - [`lavInspect`](https://rdrr.io/cran/lavaan/man/lavInspect.html) — `"gamma"` (documented as
    *N times the asymptotic covariance matrix of the sample statistics*; alias
    `"sampstat.nacov"` — this is Γ under this spec's Section 4.4 convention) vs. `"wls.v"` (the
    WLS *weight vector*, not the covariance matrix — these are different objects and the
    validation code must not conflate them).
  - [`lav_options.R` source](https://github.com/yrosseel/lavaan/blob/master/R/lav_options.R) —
    confirms `wlsm`, `wlsmv`, and `wlsmvs` are all valid estimator names within lavaan's `"DWLS"`
    estimator group (see Section 4.9).
- Kolenikov, S. `polychoric` (Stata Journal / SSC). Existing Stata command for pairwise polychoric
  correlations — reusable or rewritable for Module 2, but does **not** currently expose the joint
  asymptotic covariance matrix needed for Module 4.
- Asparouhov, T., & Muthén, B. (2009). *Exploratory structural equation modeling.* Structural
  Equation Modeling, 16(3), 397–438 ([PDF](https://www.statmodel.com/download/SEM-Asparouhov2009.pdf)).
  — reference for the future ESEM extension (Module 12).

## 3. Scope for v1

**In scope:**
- Single-group CFA and full SEM (measurement model + structural regressions among factors).
- **Indicators: ordinal and binary only** (see the scope decision below).
  Continuous and mixed continuous/ordinal indicators are an explicit, separate v2 extension
  (Section 3's compatibility target below explains why they don't belong in v1).
- DWLS point estimation (diagonal weight matrix) off polychoric correlations and thresholds.
- WLSMV as the default: robust (sandwich) standard errors and the mean- and variance-adjusted
  ("scaled and shifted") test statistic — see Section 4.9 for the precise variant targeted.
- Listwise-deleted sample statistics.

**Explicitly out of scope for v1:**
- Continuous and mixed continuous/ordinal indicators (see the scope decision below).
- Multiple-group models (planned for v2; Module 12).
- Full-information / FIML missing-data handling.
- Pairwise-present sample statistics — deferred to v2, built on the same Γ machinery once
  validated on the listwise case. (Note: lavaan's own documented default for `missing` is
  listwise; pairwise presence is a WLS-specific option, not "the" default of either package —
  don't attribute it more broadly than that when this gets referenced later.)
- Survey design corrections (clustering, weights): planned for v2 (decided by the user 2026-10-05;
  until then `vce()` is refused, D10).
- **Observed exogenous covariates in structural regressions** (e.g. `reg(F2 = F1 x)`) — v1
  supports structural regressions among factors only (Module 11). Estimating a path from an
  observed covariate to a factor requires that covariate's covariance with the factor's ordinal
  indicators, which is a polyserial-correlation problem — i.e. the same mixed-type sample-statistics
  extension already deferred below, not a special case that avoids it.
- ESEM (Module 12, described below for architectural continuity; planned for v3, after v2's
  extensions).
  Exception: single-group EFA with geomin point estimates is pulled into 0.3.0 (user,
  2026-10-26; Section 5, Module 12 release plan).
- Bifactor-ESEM, higher-order factor rotation.

**Scope decision: ordinal/binary indicators only for v1.** A full WLS sample-statistics vector for
genuinely mixed continuous/ordinal data needs to include continuous-variable variances/covariances,
not just correlations — otherwise continuous indicators are implicitly standardized with no
independently estimated residual variance, which is a real but different estimator design choice
that would need its own explicit specification (see the discussion this replaced, preserved for
context: v2 said `s` contains "only correlations and thresholds," which is correct and sufficient
for ordinal/binary data but understates what a *mixed* model's sample-statistics vector needs).
Rather than half-specify the mixed case, v1 restricts to ordinal/binary indicators, where `s` is
unambiguously thresholds + polychoric correlations. Continuous and mixed models become a v2
extension with their own Section-4-style formal specification once v1 is validated — this also
means Modules 2, 4, 5, and 9 don't have to carry unresolved mixed-type complexity into the first
implementation pass.

**Missing-data handling — decided now, not deferred.** Listwise deletion gives a constant sample
size $N$ across every cell of $s$, which makes the joint influence-function construction in
Module 4 straightforward. Pairwise presence gives each pair $(j,k)$ its own effective $N_{jk}$,
which requires sample-size scaling factors inside the sandwich estimator and changes the
asymptotic-distribution derivation. **Decision for v1: listwise deletion.**

**Compatibility target.**

> **Primary compatibility target:** lavaan `estimator = "WLSMV"` (DWLS point estimates, robust
> sandwich SEs, `test = "scaled.shifted"`), single group, listwise-deleted, ordinal/binary
> indicators only, delta parameterization (Section 4.3).
> **Secondary compatibility target:** Mplus WLSMV under equivalent identification and scaling
> conventions.
> **Ultimate reference for correctness:** the Muthén (1984) / Muthén & Satorra (1995) / Muthén,
> du Toit & Spisic (1997) derivations — if Mplus and lavaan ever disagree in a corner case, that's
> the tiebreaker, not "whichever one we happened to validate against first."

### Parameterization compatibility matrix

Pin every choice down in one place, so a future "Stata and lavaan don't match" turns out to be an
actual bug rather than two different parameterizations being compared as if they were the same
estimator:

| Choice | v1 |
|---|---|
| Indicator types | ordinal + binary only |
| Point estimator | DWLS |
| Robust SE | sandwich, $\text{ACOV}(\hat\theta)/(N-1)$ (Section 4.10) |
| Test statistic | scaled-and-shifted (Section 4.9) |
| Ordinal response parameterization | delta (total latent-response variance fixed to 1) |
| Factor metric | chosen per factor: marker-variable or unit-variance (Section 4.8) — record which |
| Missing data | listwise |
| Sample statistics | thresholds + polychoric correlations only (Section 4.2) |
| Statistic ordering | thresholds (by variable, ascending index) then correlations (column-major lower triangle) — matches lavaan, verified (Section 4.2) |
| Γ construction | empirical, $\frac{1}{N}\sum_i \hat\psi_i\hat\psi_i'$ — **required**, not preferred: lavaan does the same, and an analytic plug-in would disagree at $O(N^{-1/2})$ (Section 4.5) |
| Continuous variables | out of scope for v1 |
| Groups | single |
| Clustering / sampling weights | none |

## 4. Formal specification

Everything here should be pinned down *before* substantial Mata code gets written — this section
exists specifically because "stack thresholds and correlations into a vector and weight them" is
directionally correct but not precise enough for two implementers to build the same estimator from
it.

### 4.1 Underlying-variable model

For each ordinal (including binary) indicator $y_j$ with $c_j$ categories, posit a continuous
latent response $y_j^*$ with thresholds $\tau_{j,1} < \dots < \tau_{j,c_j-1}$ such that
$y_j = k \iff \tau_{j,k-1} < y_j^* \le \tau_{j,k}$.

### 4.2 The sample statistics vector $s$ and its index map

For v1 (ordinal/binary only), $s$ is unambiguous: all thresholds, plus all pairwise polychoric
correlations. No polyserial or Pearson entries, no continuous-variable variances — that
simplification is exactly what the Section 3 scope decision buys.

Build, as part of Module 0, a formal **sample-statistic index map** — a table with one row per
element of $s$:

```text
stat_id       (unique index into s)
stat_type     (threshold | polychoric)
var_i, var_j  (which observed variables; var_j blank for a threshold)
threshold_k   (which threshold, for stat_type = threshold)
position      (row/column position in s, Δ, Γ)
```

**Position arithmetic, closed form (Gate 1).** With $C_j = K_j - 1$ thresholds on variable
$j$, $T = \sum_j C_j$, and $P = \binom{p}{2}$, so that $s$ has length $D_s = T + P$:

$$
\text{pos}(\tau_{jc}) = \sum_{m<j} C_m + c,
\qquad
\text{pos}(\rho_{ij}) = T + (i-1)p - \tfrac{i(i-1)}{2} + (j-i) \quad (i<j).
$$

The inverse, needed only for diagnostics, is also closed form: for $r = s - T$ in the correlation
block, $i = \lceil \{(2p-1) - \sqrt{(2p-1)^2 - 8r}\}/2 \rceil$ and
$j = r - \{(i-1)p - i(i-1)/2\} + i$. The ceiling is safe because at every block boundary the
discriminant equals $(2p-2i-1)^2$ — a perfect square, so the root is exact and cannot tip the
ceiling by one. All three are verified in the derivation file over the four fixture designs plus
four shapes the fixtures do not cover, for order, bijectivity and inverse round-trip.

The map is keyed on $C$ alone — not on the data, the estimates, or the model. That is what lets
Module 0 build it before Stage 1 runs. It also means an empty interior category does **not** shrink
the map: $C_j$ as declared by the model still governs it, and the resulting degenerate row of
$\Gamma$ is Section 4.5.6's concern, caught by Module 10's rank check. A map that silently resized
when a category happened to be unobserved would make $\Gamma$ non-conformable with $\Delta$ for
reasons no error message would explain.

Every Mata matrix refers back to this same indexing object. This also directly pays off for two
later extensions: the ESEM extension (Module 12), which needs to relate rotated-solution
quantities back to the same underlying statistics, and the v2 continuous-indicator extension, which
just adds new `stat_type` values and rows to the same table rather than restructuring it.

**Ordering — fixed to match lavaan, and verified against it.** All thresholds first, grouped by
variable and in ascending threshold order within variable, then all pairwise polychoric
correlations in column-major lower-triangle order over the variable list. This ordering is
not arbitrary and must not be changed. On a four-indicator design
(two 4-category, two 3-category items) `lavInspect(fit, "gamma")` is labelled

```text
y1|t1 y1|t2 y1|t3 y2|t1 y2|t2 y2|t3 y3|t1 y3|t2 y4|t1 y4|t2   y1~~y2 y1~~y3 y1~~y4 y2~~y3 y2~~y4 y3~~y4
```

which is exactly the above — $T = 10$ thresholds then $\binom{4}{2} = 6$ correlations. Adopting it
means Module 0's cell-by-cell Γ comparison needs no permutation layer, and a permutation layer is
itself a place to introduce a silent bug. **The harness should nonetheless index by name rather
than by position**, so the comparison survives a future lavaan changing this; `wls.obs`'s names are
the authoritative source (Module 0).

**Multiple groups (0.3.0; specified, not built).** With $G$ groups, Stage 1 runs on each group's
own $(Y_g, K_g)$, so each group has its own index map, keyed on its own $C_g$ exactly as above. The
maps need not agree: an item's category count may differ between groups (what that does to a
threshold tie is Section 4.8's). $s$ is the concatenation $(s_1', \dots, s_G')'$, each block in the
single-group order above. With $D_g = T_g + P$ the length of block $g$, an element of group $g$
sits at its position in $g$'s own map plus $\sum_{h<g} D_h$; the closed forms above apply within
each block unchanged. No field is needed for the offset: it is derived from the maps. $\hat\Gamma$
is block diagonal across independent groups. How the blocks are weighted in the pooled objective,
ACOV and $T_3$ is the pooled multiple-group scaling, derived in `@sec-multigroup-scaling`.
**Decided 2026-10-09 (user), after the check 10 review:** $T = (N-G)F$ (candidate B) and
`e(V)` $= \mathrm{ACOV}/(N-G)$, one constant shared by both, with check 10's
$\varepsilon = (N-G)\,\mathrm{nrtol}$ (Module 12). lavaan's unscaled $\chi^2$ at $G > 1$ is
$\sum_g(N_g-1)F_g$, so dwlsmv does not reproduce it exactly there; it does reproduce lavaan's
`vcov()`. A single-group fit is
$G = 1$ with offset 0, and nothing in this section changes for it.

**Conditional-covariate mode (v0.4.0; specified 2026-10-25, not built).** With $q$ observed
exogenous covariates entering conditionally (`@sec-condx`), `stat_type` gains `slope` (`var_i` an
indicator, `var_j` the covariate). $s$ is thresholds, then slopes, then the **residual** correlations
given $x$, so $D_s = T + pq + P$, with
$\text{pos}(b_{jk}) = T + (k-1)p + j$ (covariate-major, items within each covariate), and every
$\text{pos}(\rho_{ij})$ above shifted by $pq$. The thresholds are at $x = 0$. This order was proved
from lavaan's own labels, not assumed (`@sec-condx-norm`; `fixtures/generate/99h_lavaan_condx.R`).
The map is still keyed on $C$, and now on $q$, never on the data.

### 4.3 Scale convention (simplified by the v1 scope decision)

With continuous indicators out of scope, this is no longer the identification minefield it would
be for a mixed model — there's no continuous-variable variance question to resolve. State it
anyway, since it's still a real choice:

**v1 convention:** all indicators are ordinal/binary; latent response variables $y_j^*$ use the
**delta parameterization** (total variance fixed to 1) — this is the parameterization selected for
the primary compatibility target (lavaan's default for WLSMV on categorical data), not an intrinsic
property of "WLSMV" in general (lavaan also supports a `theta` parameterization, where the
*residual* variance is fixed to 1 instead — a legitimate alternative, just not the v1 target).
Thresholds are parameterized relative to this fixed-scale latent variable. The model-implied moment
structure (Section 4.6) is computed on this same correlation metric.

### 4.4 Asymptotic normalization: what exactly is Γ?

**Convention used throughout this document:**

$$\sqrt{N}(\hat s - s) \xrightarrow{d} N(0, \Gamma)$$

so that $\text{Var}(\hat s) \approx \Gamma / N$. This matches lavaan's own definition of
`lavInspect(fit, "gamma")` as *N times the asymptotic covariance matrix of the sample statistics* —
so that call is the direct, confirmed cross-language validation target for Γ under this convention
(Section 5, Module 0). Every downstream formula ($W_D$, $\text{ACOV}(\hat\theta)$, `e(V)`) is
written in terms of this Γ; the module that assembles the final `e(V)` matrix must apply the
sample-size divisor explicitly — $1/(N-1)$, lavaan's convention, see Section 4.10 — and that step must be unit-tested in isolation — it's exactly the kind of thing that's
easy to get right in the math and wrong in the code.

### 4.5 Γ via observation-level influence functions — including the two-stage correction

Define Γ through the observation-level influence function of every element of $s$:

$$\Gamma = E[\psi_i \psi_i']$$

built, for every element of $s$, as

```text
observation i  →  influence vector ψ_i (one entry per element of s, each correlation's entry
                   including its threshold-propagation terms)  →  Γ = crossprod(ψ) / N
```

Getting $\psi_i$ right is the hard problem — the crossprod step that turns it into Γ is
mechanical. Don't read the crossprod step as "cross-pair covariances fall out automatically
regardless of how $\psi_i$ was built" — they fall out correctly only once $\psi_i$ correctly
includes the two-stage propagation term of Section 4.5.2.

This is the project's highest-risk derivation, so the section is long. Its parts, and what each is
for:

| Part | Covers | Status |
|---|---|---|
| **4.5.1** | The threshold entries — Stage 1, closed form, plus an exact identity worth exploiting | Closed |
| **4.5.2** | The correlation entries — Stage 2 and its threshold-propagation term | Closed |
| **4.5.3** | Sub-phase 4a's acceptance test: the binary case in closed form | Closed, binary only |
| **4.5.4** | The construction rule: empirical, never plug-in | Decided |
| **4.5.5** | Centering, as a diagnostic | Decided |
| **4.5.6** | Numerical hazards: sparse cells have two distinct failure modes | Decided |

**Read 4.5.4 before writing any code**, even while 4.5.2 is still in progress — it constrains how Γ
is assembled regardless of which closed forms are available, and retrofitting it later means
rewriting Module 4's core.

#### 4.5.1 The threshold entries (Stage 1)

Roughly a third of Γ is the threshold-by-threshold block, and unlike the correlation entries of
4.5.2 it is entirely free. Let $q_{j,c} = \Phi(\tau_{j,c})$ be the cumulative probability below
variable $j$'s $c$-th threshold, and $Z_{j,c,i}$ the indicator that observation $i$ falls below it.
The first-stage estimator is $\hat\tau_{j,c} = \Phi^{-1}(\hat q_{j,c})$ with $\hat q_{j,c}$ the
cumulative sample proportion, so by the delta method

$$\psi_{\tau_{j,c},i} = \frac{Z_{j,c,i} - q_{j,c}}{\phi(\tau_{j,c})},
\qquad
\Gamma_{(j,c),(k,d)} = \frac{P(Z_{j,c}=1,\, Z_{k,d}=1) - q_{j,c}q_{k,d}}
                            {\phi(\tau_{j,c})\,\phi(\tau_{k,d})}$$

with the joint probability being $\min(q_{j,c}, q_{j,d})$ when $j = k$ and
$\Phi_2(\tau_{j,c}, \tau_{k,d}; \rho_{jk})$ when $j \neq k$.

Two consequences worth knowing before Module 4 starts. **No entry of the threshold block ever needs
more than a bivariate normal CDF**, at any number of categories or variables — the
higher-dimensional integrals people expect are confined to the correlation entries. And **within a
variable the thresholds are always positively correlated**, whatever $\rho$ is, since the
within-variable entry reduces to $q_{j,c}(1-q_{j,d})/[\phi\phi]$ for $c < d$.

Two implementation notes fall out of the same algebra. $\psi_{\tau}$ depends on the observation
only through which category it landed in, so a $c_j$-category variable has exactly $c_j$ distinct
threshold influence vectors: build a $c_j \times (c_j-1)$ lookup table once and index it by the
observed category rather than computing per-observation. And the first-stage estimator satisfies
the ordering constraint $\tau_{j,1} < \dots < \tau_{j,c_j-1}$ automatically, with no constrained
optimization — provided no interior category is empty, which is exactly the degeneracy flagged in
4.5.6.

**An exact identity: the within-variable block needs no estimation.** Within a single variable, the
empirical and closed-form constructions are not merely close — they are the same arithmetic. The
product of two cumulative indicators on one variable collapses, $Z_{j,c}Z_{j,d} =
Z_{j,\min(c,d)}$, so the sample mean of that product is $\hat q_{j,\min(c,d)}$ *exactly*, and

$$\frac1N\sum_i \hat\psi_{\tau_{j,c},i}\,\hat\psi_{\tau_{j,d},i}
= \frac{\min(\hat q_{j,c}, \hat q_{j,d}) - \hat q_{j,c}\hat q_{j,d}}
       {\phi(\hat\tau_{j,c})\,\phi(\hat\tau_{j,d})}$$

for any sample at any $N$ — not asymptotically. Verified at $N$ = 200, 5,000 and 100,000: gaps of
8e-15, 2e-13 and 4e-12, growing only with floating-point accumulation in the crossproduct.

This is worth more to Module 4 than it looks. It is a **dependency-free, machine-precision unit
test available on day one** — no lavaan, no Monte Carlo, no tolerance to argue about. A failure is
a coding error, full stop. It also localizes where sampling error can possibly live: all of the
discrepancy between an empirical and a closed-form Γ sits in the cross-variable entries, where the
sample proportion *estimates* $\Phi_2$ rather than equalling it. That is what makes the tolerance
split in Section 7.2 principled rather than arbitrary.

#### 4.5.2 The correlation entries (Stage 2)

**This has to account for the two-stage structure of the Stage 1 estimator, not just treat each
statistic as its own independent M-estimator.** Module 1 estimates the thresholds; Module 2 then
holds them fixed while estimating each pairwise correlation. That two-module split *is* the
two-stage structure, which is why the propagation term below crosses the boundary between them
rather than living inside either. That means a pairwise correlation estimator
is a function of the (random) estimated thresholds: $\hat\rho_{jk} = \hat\rho_{jk}(\hat\tau_j,
\hat\tau_k)$, so its influence function is **not** simply the direct bivariate-likelihood score
treated as if the thresholds were fixed and known. It has to include the propagated contribution
from threshold estimation. **The expected form** — stated as a target for the derivation to
confirm, not as an already-proven result to code from directly — is:

$$\psi_{\rho_{jk}} = \psi_{\rho_{jk}}^{\text{direct}} + \frac{\partial \rho_{jk}}{\partial \tau_j} \psi_{\tau_j} + \frac{\partial \rho_{jk}}{\partial \tau_k} \psi_{\tau_k}$$

The exact sensitivity terms must actually be derived from the stacked two-stage estimating system
(thresholds' first-order conditions together with each pairwise correlation's), not assumed from
this chain-rule sketch — this is the single most mathematically consequential piece of the
project, and a provisional expression should never get treated as a settled theorem just because
it's displayed with an equals sign. Where $\psi^{\text{direct}}$ is the usual
inverse-Hessian-times-score term from the bivariate log-likelihood (thresholds held fixed), and the
other two terms propagate the uncertainty in the first-stage threshold estimates through the
second-stage correlation estimate. This propagated term is exactly what generates a meaningful
share of Γ's off-diagonal structure — omitting it doesn't just shrink Γ slightly, it's a different
(wrong) estimator, since it's dropping a real source of correlation between $\hat\tau_j$ and every
$\hat\rho_{jk}$ involving variable $j$.

**How to actually compute $\partial \rho_{jk}/\partial \tau_j$:** apply the implicit function
theorem to the stage-2 estimating equation. $\hat\rho_{jk}$ solves the bivariate-likelihood score
equation $S_\rho(\rho; \tau_j, \tau_k) = 0$ (holding thresholds fixed, as Module 2 does), so

$$\frac{\partial \rho_{jk}}{\partial \tau_j} = -\frac{\partial^2 \ln L_{\text{bi}} / \partial \rho\, \partial \tau_j}{\partial^2 \ln L_{\text{bi}} / \partial \rho^2}$$

evaluated at the estimated values — this is the standard two-step M-estimator sensitivity
argument, not something specific to this problem. Building this out uses the same ingredients as
Section 4.7 (Plackett's identity for the $\rho$-derivatives, the univariate-density-times-
conditional-CDF form for the threshold derivatives), but don't assume in advance that the full
expression collapses into one clean closed form beyond the simplest tetrachoric (single-threshold)
case — for a general polychoric pair with multiple thresholds the numerator and denominator above
are each sums of several such terms, and the exact combination needs to be worked through and
checked, not asserted. The single-threshold case *is* now solved in closed form (4.5.3), which
gives this derivation a reduction target but says nothing about its general structure.

lavaan's own `lav_muthen1984.R` implementation is a useful reference here: it explicitly separates
univariate and bivariate stages and constructs the ACOV through delta-method machinery connecting
them, which is the same two-stage structure this section is describing.

**The result.** Derived and verified in the derivation file. With $\pi_{rs}$ the
cell probabilities and $u_{rs} = \partial_\rho\pi_{rs}/\pi_{rs}$,

$$\psi_{\rho,i} = \frac{1}{I_{\rho\rho}}\left[u_{r_i s_i}
  - \sum_c I_{\rho a_c}\psi_{a_c,i} - \sum_d I_{\rho b_d}\psi_{b_d,i}\right],
\qquad
\frac{\partial\rho}{\partial a_c} = -\frac{I_{\rho a_c}}{I_{\rho\rho}}$$

where $I_{\rho\rho} = \sum_{rs}(\partial_\rho\pi_{rs})^2/\pi_{rs}$ and $I_{\rho a_c} =
\sum_{rs}(\partial_\rho\pi_{rs})(\partial_{a_c}\pi_{rs})/\pi_{rs}$. Both expectations in the
two-step formula collapse for the same reason: the term carrying a second derivative of
$\pi$ sums to $\partial(1) = 0$ over the table. A threshold couples only to the two rows it
separates, so $I_{\rho a_c}$ reduces to a sum over one variable's categories rather than the
whole table — the $K>2$ analogue of the tridiagonal sparsity in 4.5.1.

**Which information normalises $\psi_\rho$ — settled 2026-09-20, and v1 does *not* use the
expected form written above.** The two-step expansion needs $A = -E[\partial g_i/\partial\rho]$
and $G_{a_c} = E[\partial g_i/\partial a_c]$, and *any consistent estimator of those
expectations* gives a valid influence function. Writing both at once by taking the weight
matrix as an argument,

$$\hat I_{\rho\rho}(W) = \sum_{rs} W_{rs}u_{rs}^2,
\qquad
\hat I_{\rho a_c}(W) = \sum_{rs} W_{rs}\,u_{rs}\,\partial_{a_c}\pi_{rs}/\pi_{rs},$$

$W = \pi$ recovers the expression above (the *expected* information) and $W = \hat\pi = n/N$
gives the *empirical*, outer-product or BHHH, form in which both quantities are literal
per-observation crossproducts. **lavaan uses $W = \hat\pi$**, which was read from its source
rather than inferred: `muthen1984()` returns `b_inv %*% crossprod(sc) %*% t(b_inv)` with
`a22[i,i] <- sum(sc_cor[,i]^2)` and `a21[...] <- crossprod(sc_cor[,pstar], pair$dx_th_i)`.

**v1 ships $W = \hat\pi$**, by the same compatibility argument 4.5.4 makes one level up: the
two differ at $O(N^{-1/2})$, and on `case_D` at $n = 5000$ the expected form disagrees with
`lavInspect(fit,"gamma")` by 0.051 absolute and 6.56 relative against 1.2e-09 relative for
the empirical form — a difference that propagates to 2–3% on standard errors, against 7.2's
$10^{-3}$. Corroboration that this tracks the mathematics rather than one package: lavaan
and Mplus agree on `case_D`'s $z$ statistics to ~1.1e-4 relative, so Mplus is on the same
side.

Three things keep this from being a licence to substitute estimators freely. It is *not* a
correctness question — the derivation file (`@sec-normaliser`) measures the
empirical-versus-expected gap shrinking at $N^{-1/2}$ while the naive-versus-expected gap
sits at 8.2% and does not move, so the two normalisers estimate the same $\Gamma$ and
dropping the propagation term does not. The two **coincide exactly** whenever the table is
saturated, which covers every univariate threshold — hence the threshold block of $\Gamma$ is
untouched by this choice — and every binary pair, so Cases A and B cannot detect it and only
$K > 2$ separates them. And $\hat\rho$ is unaffected: the information enters Module 2 only as
a Fisher step scaler, whose fixed point is the same for any positive scaler, so Module 2
keeps the **expected** form there (its positivity is guaranteed) and reports the empirical
one. The choice changes what Module 2 *reports*, not what it *finds*. The fixed point is
not the whole story: the iteration converges only where $|1 - H/c| < 1$, and on a sparse
table the empirical form breaks that where the expected one does not (derivation
`@sec-emptycells`, 2026-09-26). The expected form is summed over cells whose *computed*
probability is positive, and the score runs over **observed** cells only. An empty cell
contributes $0\cdot\log\pi = 0$, and near $\rho = \pm1$ its computed $\pi$ is rounding noise.
Module 2 searches only **usable** points (derivation `@sec-rootfind`, decided 2026-09-29):
$$\mathcal R_P = \{\rho : \min_{\mathcal O}\tilde\pi_{rs}(\rho) > \delta_P\},$$
where $\mathcal O$ is the observed cells and $\delta_P$ the absolute error of a computed cell
probability (Module 2, below). Usability replaces `@sec-emptycells`' earlier rule "every
observed cell's computed $\pi > 0$" for the search: at a point that is positive but not
usable, some observed $\tilde u_{rs}$ has relative error above 1, so the sign of the
computed score is not determined by the data. Usable implies positive.

`dwlsmv_pairinfo()` therefore takes $W$ as an argument rather than hard-coding $\hat\pi$:
$W = \pi$ remains the unit-test oracle that ties the shipped routine to the compact
$w$-weighted form validated against numerical differentiation, and it is the candidate
fallback if the empirical form proves worse behaved under sparse cells, where $\hat\pi$ gives
an empty cell zero weight and $\pi$ retains it. That question is open — see
E-001.

**Module 2's API, settled by Gate 3.** The obligation is lighter than it looks: **no numerical second derivatives
are required.** $I_{\rho\rho}$ and $I_{\rho a_c}$ are expected informations, closed-form in
the cell probabilities and their first derivatives. Module 2 returns, per pair:

| Quantity | Shape | Used for |
|---|---|---|
| $\hat\rho_{jk}$ | scalar | the statistic vector $s$ |
| $u_{rs}$ | $K_j \times K_k$ | the direct term |
| $I_{\rho\rho}$ | scalar | the normalisation |
| $I_{\rho a_c}$, $I_{\rho b_d}$ | $K_j-1$, $K_k-1$ | the propagation weights |

Module 4 combines those with the Stage 1 threshold influence functions, so the Section 5
boundary is unchanged: Module 2 estimates, Module 4 assembles.

**$\psi_\rho$ is a lookup table.** It depends on observation $i$ only through its cell
$(r_i, s_i)$, since the threshold influence functions are themselves category lookups. A
$K_j \times K_k$ table per pair, built once. Together with 4.5.1's staircase, this means no
per-observation influence function is computed anywhere in Γ's construction.

Build in two sub-phases of increasing complexity (Module 4 in the breakdown below):

- *4a — purely binary (tetrachoric) indicators.* The bivariate normal integral for a 2×2 table has
  only 4 cells, and Plackett's identity (Section 4.7) keeps the derivatives simple. Easiest case to
  debug $\psi_i$ — including its threshold-propagation term — against, and it has a verified
  closed-form target: see 4.5.3.
- *4b — purely ordinal (polychoric) indicators, ≥3 categories.* Same structure, more
  categories/thresholds per pair, so more propagation terms per correlation.

(Sub-phase 4c, mixed continuous/ordinal, moves to the v2 continuous-indicator extension along with
the rest of the mixed-type machinery — see Section 3.)

#### 4.5.3 Sub-phase 4a's acceptance test: the binary case in closed form

A $2\times2$ table is **saturated** — 3 free cell probabilities, 3 parameters
$(\tau_1, \tau_2, \rho)$. Two things follow that hold for binary and for binary only. The
two-stage estimator coincides *exactly* with joint ML (the saturated MLE reproduces the observed
table, hence its margins, hence $\Phi(\hat\tau_j) = \hat q_j$), so $\hat\rho$ is a direct inversion
of $\Phi_2(\hat\tau_1,\hat\tau_2;\hat\rho) = \hat p_{00}$ with no optimization. And $\hat\rho$ is
therefore a smooth function of the trivariate sample mean $(\hat q_1, \hat q_2, \hat p_{00})$, so
the delta method delivers its influence function in closed form, threshold propagation included:

$$\psi_{\rho,i} = \frac{1}{\phi_2(\tau_1,\tau_2;\rho)}
\Big[(Z_{1i}Z_{2i} - p_{00})
 - \Phi\!\Big(\tfrac{\tau_2-\rho\tau_1}{\sqrt{1-\rho^2}}\Big)(Z_{1i} - q_1)
 - \Phi\!\Big(\tfrac{\tau_1-\rho\tau_2}{\sqrt{1-\rho^2}}\Big)(Z_{2i} - q_2)\Big]$$

The shape is the one 4.5.2 predicts: a direct term in the joint cell, minus two first-stage
correction terms, each weighted by a conditional normal probability — the weight on the $q_1$
correction is $P(y_2^* \le \tau_2 \mid y_1^* = \tau_1)$, the rate at which a shift in the first
threshold propagates into the joint cell.

**This is genuine independent confirmation for Gate 3, not a shortcut around it.** It was derived
by applying the implicit function theorem to the saturated inversion, which is a *different route*
from the stacked-estimating-equations argument 4.5.2 prescribes. Two routes reaching the same
expression is exactly the corroboration Gate 3 is asking for. Verified numerically: the three
derivatives match numerical differentiation to 8–9 significant digits, and the implied
$N\cdot\text{Var}(\hat\rho)$ of 1.96728 matches a Monte Carlo value of 1.95368 with a Monte Carlo
standard error of 0.0714 ($z = 0.19$).

**Do not let this leak into 4b.** Over-identification for $K > 2$ means two-stage and joint ML
genuinely differ there, so this closed form is a binary-only oracle. Confirmed on a $4\times3$
table: the two estimators separate, by roughly 70× the $2\times2$ case's optimizer-tolerance
residual on a table built and perturbed identically. The general derivation must *reduce* to the
boxed expression at $c_1 = c_2 = 2$ — that is 4b's acceptance test — but must not be assumed to
share its structure.

#### 4.5.4 The construction rule: empirical, never plug-in

**The influence-function route is a compatibility requirement, not just a convenience.** It would
be tempting, once the algebra of 4.5.1–4.5.3 is done, to skip the per-observation $\psi_i$ for the
parts of Γ that have closed forms and plug estimates into the formula instead. Don't.
**lavaan builds $\hat\Gamma$ empirically, by averaging estimated influence functions** — verified
by comparing its `gamma` against a closed-form Γ evaluated at lavaan's own estimates on a
four-indicator design at $n = 5000$: the maximum absolute discrepancy was 0.058, against a
simulated reference distribution for that discrepancy with median 0.057 and a 90% interval of
[0.039, 0.087], i.e. the 55th percentile. An analytic plug-in would have agreed to machine
precision. Since the two constructions differ at order $N^{-1/2}$, a Stata implementation using
the plug-in would produce standard errors that disagree with the primary compatibility target by a
visible amount at realistic $N$, and that would be reported as a bug.

**Production code uses $\hat\Gamma = \frac{1}{N}\sum_i \hat\psi_i\hat\psi_i'$. The closed forms in
4.5.1 and 4.5.3 are unit-test oracles, never the production path.**

Exploit the symmetry of Γ and compute only its lower triangle.

#### 4.5.5 Centering, as a diagnostic rather than a required correction

In exact arithmetic, a correctly constructed $\psi_i$ should already sum to (numerically) zero
across observations: both the threshold and pairwise-correlation stage-1/stage-2 estimating
equations are satisfied exactly at $\hat\tau$, $\hat\rho$ by construction, so $\sum_i \psi_i
\approx 0$ should hold automatically rather than needing to be enforced by explicitly centering
($\hat\Gamma = \frac{1}{N}\sum_i (\psi_i - \bar\psi)(\psi_i - \bar\psi)'$ vs. the uncentered
$\frac{1}{N}\sum_i \psi_i \psi_i'$ should give numerically indistinguishable results). Precisely
because of that, implement the check rather than skip it: compute $\bar\psi$ and treat it as a
diagnostic — if it isn't close to zero (beyond ordinary floating-point roundoff), that's a signal
something in $\psi_i$'s construction is wrong (a missing propagation term, a derivative error), not
a reason to quietly switch to the centered formula and move on. Decide explicitly which formula
production code uses, but don't treat centering as a fix for a nonzero mean that should be treated
as a bug report instead.

**Split that diagnostic by block — it is far more sensitive that way.** The two halves of $\psi$
are not equally exact, and a single tolerance covering both is set by the looser one, which hides
real bugs in the tighter one. Because $\hat q_{j,c}$ is literally the sample proportion, the
threshold entries of $\hat\psi$ have sample mean **exactly** zero, to machine precision, with no
centring step and no dependence on any optimizer — verified at 1e-16 and below. The correlation
entries are only approximately zero, since $\hat\rho_{jk}$ solves its score equation numerically,
so their mean is zero only to Module 2's stopping rule, $|\bar S(\hat\rho)| \le \tau_S(\hat\rho)$
(Module 2, below). So: assert machine precision on the threshold entries, and tie the correlation
entries' tolerance explicitly to that rule. A threshold-entry mean of 1e-8 would pass a single
blanket check and is already a bug.

**The one exception: a one-empty-cell pair repaired in Module 2 (decided 2026-10-06).** Its $\hat\rho^+$
solves the score on the repaired table, while ψ is scored on the real observations. So its
correlation entry has $\bar\psi_\rho = \bar S(\hat\rho^+)/\hat I_{\rho\rho} \ne 0$, at $O(1)$, by
construction (`@sec-o6`, *The v1 decision*). That entry is excluded from the bug-report reading
and reported as the repaired pair's value. Its threshold entries keep the machine-precision
assertion. Every other correlation entry keeps the τ_S tie.

#### 4.5.6 Numerical hazards: sparse cells have two distinct failure modes

The densities entering both the direct score and the propagation term of 4.5.2 can get very small
for cells with low counts, which risks overflow or division-by-zero in the sensitivity ratio. Tie
this into the existing fail/warn-by-default sparse-cell handling (Module 2) rather than letting it
surface as a silent NaN/Inf deep inside Γ's construction.

But Module 2's sparse-cell handling frames this as one problem. It is two — three, counting the
near-boundary case — and they fail differently and break different modules.

| Mode | Trigger | Symptom | What breaks |
|---|---|---|---|
| **Boundary** | First or last category empty; $\hat q \in \{0,1\}$ | $\hat\tau = \pm\infty$, $\phi(\hat\tau) = 0$ | Loud, easy to detect, kills that threshold outright |
| **Interior** | An interior category empty; $\hat q_{j,c} = \hat q_{j,c+1}$ | Two thresholds collide at a *finite* value, $\phi$ finite, nothing overflows | Two $\psi$ columns become identical, so $\hat\Gamma$ is rank-deficient — **silently** |
| **Near-boundary** | Category present but very sparse | $\phi(\hat\tau)$ small but nonzero | $\psi$ and hence $\hat\Gamma$ inflated, no error raised at all |

The interior mode is the one to design against, because its consequence is asymmetric and
counterintuitive: **DWLS point estimation is unaffected**, since $W_D = \text{diag}(\Gamma)^{-1}$
never inverts the full matrix, so Module 7 converges and returns a plausible $\hat\theta$. It is
Modules 8 and 9 that fail, because the robust sandwich and the scaled-and-shifted statistic both
need Γ in full. In other words a Module 1/2 data condition surfaces as wrong standard errors and a
wrong test statistic with correct point estimates — precisely the "runs without error but is
quietly wrong" failure mode Section 6's risk matrix exists to catch, and one that neither Module 7's
admissibility check nor a comparison of point estimates against lavaan would reveal.

**Therefore:** Module 1's category validation must distinguish empty *interior* categories from
empty boundary categories rather than treating "degenerate category" as one class — and v1's
Module 1 aborts on either, so the interior route above cannot reach Γ.

**Exact rank deficiency of $\hat\Gamma$ in v1 has causes that can be named, and Module 10 detects
each from its cause rather than by a singular-value cutoff** (derivation `@sec-gamma-rank`, added
2026-09-24; every statement below is asserted there).

*Response pattern.* $n_{\text{pat}}$ is defined on the matrix $Y$ exactly as
`dwls_sample_stats()` receives it: after listwise deletion, restricted to the $p$ indicators whose
thresholds and polychorics make up $s$, in $K$'s column order, coded $0,\dots,K_j-1$. A pattern is
one row of that matrix, and $n_{\text{pat}}$ is the number of distinct rows. Nothing else counts:
variables outside the $p$ indicators, observations removed by listwise deletion, and declared
categories that are never observed. It is not called $P$, because $P$ is already the number of
pairs (Section 4.2, `map.P`). This is the right object because $\hat\psi_i$ is a function of
exactly that row: $\psi_\tau$ is looked up by category and $\psi_\rho$ by observed cell.
`dwlsmv_npatterns(Y, K)` (Module 10) counts it. It takes the same $(Y, K)$ as
`dwls_sample_stats()` and applies the same validation, aborting with rc 3498 on a missing value, a
non-integer code or an out-of-range code rather than counting a pattern Stage 1 would reject. A
`.` counted as a code would inflate $n_{\text{pat}}$ and make the check below pass wrongly.

**This definition and the bounds below are v1 only.** A continuous indicator makes $\psi$ depend on
a continuous value, so the pattern set is no longer finite. The v2 mixed-type extension must
re-derive them, not reuse them.

| Cause | Exact statement | Detected by |
|---|---|---|
| Too few patterns | $\operatorname{rank}\hat\Gamma \le \min(D,\, n_{\text{pat}}-1)$ when every pair's score vanishes at $\hat\rho$ | $n_{\text{pat}} - 1 < D$, by counting |
| Threshold aliasing | $\operatorname{rank}\hat\Gamma_{\tau\tau} = \operatorname{rank}[\mathbf 1 \mid Z] - 1$, where $Z$ holds the cumulative indicators $\mathbb 1\{Y_{ij} < c\}$ | $\operatorname{rank}[\mathbf 1 \mid Z] < T + 1$ |

The interior-empty mode above is one case of threshold aliasing: two nested columns of one variable
coincide. The general case needs no empty category. In $\mathbb 1\{Y_j < c\} = \mathbb 1\{Y_k < d\}$
on every row, every category can be non-empty, $\hat\rho_{jk}$ interior and $n_{\text{pat}} - 1 \gg D$,
and $\hat\Gamma$ is still exactly singular. **This route does reach Γ in v1, silently**, which was
measured in production Stata on 2026-09-24 as rank 15 of 16.

*Where the pattern count's "−1" comes from, and when it is lost.* Production builds $\hat\Gamma$
uncentred, so strictly $\operatorname{rank}\hat\Gamma \le \min(D, n_{\text{pat}})$. The
$(n_{\text{pat}})$-th eigenvalue of the scale-free $\hat\Gamma$ is bounded by $\sum_k z_k^2$, with
$z_k = |\bar S_k|/\sqrt{J_k}$ per pair (`@sec-gamma-budgets`). For a pair stopped by Module 2's
rule, $|\bar S| \le \tau_S$, this is at most $4\tau_S^2/J$. On every measured table that is far below
the eigensolver's floor, but the bound depends on the table's smallest observed cell probability, so
it is measured, not a priori (`@sec-gamma-rank`; re-derived 2026-10-02 for option B). A pair whose
$\hat\rho$ does not satisfy $|\bar S| \le \tau_S$ on recomputation breaks this premise: statuses 2,
3, 4 and 8, and status 7 unless its cap hit while probing a root already found. `converged = 0` is
the wider set: statuses 5 and 6, and status 7 from the endpoint check's cap, stopped on the rule.
Under the v1 policy (Module 2, *Which pairs enter `s` and Γ̂*) only status 1 and status 6 by the
endpoint check reach Γ̂; every pair breaking the premise aborts Stage 1.
That is also the case in which
$\hat\Gamma$ is *wrong while well conditioned*, which a spectrum cannot detect. Module 10 therefore
reports the centring residual $\bar\psi_\rho$ alongside the spectral quantities. $\bar\psi$ detects
a failed score equation, never a missing propagation term (Section 4.5.5).

A category of count 1 and a near-duplicate pair are **not** exact degeneracies. They are
ill-conditioning. v1 sets no numerical cutoff for them (decided 2026-10-05, below).

Module 10's rank/condition check on Γ is **not** a precondition of Modules 8–9 (decided
2026-09-23, superseding an earlier version of this paragraph). "Need Γ in full" is not "invert
Γ": the sandwich of Section 4.10 inverts only $\Delta' W_D \Delta$, and $T_3$ takes traces of
$U\Gamma$, so a singular positive-semidefinite Γ causes no numerical failure in either.
**Handling in v1 (decided by the user 2026-10-05): no Γ̂ cutoff.** The scale-free and raw
ratios are reported as values. The mechanisms that act are Stage 1's refusals: the bound-pair
policy, the empty-category abort, and assembly when I_ρρ underflows. Module 3's Δ rank check at
identification also acts. Module 10's structural rows, centring included, are reported in the
notes block and do not refuse. The evidence is the threshold-aliasing Monte Carlo, under both a
misspecified and a correctly specified design, and the paired run re-done after the near-±1 fix.
In none did the ratio identify an inference failure, down to a median scale-free ratio of about
5e-3. Severe near-aliasing stays open and is not closed by inventing a cutoff
(E-002). **No centring refusal either (decided 2026-10-05):** under the default pair
policy the pair status refuses first, and under `keepbound` a kept pair flags by construction. A
repaired one-empty-cell pair is excluded from the row's $z$ and reported separately (§4.5.5; decided
2026-10-06). The bread's numerical-rank rule is unaffected.

### 4.6 Model-implied moments $\sigma(\theta)$

For a measurement model with structural regressions among factors:

$$y^* = \Lambda \eta + \epsilon, \qquad \eta = B\eta + \zeta$$

$$\Phi_\eta = (I - B)^{-1} \Psi (I - B)^{-\top}$$

$$\Sigma_{y^*} = \Lambda \Phi_\eta \Lambda' + \Theta$$

Because all indicators are ordinal/binary with the delta parameterization (Section 4.3), introduce
the standardized correlation matrix explicitly rather than leaving the scale transformation
implicit:

$$R_{y^*} = D^{-1/2} \Sigma_{y^*} D^{-1/2}, \qquad D = \text{diag}(\Sigma_{y^*})$$

**Under the delta parameterization $\Theta$ is derived, not free, and the standardization above
is a no-op.** Fixing $\operatorname{Var}(y_j^*) = 1$ means
$\theta_{jj} = 1 - (\Lambda\Phi_\eta\Lambda')_{jj}$, so $\Sigma_{y^*}$ already has unit diagonal
and $R_{y^*} = \Sigma_{y^*}$ identically — verified to machine zero. Leaving $\Theta$ free and
standardizing afterwards is not a harmless alternative: scaling $\lambda_j$ by $c$ and
$\theta_{jj}$ by $c^2$ leaves $R_{y^*}$ unchanged, so correlations alone cannot separate them.
Implementations should compute $R_{y^*}$ as $\Sigma_{y^*}$ and treat a request to free a residual
variance as an error, not an option.

**The model-implied polychoric correlation for pair $(j,k)$ is simply the corresponding entry of
$R_{y^*}$** — i.e. $\sigma_{\rho_{jk}}(\theta) = \text{Corr}(y_j^*, y_k^*)$. This is a point worth
stating explicitly because it's easy to conflate: thresholds determine the *marginal category
probabilities* for each indicator, but they do not enter the population polychoric correlation
itself — that's purely a function of the latent-response correlation. Thresholds are separate,
directly-parameterized elements of $\sigma(\theta)$, not something the correlation mapping needs to
account for. The full pipeline:

```text
θ  →  Σ_y*  →  R_y*  →  model-implied polychoric correlations (= entries of R_y*)
                          + threshold parameters (directly from θ)
      =  σ(θ), indexed exactly as s is (Section 4.2)
```

This mapping is what Module 6 (derivatives) and Module 7 (the optimizer) are built against. Stated
as an explicit partitioned vector, matching $s$'s index map (Section 4.2) directly:

$$\sigma(\theta) = \begin{bmatrix} \tau(\theta) \\ \text{vech}_{\text{off}}\big(R_{y^*}(\theta)\big) \end{bmatrix}, \qquad \Delta = \begin{bmatrix} \partial \tau / \partial \theta' \\ \partial\, \text{vech}_{\text{off}}(R_{y^*}) / \partial \theta' \end{bmatrix}$$

where $\text{vech}_{\text{off}}$ stacks the strictly-lower-triangular (off-diagonal) entries of
$R_{y^*}$ in the same order Section 4.2 fixed for the correlation entries of $s$. This gives
Module 6 an unambiguous target: $\Delta$ is just the derivative of this same partitioned vector.

### 4.7 Bivariate-ML derivative identity (for Module 2 performance work)

Plackett's formula: $\partial \Phi_2 / \partial \rho = \phi_2(\tau_1, \tau_2; \rho)$ — the
derivative of the bivariate normal CDF with respect to $\rho$ is the bivariate normal *density*
evaluated at the thresholds. Derivatives with respect to the thresholds involve a univariate
density times a conditional normal CDF term, not the same simplification — needed both for an
analytic-derivative version of Module 2's correlation estimation, and for the $\partial \rho_{jk} /
\partial \tau_j$ terms in Section 4.5's influence function.

### 4.8 Identification and parameter mapping

Two genuinely separate identification questions have to be resolved into a concrete free/fixed
parameter table — don't let them get tangled together, since they're independent choices:

**(a) Latent factor metric** — how each factor $\eta$'s scale is fixed:
- marker-variable ($\lambda_{1f} = 1$ for a chosen indicator), or
- unit-variance / standardized-latent-variable, which fixes $\psi_{ff} = 1$ —
  **the disturbance variance, not the total variance.** See below; the two are
  the same number for an exogenous factor and different ones for an endogenous
  factor, so the distinction has to be written into the table rather than left
  to the phrase "unit variance".

**(b) Ordinal response-variable scale** — already fixed at the delta parameterization for v1
(Section 4.3), i.e. $\text{Var}(y_j^*) = 1$ for every ordinal indicator; not a per-factor choice,
a global one for this spec.

These are different kinds of identification (one about the latent factors, one about the observed
ordinal indicators' underlying scale) and belong in separate rows of the parameter table, not
merged into one "identification" concept:

```text
factor metric:               marker-variable | unit-variance   (per factor)
ordinal response scale:      delta (fixed for v1)               (global)
residual variances (Theta):  DERIVED, never free                (per indicator)
thresholds:                   free | tied across groups         (per threshold; ties only with groups)
structural paths (B):        free | fixed-to-zero               (per directed factor pair; no cycles in v1)
exogenous-factor covariances: free | fixed-to-zero               (per pair of exogenous factors)
disturbance variances (Ψ):    free | fixed-to-1 by the metric    (per endogenous factor)
disturbance covariances (Ψ):  free | fixed-to-zero               (per pair of endogenous factors)
total factor variances (Φ_η): DERIVED, never free                (per factor)
latent-response intercepts ν: fixed-at-0 | free | tied          (per indicator; groups g ≥ 2 only)
response scales ω:            fixed-at-1 | free | tied          (per indicator; groups g ≥ 2 only)
factor intercepts α:          fixed-at-0 | free | tied          (per factor; groups g ≥ 2 only)
```

The last three rows exist only with multiple groups (0.3.0; *The canonical model specification*,
**Multiple groups**, below). Under groups the delta convention is unchanged: the standardized
$y^*_{gj}$ has unit variance in every group, and Θ stays derived. $\omega_{gj}$ is the scale
linking the invariant (unscaled) thresholds to the standardized ones, not a variance of $y^*$
(derivation `@sec-invariance-id`). $\alpha$ is the factor intercept in
$\eta = \alpha + B\eta + \zeta$, never the mean: the mean $\kappa_g = (I-B)^{-1}\alpha_g$ is derived,
as $\Phi$ is from $\Psi$. With structural paths refused under groups in 0.3.0, $\kappa_g = \alpha_g$
in every admitted model.

**What "unit variance" fixes for an endogenous factor.** The quantity set to 1
is $\psi_{ff}$, the diagonal of $\Psi$ — the variance of the *disturbance*
$\zeta_f$. It is not $\phi_{ff}$, the diagonal of
$\Phi_\eta = (I-B)^{-1}\Psi(I-B)^{-\top}$. For an exogenous factor the two are
the same number and the choice is invisible; for a factor that is predicted
they differ by the explained variance, and the difference rescales every
loading on that factor. Module 3 must therefore write the constraint into the
$\Psi$ row of the table, never the $\Phi_\eta$ row — and $\Phi_\eta$ is derived
throughout, which the last row above says explicitly. A user request to fix
an endogenous factor's *total* variance is not supported in v1 and is not the
same request.

This is a measured convention, not a reading of the literature. On Case E
(`fixtures/data/case_E.csv`, two factors, $F_2 \leftarrow F_1$) under lavaan
0.7.2 with `std.lv = TRUE`, the parameter table carries `f2 ~~ f2` as **fixed**
at 1, and the implied $\operatorname{Var}(F_2)$ is $1.2800 = \beta^2 + 1$ at
$\hat\beta = 0.5291$ — so the total variance is not 1 and the disturbance
variance is. In the same fit $\operatorname{Var}(F_1) = 1$ exactly, $F_1$ being
exogenous, and in the correlated two-factor fit of the same data both factor
variances are 1: that is why no design in this project before Case E could have
settled the question, and why Case E exists.

Stata's own `sem` reaches the same convention and makes it syntactically
explicit, which matters because Module 11's surface syntax is meant to accept a
user's existing `sem` model statement. `var(F2@1)` is *rejected* there once
`F2` is a latent dependent variable — `invalid specification of variance of
'F2'; 'F2' is a latent dependent variable`, rc 498 — and the unit-variance
metric must be written `var(e.F2@1)`, on the disturbance. Fitted, it likewise
returns $\operatorname{Var}(F_2) = \beta^2 + 1$.

Effects coding is a third metric method and is out of v1 scope; it does not
make the convention ambiguous, because it constrains the *loadings* on a factor
to average 1 and leaves that factor's variance free, so the
total-versus-disturbance question does not arise under it at all.

**The equivalent-models identity, and what it is for.** A correlated
two-factor model and an $F_2 \leftarrow F_1$ model on the same indicators are
reparameterizations of each other: each spends one parameter on the link
between the factors and each maps onto the other without changing
$\sigma(\theta)$. Fitted to the same data they must return the same number of
free parameters, the same $df$, the same fitted statistics and the same $T_3$.
This holds under either factor metric and regardless of how the question above
was resolved, which is what makes it an independent check on Module 7 rather
than a restatement. It is asserted at Section 7.2's test-statistic tier;
what caps the agreement in practice is that the two are separately converged
optimizations, measured at ~1e-09 relative.

**A parameterization detail worth stating explicitly, not just implying:** once any factor is
predicted by another ($\eta = B\eta + \zeta$, Section 4.6), its covariance with everything else is
*implied* by $\Phi_\eta = (I - B)^{-1}\Psi(I-B)^{-\top}$, not independently free. The parameter
table above reflects that: only exogenous-factor covariances and endogenous-factor disturbance
(co)variances are free parameters; a covariance between an exogenous and an endogenous factor, or
between two endogenous factors, is a derived quantity, never a row Module 3 should let a user set
directly. Getting this distinction wrong — e.g. letting $\text{Cov}(F_1, F_2)$ be independently
free in a model where $F_1 \to F_2$ is also specified — overparameterizes the model and Module 3's
admissibility/identification checks should catch it rather than let optimization silently absorb
the redundancy.

**v1 is recursive: B may contain no cycle (decided 2026-09-23).** `dwls_identify()` rejects any
specification whose `beta` directives — free, or fixed at a nonzero value; a fixed zero is no path
— lead from a factor back to itself, whether through a reciprocal pair or a longer loop. The
relevant matrix is $(I-B)$, never $B$: a recursive $B$ is nilpotent and always singular, yet $(I-B)$
is unit triangular and always invertible. A cyclic $B$ with $(I-B)$ nonsingular is computable, and
refusing it is a scope decision rather than a numerical one, for two reasons:
- **The rank check cannot judge it.** At the default starts ($B = 0$, $\Psi = I$),
  $\mathrm{d}\Phi = \mathrm{d}B + \mathrm{d}B'$, so a reciprocal pair's two Δ columns are
  identical. The rank check would therefore report *every* such model as overparameterized,
  identified or not, and give the wrong reason.
- **A feedback loop needs a condition nothing checks.** It needs a stability condition — spectral
  radius of $B$ below 1 — that neither the rank check nor the admissibility conditions of Section
  5, Module 7 test.

Supporting non-recursive models would need nonzero default starts for cyclic paths, a stability
row in the admissibility check, and a fixture with an identified non-recursive model; it is a v1b
candidate, not a v1a feature. The check is a transitive closure of $B$'s directed graph;
`tests/unit/test_identify.do` rejects 16–18 and check 11.

**The concrete test, a numerical rank check on Δ at the starting values.** Build
Δ before optimizing and check $\operatorname{rank}(\Delta)$ equals the number of free parameters.
Freeing $\psi_{12}$ alongside $\beta$ — the exact error above — drops the rank by exactly one, and
the null-space direction is $\pm 1/\sqrt 2$ on $(\beta, \psi_{12})$ with numerical zero elsewhere,
naming the offending pair. The singular-value ratio is machine zero at $\psi_{12} = 0$ and about
$10^{-12}$ elsewhere, against a threshold of $10^{-8}$ and genuine singular values of order
$10^{-1}$, so this is a reliable detector rather than a knife-edge numerical judgment. A
deficiency of $k$ means $k$ redundant parameters and the null space says which.

**The deficiency holds at every $\psi_{12}$; what "at the starting values" buys is a clean null
direction.** Under the disturbance-variance convention settled above,
$\Phi_{12} = \beta + \psi_{12}$ and $\Phi_{22} = \beta^2 + 2\beta\psi_{12} + \psi_{22}$, and with
$\psi_{22}$ fixed the map $(\beta, \psi_{12}) \mapsto (\Phi_{12}, \Phi_{22})$ has Jacobian
determinant $-2\psi_{12}$. At $\psi_{12} = 0$ it is singular: $\mathrm{d}\beta = -\mathrm{d}\psi_{12}$
holds both entries of $\Phi$ fixed, the null direction is exactly $\pm 1/\sqrt2$ on
$(\beta, \psi_{12})$, and the two Δ columns are the same arithmetic, so the ratio is machine zero.
At $\psi_{12} \ne 0$ the pair can move $\Phi_{22}$, which re-opens the factor-scale indeterminacy
on $F_2$ that the metric constraint existed to close: **the deficiency is still exactly one**, and
the null direction now rescales $F_2$'s loadings with $(\beta, \psi_{12})$ absorbing it. The
model is unidentified everywhere and the rank check sees it everywhere. Since $0$ is a disturbance
covariance's starting value, the check as specified also gets the null direction that names only
the offending pair. Measured on case E: ratio $6.0\times10^{-17}$ at the default starts and
$2.1\times10^{-17}$ at $\beta = 0.4$ with $\psi_{12} = 0$; $4.8\times10^{-13}$ and
$8.9\times10^{-13}$ at $\psi_{12} = 0.1$ and $0.3$, with the null vector within $4.3\times10^{-12}$
of the analytic loading rescaling (`tests/unit/test_delta.do` checks 6–8). The smallest genuinely
nonzero singular value across four identified configurations is $0.226$. (Quote the *ratio*, not
the smallest singular value — $1.2\times10^{-16}$ at the default starts — since the check
thresholds on $\sigma_{\min}/\sigma_{\max}$ and the two differ by the largest singular value,
about 2.1 here.)

Until 2026-09-23 this paragraph said a Δ evaluated at a nonzero $\psi_{12}$ comes back full rank.
Nothing had measured it, and it is false; the derivation's Gate 2 Finding 4, rewritten under the
settled metric that day, derives the rotating null direction and shows the deficiency at
$\psi_{12} = 0$, $0.1$ and $0.3$.

**A second point before aborting (decided 2026-10-04).** The starting values are a symmetric point
(every free loading 0.7, $B = 0$, disturbance covariances 0), and a model of full rank almost
everywhere can lose rank there; Finding 4's model is deficient at every point. The derivation's
Finding 5 (`@sec-rank-second`) derives two such models: a random intercept with $\psi_{21}$ free
loses one rank at equal loadings (`tools/probe_rank_starts.do`), and two factors of two indicators
with $\psi_{21}$ free lose two at $\psi_{21} = 0$ whatever the loadings. The rule:
- Check rank(Δ) at the starting values $\theta^{(0)}$. Full rank passes, as before.
- Otherwise check at $\theta^{(2)}$: $\theta^{(2)}_i = \theta^{(0)}_i + 0.1(\{i\varphi\} - \tfrac12)$
  for every free λ, β and ψ entry, $i$ its position in θ, $\{\cdot\}$ the fractional part,
  $\varphi = (\sqrt5-1)/2$; thresholds unchanged. The shifts are distinct, at most 0.05, and
  deterministic (no random-number state); the scale and sequence were chosen without measuring.
  Perturbing loadings alone does not suffice (the two-by-two model stays deficient by two).
- Abort only if $\theta^{(2)}$ is also deficient, with the same criterion, and report the null
  directions from $\theta^{(0)}$, which name the offending parameters cleanly. A model that passes
  only at $\theta^{(2)}$ gets no note: that describes the check, not the model.
- **When $\theta^{(2)}$ is what passes, it becomes the optimizer's start** (`dwls_identify()` writes
  it into `par.value` for the free rows; no API or struct change). At $\theta^{(0)}$ the matrix
  $\Delta'W\Delta$ is singular, and a fit of the random-intercept model started there did not
  converge ($\psi_{21} \to -25.5$), while from $\theta^{(2)}$ it converged to the solution reached
  from distinct user starts (`tools/probe_rank_starts.do`, measured 2026-10-04, after the second
  point was first decided with $\theta^{(0)}$ kept as the start). A model full rank at
  $\theta^{(0)}$ keeps $\theta^{(0)}$. A model whose estimate lands near a point of lost rank is
  Module 10's to catch, at $\hat\theta$, through the bread's rank criterion.

No fixed point can be shown to avoid every model's exceptional set; $\theta^{(2)}$ is chosen
without reference to any model.

**Under multiple groups (O7, decided 2026-10-09): the rule as written, read literally.** Every
free λ, β and ψ entry is perturbed, a later group's free ψ included; ν, ω, α and thresholds are
not. $i$ is the element's θ position, so a tie (one θ element) gets one shift. Adding ν, ω and α
gains nothing on the known singular start (@sec-singular-guards Result G4 (ii)); perturbing
thresholds would change the single-group rule G0 keeps. A start deficient on I5's set (equal
sample thresholds) therefore aborts, and the pre-fit guard is what speaks to the population.
The derivation's `inv-start` chunk measured the level models' rank at $\theta^{(0)}$ only;
the chunk owed before Module 3's group code landed 2026-10-10 (`@sec-mg-start2`, `mg2-start2`):
full rank at $\theta^{(2)}$ in all 20 design × level cells, design E's four levels included; with
every binary threshold at 0, only the all-binary intercepts cell is deficient, by one, at both
points, so the abort stands there.
The rank refusal says so (user, 2026-11-01): when a free ν, ω or α sits at its reference value it lists them, says the deficiency may hold only at those values and that a different start does not identify the model if the true value is the reference, in place of the "not an artefact" sentence.

Also resolve, per Section 4.8's original scope: which loadings are free vs. fixed (fixed
cross-loadings for standard CFA — this free/fixed table format is exactly what Module 12/ESEM
reuses with cross-loadings freed), residual covariances among indicators (permitted or not).

The pipeline this module owns:

```text
user model specification            sem-style surface syntax
       ↓                            parse, validate, normalize   — Module 11
canonical model specification       Module 3's ONLY input
       ↓                            identify, check, order       — Module 3
free/fixed parameter table (struct dwlsmv_parmap)
       ↓
θ vector (stable, documented ordering — same principle as Section 4.2's index map for s)
       ↓
σ(θ)  (Section 4.6)
       ↓
Δ = ∂σ(θ)/∂θ'
```

#### The canonical model specification

**`dwls_identify()`'s `modelspec` argument is not the user's text.** Module 11 owns the surface
syntax — sem-style, so that a user's existing `sem` model statement ports across — and translates
it into the canonical form below. Module 3 parses only the canonical form.

The seam earns its keep three ways. Stage 2 is not built against a parser whose surface grammar is
still moving. Malformed *user* input is diagnosed at the ado level, where `syntax` and ordinary
string handling do the work, rather than inside Mata, where `_error()`'s 100-character limit fails
destructively. And Modules 3, 5 and 6 become testable before Module 11 exists at all, because a
canonical specification is writable by hand in a do-file.

```text
factor  <F> [<F> ...]         declare factors; declaration order fixes the row
                              and column order of B and Ψ
ind     <y> [<y> ...]         declare indicators, in the order of the observed
                              variables behind the Section 4.2 index map

lambda  <F> <y>      free [<start>] | fixed <value>
beta    <to> <from>  free [<start>] | fixed <value>
psi     <F> <F>      free [<start>] | fixed <value>

metric  marker <F> <y>  |  unitvariance <F>
```

One directive per line; `#` begins a comment; blank lines are ignored. **Order within the file is
irrelevant** — θ's ordering is fixed below as a function of the model, not of the text. (Under
multiple groups a directive's section is part of what it says; order within and across sections
stays irrelevant. **Multiple groups**, below.)

**The form is a literal serialization of the parameter table, with exactly one exception.** Every
directive but `metric` names a single matrix element and states its free/fixed status outright, so
a canonical specification and a `dwlsmv_parmap` carry the same facts in the same shape. (`tie`,
multiple groups only, names several elements that other lines already declare and adds no row: it
identifies rows with one another and expands nothing, so `metric` stays the only exception.) Four
rules follow, and each of them is a trap this section already paid for:

- **`psi` is Ψ, always.** $\Phi_\eta$ has no directive and cannot acquire one: it is derived, per
  the table above. A canonical form offering a `phi` line would re-open the exact ambiguity this
  section exists to close, one layer further down.
- **Threshold parameters are not written.** Module 3 generates one free τ row per threshold from
  the index map, which knows $K$. A `tau` directive is rejected rather than accepted and ignored,
  with or without groups. Constrained thresholds are measurement invariance and exist only as
  cross-group ties: the `tau` operand of a `tie` line (**Multiple groups**, rule G3) is the only
  place a threshold is named.
- **`metric` is the only structural directive**, and it is structural precisely because what it
  fixes depends on the model graph rather than on the line itself. `unitvariance F` fixes
  $\psi_{ff} = 1$ — the *disturbance* variance — whether or not $F$ is endogenous. Module 11 does
  **not** pre-expand it: the endogenous/exogenous reading is only visible once the whole table
  exists, and expanding it here keeps this section's hardest-won result asserted in Module 3's own
  tests instead of in the command layer, far from the derivation behind it. `marker F y` fixes
  $\lambda_{Fy} = 1$ and names the marker explicitly, so Module 11's *first-listed-indicator*
  convention stays a surface convention.
- **A `metric` directive colliding with an explicit `lambda` or `psi` line on the same element is
  an error.** Two constraints on one parameter is the overparameterization this module exists to
  catch, arriving early and cheaply.

**θ's ordering: thresholds first**, in Section 4.2's index-map order, then the structural
parameters grouped `lambda`, `beta`, `psi`, each ascending by (row, col). Putting τ first makes
Section 5's block-diagonal Δ *literal* — the leading block is an identity of order $T$, needing no
computation and no permutation to line up with $s$. The ordering is a function of the model in the
same way Section 4.2's map is a function of $K$ alone; two callers building θ from one canonical
specification cannot disagree about it. The identity leading block holds for single-group and
configural texts; with ties or ν, ω rows it carries selection and scale entries instead
(**Multiple groups**, G6).

**ψ parameter order** *(settled 2026-09-28)*. Factor variance and covariance parameters use the
**lower triangle** of the symmetric Ψ. For each factor row $r$, list the columns $c \le r$; order all
ψ parameters lexicographically by $(r, c)$, with factors in their declaration order on the `factor`
line — that is, **row-major** over the lower triangle. For factors F1, F2, F3 the order is
ψ F1 F1, ψ F2 F1, ψ F2 F2, ψ F3 F1, ψ F3 F2, ψ F3 F3. An off-diagonal covariance is therefore
represented canonically as ψ Fi Fj with $i > j$, and its parameter-table name is `psi[Fi,Fj]`.

- **Both spellings are accepted; only the lower one orders.** `psi Fj Fi` and `psi Fi Fj` denote the
  same symmetric covariance, and `dwls_identify()` stores either as the lower-triangle element. The
  upper-triangular spelling does not participate in θ ordering. Writing both spellings is *the same
  element specified twice*, an error (rc 3498).
- **Row-major here, column-major in $s$, on purpose.** Section 4.2 orders the correlations in $s$
  column-major over their lower triangle, lavaan's order; θ's ψ block is row-major. The two orders
  never meet: Δ's columns are indexed by the parameter table and its rows by the index map. Do not
  "harmonize" them.
- **Module 3 is the sole owner of θ's order.** `dwls_identify()` builds the parameter table in θ
  order; `dwls_model()`, `dwls_delta()` and everything after them index the table (`ptype`, `row`,
  `col`; under groups the key gains the group, G6) and never reconstruct the order.

Worked example — Case E, six 3-category indicators, $F_2 \leftarrow F_1$, unit-variance metric:

```text
factor  F1 F2
ind     y1 y2 y3 y4 y5 y6
lambda  F1 y1 free
lambda  F1 y2 free
lambda  F1 y3 free
lambda  F2 y4 free
lambda  F2 y5 free
lambda  F2 y6 free
beta    F2 F1 free
metric  unitvariance F1
metric  unitvariance F2
```

**The error this section is built around is expressible in the form**, which it has to be for the
rank check to have anything to catch. Adding

```text
psi     F2 F1 free
```

to the specification above frees $\psi_{12}$ alongside $\beta$, and
$\operatorname{rank}(\Delta)$ drops by exactly one. A canonical form whose defaults quietly
refused to emit that line would leave the check unable to see the case it was written for.

#### Multiple groups (0.3.0)

**Module 3's group path built 2026-10-10** (`dwls_identify()` → `dwlsmv_id_groups()`; Modules 5/6
for groups the same day): every G0–G7 rule below, with O2's β refusal standing until its derivation
lands (scope (B)). The level generator (G8) is built (2026-10-12, `dwlsmv_levelgen`, Module 12). Module 11's group surface is built (2026-10-12: the `invariance()` flow, Stage 1 per group
in the command, posting; `test_invariance_flow.do`, `test_invariance_post.do`; the pooled fit and inference since 2026-10-10). Decisions and their evidence are the
2026-10-08 rows of E-003; the identification results cited (I0–I6, the level table)
are the derivation's `@sec-invariance-id`. Points marked **open** are not settled and are collected
at the end of this subsection; the text beside each states the reading a builder would otherwise
assume, so that it can be checked, not adopted.

```text
groups  <G>                                    G >= 2; once
group   <g>                                    opens section g; g = 1..G, each once

nu      <y>   free [<start>] | fixed <value>   groups g >= 2
omega   <y>   free [<start>] | fixed <value>   groups g >= 2
alpha   <F>   free [<start>] | fixed <value>   groups g >= 2

tie     <elem> <elem> [<elem> ...]
<elem>  := lambda <g> <F> <y> | beta <g> <to> <from> | psi <g> <F> <F>
         | tau <g> <y> <k>    | nu <g> <y> | omega <g> <y> | alpha <g> <F>
```

An operand is the element's type, its group, and its names in that type's directive order.
`tau g y k` is threshold $k$ of indicator $y$ in group $g$. For `nu`, `omega` and `alpha` the
status keyword is the line's third token, not the fourth as for `lambda`, `beta` and `psi`.

**G0. Single group is unchanged.** A text with no `groups` line is single-group and parses exactly
as before: the same texts are accepted, with the same table, θ, order and starts, and every refusal
is the same rc 3498 for the same reason. `group`, `tie`, `nu`, `omega` and `alpha` are unrecognized
directives there, and `tau` stays rejected. Tests must show it three ways: every existing accept or
reject test passes unedited; for each corpus text T, a `groups 2` text whose sections are T and T
with an explicit `psi F F fixed 1` reproduces T's own τ, λ, β and ψ rows as group 1's, field for
field; and a text without `groups` carrying any of the five tokens gives `_rc == 3498` (killing a
mutant that reads a missing `groups` line as `groups 1`, or accepts `nu` without `groups`).
**One field differs by decision (2026-10-10, agreed by the user):** the threshold *start*. A
single group starts every τ at 0 (unchanged); under groups τ starts come from Stage 1 (O5, the
O14 C3 note). So the second test compares every field except τ starts, and asserts that group 1's
τ starts equal its sample thresholds.

**G1. Sections.**
- `groups G` appears once, anywhere at top level, with $G$ an integer of at least 2. `groups 1` is
  refused: it would be a second spelling of single-group.
- A section runs from a `group` line to the next `group`, `groups` or `tie` line. Every other
  directive lies in a section; `groups` and `tie` lines belong to none.
- Sections may appear in any textual order: θ depends on the model, not the text. (The generator's
  *emission* order is a different matter, normative under G8.)
- Each section is a complete single-group text under the existing rules, with its own `factor` and
  `ind` lines. Every group declares **the same factors and the same indicators, in the same order as
  group 1** (decided 2026-10-08; per-group indicator order is an `e()` reporting matter). A later
  group's factor carries an explicit `psi F F` row and no `metric` line (G5).

**G2. `tie`.**
- A tie names at least two elements, all of one type, in pairwise distinct groups, each present in
  the table and free.
- Its operands share one θ element. It adds no row, changes no row's status, and carries no start
  value of its own.
- An element lies in at most one tie, so an equality class over $G$ groups is one line.
- `psi` operands are accepted in either spelling and normalized to the lower triangle; both
  spellings of one element count as that element twice.
- Module 3 never generates ties. Levels belong to the generator (G8).
- Nothing longitudinal is added; the group token is the only place the dimension is named.
- **O3 (decided 2026-10-09): operands name the same element in each group** — identical names
  in every key position (`lambda`: F and y; `tau`: y and k), differing only in the group. A
  tie across different rows or columns is refused (G7). Longitudinal ties stay deferred, and
  would not use `tie` across independent groups.
- **O4 (decided 2026-10-09): only `tau` and `lambda` operands are admitted**, the types G8's
  levels emit. The grammar above still parses all seven element spellings, so that a known but
  unadmitted type (`beta`, `psi`, `nu`, `omega`, `alpha`) gets its own refusal, distinct from an
  unrecognized token. The whitelist is applied by Module 3 to **every** `tie` line, whatever
  produced it; it does not rely on the generator emitting only τ and λ.
- **O5 (decided 2026-10-09).** (a, amended 2026-11-01, user) A tied threshold starts at the sample threshold of its **lowest-numbered member group**; when group 1 is a member, **group 1's sample threshold**
  $s_{1\tau}$: with ν = 0, ω = 1, κ = 0 in the reference group, $s_{1\tau}$ estimates the tied τ
  itself, while a later group's $s_{g\tau}$ estimates $\omega_{gj}(\tau-\nu_{gj}-\lambda_{gj}\kappa_g)$
  (C11-1). A tie excluding group 1 is admitted: at the starts (ν = 0, ω = 1, α = 0) its lowest member's $s_{g	au}$ reads τ too, and group 1's threshold is not used, since it may belong to a different category split. Other members' sample thresholds are not used and never compared. K must match among members only (G7); a non-member's K does not enter. (b) Tied `lambda`
  operands whose start values differ are **refused** (G7), never resolved by a winner.

**G3. Thresholds.** Threshold rows are still never written: Module 3 generates exactly one free τ
row per threshold per group from that group's map. A standalone `tau` directive is refused at
rc 3498, with or without `groups`; only a `tie` operand names a threshold. Because a threshold is
always present and free, a `tau` operand can fail only by being absent ($k$ outside $1..C_{gj}$).
**No silent partial tie (decided 2026-10-08):** a `tie` of thresholds of an indicator whose
category count differs between the operands' groups is an incompatible restriction and an error
naming the indicator and the groups — at Module 11 (rc 198) first, with Module 3's `tie` rule
(rc 3498) as backstop. Since the thresholds, loadings and intercepts levels all tie τ, the chain
past configural fails for such data, and the message says so. Partial threshold invariance is not
in 0.3.0.
**Order against the binary scalar refusal, F4 (decided 2026-10-12, user, option (a); built
2026-10-12, test_invariance_flow.do check 3, swap mutant red):** in Module 11's `invariance()` flow, past configural, the equal-K check runs first, on
every τ-tied indicator (asserted, not assumed), and refuses unequal K at rc 198 naming the
indicator and the groups. Only then does the binary scalar refusal (T0, T2) run, on the now
common K: a factor is **all-binary** iff every indicator of it has K = 2, so "in every group"
and "in any group" coincide by construction and the refusal needs no group qualifier. Test
(user): an item binary in group 1 and three-category in group 2 under `invariance(scalar)` gets
the unequal-K refusal at rc 198, never the binary one; a mutant reversing the order fails it.
Declined: (b) binary check first, any group (its message would suggest levels the data cannot
fit); (c) binary check first, every group (same outcome as (a), on unvalidated K).
**D9 under groups: supports compared, option (b) (decided 2026-10-12, user; built the same day, `dwlsmv_m11_gsupport`, `dwlsmv_m11_gcheck`, `test_invariance_flow.do`).**
Module 11 recodes each indicator on each group's own observed support (D9 per group), so K is
G × p. At configural nothing is tied and differing supports are fitted as they are. Past
configural the equal-K check compares each τ-tied indicator's **set of observed original values**
across groups, not only its count: `{1,2,3}` in one group and `{1,2,4}` in another have K = 3 in
both but their second thresholds bound different categories, so a tie would be silent and wrong
(G3's reason). Such an indicator is refused at rc 198 before any fit, in the same check and
before the binary scalar refusal (F4). Its message line is in the unequal-K
wording below. Tests (user): configural is not refused for support differences
alone; a τ-tied level refuses before fitting; the refusal fires on equal counts with different
supports, not only on different K; equal supports go through the existing path unchanged. No
committed fixture has an equal-count, different-support item (`mg_E.csv`, `mg_beta.csv` checked
2026-10-12: identical supports everywhere), so the test builds its own data. Declined: (a) per
group, counts compared (silent mis-tie); (c) pooled support (refuses or falls back at
configural).
**The unequal-K wording (decided 2026-10-12, user, option (a); amended the same day for D9 (b)):** one line per offending indicator,
"y4 has 2 categories in group 1 (grp = 0) and 3 in group 2 (grp = 1).", naming every group by
position and value; then once, "Every invariance level past configural ties each indicator's
thresholds across groups, which needs the same categories in every group. No model
was fitted. invariance(configural) is not affected; partial threshold invariance is not
available in this version." The same text for every level past configural. **Equal count, different values (user, option (a)):** such an indicator gets its own line, "y4 takes the values 1 2 3 in group 1 (grp = 0) and 1 2 4 in group 2 (grp = 1).", original codes as observed, every group named; an indicator whose counts differ keeps the count line; one block may mix both. The closing sentence is shared (hence "the same categories", amended from "the same number of categories"). Declined: (b) values for every item (replaces approved wording); (c) an extra sentence after the unchanged closing.

**G4. ν, ω, α.**

| Directive | Element | Groups | Absent means | Start default |
|---|---|---|---|---|
| `nu y` | ptype 5, ν (y, 1) | g ≥ 2 | no row; ν = 0 (O1) | 0 |
| `omega y` | ptype 6, ω (y, 1) | g ≥ 2 | no row; ω = 1 (O1) | 1 |
| `alpha F` | ptype 7, α (F, 1) | g ≥ 2 | no row; α = 0 (O1) | 0 |

**O1 (decided 2026-10-09): no row at the reference value.** A ν, ω or α that stays at its
reference value is not written; the line appears only when the parameter is free or fixed at
another value. Writing every one as a literal fixed row was declined.

**Q-2 (decided 2026-10-09): ω multiplies the loading at every level.** Group $g$'s implied
statistics are $\omega_{gj}(\tau_{gjc}-\nu_{gj}-\lambda_{gj}\kappa_g)$ and
$\rho_{g,ij}=\omega_{gi}\omega_{gj}(\Lambda_g\Phi_g\Lambda_g')_{ij}$ whatever the level, so Modules
5 and 6 need no level or tie information (`@sec-inv-check11` C11-1, $m_{gj}=\omega_{gj}$). Rank and
df do not discriminate this from ω multiplying only at the loadings and intercepts levels
(`inv_model`'s reading); both pass the same checks, which establish compatibility, not
preference. **Output convention:** `e(b)` is θ̂ as stored; no derived row is posted. **Corrected
2026-10-09:** this θ is lavaan's own parameterization at every level (lavaan 0.7.2, delta:
σ at lavaan's estimates with lavaan's λ reproduces its implied moments to ≤ 5.6e-16, and with
λ/ω misses by 0.02–0.31 at the thresholds level; `@sec-mg-sigma-delta` `mg2-lavaan`), so
group-$g$ loadings are directly comparable with lavaan's and no `nlcom` recipe is needed. Only
the derivation's helper `inv_model` uses $u=\omega\lambda$ at the thresholds level. Θ_g in this
section means the standardized residual variance $\Theta^*_{g,jj}=1-\omega_{gj}^2(\Lambda_g\Phi_g\Lambda_g')_{jj}$
(C11-1); the unscaled $\Theta^\circ_g$ is named as such where used.

- Group 1 is the reference: ν = 0, ω = 1, α = 0, ψ fixed by the metric. **O9 (decided
  2026-10-09): group 1 refuses every `nu`, `omega` and `alpha` line**, whatever its status —
  free, fixed at another value, or fixed at the reference value — so the reference group has one
  canonical spelling. Under O1 the generator never writes these lines in group 1. The refusal is
  G7's row; its `errprintf()` names group 1, the directive (`nu`, `omega` or `alpha`), its
  indicator or factor and its status as written, and says the line is refused because group 1 is
  the reference, not because of its value. It is thereby distinguishable from the ω-positivity
  refusal and from any error about a value or a free parameter in groups $g \ge 2$.
- The start defaults are the values at which `inv-start` measured the level models' rank
  (2026-10-08). A tied threshold starts at its lowest member group's sample threshold, group 1's when a member (O5, amended 2026-11-01; G2), the
  start `inv-start` and `@sec-singular-guards` Result G4 assumed.
- A `fixed` value or a start for ω must be positive. **A free ω going non-positive is Module 7's
  admissibility (O12, decided 2026-10-09):** ω̂_gj ≤ 0 is a new admissibility condition, derived
  and reported as a note naming the indicator, the group and the value, like a Heywood case. The
  fit is not constrained (no bounds, no reparameterization), not aborted, and not marked
  non-converged; θ̂ is posted as is. ω's sign is not a reflection (Q-2: it enters the thresholds
  and the correlations), so nothing canonicalizes it.
- Modules 5 and 6 read an absent ν, ω or α as its reference value.
- **Absent-row default (decided 2026-10-09): per type, one rule for every consumer.** A row absent
  from a table reads as ω = 1 for ptype 6 and as 0 for every other ptype (ν and α included).
  Modules 5 and 6 and checks 8 and 11 all call the same derived rule; it is not stored. This
  applies whenever a row is absent; under O1 (no row at the reference value) that is every
  reference-valued ν, ω and α. Before it, checks 8 and 11 read every absent row as 0, which would refuse a
  pair with ω absent in $H_0$ and `omega y fixed 1` in $H_1$.

**G5. `metric` and scope under groups.**
- `metric` appears only in group 1's section, and only as `unitvariance` (unit-variance metric
  only, decided 2026-10-08). It expands exactly as today, to $\psi_{FF}$ fixed at 1 in group 1.
- A later group's factor variance is its `psi F F` row, fixed at 1 or free (G8); any other fixed
  value is refused at rc 3498 (decided 2026-11-02, user; single-group text is unchanged). A `metric`
  line in group $g \ge 2$ is refused.
- `metric marker` and structural paths (a β free or fixed nonzero, the recursion rule's notion) are
  refused with invariance levels in 0.3.0 (decided 2026-10-08). Module 11 refuses first, at rc 198,
  with Module 3 at rc 3498 as backstop. **At the surface (decided 2026-10-12, user, option (b);
  not built):** `invariance()` requires an explicit `metric(unitvariance)`; D2's default is not
  changed by another option. Absent `metric()` or `metric(marker)` with `invariance()`, rc 198:
  "invariance() requires metric(unitvariance) in this version; specify metric(unitvariance). /
  The invariance levels are defined with each factor's variance fixed at 1 in / group 1; the
  marker metric is not yet available with them. No model was fitted." A loading fixed at 1 in
  the paths (D2's user normalisation), rc 198, per factor: "Factor F has a loading fixed at 1
  (F -> y1@1). With invariance(), each factor's / scale is set by metric(unitvariance) instead;
  remove the fixed loading and / specify metric(unitvariance). No model was fitted." Declined:
  (a) `invariance()` implies unit variance (an absent option's meaning would change if the marker
  metric is admitted later); (c) marker default refused literally (names an option never typed).
  **O2 (decided 2026-10-09): Module 3 refuses both in every
  `groups` text**, configural included, whether or not it carries a `tie`, ν, ω or α. The β
  refusal is decided on the assembled table after the whole text is parsed: it fires when any
  group's table holds a structural path, whichever group, section order or line order put it
  there, never on the parser meeting a particular sequence of directives. The marker refusal is
  the `metric marker` line's (this bullet's first sentence). A hand-written marker pattern (a
  loading fixed at 1 with `psi F F free`) has no marker meaning and no detector (corrected 2026-11-01, user): in group 1 it is refused by the rule that a group-1 factor needs a metric, or as a collision if `metric unitvariance` is added (G7); in a later group, which needs no metric line, it is ordinary text and passes to the rank check. **Tests
  (user):** a β row in group 2 only; the same text with its sections in reverse textual order;
  and with the β line placed before and after the factor's `lambda` lines — each `_rc == 3498`
  with the structural-path message; a configural text without ties carrying a β is refused too.

**G6. θ order.**
- Blocks come in the order τ, λ, β, ψ, ν, ω, α (ptypes 4, 1, 3, 2, 5, 6, 7), so a single-group θ
  is unchanged.
- ν, ω and α are each ordered by (group, row) (decided 2026-10-08).
- **Within the τ, λ, β and ψ blocks, elements are ordered by (group, row, col)** (O6, decided
  2026-10-09): all of group 1's elements of the block, then group 2's, and so on, each group's
  in the single-group order. This matches ν, ω, α's (group, row) and the group-major stacked
  $s$, so the configural threshold submatrix of Δ is the identity, not a permutation. **Test
  (user):** a configural text with $G \ge 2$ and several thresholds per group (design E
  qualifies) asserts each element's actual θ position (`fidx`) against a hand-written
  expectation, and asserts that Δ restricted to the threshold rows (the selection below) and
  the threshold columns is exactly the identity; a parameter count alone does not pass. A
  mutant ordering (row, col, group) must fail.
- **Table row order (user, 2026-10-10, option (a)):** group-major within each block. Blocks run
  τ, λ, β, ψ, ν, ω, α; inside a block, rows are ordered (group, row, col), fixed rows interleaved
  in place, a tied member in its own group's position. So the free rows in table order are θ's
  order, as for a single group, and `e(b)` follows the table. Declined: (b) all of group 1's rows,
  then group 2's (it mirrors the stacked s but separates table order from θ order).
- A tie contributes one θ element, placed where its lowest-group member sits; the other members
  take the same `fidx`. `npar` counts θ elements and `nall` counts rows.
- The key `(ptype, row, col)` gains a group component: **`struct dwlsmv_parmap` field
  `real colvector grp`**, one entry per table row, the row's group 1..G, a tied member keeping its
  own group (user, 2026-10-10, option (a) of three; packing the group into `row`/`col` and a
  vector of per-group tables were declined). ptypes 5, 6, 7 are ν (y, 1), ω (y, 1), α (F, 1).
  Module 3 fills `grp` (all 1 for a single group); Modules 5/6 read it only when more than one
  map is passed, and refuse a key that is not one integer in 1..G per row with every group
  present (rc 3498). The single-group path never reads it and refuses ptypes 5–7.
- **Stacked statistics are ordered group by group** (decided 2026-10-09): $s = (s_1', \dots,
  s_G')'$, each $s_g$ in its own map's order (@sec-multigroup-scaling MG1). Hence the threshold
  invariant is stated as a selection, not a leading block: for single-group and configural texts,
  Δ restricted to the threshold rows and the threshold columns is the identity; with several
  groups those rows are $1..T_1$, $D_1+1..D_1+T_2$, …, a selection. With ties, ν or ω it carries
  selection and scale entries (the 2026-10-08 scoping of the τ̂ = s_τ invariant).
- The Δ rank check runs once, on the stacked table.

**G7. Validation.** Every row is rc 3498, prefixed `dwls_identify: `; each `_error()` text is under
100 characters and `errprintf()` carries the detail. The existing per-group rules still apply (a
factor needs loadings, a group-1 factor needs a metric, no duplicates).

| Rejection | `_error()` message |
|---|---|
| `groups` not one integer ≥ 2 | `'groups' needs one integer, 2 or more` |
| `groups` twice | `more than one groups directive` |
| `group g` outside 1..G | `group number outside 1..G` |
| the same group section twice | `a group section appears twice` |
| numbering gap or missing section | `group sections are not numbered 1..G` |
| directive outside a section | `a directive lies outside any group section` |
| factors differ from group 1 | `a group declares different factors from group 1` |
| indicators differ from group 1 (names or order) | `a group declares different indicators from group 1` |
| fewer than 2 operands | `a tie needs at least two elements` |
| unknown or incomplete operand | `unrecognized element type in a tie`, `an element in a tie is incomplete` |
| known type other than τ, λ (O4) | `a tie may join only tau or lambda elements`; `errprintf()` names the type |
| operands name different elements (O3) | `a tie joins different elements across groups`; `errprintf()` names both operands |
| different ptypes | `a tie joins elements of different types` |
| τ operands with different K (G3) | `a tie joins thresholds of indicators with different K`; `errprintf()` names the indicator and the groups |
| absent element | `a tie names an element that is not in the table` |
| fixed element | `a tie names a fixed element` |
| two operands in one group | `a tie names two elements of one group` |
| duplicate or overlapping ties | `an element appears in a tie more than once` |
| differing starts, compared only for elements the start rule requires to share a start (decided 2026-10-09; τ never, its start is its lowest member group's data; under O4 and O5 this is tied λ) | `tied elements carry different start values` |
| ν/ω/α in group 1, any status (O9, decided 2026-10-09) | `nu, omega and alpha are not allowed in group 1`; `errprintf()` names group 1, the directive, its indicator or factor, its status as written, and that group 1 is the reference |
| ω fixed or start ≤ 0 | `omega must be positive` |
| `metric` outside group 1 | `metric appears in group 1 only` |
| marker metric, any `groups` text (O2) | `marker metric is not supported with groups` |
| structural path in any group, any `groups` text (O2) | `structural paths are not supported with groups` |
| map count ≠ G, including more than one map for a text with no `groups` line (user, 2026-10-10) | `index maps do not match the number of groups` |
| later-group factor without ψ diagonal | `a group other than 1 has a factor with no psi diagonal` |

**Precedence where two rows fire (user, 2026-10-10, option (a)): the most specific message
first**, fixed and asserted (one message per case, never "either"): different types before
different elements (a λ tied to a τ); two elements of one group before an element appearing more
than once (an element named twice in one tie); `metric` in group 1 only before the missing-ψ
diagonal and before the existing metric/ψ collision (a group-2 `metric` line). Declined: (b) an
unspecified order.

Error-path tests assert `_rc == 3498`. The tie rows need a mutant each at least for the
distinct-group check, the fixed check and the K-mismatch check. The group-1 row (O9) is tested
for each of `nu`, `omega` and `alpha` in three statuses (fixed at the reference value, fixed
elsewhere, free), each `_rc == 3498`, with the `errprintf()` text asserted to name the
directive and group 1; a mutant that accepts `fixed` at the reference value must fail.
The O3 row is tested on each key position independently: `lambda` with the same factor and a
different indicator, and with a different factor and the same indicator; `tau` with the same
indicator and a different $k$, and with a different indicator and the same $k$ (a mutant
comparing only one position must fail). The O4 row is tested with at least one admitted-syntax
type (e.g. `psi`), asserted distinct from the unrecognized-token refusal, on a hand-written
`tie` line, so the whitelist is shown to act on text the generator did not produce.
O5 (user): the differing-starts row keeps its mutant (a parser that accepts tied λ with
different starts must fail); and a tied threshold including group 1 is asserted to start at group 1's sample
threshold on data whose later groups' sample thresholds differ from group 1's, so a mutant
taking group 2's, or any mean, fails.

**G8. What the generator writes (normative emission).** Levels are not a Module 3 concept: the
level generator emits a literal text plus `tie` lines. Its **emission order is normative**
(decided 2026-10-08): check 8's exemption regenerates each level's text from the posted
configural specification and the level name and compares exactly. The order is the one in the
examples below: `groups`; then each section in group order, each laid out `group`, `factor`, `ind`,
`lambda` lines, then `metric` (group 1) or `psi` (later groups), `nu`, `omega`, `alpha`; then the
`tie` lines, τ ties first, then λ. **Detail fixed by the build (2026-10-12, `dwlsmv_levelgen`;
within the order above, nothing reordered):** one `factor` and one `ind` line; `lambda` lines in
group 1's input order; group 1's `metric unitvariance F`, or a later group's `psi F F`, in factor
order; then the off-diagonal `psi` lines in group 1's input order; `nu` and `omega` in indicator
order, `alpha` in factor order; τ ties by indicator then threshold, λ ties in group 1's `lambda`
order, one line per element naming all G groups; tokens single-spaced, each line ended by a
newline, copied lines keeping their tokens as written. In each group $g \ge 2$, from
`@sec-invariance-id`'s level table:

| Level | `psi F F` | ν | ω | α | Ties |
|---|---|---|---|---|---|
| configural | fixed 1 | — | — | — | none |
| thresholds | fixed 1 | all items free | polytomous items free | — | τ |
| loadings | free iff F has a polytomous item, else fixed 1 | all free | all free | — | τ, λ |
| intercepts | free (every factor) | — | all free | free | τ, λ |

**Level names (decided 2026-10-12, user, L2 (a)):** code and `e()` use the surface names throughout:
the generator's level argument and `e(invariance)` hold `configural`, `thresholds`, `metric` or
`scalar`, spelled out in full. The table's **loadings** is `metric` and **intercepts** is `scalar`;
the derivation keeps its own terms.

**What else the generator copies (Q-5, decided 2026-10-09).** Every user-written λ and ψ row of
the configural text (cross-loadings and factor covariances included, free or fixed) is copied
into each later group with the same status, value and start; free rows become that group's own
parameters, tied only by the level's τ and λ ties. A fixed row is never tied (G2). This is
`@sec-inv-check11`'s assumption (G1), so Result C11-4 holds as stated. **Tests (user):** one
generated chain whose configural text has a fixed nonzero cross-loading, and one with a fixed
nonzero factor covariance; at every level each such row appears in every group with its status
and value unchanged, and the generated ties name exactly the eligible free τ and λ elements,
never the fixed row.

The loadings row's "free iff" applies at the loadings level only; at the intercepts level the
factor variance is free for every factor (decided 2026-10-08). **The generator lives in Module 12
(decided 2026-10-12, user, option (b); not built)**: a Mata function from (configural canonical
text, category counts K, level name) to the level's canonical text, called by Module 11's
`invariance()` and by check 8's exemption. The exemption reads K and the configural text from
stored results only: never the current dataset, never a reparse of `invariance()` (user). **K is posted as `e(K)` (decided 2026-10-12, user, option (ii); built
2026-10-12, `test_invariance_post.do` check 2):** a G × p matrix, columns in `e(indicators)` order and named by indicator (the item-order
contract), rows the groups in O13 order (sorted group values, group 1 the reference) and named by
group value; at G = 1 one row. Asserted when posted against the threshold count in `e(b)`; the
DIFFTEST refit asserts exact equality with its own Stage 1's K. Declined: (i) deriving K from
`e(b)`'s threshold names (a display convention made into a data contract); (iii) a K line in the
canonical text (K describes the data, not the parameter table).
**The generator's input contract (decided 2026-10-12, user, L1 option (b); not built):** it takes
the full multi-group configural canonical text and **validates** that every later group mirrors
group 1 (same factors, indicators and order; the same λ and ψ rows with the same status, value and
start, per Q-5), the only allowed difference being G8's: `metric unitvariance` in group 1 and
`psi F F fixed 1` in later groups. Anything else (e.g. a cross-loading freed in group 2 only) is
refused: such a chain lies outside the derivation, so nesting would not hold by construction. It
then writes every level in the normative order, **configural included**: a configural fit posts the
generator's configural output as `e(canon)`, never the text Module 11 assembled, so the exemption's
regeneration reproduces it exactly. Declined: (a) a group-1 template input; (c) no validation.
**What each fit posts for the exemption (decided 2026-10-12, user, option (a); built 2026-10-12, `test_invariance_post.do` check 1):**
`e(canon)`, the fit's own canonical text (the refit's input; the name open since 2026-10-08);
`e(canon_configural)`, the configural text the level was generated from; `e(invariance)`, the level.
For a configural fit `e(canon_configural)` equals `e(canon)`. For every higher level it is the
**original configural text, never an intermediate level**: the generator maps configural to a
level directly and is never chained level to level, so each level is reproducible from the same
base whatever order models were fitted or restored in (user). A hand-written model posts
`e(invariance)` and `e(canon_configural)` empty and never qualifies. The exemption passes a pair
only if (1) both post the same `e(canon_configural)`, (2) the same `e(K)`, (3) regenerating each
fit's level from that base reproduces its `e(canon)` exactly, and (4) check 6's refit reproduces
each fit. **Where condition (4) is enforced (decided 2026-10-13, user, option (ii)):**
`dwlsmv_dtgenerated` checks (1)–(3) itself; (4) is met by check 6, which the front end runs on
every model before checks 8 and 9 (the contract's order, B before C). The exemption does not
call `dwlsmv_dtrefit` itself. Declined: (i) the exemption takes stored and refit results and calls
`dwlsmv_dtrefit` (a second enforcement path applying check 6's tolerance at $G>1$ before step 4
measures it there). Hardening the internals against a caller that bypasses the front end is a
separate decision, not taken. Declined: (b) a hash of the base text (the exemption must regenerate, so it needs the
text); (c) nothing new, in-call chains only. Not established: that Stata 16's `estimates
store`/`restore` keep long `e()` locals intact (a milestone probe). **Probe written 2026-10-13,
`tools/probe_long_elocal.do`:** on StataMP 18, `e(canon)` (815 bytes) and `e(canon_configural)`
of an `invariance(scalar)` fit, and a synthetic 200,000-byte local with 10,000 newlines, are
identical (length, newlines, exact checksum, regeneration) after `estimates store`/`restore` and
after `estimates save`/`use`. **Stata 16.1 run 2026-10-13 (milestone machine, 4 processors):
identical to the StataMP 18 record except the environment line**
(`tools/probe_long_elocal_stata16.txt`), so this point is established on the floor for
texts up to 200,000 bytes.

**Worked examples** (design E of `@sec-invariance-id`): one factor F; y1–y3 with 4 categories,
y4–y5 binary; two groups. So $p = 5$, $p^* = 10$, $T = 11$, $D = 2(11 + 10) = 42$, $|\mathcal P| = 3$.

Configural:

```text
groups  2
group   1
factor  F
ind     y1 y2 y3 y4 y5
lambda  F y1 free
lambda  F y2 free
lambda  F y3 free
lambda  F y4 free
lambda  F y5 free
metric  unitvariance F
group   2
factor  F
ind     y1 y2 y3 y4 y5
lambda  F y1 free
lambda  F y2 free
lambda  F y3 free
lambda  F y4 free
lambda  F y5 free
psi     F F fixed 1
```

θ has 32 elements (22 τ, 10 λ; both ψ rows fixed), `nall` = 34. df = 42 − 32 = **10** =
$G(p^* - p)$, Result I4.

Intercepts level: group 1 as above; group 2 and the ties (abbreviated here; the generator writes
every line in full):

```text
group   2
factor  F
ind     y1 y2 y3 y4 y5
lambda  F y1 free                     # ... through y5
psi     F F free
omega   y1 free                       # ... through y5
alpha   F free
tie     tau 1 y1 1  tau 2 y1 1        # 11 lines: y1-y3 k = 1..3, y4 1, y5 1
tie     lambda 1 F y1  lambda 2 F y1  # 5 lines, y1..y5
```

θ: τ 1–11, λ 12–16, group-2 ψ 17, ω 18–22, α 23; 23 elements, `nall` = 40. df = 42 − 23 =
**19** = df_conf + (G − 1)(T − 2F), Result I4. The chain is df 10, 13, 15, 19 with 32, 29, 27, 23
free parameters, verified against lavaan 0.7.2 and Result I4 (design E, `inv-algebra`,
`inv-lavaan`). The `nall` counts follow O1 (no row at the reference value); θ and df do not depend on it.

**Open points (not settled; decided only by the user, recorded in the development record).**
- **O1 (decided 2026-10-09: no row at the reference value).** Absent ν, ω, α: no row at the
  reference value, or every one written as a literal row? See G4.
- **O2 (decided 2026-10-09: every `groups` text; no detector for a hand-written marker pattern; see G5).**
  Module 3's marker/structural refusal: every `groups` text, or only one with a `tie`, ν, ω or α?
  See G5. **0.3.0 scope (B), decided 2026-10-10 (user):** cross-sectional carry-forward. β is to be
  allowed in *user-written* `groups` texts carrying metric- or scalar-level ties, never in
  generated chains, and a user can obtain a generated level's canonical text to build from. **The
  refusal stays as written until its derivation lands** (identification and rank under groups with
  $B\ne0$ at the metric and scalar levels; a reference-tool check; DIFFTEST validity), per Module
  12's *Architecture principles* 3. Not established: whether metric-level ties alone identify
  cross-group β comparisons under delta with ω free. **A later group's endogenous ψ, decided
  2026-10-10 (user, task R's R-Q1, draft recommendation):** ψ is the disturbance variance; generated
  and carried-forward texts follow the B = 0 rule unchanged (ψ_gff free iff factor f has a
  polytomous item, else fixed at 1, endogenous or not); for user-written texts Module 3 imposes no
  rule and the rank check at θ⁽⁰⁾/θ⁽²⁾ is the only arbiter, the help stating the B = 0 rule. G5's later-group rule still applies: a diagonal `psi F F` is fixed at 1 or free, endogenous or not (user, 2026-11-02).
  Evidence: draft R's `mgb-rank` grid and H-B3 controls, run on lavaan 0.7.2 (decisions
  2026-10-10, E-005); the derivation has not landed. **What a tied β compares, decided 2026-10-10
  (user, R-Q2 option (d)):** O4 is not lifted for β in 0.3.0; β may be free per group (untied) in
  user-written texts, and a cross-group β comparison is a Wald test on `e(b)`, `e(V)` (groups
  independent; `e(V)` = ACOV/(N − G)). The help states draft R's Result B6: at the metric level a
  β on a path with an all-binary factor is on a convention scale, so a raw difference compares a
  standardized association. Evidence: draft R's `mgb-pop`. Consequences: a tied β cannot be
  written, so the tied-β escape at θ⁽²⁾ (R-Q4) is not reachable from a 0.3.0 text, check 8 needs
  no closure rule for β ties yet, and task S has no tied-β pair to validate. **Which texts admit β,
  decided 2026-10-10 (user, R-Q3 option (a)); scope (B) revised on evidence:** O2's refusal is
  lifted for user-written `groups` texts with untied β at any level; Module 3 checks no level
  condition and the rank check is the gate (principle 1). Ground: draft R's `mgb-rank` found
  configural and thresholds-level β tables locally identified with Result B4's df in every cell,
  so the tie condition is not an established identification requirement. The lift takes effect
  when the derivation lands (principle 3). **Principle 3's DIFFTEST validity for this lift,
  decided 2026-10-10 (user, option (b), before any task S result was analyzed):** task S's X0
  gate holds; its agreement item succeeds in every positive-control-passing design; all four
  primary weighting cells read `REPRODUCED: A`; and no size cell is proven above or below the
  band (UNDETERMINED does not block and is stated in the help). E-004
  §13. **Read 2026-10-10: not met** (agreement fails in design M at near-zero statistics, against
  lavaan values rounded to 1e-10; E-004); the refusal stays and a new
  preregistered study, its criterion fixed before it runs, decides the lift (user, option (a)). **Follow-up
  read 2026-10-10: not met** (one replication with ε_A 1.44e-8 against TOL 1e-8 at a statistic of
  1.08e-7, cell MIXED; every other condition held; E-004);
  the refusal stays. **Release scope after the follow-up, decided 2026-10-10 (user, option (e), with
  (ii)):** (1) β fitting and inference are gated by their own validation evidence, not by DIFFTEST:
  the lift of O2 for untied β takes effect once a **preregistered** re-run of `mgb-infer`'s H-I1 (T3)
  and H-I2 (`e(V)`) on fresh designs passes (the existing run is diagnostic, after H-L5's failure),
  together with the build's own Stata tests; (2) DIFFTEST for nested pairs involving β stays
  refused (check 8) until its own gate is met; (3) the follow-up remains NOT MET, with no
  reinterpretation or retrospective pass; (4) the help documents the limitation and its precise
  scope; (5) a later successful DIFFTEST validation lifts the remaining restriction independently;
  (6) a cross-group Wald comparison of β is available, and its finite-sample size is stated as not
  established. A floating-point allowance for near-zero statistics (option (d)) is pursued as a
  separate numerical-methods investigation; the release does not wait for it. **Gate (1) met 2026-10-10:**
  `mgb-infer-prereg` passed (84/84 fits usable; T3 ≤ 2.8e-13, `e(V)` ≤ 1.3e-11 relative; the
  mutant control behaved as preregistered). The lift for untied β takes effect in the commit that
  builds it, with the build's own Stata tests; β DIFFTEST stays refused. **Admitted β patterns, decided
  2026-10-10 (user, option (a)):** a β free in every group only. A β fixed at a nonzero value, or
  free in some groups and absent (= 0) in others, stays refused (Module 3, rc 3498) until a
  preregistered extension validates it; a β fixed at 0 is an absent path, as before. The β DIFFTEST
  refusal is a guard in check 8 (`dwlsmv_dtnest`, rc 3498, the Mata convention; "rc 198" above
  was the surface code and is corrected here): a pair is refused when either table has groups and a
  β row. Check 8 is G = 1 only today; the guard holds when G > 1 DIFFTEST is built. **Built 2026-10-10:**
  `dwlsmv_id_groups` step 4 (test_identify_groups n24a–l, mutants M48, M64–M66); the check 8
  guard (`test_dtnest_mgbeta.do`, mutants G1–G5); end-to-end against lavaan 0.7.2 on design M
  with β free, four levels (`99_lavaan_mg_beta.R`, `test_fit_groups_beta.do`; its check 4 carries
  an attribution rule for a and b, chosen after the first run failed at scalar b by 1.2e-6 and
  fixed before the rerun). **Widened 2026-10-11 (the user's sequence, after `mgb-fixed-prereg`
  passed: 168/168 fits, T3 ≤ 8.7e-14, `e(V)` ≤ 9.1e-12, negative control as preregistered):** a β
  fixed at the *same* nonzero value in every group is admitted as well. Still refused (rc 3498): a
  fixed value that differs across groups, a β fixed in some groups and free or absent in others, and
  a β free in some groups only. The existing result, criteria and the β DIFFTEST guard are unchanged
  (the guard already covers a fixed nonzero β). Built: step 4 (n24f, n24i, n24j, n24m–o; mutants
  M48, M67, M68); end-to-end against lavaan on design M with β fixed at 0.4
  (`99b_lavaan_mg_betafix.R`, `test_fit_groups_betafix.do`, green on its first run, a β = 0.3
  text red). **Widened again 2026-10-11 (after `mgb-absent-prereg` passed in all three scopes:
  240/240 fits from 21 datasets, T3 ≤ 2.0e-13, `e(V)` ≤ 1.0e-11; agreement on those datasets and
  specifications, not a population failure rate):** with **exactly two groups**, a β free in one
  group and absent (or fixed at 0) in the other is admitted, in either orientation and on any
  path, including g2only on one path and g1only on another. With more than two groups a β must
  still be free in every group (rc 3498). Still refused at any G: a fixed value differing across
  groups, a β fixed nonzero in some groups and absent or free in others. The β DIFFTEST guard is
  unchanged and covers these tables. Built: step 4 (n24a–e, h accepted; n24p–r refused at G = 3,
  n24s accepted; mutants M69, M70); end-to-end against lavaan on design M, g2only
  (`99c_lavaan_mg_betag2.R`, `test_fit_groups_betag2.do`, green on its first run; the text
  without β red). Help passages owed with Module 11: β DIFFTEST unavailable and why; Wald
  size not established. Conditions carried with it: (1) **help wording
  (refined by the user):** comparability of a raw β difference is a property of the path's
  factors, not of the level. A raw difference compares slopes only where *both* factors of the
  path have a group-2 scale identified by data (Result B6): a factor with a polytomous item at the
  metric level or above, or, at the scalar level, a factor off its I5 set. The help must not say
  the scalar level alone guarantees it; elsewhere (every configural model; a convention-scaled
  all-binary factor) the difference compares associations on convention scales. (2) **Not
  established**, the rank check being the only gate there: G > 2, cross-loadings, non-recursive
  B, ψ_fh free beside β, single-indicator factors. (3) **Validation scope:** every β model, at
  any level, is outside the singular-case guard's calibration and MG7's ĥ result. **R-Q4
  deferred** with the β-tie lift: a model identified only off β = 0 needs a tied β, which no
  0.3.0 text can write; untied, `mgb-start` shows the same deficiency at θ⁽⁰⁾ and θ⁽²⁾.
  **R-Q5, decided 2026-10-10 (user):** α of an endogenous factor is posted as stored (an
  intercept; Q-2's convention), κ_g = (I − B_g)⁻¹α_g is not posted; the help says α is not a
  mean and gives the `nlcom` expression for a group's factor means (H-L4: lavaan's `~1` is α,
  its `mean.lv` is κ). Not established: Stata `sem`'s labelling. **R-Q6, decided 2026-10-10
  (user, option (a)):** help text only. Near-equal binary thresholds on a factor at the scalar
  level can inflate the SEs of structural paths into or out of that factor (Result B5: the
  near-null direction loads on the outcome's β, ψ and α, which the guard's κ-based harm criterion
  does not see). A β-based post-fit diagnostic is recorded as a possible later derivation, not
  0.3.0 work; no threshold for "near" is established. Generated chains,
  the level table, the singular-case guard's calibration and check 11's C11-4 stay B = 0.
- **O3 (decided 2026-10-09: identical elements).** Tie operands across different rows or
  columns (occasion-ready), or identical (row, col)? See G2.
- **O4 (decided 2026-10-09: τ and λ only, enforced on every tie).** All seven element types in
  `tie`, or only τ and λ, the ones levels use? See G2. β stays out of `tie` in 0.3.0 (O2, R-Q2,
  2026-10-10); draft R section 6 lists what lifting it would require.
- **O5 (decided 2026-10-09: reject differing λ starts; a tied τ starts at its lowest member group's sample
  threshold).** Operands with different starts: reject, or the lowest group's wins? The start of
  a tied threshold. See G2.
- **O6 (decided 2026-10-09: (group, row, col)).** θ order across groups within the τ, λ, β, ψ
  blocks. See G6.
- **O7 (decided 2026-10-09: no; later-group ψ yes, by the literal rule).** Whether ν, ω, α (and
  later-group ψ) join the second-point perturbation. See *A second point before aborting*.
- **O9 (decided 2026-10-09: refuse, any status).** ν, ω, α in group 1: refuse, or accept
  `fixed` at the reference value? See G4 and G7.
- **O12 (decided 2026-10-09: yes, an admissibility condition and note; unconstrained).** A free
  ω going non-positive: Module 7's admissibility? See G4.
- **O13 (decided 2026-10-09).** Mapping the group variable's values to 1..G, and choosing the
  reference group (Module 11's). **Mapping:** the sorted distinct values of the group variable
  in the estimation sample map to 1..G (D9's rule for categories), independent of data order, so
  a re-sort cannot change which group is group 1; value labels are display only. **Reference:**
  always group 1, the lowest value; no option to choose another (a user recodes). **Values (amended 2026-10-10, user, option (a)):** nonnegative integers only, as sem's
  `group()`, because names carry the value (`2.grp#c.F`); a non-integer or negative value is refused
  at rc 198 and a string variable at rc 109, each with an `errprintf` hint (`egen newgrp =
  group(grp)`, `encode`) (`tools/probe_group_values.do`, StataMP 18). Was: any nonmissing numeric
  value. **Missing:** an observation with
  a missing group value is excluded listwise and is outside `e(sample)`, like any other missing
  value. Not established: whether T and df are
  invariant to which group is the reference (no option depends on it).
- **O14 (decided 2026-10-09: vectors throughout).** The per-group representation is explicit
  across the pipeline: `dwls_identify()`, `dwls_model()` and `dwls_delta()` take
  `struct dwlsmv_statmap vector maps`, and `dwls_fit()` takes `struct dwlsmv_stage1 vector`
  (the 2026-10-08 "array of `dwlsmv_stage1`" decision, reconciled with all four). A single group
  is a 1-vector; there is no second way to pass it. One map built from the stacked K is refused as
  a design: it would create cross-group ρ pairs. **Built 2026-10-10** (§5.1 changed in the same
  commit): the four signatures take vectors; `struct dwlsmv_fit`'s `stage1` field is
  `struct dwlsmv_stage1 vector` (user, 2026-10-10, option (a)). All five take groups since 2026-10-10 (Modules 3, 5, 6; the pooled
  fit and inference in commit 4, `test_fit_groups.do`); a group count that does not match the table
  or fit is refused at rc 3498 (`test_o14.do`). Internal helpers (`dwlsmv_fdwls`, `dwlsmv_minimize`, `dwlsmv_corrvec`,
  `dwlsmv_delta_cols`, `dwlsmv_fit_from`, …) keep a single map or Stage 1: they do one group's
  arithmetic, and the public functions own the loop over groups. **Amended 2026-10-10 (user):**
  the objective's evaluators (`dwlsmv_fdwls`, `dwlsmv_minimize`, `dwlsmv_fit_from`) take
  `statmap vector maps` with the stacked s and the weighting-(a) vector w = π_g/diag(Γ_g), since
  the pooled objective is one function of θ; per-group helpers (`dwlsmv_corrvec`,
  `dwlsmv_delta_cols`, …) keep a single map, and the public functions own Stage 1's per-group work.
  A scalar passes where a vector is declared, so G = 1 callers are unchanged. **The contracts between the
  building commits (user, 2026-10-10; build order O14 → Modules 5/6 → Module 3):**
  - *C1, O14:* the public API above. `dwls_model()` returns the stacked
    $\sigma = (\sigma_1', \dots, \sigma_G')'$, group-major, each $\sigma_g$ in its own map's
    order (G6); `dwls_delta()` the stacked $(\sum_g D_g) \times$ `npar` Δ.
  - *C2, Modules 5/6:* they take a parameter table (with the group key and ptypes 1–7) and the
    maps, and produce the stacked representation; they read no grammar, so they are tested
    from hand-built tables (`@sec-mg-sigma-delta`), and apply the absent-row default (G4).
    **Built 2026-10-10** (`dwlsmv_mg_sub`, `_mg_parts`, `_mg_sigma`, `_mg_dcols`,
    `_mg_check`; `test_mg_model.do`): σ_g by MS1, Δ's columns by MS2 (τ, ν, ω, α closed
    form; λ, β, ψ through `dwlsmv_delta_cols` on the group's rows, scaled by ω_jω_k), ties
    summed into one column (MS3). Not yet reachable from a `groups` text: Module 3 (C3) is next.
  - *C3, Module 3:* it builds that table from the canonical text and runs the rank check once,
    on the complete stacked Δ at θ⁽⁰⁾ (and θ⁽²⁾), never separately by group (G6).
  - **Threshold starts for the group rank check (user, 2026-10-10, option (a)).** Measured the
    same day (scratch, not preregistered; E-005): with every τ started at 0, as
    for a single group, Δ is rank-deficient at θ⁽⁰⁾ and θ⁽²⁾ in 14 of the 20 design × level
    cells of `@sec-mg-start2`; with the data-free start Φ⁻¹(k/K) it is still deficient in the
    all-binary scalar cell, which Result I5 says data usually identify. So at G > 1 `dwls_identify()` takes the Stage 1
    results as an optional third argument, `struct dwlsmv_stage1 vector st` (preferred over bare
    threshold vectors, user): an untied τ starts at its own group's sample threshold, a tied τ
    at its lowest member group's (O5), and C3's single rank check runs at that θ⁽⁰⁾. Stage 1 still sees nothing
    model-implied; Module 3 reads Stage 1's output. At G = 1 the argument is optional and not
    read: single-group starts and checks are unchanged. The signature lands in §5.1 with the
    building commit. **Approved with it (user, 2026-10-10):** at G > 1 the argument is
    required (refused at rc 3498 if absent, since a data-free θ⁽⁰⁾ can refuse identified
    models), and `st` must have one entry per map with `st[g].map.K` equal to `maps[g].K`.
    Declined: (b) the caller runs the rank check (contradicts C3); (c) a data-free start
    (falsely refuses every all-binary scalar model).
- **O15 (provisional; registered 2026-10-09, outside the dependency order).** Per-group
  admissibility: is Module 7's admissibility (derived Θ, Φ positive definite, …) assessed per
  group's Θ*_g (standardized) and Φ_g, and what is reported or done when one group fails?
  **Decided 2026-10-31 (user):** yes, per group, on group g's own scale with ω entering, by
  Module 7's routine itself on a transformed table (loadings ω_jΛ_j·, thresholds σ_gjc, MS1), so
  row 1 reads Θ*_gj = 1 − ω²_gj(ΛΦΛ′)_jj, row 4 R*_g, row 6 the order of σ_gj· (a negative ω
  reverses it), rows 2, 3 and 5 Φ_g and (I − B_g); plus O12's row 7, ω_gj ≤ 0. A failing group
  is reported as Module 7's notes under that group's heading, only where it fails; nothing is
  refused and nothing is stored (no `e()` change). D11's admissibility trigger is not built at
  G > 1. `dwlsmv_mg_admissible()`, `test_mg_admissible.do`.
- **O16 (provisional; registered 2026-10-09, outside the dependency order).** Sign
  canonicalization with `tie`: `dwlsmv_signcanon()` reverses a factor whose first free loading
  is negative; reversing it in one group alone breaks a tie on that loading. Whether
  canonicalization acts on the pooled θ, only on factors with no tied loading, or otherwise, is
  open; the result must not depend on the order of groups or ties. Not before O3/O4 and O6.
  **Decided 2026-10-10 (user, option (a)):** each factor is reflected in every group at once,
  keyed on its first free loading in θ order (group 1's when tied; θ order fixed by O6), so a
  tie is never broken and the result does not depend on group or tie order. Reflection flips λ,
  the β rows and columns, off-diagonal ψ and α of that factor; τ, ν, ω, Θ are unchanged. In a
  configural fit a later group may show a negative first loading. Tests: a metamorphic reflection
  (all groups at once leaves σ, F̂, T_raw, T3, a, b, df and `fidx` unchanged and canonicalizes to
  the same θ̂), and a negative control (one group alone on a tied table is not a symmetry). Not
  established: interaction with D11's marker transform at G > 1 (not built).
- **O17 (scope decided 2026-10-10, user: in 0.3.0; Q-A frozen 2026-10-10, user).** D7's
  `standardized` applies at $G > 1$ in 0.3.0. **Q-A, frozen: convention P.** Each group is
  standardized on its own scale with ω entering: the loading is $\omega_{gj}\lambda_{gjf}s_{gf}$,
  the Θ row is $\Theta^*_{gj} = 1 - \omega^2_{gj}V_{gj}$ (the spec's $\Theta_g$, unchanged), so
  $\lambda^{*2} + \Theta^* = 1$ for one factor without cross-loading. Evidence, recorded as prior
  information (seen before this freeze): cloud draft Q run on lavaan 0.7.2, design A–E × four
  levels; P matches lavaan's standardized λ, ψ and Θ to ≤ 3.3e-16 in every cell, L differs by up
  to 0.25 where ω is free. Draft Q **predicted the opposite** (lavaan = L, P only where ω = 1); that
  prediction failed and is kept as a failure when the derivation lands, its equivalence table
  corrected to the measurement. **Q-B, frozen 2026-10-10 (user): option (b) for all three rows.**
  Standardized τ_gjk is ω_gj·τ_gjk, ν_gj is ω_gj·ν_gj, and the ω row is the constant 1, displayed
  `(constrained)`; α_gf/s_gf is not a choice. Evidence, chosen after measuring the Q run (not
  preregistered) and seen before this freeze: (b) matches lavaan 0.7.2 `std.all` to ≤ 5.6e-16 in
  every cell and group; (a) and (c) miss by up to 0.29 (decisions row 2026-10-10, E-005). Draft Q predicted
  (a); kept as a failure. Still open, each a separate question: tied-row display
  (Q-C) and constant rows (Q-D), as presentation details; and the delta-method SEs. **SE gaps explained 2026-10-10** (diagnosis after the Q run, not
  preregistered): the run compared lavaan with the L Jacobian; under the frozen conventions the SEs
  match lavaan 0.7.2 to at most 2.8e-8 relative in every cell, and the ω row (constant 1) has SE 0
  as lavaan. No formula change is indicated; the derivation carries this when it lands. **Cross-software questions are not part of this decision:**
  Mplus STDYX at $g \ge 2$ is unread and has no tier-3 fixture; why draft Q misread lavaan is not
  established; agreement with lavaan 0.7.2 is a measured fact, not a claim about other versions.
  A derivation with executable verification lands before any code. **Landed 2026-10-10:**
  `@sec-standardize-mg` (draft Q, corrected to P and (b)); its failed preregistered assertions are
  kept as a record that asserts they reproduce; `std-mg-frozen` (written after measuring) verifies
  the frozen map against lavaan 0.7.2 `std.all`, points and SEs, designs A–E × four levels. **Q-C
  and Q-D decided 2026-11-01 (user), both (a).** Q-C: every group's row is shown with its own value
  and SE, with one table note past configural. Q-D: constant rows stay, as `(constrained)`, or Stata's
  `0 (omitted)` for a constant 0. **Built 2026-11-01** (`dwlsmv_m11_std_g`): from `e()` alone, and it
  fails closed at rc 3498 unless every row has one group and one kind, the rows are `e(K)`'s complete
  table, `e(canon)`'s ties agree with `e(b)`, `e(V)` and `e(rank)`, and `e(V)` fits `e(b)` (user). A β
  row is an internal error, since G5 refuses structural paths under `invariance()`. A factor variance
  at or below 0 is rc 3498 with its own message (SG7 assumes Φ_ff > 0).
- The `e()` naming (D5, D6, D8d) with the group in row names, tied members repeated in `e(b)`,
  lands with the building commit. **Form of the group in names, decided 2026-10-10 (user, option
  (a)): sem's own `group()` form**, with the group variable's values and the first group as base
  (`bn`): `y1:1bn.grp#c.F`, `G:2.grp#c.F`, `y1:1bn.grp` (intercept), `/var(e.y1)#1bn.grp`,
  `/mean(F)#2.grp` (`tools/probe_sem_group_names.do`, StataMP 18; Stata 16 not yet run). D6's
  runtime check against sem's stripe extends to G > 1. Module 11's parameter references and `e()`
  use the same names (user's goal: estout compatibility). Declined: (b) the group in the equation
  name. **Thresholds at G > 1, decided 2026-10-10 (user, option (a)): gsem's form** `/y1:1bn.grp#c.cut1`
  (`tools/probe_group_names_limits.do`, StataMP 18); declined (b) `/y1:cut1#2.grp`. At G = 1 every name
  stays D6's, with no group tag, as sem and gsem post without `group()`. The 32-character limit
  applies per name part (same probe), so group tags need no abbreviation rule. **ν, ω, α, decided 2026-10-10 (user, option (a)):** α `/mean(F)#2.grp` and ν `y1:2.grp` (sem's
  names, covered by D6's check); ω `/scale(y1)#2.grp` (our name, sem's variance-block pattern).
  Group 1's fixed rows (ω = 1, ν = α = 0) are posted per D5. `/scale(y1)#…` posts, resolves in
  `_b[]` and `test`, and displays as its own block (`tools/probe_scale_name.do`, StataMP 18).
  Declined (b) `/y1:2.grp#c.scale`. These rows exist only at G > 1. Group values: O13 (nonnegative integers; strings refused). **Testing (user):** a unit test asserts `e(b)`/`e(V)` names are
  unique, conformable and follow the rule; an optional `tools/` diagnostic runs `esttab` only when
  estout is installed. The suite never depends on estout.

On is the draft's question Qn, the numbering the 2026-10-08 decisions row (E-003) cites; the gaps (Q8,
Q10, Q11) are the points decided there.


#### Conditional covariates (v0.4.0; specified 2026-10-25, not built)

The derivations are `@sec-condx` (one binary $x$, single group) and `@sec-condx-norm`. Several
covariates, continuous $x$, and multiple groups are specified only once their derivations land
(E-006, 2026-10-24).

```text
covariate <x> [<x> ...]       declare observed exogenous covariates; declaration
                              order fixes covariate order in s and in θ
gamma   <F> <x>  free [<start>] | fixed <value>     x → η   (ptype 8)
kappa   <y> <x>  free [<start>] | fixed <value>     x → y*  (ptype 9)
```

- **Meaning.** $\gamma_{fk}$ is the effect of a one-unit difference in the observed $x_k$ on $\eta_f$,
  in units of $\eta_f$'s **disturbance** SD, given the covariates in the model. For a binary $x_k$ it
  is the within-group standardized treatment contrast. $\kappa_{jk}$ is a direct effect of $x_k$ on
  $y_j^*$ at fixed η, i.e. DIF in threshold location: **$\kappa = -(\text{threshold shift})$**
  (`@sec-condx`: $\kappa_3 = -0.300$ for DIF of $+0.3$). These are different parameters with
  different `ptype`s. Reusing β's ptype for γ would blur $\eta\leftarrow\eta$ with $\eta\leftarrow x$.
- **Reduced form.** $b = \Lambda(I-B)^{-1}\Gamma_x + K$ and
  $\rho = \operatorname{offdiag}\{\Lambda(I-B)^{-1}\Psi(I-B)^{-\top}\Lambda'\}$, with
  $\operatorname{Var}(y^*\mid x) = 1$. Θ stays derived. A factor that any $x$ predicts is
  endogenous, so `unitvariance` fixes its disturbance ψ, as everywhere else.
- **θ order.** Blocks τ, λ, β, γ, κ, ψ, ν, ω, α (ptypes 4, 1, 3, 8, 9, 2, 5, 6, 7). γ is ordered by
  (factor, covariate) and κ by (indicator, covariate), each in declaration order. A model without
  covariates has empty γ and κ blocks, so every existing θ is unchanged. The threshold block of Δ
  stays the identity, with or without κ (`@sec-condx` (d)), at a single group.
- **Identification.** Each factor needs at least one indicator with no κ. Otherwise γ and the κ's are
  confounded along one direction (`@sec-condx` (e)). `dwls_identify()` checks this by name and
  refuses at rc 3498, rather than leaving it to the rank check.
- **Coding.** $x$ is used as given. A covariate taking exactly two values must take **0 and 1**, and
  anything else is refused at rc 3498 with guidance. Otherwise γ's unit and the reported thresholds
  silently change: ±1 coding halves γ, and centring shifts every threshold by $-b_j\bar x$. An
  ordinal predictor enters as $K-1$ dummies, each 0/1, once several covariates are specified.
- **Missing $x$:** listwise, jointly with the indicators' rule (one complete-case $N$).
- **Stage 1's input becomes $(Y, K, X)$.** $X$ is data, and nothing model-implied enters. The §5.1
  signature change is specified with the Stage 1 build, in the same commit as its code.
- **Stage 1 covariance: the Hessian (M-estimator) form** (decided 2026-10-25, the user). The item
  and pair blocks of the two-step influence function are estimated as derivatives of the mean
  score, not as outer products (BHHH). Under a correctly specified Stage 1 the two forms estimate
  the same Γ and differ by $O(N^{-1/2})$. Under Stage 1 misspecification only the Hessian form is
  the estimator's influence function (`@sec-condx-norm` Result N1). Production uses analytic
  second derivatives, derived and verified in `@sec-condx-hess` (closed forms for $H_j$,
  $U_\rho$, $U_\theta$). An item's $P$ uses the upper-tail form $\Phi(-\ell)-\Phi(-\upsilon)$
  when $\ell>0$. **Caveat:** the closed forms inherit $P$'s accuracy. A pair's $P$ is a
  four-corner sum of $\Phi_2$ values, so it cancels in the joint tails. Whether the existing
  $\Phi_2$ routine is adequate for the range the covariate shift reaches is not yet
  established. A separate $\Phi_2$ tail-validation study will settle it. **Limitations:** the
  Hessian form makes the SEs right for the *pseudo-true* target, not the estimand robust (γ can be
  badly biased under DIF, `@sec-condx` (a)). Γ̂ differs from lavaan's and Mplus's BHHH form, by
  3.3% on the reference data, so their SEs will differ slightly. v1's `psi_rho` keeps BHHH; the
  asymmetry is known and open (E-007).
- **References.** lavaan 0.7.2's `conditional.x = TRUE` is the Stage 1 reference for $s$ (an equality
  target) and, through the oracle's BHHH route, for Γ̂. It is not the authority for the covariance
  estimator. Stata's Γ̂ is tested against a tier-1 fixture written from the R oracle's Hessian form.

### 4.9 The exact WLSMV variant being targeted

**lavaan does have a `"WLSMVS"` estimator**, which is easy to assume otherwise and wrong.
Checked directly against lavaan's source (`lav_options.R`, `lav_options_estimatorgroup()`):

```r
dwls = ,
wlsm = ,
wlsmv = ,
wlsmvs = "DWLS",
```

`wlsm`, `wlsmv`, and `wlsmvs` are all valid lavaan `estimator` names, all classified within the
`"DWLS"` estimator group (i.e. all three use DWLS point estimates; they differ in which robust
`test=` statistic they default to). Per lavaan's documented `test` options
([`lavOptions`](https://rdrr.io/cran/lavaan/man/lavOptions.html)): `"scaled.shifted"` is the
mean-*and*-variance-adjusted correction via a scale-and-shift transformation, and
`"mean.var.adjusted"` is a distinct Satterthwaite-style adjustment. The specific per-estimator
default mapping (`wlsmv` → `scaled.shifted`, `wlsmvs` → `mean.var.adjusted`) is documented lavaan
behavior but was not confirmed by direct inspection of the relevant source function in this round
of review (`lav_options_est_dwls()` or equivalent was not located) — worth a direct source check
before this specific mapping goes into code, even though it's consistent with lavaan's
documentation and with general knowledge of the package.

**v1 target: lavaan `estimator = "WLSMV"` with `test = "scaled.shifted"`** — the scaled-and-shifted
WLSMV variant selected as the Mplus compatibility target. Keep that phrasing deliberate: this is
the *target*, not a claim of demonstrated equivalence — "matching Mplus's WLSMV" shouldn't be
asserted as fact until Module 9's output is actually checked numerically against Mplus (Section 7),
at which point it can be documented as "verified numerical equivalence to Mplus under [tested
configuration]," not before. **Gate 5: discharged.** The statistic is derived and verified in the companion
derivation file. With $M = U\Gamma$ and $D$ the model degrees of freedom,

$$T_3 = aT + b, \qquad a = \sqrt{\frac{D}{\operatorname{Tr}(M^2)}}, \qquad
b = D - \sqrt{\frac{D\operatorname{Tr}(M)^2}{\operatorname{Tr}(M^2)}}$$

chosen so $E(T_3) = D$ and $\operatorname{Var}(T_3) = 2D$, with $df = D$ unchanged. Verified: the
moment equations hold to machine zero (they are solved exactly); $\operatorname{rank}(U) = D$ and
$M$ has exactly $D$ nonzero, real, non-negative eigenvalues; and in the degenerate case
$W_D = \Gamma^{-1}$, $M$ is exactly idempotent and $a \to 1$, $b \to 0$, so every correction
collapses to the identity — which is the strongest available structural check, exercising $U$, both
traces, $a$, $b$ and $\hat D$ at once against a known answer.

**Neither $U$ nor $M$ is formed in production (decided 2026-09-24).** Both traces expand into
$q \times q$ quantities built from $G = \Gamma W_D\Delta$, $B = (W_D\Delta)'G$, $C = G'W_DG$ and the
bread $E = \Delta'W_D\Delta$: $\operatorname{Tr}(M) = w'\operatorname{diag}(\Gamma) -
\operatorname{Tr}(E^{-1}B)$ and $\operatorname{Tr}(M^2) = w'(\Gamma\circ\Gamma)w -
2\operatorname{Tr}(E^{-1}C) + \operatorname{Tr}((E^{-1}B)^2)$, $w = \operatorname{diag}(W_D)$. No
$O(p^{*3})$ step remains, and it is lavaan 0.7.2's own default route, reproduced from its objects to
≤ 5.7e-15 on all five tier-2 fits (derivation `@sec-gate5-traces`, `@sec-gate5-lavaan`). $U$
remains the definition and the unit-test oracle. $b = D - a\operatorname{Tr}(M)$ is a cancellation
when the correction is small, so it is asserted with an absolute tolerance.

**$T = (N-1)\,F_{DWLS}(\hat\theta)$ — lavaan's scaling, not Mplus's $N F$.** A compatibility
convention, decided 2026-09-23: follow lavaan, and where lavaan and Mplus agree, follow both.
Here they do not. lavaan 0.7.2 reports $(N-1)F$ as its unscaled statistic (read from
`lav_model_objective()` and `lav_model_test()`); Mplus uses $N F$. The two share a limit and give
$T_3$ values at most $1/(N-1)$ apart relative — $2\times10^{-4}$ at $N = 5000$. That is inside
the $10^{-3}$ test-statistic tier Section 7.2 carried until 2026-09-24, which is why that tier never saw
it; the current tier rejects $N F$ against lavaan by 14–18× on all five fits. Module 7 stores $T$ as `T_raw`. The
tier-3 comparison against Mplus sits about $1/N$ away for this reason and no other: with our $a$
and $b$, $aNF + b$ reproduces Mplus's statistic to $9.5\times10^{-6}$ and our $a(N-1)F + b$ misses
by $1.9\times10^{-4}$. `tests/unit/test_fit.do` check 5 tells the conventions apart on $T$, and
`tests/unit/test_inference.do` check 4 on $T_3$, against lavaan and against Mplus at
$3\times10^{-5}$.

**`dwls_inference()` (Modules 8–9, decided 2026-09-24).** It fills `acov`, `shift_a`, `shift_b`,
`T3` $= a\,T_{\text{raw}} + b$ and `T3_df` $= D$ together, or raises an error and fills nothing: every
quantity is computed before any field is assigned. It calls `dwlsmv_acov()`, `dwlsmv_traces()`
(the expansion above) and `dwlsmv_scaleshift()`. It shares no intermediates with
`dwlsmv_acov()`, whose arithmetic stays as tested. That routine never forms $G$ or $B$, and
sharing would save about 1.3% of Stage 1. $\operatorname{Tr}(M)$ and $\operatorname{Tr}(M^2)$ are
not stored: they are $(D-b)/a$ and $D/a^2$. Two cases the formulas above leave undefined are
settled as follows:

- **$D = 0$, a just-identified model.** $a = \sqrt{0/0}$, so no scale or shift exists. `acov` is
  filled, $T_3 = T_{\text{raw}}$, `T3_df` $= 0$, and `shift_a` and `shift_b` are missing. This is
  lavaan 0.7.2's behaviour, measured on a three-indicator model of case D: `scaling.factor` NA,
  no shift, scaled statistic equal to the unscaled one, SEs reported. It is a complete state:
  `rows(acov) > 0` distinguishes it from an unfilled fit.
- **An unconverged fit** is processed like any other. Module 7 reports non-convergence and never
  enforces it, and `converged` remains the record. Results at a non-optimum describe the last
  iterate, and whether to show them is the caller's decision.

With $D > 0$, a non-positive $\operatorname{Tr}(M)$ or $\operatorname{Tr}(M^2)$ can arise only
from a degenerate Γ, and it aborts.

**A notational trap, recorded because it is silent.** Muthén et al.'s $W$ is the *inverse* of
this spec's $W_D$: their fit function is $(s-\sigma)'W^{-1}(s-\sigma)$ and their eq. (14) is
$U = W^{-1} - W^{-1}\Delta(\Delta'W^{-1}\Delta)^{-1}\Delta'W^{-1}$. Substituting
$\operatorname{diag}(\Gamma)^{-1}$ where the paper writes $W$ inverts the weight matrix and
produces a wrong $U$ with no error raised. The derivation write-up (a short appendix, separate
from this spec) must define, without ambiguity: the residual vector, $W_D$, $\Delta$, $\Gamma$, the
construction of the residual-weight matrix $U$, its eigenvalues, the unscaled statistic, the
scaling factor, the shift parameter, the adjusted statistic, and the adjusted degrees of freedom.

**A concrete anchor for that derivation** (this doesn't replace the derivation — it's the starting
point the appendix should verify, extend, and cite Muthén, du Toit & Spisic (1997) for): the
residual-weight matrix takes the standard GMM/minimum-distance form

$$U = W_D - W_D \Delta (\Delta' W_D \Delta)^{-1} \Delta' W_D$$

— the $W_D$-metric analogue of an annihilator/residual-maker matrix, projecting out the column
space of $\Delta$. The scaling factor $a$ and shift parameter $b$ of the scaled-and-shifted
correction are functions of $\text{tr}(U\Gamma)$ and $\text{tr}(U\Gamma U\Gamma)$ (matching the
first two moments of the naive statistic's asymptotic distribution — a weighted sum of $\chi^2_1$
variables with weights equal to $U\Gamma$'s eigenvalues — to those of a scaled-and-shifted
$\chi^2$). Worth being explicit about *why* the naive statistic isn't already chi-square
distributed: if $W_D$ equaled $\Gamma^{-1}$ exactly, $U\Gamma$ would be idempotent and its
eigenvalues would be exactly 0 or 1, giving a clean chi-square with integer df and no correction
needed — the standard **full WLS case using $W = \Gamma^{-1}$** (deliberately not calling this
"full-information WLS" — that term means something different and unrelated in the SEM literature,
contrasting limited-information two-stage estimators like this one against full-information joint
likelihood; using it here would conflate two different axes of the estimator taxonomy). Because
$W_D$ is only $\text{diag}(\Gamma)^{-1}$
rather than the full inverse, $U\Gamma$ is *not* idempotent, its eigenvalues are generally
fractional, and that's the precise reason the mean-and-variance (scaled-and-shift) correction is
needed at all. The $a$ and $b$ formulas are given above and verified in the derivation file. Two findings from
that verification that bear on implementation. **$T_1$ is not an option**: matching only the mean
leaves the variance badly wrong (153 against a target of 20 under an extreme design), and rejection
climbs to 16%. **$T_2$ rounds its estimated degrees of freedom to an integer**, and the rounding is
lossy precisely when the correction is working hardest — at extreme eigenvalue spread $\hat D$
rounds to 1 and the statistic turns conservative. $T_3$ avoids both by keeping $df = D$. It is not
exact either — matching two moments is not matching a distribution, and it rejects at about 5.6%
rather than 5.0% under the most extreme design tested. That residual is a property of the
approximation, not a bug, and the test suite should not chase it.

### 4.10 Robust standard errors

$$\text{ACOV}(\hat\theta) = (\Delta' W_D \Delta)^{-1} \, \Delta' W_D \, \Gamma \, W_D \Delta \, (\Delta' W_D \Delta)^{-1}$$

**What ACOV means.** The estimated asymptotic covariance of $\sqrt N(\hat\theta - \theta_0)$, built
from Section 4.4's Γ exactly as stored — $\hat\Gamma = N^{-1}\sum_i \hat\psi_i\hat\psi_i'$, with no
$N/(N-1)$ factor — and the expected-information bread $\Delta' W_D \Delta$ at the
sign-canonicalized $\hat\theta$. $W_D = \operatorname{diag}(\Gamma)^{-1}$ is re-formed, never
stored. `dwlsmv_acov(Delta, Gamma)` returns this quantity and nothing else.

**What e(V) is (decided 2026-09-23):**

$$\boxed{e(V) = \operatorname{ACOV}(\hat\theta)/(N-1)}$$

and `dwlsmv_ev(ACOV, N)` implements literally that transformation. It takes $N$, not a divisor, so
the convention lives in one line. It is called where `e(V)` is posted (Module 11), not by
`dwls_inference()`: the fit stores ACOV, never `e(V)`. Unit-test the step in isolation.

**$N-1$, not $N$, is a compatibility convention and not a result** — the same one as
$T = (N-1)F$ in Section 4.9, so the two inference outputs share it. lavaan 0.7.2 divides by
`ntotal - ngroups` in `lav_model_vcov()` for every estimator but normal-likelihood ML, having
built Γ N-scaled with `gamma.n.minus.one = FALSE`. The fixtures pin only the product: "Γ, then
$/(N-1)$" and "Γ·N/(N−1), then $/N$" are the same number, and which lavaan does is read from its
code. $/N$ misses lavaan's SEs by exactly $\sqrt{N/(N-1)} - 1$ — 1.0e-4 at $N = 5000$, inside
the SE tier Section 7.2 carried until 2026-09-24 (the current $10^{-5}$ tier sees it) — so the test
that pins this must require $/(N-1)$ and reject $/N$, as
`test_fit.do` check 5 does for $T$; `test_acov.do` checks 2, 4 and 6 do (check 2b shows the
fixture can tell the divisors apart). Mplus, like its $T$, divides by $N$: measured on case D's
threshold SEs, which by the derivation's block identity depend on Γ alone (`test_acov.do` check
6). As for lavaan, that pins the product of Mplus's Γ normalization and its divisor.

## 5. Module breakdown and build order

Twelve modules (Module 12 is the future ESEM extension). Each is validated against Section 4's
specification, and against lavaan/Mplus output on identical simulated data where applicable,
before moving to the next.

**The conceptual spine underneath all twelve modules, worth keeping explicit rather than letting
it stay implicit in the module list:**

```text
observed data
    ↓
Stage 1 estimator  (Module 1: thresholds; Module 2: pairwise correlations)
    ↓
ŝ, Γ̂              (Module 4 estimates Γ's asymptotic covariance)
    ↓
Stage 2 minimum-distance estimator  (Modules 3, 5–7: identify, model, differentiate, fit)
    ↓
θ̂
    ↓
robust inference  (Modules 8–9: sandwich SE, scaled-and-shifted test statistic)
```

**Modules 1–2 estimate the limited-information sample statistics — Module 1 the univariate half,
Module 2 the pairwise half; Module 4 estimates their asymptotic covariance; Modules 5–7 fit the
structural model to those statistics.** Nothing model-implied (Module 5's $\sigma(\theta)$, Module
3's identification choices) should leak into Modules 1–2's or Module 4's code — Stage 1 only ever
sees the raw data and produces $s$/Γ, and has no knowledge of what structural model will eventually
be fit to them. Keeping that boundary firm is
what makes the paired Monte Carlo validation (Section 7.1) and the API separation (Section 5.1)
actually mean something, rather than being an organizational fiction that later modules quietly
violate.

**Module 0 — scaffolding, test harness, and the index map.** Project structure, the §4.2 index
map, the reference fixtures, and a cross-language harness that diffs Stata against lavaan and
Mplus at every intermediate stage rather than at final estimates only: $s$, Γ, diag(Γ), $W_D$, Δ,
the sandwich bread and meat, ACOV(θ̂), and the test-statistic components. How the lavaan reference
is fitted, and which objects it is compared against, is specified in §7.

The harness is deliberately staged rather than built wholesale — the stages past $W_D$ diff
objects that Modules 5–9 have not produced yet. Build status and what is deliberately deferred are
in the development record, with the reasoning.

Settled here: the §4.3 scale convention and the missing-data rule (listwise).

**The canonical reference cases**, which give Module 4 a controlled staircase rather than a jump
straight to multi-factor models. Each exists because it exercises something the previous one
cannot, and each has its own acceptance test:

| Case | Design | What only it provides |
|---|---|---|
| A | two binary indicators | The *only* fully closed-form oracle. A 2×2 table is saturated, so two-stage equals joint ML exactly, ρ̂ is a direct inversion needing no optimizer, and ψ_ρ has the closed form of §4.5.3 |
| B | three binary indicators | The smallest genuinely joint Γ, and cheap enough for the exhaustive Monte Carlo of §7.3 |
| C | two 3-category ordinal indicators | Multiple thresholds per variable, a genuine polychoric, and a non-degenerate propagation term. **Over-identified**, so two-stage no longer equals joint ML and Case A's oracle does not extend — C's acceptance test is that the general ψ_ρ *reduces* to A's when both variables collapse to binary, not that it shares its form |
| D | four ordinal indicators, unequal category counts | Structure that no single pair can exhibit: several different correlations at once, Γ ≻ 0 with all of them present, blocks that must repeat wherever the design repeats, and the first case where the §4.2 index map does real work rather than being trivially satisfiable |

Only after A–D pass should Module 4 move to the one-factor → two-factor → structural progression
of §7.4.

**The boundary between Modules 1 and 2 is univariate versus bivariate.** Module 1 owns
everything computed from one variable's marginal, Module 2 everything computed from a pair.
Three things follow from that seam and are relied on elsewhere: §4.5.6 assigns the
interior-versus-boundary category distinction to "Module 1's category validation", which only
makes sense if Module 1 forms $\hat q$; §4.5.1's threshold influence functions are univariate and
need nothing from Module 2; and the two halves differ in character — one closed-form and
vectorizing, the other needing an optimizer per pair — so the validation splits there too.
(An earlier split by *preprocessing* versus *sample statistics* is recorded as superseded in
the development record.)

**Module 1 — univariate sample statistics and data validation.** Variable typing (ordinal/binary
only for v1 — flag and reject continuous indicators rather than silently mishandling them),
listwise deletion, category validation, and the univariate half of Stage 1: thresholds from the
empirical cumulative proportions, $\hat\tau_{j,c} = \Phi^{-1}(\hat q_{j,c})$ through `invnormal()`,
plus the Section 4.5.1 influence-function staircase they imply. This vectorizes cleanly and needs
no optimizer, no starting value and no iteration; a Module 1 that reaches for one is wrong rather
than merely slow.

Category validation must **distinguish empty *interior* categories from empty *boundary* ones**
and carry that distinction outward rather than collapsing it — see Section 4.5.6, where the two
modes break different modules and only one of them is loud. Classifying is Module 1's job;
*deciding* what to do about it is not. Whether a degenerate category collapses, warns or aborts,
and whether a continuity correction applies, depends on Module 2's estimator and on Module 10's Γ
rank check, so it is settled there. Module 1's own routines abort on a degenerate category, which
is this spec's v1 default (fail or warn, never silently repair) and not a resolution of that
question.

**Module 2 — pairwise sample statistics.** Mata routines to estimate pairwise polychoric
correlations via bivariate ML, following the standard two-stage approach (Muthén, 1984): Module 1's
thresholds are **held fixed** while each pairwise correlation is estimated. That the thresholds are
estimated first and treated as fixed here is exactly what creates Section 4.5.2's propagation term,
so the module boundary and the two-stage correction describe the same fact. Module 2 also returns
the information quantities Gate 3 settled (Section 4.5.2): `rho_hat`, `u[r,s]`, `I_rho_rho`,
`I_rho_a`, `I_rho_b` per pair, evaluated at $\hat\rho$ by `dwlsmv_pairinfo()`; every observed cell
is usable there, so that call cannot abort on an observed cell.

*The estimator is unchanged; the algorithm is a bracketed root-finder on the analytic score*
(derivation `@sec-rootfind`, "option B", decided 2026-09-29). $\hat\rho$ is the root of the
observed-cells score $S(\rho) = \sum_{\mathcal O} n_{rs}u_{rs}$ of `@sec-emptycells`; write
$\bar S = S/N$. It replaces the earlier `optimize()` d0 search and its Fisher-scoring polish,
whose d0 step collapsed on high-ρ pairs that sat 8.6e-10 to 5.2e-8 from the root
(E-008, 2026-09-28). The evaluator is the polish's: `dwlsmv_cellprob_raw()`,
`dwlsmv_dcellprob_rho()`, and the observed-cells mask applied **before** any division (in Mata,
$x/0$ is missing and $0 \times .$ is still missing). The correctness order that governed the old
prototype still holds in substance: the Plackett derivatives are validated against numerical
differentiation and the 2×2 closed form (Section 4.5.3) before the root-finder relies on them.

- **Usable points only.** Every point Module 2 evaluates for a decision is tested against
  $\mathcal R_P$ (Section 4.5.2); an unusable point is treated as infeasible, never used for a
  sign. $\rho = 0$ is usable whenever $N < 1.5\times10^{7}$, since $\pi_{rs}(0) = p_rq_s \ge N^{-2}$.
  If it is not, the pair has no starting point and **the fit aborts with rc 3498** (decided
  2026-09-30); no status describes it.
- **The budget.** $\delta_P = 4\delta_\Phi + 3\varepsilon = 4.42\times10^{-15}$, with
  $\delta_\Phi = 4.23\cdot2^{-52} = 9.39\times10^{-16}$ per `binormal()` corner, **chosen after
  measuring**: the largest error of StataMP 18's `binormal()` against a 160-bit reference over
  $|\rho| \le B$ and thresholds in $[-2.5, 2.5]$ (`tools/measure_binormal.do`), reproduced by
  StataMP 16.1 on 2026-09-30. Thresholds outside $[-2.5, 2.5]$ are unmeasured; *whether to extend
  the range is open.*
- **The bound.** $B = 0.999$ is a declared hard search limit on $|\rho|$ (decided 2026-09-29),
  a parameter of the code, replacing the implicit cap of $\rho = 0.999\tanh z$.
- **The walk.** From $x_0 = 0$, with $d = \operatorname{sign}\bar S(0)$, trial points
  $x' = d\min(B,\ 1 - (1 - dx)/2)$ halve the distance to $\pm1$, capped at $\pm B$. A usable $x'$
  with $d\,\bar S(x') < 0$ (strict: an exact $0$ does not close a bracket) closes the bracket;
  a usable $x'$ otherwise becomes the inner point. The walk reaches $B$ in ten trials.
- **The frontier.** If $x'$ is unusable and is not $\pm B$, $dB$ is evaluated first; if it is
  usable the search stops (*gap in the usable region*). Otherwise Module 2 bisects between the
  inner point and $x'$ on usability, evaluating $\bar S$ at each usable midpoint, until a sign
  change closes the bracket or the interval is four ulps wide.
- **The bracket** $[\ell, h]$ has evaluated, usable ends with $\bar S(\ell) > 0 > \bar S(h)$,
  a maximum's orientation, so a minimum is never returned. The walk guarantees both properties;
  where the internal root finder is handed a bracket (internal testing only, §5.1), an
  unusable end or ends not of sign $(+,-)$ is a broken invariant and **aborts with rc 3498**,
  not a status (decided 2026-09-30).
- **The step: safeguarded Fisher scoring.** From the end with smaller $|\bar S|$, the trial is
  $x^+ = x + \bar S(x)/c(x)$, $c$ the masked expected information of Section 4.5.2. It is replaced
  by $(\ell + h)/2$ when $x^+ \notin (\ell, h)$ or $|x^+ - x|$ exceeds half the previous step. The
  trial replaces the end whose computed sign it shares, so every iteration shrinks a
  sign-certified bracket and the loop terminates whatever the Fisher step does. Once
  $|\bar S(x)| \le \tau_S$, the trial is instead a probe of the far side at $x \pm \tau_\rho/2$.
  No second derivative of $\pi$ is used anywhere.
- **Both stopping conditions** (decided 2026-09-29). Stop with $\hat\rho = x$ when
  $|\bar S(x)| \le \tau_S(x)$ **and** $h - \ell \le \tau_\rho(x)$. $\tau_S$ is the score's own
  first-order evaluation error, computed per evaluation from the same cell quantities:
  $$\tau_S = \sum_{\mathcal O}\hat\pi_{rs}\,\frac{\delta_{D,rs} + |\tilde u_{rs}|\,\delta_P}{\tilde\pi_{rs}}
    + (|\mathcal O| + 2)\,\varepsilon\sum_{\mathcal O}\hat\pi_{rs}|\tilde u_{rs}|,
  \qquad \delta_{D,rs} = \varepsilon\sum_{\text{corners}}|\phi_2|\,(|z| + 7),$$
  as computed. This $\delta_D$ is **not a per-corner bound**: near $|\rho|\to1$ with
  $a\approx\pm b$, the density's exponent and $1-\rho\cdot\rho$ cancel. The corrected model
  is in derivation `@sec-deltad` (decided 2026-10-25: the stated model is corrected, and the
  code is unchanged). $\tau_S$ as computed bounds the score's density error with a
  measured margin ($\le 0.0042\,\tau_S$), not a derived one.
  and $\tau_\rho = 16\,\tau_S/c + 8\,\mathrm{ulp}(x)$. The factor 16 is derived from the premise
  $H/c \ge 1/2$, which is **measured, not proved** (0.983–1.055 on the verification's pairs); if
  it fails the loop only slows. The width condition certifies a sign change; the score condition
  rejects a sign change that is not a root. A bracket that shrinks while $|\bar S|$ stays large (a
  jump) is not tested for separately: the loop runs until the evaluation cap ends it, and the
  pair is returned as *evaluation cap* (or *unresolvable point inside the bracket*, if it lands on
  an unusable point first). There is no separate "stalled bracket" status (decided 2026-09-30).
- **The endpoint check** (decided 2026-09-29; frontier width four ulps, decided 2026-09-30).
  After the stopping rule holds, $\bar S$ is evaluated at the **last usable evaluation on each
  side**: $\pm B$ if usable, otherwise the usable end of a frontier bisection towards $\pm B$,
  stopped at four ulps. It is never evaluated at $\pm1$ or at an unusable point. The root is
  accepted only if $\bar S > \tau_S$ at the lower end and $\bar S < -\tau_S$ at the upper end. The
  check costs 54–105 evaluations per pair on the verification's designs, more than the root
  (7–17); that cost was accepted. It detects an odd number of extra sign changes on one side, not
  an even number. It has **its own cap, 54 evaluations per side, derived** from the four-ulp
  width (`@sec-rootfind`, *Uniqueness*): one evaluation at $\pm B$ plus at most 53 bisections,
  under the premise, measured not proved, that no unusable point lies in $[-\tfrac12, \tfrac12]$.
  Reaching it returns *evaluation cap* (decided 2026-09-30).
- **Uniqueness is unproved.** The log-likelihood is not concave in $\rho$ in general, and no
  counterexample to a single sign change was found (derivation Levels R3, R4). There is no global
  scan in production (decided 2026-09-29); multiple roots are tested for during development.
  An accepted root (status 1; status 6 by the endpoint check; status 7 by the endpoint-check cap) is
  a centred local maximum: τ_S gives centring, the maximum-oriented bracket the local maximum.
  Treating it as *the* $\hat\rho$ rests on the premise, measured not proved, of a single sign change
  on the usable part of $(-B, B)$. Status 5 contradicts that premise for its own pair, and a
  zero-probe status 6 has no derived maximum property; neither is covered by it, and both are
  decided with the non-root statuses (decided 2026-10-02). For Γ̂, the pairs covered by this
  contract take status 1's policy. A zero-probe status 6 is status 6 at three evaluations; the
  endpoint check evaluates at least once per side, so any other status 6 has five or more. Every
  status 7 is treated with the non-root statuses until the root finder itself reports which cap
  fired, a change decided together with what the cap status means (decided 2026-10-02 and
  2026-10-03: the evaluation count separates the two exits only under a measured premise, and a
  per-pair check would replay the root finder's opening).
- **The evaluation cap** is 40 evaluations for the walk and the refinement together, **chosen
  after measuring** (never reached). The endpoint check has its own cap (above). The walk's
  frontier bisection counts toward the 40 and can alone cost about 64, so a pair whose walk meets
  the frontier may reach the cap before any Fisher step; no verified pair does. *The value 40,
  what the status implies, and the fate of the old 25-step polish limit, are open* (settled
  together).
- **`iterations`** counts evaluations of the score and cell probabilities, not loop passes, and
  includes the endpoint check's (decided 2026-09-30), so that the caps bound all the work done.

*`status`, `converged` and the distinct statuses.* The pair's **`status`** is the stored outcome;
**`converged`** is kept for compatibility and is **set from `status`** (1 exactly when `status`
is *root*), never maintained independently (decided 2026-09-30). `converged = 1` only when both
stopping conditions and the endpoint check hold. Every other exit is `converged = 0`, with
$\hat\rho$ the last usable point, and carries one of these **distinct statuses** (decided
2026-09-29), so that the open decisions on boundary pairs and on Γ̂ can tell them apart without
refitting. Statuses name outcomes a user can act on, not every internal exit; broken invariants
abort instead.

| Status | Name | Arises when |
|---:|---|---|
| 1 | root | Stopping rule and endpoint verification both hold (`converged = 1`) |
| 2 | no sign change up to the bound | Walk reaches $\pm B$ while retaining the sign of $\bar S(0)$; $\hat\rho = \pm B$ |
| 3 | no sign change up to the usable end | Frontier bisection reaches four ulps while retaining the sign |
| 4 | gap in the usable region | A trial point is unusable but $\pm B$ is usable |
| 5 | more than one sign change | An endpoint has a sign inconsistent with the required maximum-oriented bracket and $\lvert\bar S\rvert > \tau_S$ |
| 6 | ambiguous endpoint sign | $\lvert\bar S\rvert \le \tau_S$ at an endpoint, or the probes of an exact zero at $\rho = 0$ do not certify a sign |
| 7 | evaluation cap | Root-search (40) or endpoint-check (54 per side) evaluation cap reached |
| 8 | unresolvable point inside the bracket | A Fisher-loop trial point is unusable |

Aborts (rc 3498), never statuses: $\rho = 0$ unusable; a supplied bracket with an unusable end or
not of sign $(+,-)$.

**$\bar S(0) = 0$ exactly** (decided 2026-09-30; `@sec-rootfind`, *Uniqueness*). The walk has no
direction, and the endpoint check alone would not certify a crossing, so both sides are probed at
$\pm\tau_\rho(0)/2$. A certified $(+,-)$ pair closes the bracket at $0$ and the endpoint check
runs; a certified common sign $d$ starts the walk from the probe on side $d$; anything else is
*ambiguous endpoint sign*. The probes count as evaluations.

"No sign change up to the bound" covers both a likelihood rising to $\pm1$ (including the
empty-cell boundary estimate, e.g. a 2×2 table with one empty cell) and a root between $B$ and
$\pm1$; the two cannot be told apart without evaluating beyond $B$, and nothing is evaluated beyond
$B$ to classify. A one-empty-cell pair (called *Kind E* before 2026-10-06), a binary 2×2 table with one empty cell, is identified from the observed
table where the Γ̂ policy is applied, not inferred from a 2×2 table at status 2; its lack of an
interior maximum is to be derived before use. The remaining status-2 pairs are analysed
separately and are not assumed to have a root beyond $B$ (decided 2026-10-02).

**Which pairs enter `s` and Γ̂ (v1 default, decided 2026-10-03).** Status 1 and the pairs under
the accepted-root contract above (status 6 by the endpoint check) enter. Every other pair makes
Stage 1 abort with rc 3498, naming the pair, before `s` is built: status 5 and zero-probe status
6 (whether $\hat\rho$ is the relevant maximum is unresolved); a perfect-association pair (two empty cells on a
diagonal), and a repaired one-empty-cell pair whose $|\hat\rho^+| \ge B$ (see *One-empty-cell repair* below); and
the remaining status-2 pairs and statuses 3, 4, 7 and 8 (v1 establishes
neither an interior estimator nor that truncation at $B$ preserves the DWLS inference; a
pre-registered paired experiment did not support inclusion). This differs deliberately from
lavaan 0.7.2 and Mplus 8.10, which fit such data with a warning; the user documentation says so,
with these reasons. **An opt-in to include the remaining status-2 pairs** (O1, as lavaan does at
its own 0.999 cap) is adopted in principle; its name, message and the stored marker of the pairs
it included are decided with reporting (E-009, 2026-10-03). It does not admit Kind
E or statuses 3–8, and it cannot prevent the $I_{\rho\rho}$-underflow abort.

**One-empty-cell repair (O6; decided by the user 2026-10-06, by default, not an option).** keepbound and
the repair are separate branches: keepbound handles a score with no sign change up to $B$; an empty cell is
a data-support problem with a specified repair.

| Situation | v1 behaviour |
|---|---|
| Ordinary pair, root found | the root enters |
| Status 2 other than an empty-cell 2×2 pair | aborts; enters at $\pm B$ only under `keep_bound` |
| Binary 2×2 table, **exactly one** empty cell | repaired: lavaan 0.7.2's keep-margins rule, $c = 0.5$ (`@sec-o6`); $\hat\rho^+$ enters `s` with `bound` = 0; Γ̂ by construction (a), uncentred |
| Repaired, but $|\hat\rho^+| \ge B$ | aborts (O3); nothing is truncated |
| Binary 2×2 table, two empty cells on a diagonal (perfect-association pair) | aborts; not repaired, as in lavaan |
| B-centred, statuses 3–8 | abort (O3), as above |

Repair is applied where the pair's table is in hand. $n^+_{rs} = n_{rs} + c\,\sigma_e\sigma_{rs}$
keeps both margins, so Module 1's thresholds are unchanged and no other entry of `s` moves.
$\hat\rho^+$ is the root of the score on $n^+$. It is unique, and equals the saturated inversion
$\Phi_2(\hat a,\hat b;\hat\rho^+) = (n_{11} + c\,\sigma_e)/N$. Its ψ is the ordinary one
(@sec-normaliser), with the raw weights $W = n/N$, evaluated at $\hat\rho^+$. So Γ̂ stays
$N^{-1}\sum_i\psi_i\psi_i'$ over the real observations, which is lavaan's and Mplus's
construction (E-010, 2026-10-06). $c$ is fixed at 0.5: it is not user-settable. At
$c = 1$ a count of 1 would become a new empty cell.

The basis is the convention both reference programs use, not a statistical result. A
pre-registered Monte Carlo could not separate (a) from the repaired-data construction (b), and
rejected the hybrid. (b) is not defined outside its pair's own block. Whether a pair was repaired
is derivable from $(Y, K)$, so no field is added. Larger tables with empty cells are not one-empty-cell pairs
and are not repaired: lavaan's rule there sets empty cells to 0, which changes nothing.

**The constructions compared (definitions in `@sec-o6`).** All three take $\hat\rho^+$ from the
repaired table. They differ only in the pair's ψ and in what Γ̂ sums over.
- **(a), adopted:** ψ with $W = n/N$ at $\hat\rho^+$, summed over the real observations. This is
  the production construction, and lavaan's and Mplus's. Its $\bar\psi_\rho \ne 0$.
- **(b), repaired-data:** ψ with $W = n^+/N$, weighted by the repaired counts. That treats the
  pseudo-counts as data, two of them negative. It is exactly centred, but it defines only the
  pair's $(\tau_j, \tau_k, \rho_{jk})$ block, $\nabla h'V(\hat\pi^+)\nabla h$. The entries against
  any other statistic have no pseudo-observation to come from, so (b) is not a full Γ̂. The Monte
  Carlo therefore tested this block spliced into (a)'s Γ̂ and into the hybrid's. Both versions were
  not PSD on every repaired replicate. That follows from the missing derivation and is not a
  measured property of (b).
- **Hybrid, rejected:** (b)'s ψ summed over the real observations only. It is not close to (b).
  In the Monte Carlo its SEs were smaller and T3 over-rejected.

*Reporting (decided 2026-10-03).* The opt-in is `dwls_sample_stats()`'s optional `keep_bound`
(§5.1); the pairs it includes are marked in `stage1.bound`. Its user-facing option is
`keepbound` (Module 11 D8c). Stage 1 classifies every pair, then
reports, then aborts or continues; it never stops at the first offending pair.
- **The abort message** (`errprintf()`; `_error()` stays short, rc 3498) lists every offending
  pair by its positional label (`y3~~y5`), each with ρ̂ and its status in words. A 2×2 pair
  with two empty cells on a diagonal says "2×2 table with two empty cells: no interior maximum,
  not repaired" (`@sec-kindE`; an empty row or column never reaches Module 2). A repaired pair
  with $|\hat\rho^+| \ge B$ says that its repaired correlation reaches the bound. It then
  states once that a correlation at or near the bound has several possible readings —
  near-duplicate (redundant) items, strong skew collapsing the table, too few observations — that
  the data cannot tell them apart, and that this is not a diagnosis. Only when every offending pair
  is a status-2 pair other than an empty-cell 2×2 pair does it name the opt-in. Module 11 catches rc 3498 and
  prints the mapping from positional labels to variable names; Stage 1 is not given names.
- **Under `keep_bound`**, Stage 1 prints, on every fit, the included pairs and that they enter
  `s` and Γ̂ at $\pm B$ as lavaan does, with an effect on SEs and T3 that is not established (a
  pre-registered experiment measured it as small in one design and did not support inclusion).
- **Repaired pairs (decided 2026-10-06)** are listed on every fit that has one, by positional
  label with raw table and $\hat\rho^+$. The note states:
  - the rule: $c = 0.5$ on the empty cell and its diagonal partner, $-0.5$ on the other two,
    margins kept;
  - that lavaan and Mplus apply the same repair and the same Γ̂;
  - that its effect on SEs and T3 at small $N$ is not established.

  Repair is never silent. Stage 1 prints the repaired pairs on every call, by positional label
  with ρ̂⁺ and the raw table, as informational output (`printf`, not `errprintf`; decided by
  the user 2026-10-06), so a direct Mata caller sees it. Module 11 issues the note by variable
  name in the notes block, as for the near-bound note. Membership is never stored: the three
  consumers derive it from the same Y (`dwlsmv_emptycell_pairs()`). If a repaired pair fails,
  Stage 1 aborts, so whenever a `stage1` struct exists, "has exactly one empty cell" and "was
  repaired" are the same set.
- Module 10's `dwlsmv_diagadvice()` style applies: the fact first, a reading second, no advice
  presented as a diagnosis.
- **The near-bound note (decided 2026-10-03; its threshold chosen after measuring).** Every
  pair that enters `s` with `bound` = 0 and $|\hat\rho| \ge 0.99$ is listed, by positional
  label with ρ̂, once per fit. The rule is closeness to the bound, a substantive cutoff
  ("near-duplicate"), not one tuned to the Monte Carlo; it applies to status 1 and status 6
  alike, so it is computed from `stage1.s` and `stage1.bound` and needs no stored quantity.
  The note states the fact, then the same readings as the abort message (near-duplicate items,
  strong skew, too few observations), that the data cannot tell them apart, and that it is not a
  diagnosis. **Module 11 issues it**, with the label-to-name mapping; Stage 1 prints nothing for
  it. A z-type distance $(B - |\hat\rho|)/\mathrm{SE}$ and the score at $\pm B$ are not used: the
  first states a different fact (how far $B$ is in SE units, large at small $N$), and the second
  is defined only when $\pm B$ is evaluable, which depends on the table's sparsity rather than on
  nearness. Evidence: E-009, 2026-10-03.
- Tests assert content, rc 3498 on every abort and on contract violations (`keep_bound` not 0
  or 1), `bound` all zero by default, and `bound` marking exactly the included pairs.

Sparse cells (rare categories, small $N$) will drive some pairwise correlations toward the ±1
boundary or fail to converge. **v1 default behavior: fail or warn, do not silently repair.** The
one exception is the one-empty-cell pair, repaired by default and reported on every
fit (decided 2026-10-06; *One-empty-cell repair* above). Any other continuity or ridge correction stays
unimplemented. Patching a table changes the sample statistic being estimated. An earlier version
of this paragraph said that Γ must therefore be computed from the corrected table. For one-empty-cell pairs, v1
deliberately follows lavaan and Mplus instead: it estimates $\rho$ from the repaired table and
builds Γ̂ from scores on the real observations. So ψ̄_ρ ≠ 0 for a repaired pair (§4.5.5). The
alternative, (b), is not defined beyond its own block, and could not be separated from (a)
empirically (`@sec-o6`; E-011, E-009 and E-010).

**Module 3 — identification and the parameter map.** Owns the Section 4.8 pipeline from the
*canonical* model specification onward, keeping factor-metric identification and
ordinal-response-scale identification as explicitly separate concerns: canonical specification →
free/fixed parameter table → $\theta$ vector with documented ordering.

**It does not parse the user's syntax.** Module 11 does that, and hands Module 3 the canonical form
of Section 4.8. The module was once called "model parser and identification" and the rename is the
point: the parsing left here is a tokenizer over a form this project controls, and the validation
left here is *identification* validation — that Θ is not free, that $\Phi_\eta$ was not set
directly, that each factor carries exactly one metric constraint (under groups, in group 1; a later
group's factor carries an explicit ψ row, Section 4.8 G5), and the Δ rank check. Surface
legality is Module 11's.

*Design note for a v2 extension (Module 12).* Modification indices ask what would happen to the
fit if a currently **fixed** parameter were freed, so the free/fixed table should be able to
answer "which fixed parameters are candidates, and where would each sit in $\theta$" — not just
"which are free." This costs almost nothing to allow for now and is awkward to retrofit once the
table's shape is baked into Modules 5–7. Do not implement modification indices in v1; do keep
the table capable of expressing the question.

**Module 4 — the joint asymptotic covariance matrix (Γ).** The core novel contribution and
highest-risk module in the project. Implement the observation-level influence-function
architecture from Section 4.5 — **including the two-stage threshold-propagation term**, not just
the direct pairwise score — in sub-phases 4a (binary/tetrachoric) → 4b (ordinal/polychoric).

**Validation ladder — built, and now living in the suite.** The ordering principle is worth
keeping even though the ladder is climbed: cheapest and most decisive first, and **do not skip to
the lavaan comparison because it feels authoritative** — it is the *least* diagnostic of these,
since it can agree for the wrong reason. The rungs, and where each is now asserted:

| Rung | Asserted in |
|---|---|
| The within-variable identity at machine precision, no dependencies | `test_gamma_tau.do` check 1 |
| The whole threshold block against the observed joint proportions — cross-variable entries included, which is what exercises the gather across variables | `test_gamma_tau.do` check 2 |
| Structural invariants at p > 2 (Case D): Γ ≻ 0, symmetry | `test_gamma.do` check 1 |
| Case A's closed-form ψ_ρ (§4.5.3) — sub-phase 4a's acceptance test and 4b's reduction target | `test_gamma.do` check 2 |
| Exact moments wherever two statistics share a pair | `test_gamma.do` check 3 |
| Full matrix against `lavInspect(fit, "gamma")`, split by matched-ρ̂ and end-to-end (§7.2) | `test_gamma.do` 5a/5b |
| The pointwise expansion — the only rung that separates the correct influence function from the naive one | `test_gamma_expansion.do` |

Two rungs from the original plan were not built as stated. Olsson's (1979) closed-form single-pair
variances were superseded by the exact moment identities, which are stronger: they hold in the
sample rather than asymptotically. The §7.3 Monte Carlo check against N·Var(ŝ) remains open, and
is deliberately not a substitute for the expansion test — at feasible replication counts its
standard error is larger than the effect it would need to detect.

**Module 5 — model-implied moments $\sigma(\theta)$.** Implement the Section 4.6 mapping: $\theta$
(Module 3) → $\Sigma_{y^*}$ → $R_{y^*}$ → model-implied polychoric correlations (entries of
$R_{y^*}$) + threshold parameters, matching $s$'s index map exactly.

**Module 6 — derivatives (Δ).** $\Delta = \partial \sigma(\theta) / \partial \theta'$. Start
numerically (finite differences against Module 5), then decide whether analytic derivatives are
worth the implementation cost.

*Δ is exactly block diagonal (verified).* $\partial\tau/\partial(\text{threshold params})$ is
the identity and both off-diagonal blocks are exactly zero, so **the threshold block needs no
computation** — it is known in advance from the index map. Only the correlation block has to be
differentiated. Same shape as Gate 3's lookup-table result: a large part of what looks expensive
is structural. Validate against the numerical version before trusting an analytic
one.

*Design note for a v2 extension (Module 12).* The same modification-index consideration as Module
3: Sörbom's method needs $\Delta$ columns for parameters that are fixed in the current model.
Write the derivative routine to take a parameter *set* rather than assuming "the free ones," so a
column can be requested for a fixed parameter without restructuring it. Again: not a v1 feature,
just a v1 shape.

**Module 7 — DWLS optimizer.** Minimize $F_{DWLS}(\theta) = (s - \sigma(\theta))' W_D (s -
\sigma(\theta))$, $W_D = \text{diag}(\Gamma)^{-1}$, via `moptimize()`/`optimize()`. Use the ULS
solution ($W_D = I$) as the starting value — ULS is *not* closed-form in general; it's a separate,
typically better-behaved nonlinear optimization run first to get good starting values, not a
shortcut around optimization. Validate point estimates against lavaan `estimator = "DWLS"`.
`dwls_fit()` is ULS from the table's starts, then `dwlsmv_fit_from(st, par, start)` — the
DWLS run and everything after it — from the ULS point. **Optimization recovery (user, 2026-11-01).**
Only if that DWLS run does not converge, `dwls_fit()` runs it once more from the table's starts
(Module 3's). A converged first run is untouched. Objective and tolerances do not change. A note
says a retry happened; `e(converged)` and `e(iterations)` are those of the run reported, which is
the first run if the retry also fails. No stored result is added. **Marker fits are exempt (user,
2026-11-01):** a `metric(marker)` first fit and D11's unit-variance fit use `dwlsmv_fit_noretry()`, so
D11's trigger reads the first run and D11 runs as specified; a generic retry that converged would
otherwise erase its trigger (E-012). Evidence:
E-013. One caller skips ULS: D11's recovery
calls `dwlsmv_fit_from()` directly for its final fit, from transformed starts with
thresholds at s_τ (decided 2026-10-07, the user's; D11 Procedure step 4). `dwls_fit()`'s
signature and behaviour are unchanged.

Before accepting a fitted $\hat\theta$, run an **admissibility check**: no negative residual
variances, implied covariance matrices positive definite, factor correlations within $[-1, 1]$,
$(I - B)$ nonsingular, thresholds still ordered; under multiple groups (0.3.0, §4.8 O12) also
every free scale factor $\hat\omega_{gj} > 0$. Surface these as diagnostics rather than silently
accepting an inadmissible solution.

*How it is surfaced (decided 2026-09-23).* Every condition is a function of $(\theta, $ the
parameter table, the index map$)$, so the verdict is **derived, not stored**: `dwlsmv_admissible()`
returns one flag and worst value per condition, and `dwlsmv_admissadvice()` turns them into notes
naming the indicator or factor pair and the value. `dwls_fit()` neither calls them nor prints;
Module 11 calls both on the fitted $\hat\theta$ and shows the notes with the estimates. Point
estimates are still reported — the DWLS minimum is what it is — and the note says the robust
standard errors and $T_3$ assume an interior solution. `struct dwlsmv_fit` gains no field.

*Factor signs (decided 2026-09-23).* The optimizer runs on the unrestricted parameterization —
no sign constraint and no bound, which would put an artificial boundary into a smooth problem.
After convergence, and before $\Delta$ or anything else is derived from $\hat\theta$, `dwls_fit()`
canonicalizes: each factor whose reflection is a symmetry of the parameter table — no fixed
nonzero loading, structural path or covariance touching it, so marker factors never qualify — is
reversed if its first free loading is negative. Reversal is $\Lambda \to \Lambda D$,
$\Psi \to D\Psi D$, $B \to DBD$: the factor's loadings, every $\beta$ with exactly one end at it,
and every off-diagonal $\psi$ with exactly one end at it change sign; a path or covariance joining
two reversed factors does not. $\sigma(\theta)$, $F$ and $T$ are invariant, and because the
reported $\hat\theta$ is already canonical, Modules 8–9's $\Delta$, ACOV and statistics are
consistent with it without any separate transformation. ESEM (Module 12) may want a
sum-of-loadings anchor instead, since its first loading can be a near-zero cross-loading.

*The optimizer (measured 2026-09-23).* `optimize()`'s `"nr"` over a `"d2"` evaluator whose Hessian
is the Gauss–Newton $2\Delta' W \Delta$, with the analytic gradient, exact given Δ, $-2\,r'W\Delta$: a damped
Gauss–Newton. Against `"d1"` with `"nr"` (Hessian by differencing the gradient) and with `"bfgs"`
it is 325× and 70× faster on a 40-indicator model and converges tighter; `tools/time_fit.do`.

*Post-convergence polish (adopted 2026-10-23, user; E-014; `@sec-polish`).*
- **Why.** `optimize()`'s `"nr"` extends a good Newton step by stepping forward while F does not
  increase. Near $\hat\theta$ that comparison can be decided at the ulp level, so where `optimize()` stops
  depends on the path, by up to about $10^{-7}$ in check 6's units (E-015).
- **What it does.** After a converged DWLS run, `dwlsmv_fit_from()` calls `dwlsmv_polish()`, which takes
  pure Gauss–Newton steps $\theta\leftarrow\theta+P^{-1}\Delta'Wr$, with no step extension, no halving
  and no comparison of F. Near $\hat\theta$ this is a contraction at rate $c=\rho(P^{-1}R)$.
- **Where it sits in the stopping sequence:**
  1. `optimize()` runs to its own convergence or its iteration cap;
  2. only if it converged, the polish runs;
  3. `dwlsmv_signcanon()`;
  4. F, Δ and everything downstream are evaluated at the resulting θ.

  A run stopped early is therefore either stopped inside (1), and never polished, or stopped inside (2).
- **When it stops.** After a step whose relative size $\max_j\lvert d_j\rvert/(1+\lvert\theta_j\rvert)$
  is at most $10^{-11}$, or after 20 steps. The adopted rule also stops once a step at most $10^{-11}$
  fails to halve its predecessor; with both thresholds equal, that test never fires first. These
  constants were **chosen after measuring** (polish P0, pilot files).
- **When it is abandoned.** If a step cannot be computed, the first step exceeds $10^{-6}$, or a step
  grows while above $10^{-11}$, the polish is abandoned and `optimize()`'s θ is kept unchanged.
- **What it does not change.** `converged` and `iterations` remain `optimize()`'s. ULS is not polished,
  and neither is an unconverged run.
- **Validation** (E-014). 76 units, StataMP 16.1, processor counts 1, 2 and 4:
  - the targeted $3.9\times10^{-8}$ cross-count discrepancy fell to $5.2\times10^{-12}$;
  - every unit's cross-count residue was at most $9.7\times10^{-12}$;
  - θ̂ moved at most $1.1\times10^{-7}$.
- **What that validation does not establish.** It does not show that cross-count variation is reduced in
  general: 51 of 76 units were higher with the polish. Nor does it change or justify changing
  `dwlsmv_dt_reftol` or DIFFTEST's tolerance.
- **Stata 18 validation: FAIL-S1** (E-014, 2026-10-28). The same 76 units
  under StataMP 18, native x86, counts 1, 2, 4 and 6:
  - the targeted discrepancy fell to $5.4\times10^{-12}$;
  - one unit (Dp20261911 metric) reached $1.4\times10^{-11}$, above the $10^{-11}$ criterion, where
    the unpolished fit was at $5.8\times10^{-12}$.

  Inspection (descriptive, after the result): from bitwise-identical starting points, polished fits
  differ across counts by up to $1.7\times10^{-12}$. Corrected 2026-11-02: this was first read as the
  polish's own arithmetic depending on the count. On fixed bit-exact inputs the polish's operations were
  bitwise identical across counts while Stage 1's $w$ moved (polish-step-probe-s18: 4 units, one
  trajectory, movement below its floors not excluded), and s18-p6-diag accepts that the count reaches
  $\hat\theta$ only through $w$, for the 4-against-6 pair on its cases. The route at other pairs is not
  established.
- **Status (user, 2026-10-28).** The polish stays in production with its constants, as a targeted
  correction of the path-dependent stop. **Its cross-count stability is not validated.** It is not
  re-run; any revised method needs its own frozen validation.

**Module 8 — robust standard errors.** Implement $\text{ACOV}(\hat\theta)$ per Section 4.10,
including the explicit $/(N-1)$ scaling as its own unit-tested step (`dwlsmv_acov()` and
`dwlsmv_ev()`, internal; `dwls_inference()` lands with Module 9 so it is never partially filled).

**Module 9 — the WLSMV test statistic.** Hard gate (Section 4.9): implement the scaled-and-shifted
correction only after its derivation appendix exists and has been checked algebraically against
Muthén, du Toit & Spisic (1997). Validate the test statistic and $df$ adjustment against lavaan
`test = "scaled.shifted"` output specifically (not "WLSMV" loosely), across several model sizes,
with the tight tolerances from Section 7. Built 2026-09-24 as `dwlsmv_traces()`,
`dwlsmv_scaleshift()` and the locked `dwls_inference()` (Section 4.9); `tests/unit/test_inference.do`.

**Module 10 — diagnostics.** For Γ: symmetry, positive-semidefiniteness, condition number,
minimum/maximum eigenvalue. Same for $\Delta' W_D \Delta$. Rank checks where relevant.

**Γ symmetry is a computational input contract, not only a diagnostic** (decided 2026-09-25).
`dwlsmv_traces()` and `dwlsmv_acov()` both rely on $\Gamma' = \Gamma$, and both abort with rc
3498 when $a(\Gamma) > 4D\varepsilon$. The measure is
$a(A) = \max_{ij}|A_{ij}-A_{ji}|/\sqrt{A_{ii}A_{jj}}$. Measure and budget are defined once, in
`dwlsmv_asym()`, which returns both. The guards and the diagnostic call it, and none of them
carries its own tolerance. Derivation `@sec-gamma-sym`.

**Every diagnostic row is one of three kinds, never mixed** (decided 2026-09-25):

| Kind | Γ̂ | $\Delta'W_D\Delta$ |
|---|---|---|
| **Structural deficiency**: exact, detected from its cause, no cutoff | too few patterns, $n_{\text{pat}}-1 < D$; threshold aliasing, detected pairwise as identical, complementary or constant columns of $[\mathbf 1\mid Z]$ (§4.5.6). General linear relations among those columns are a recorded v1 limit | `npar` $> D$; an identically zero Δ column |
| **Numerical rank deficiency**: against a stated numerical criterion | none yet | Module 3's criterion, singular-value ratio $< 10^{-8}$, applied to $W_D^{1/2}\Delta$ at $\hat\theta$. It is sized against numerical Δ's noise floor and is a numerical-rank test, not an algebraic one. **Settled 2026-09-26** on `@sec-bread-rank`, which shows the row weights cancel from the noise bound and measures the margins. It is defined once, in `dwlsmv_numrank()`, which `dwls_identify()` also calls. It judges only the non-zero columns, and is not evaluated when npar > D |
| **Near-singularity ratio**: reported, **no flag in v1** (decided 2026-10-05: no Γ̂ cutoff; see §4.5.6's handling paragraph) | $\sigma_{r^*}/\sigma_1$ of the scale-free $W_D^{1/2}\hat\Gamma W_D^{1/2}$, where $r^* = \min(D - d_{\text{alias}},\ n_{\text{pat}} - 1)$ (below). The raw ratio is also reported, value only, for continuity with lavaan comparisons; any future flag attaches to the scale-free row (decided 2026-09-26: the raw ratio does not separate the small-count failure from healthy P1 data, and the count row does) | $\sigma_{\min}/\sigma_{\max}$ of $W_D^{1/2}\Delta$ |

**The Γ̂ ratio is taken within the structurally supported rank** (decided 2026-09-25). Past $r^*$,
Γ̂'s singular values are zero by construction: by Theorem 1 of `@sec-gamma-rank` when
$n_{\text{pat}} - 1 < D$, and by Theorem 3 when thresholds alias. A null vector $v$ of
Γ̂$_{\tau\tau}$ gives the null vector $(v, 0)$ of Γ̂, because Γ̂ is PSD. So when an exact deficiency is
present, the ratio rows report the smallest singular value that is **not** structurally zero. The
deficiency itself belongs to the structural rows (too few patterns, threshold aliasing), and the
ratio rows never report it a second time. $d_{\text{alias}} = T + 1 - \operatorname{rank}[\mathbf 1\mid Z]$.
The reference generator computes it from the derivation's own rank. Production counts it from
`dwlsmv_thralias()`'s pairwise relations **by components, not by relations**:
- each constant column counts 1;
- the non-constant columns joined by identical or complementary relations form connected
  components, and each component of $m$ columns counts $m - 1$, since together with $\mathbf 1$
  it spans one dimension.

Three identical columns are therefore three relations but a deficiency of two. This pairwise count
is a lower bound on the true $d_{\text{alias}}$, and exact whenever every relation is pairwise, which
is the v1 limit already recorded. Without an exact deficiency, $r^* = D$ and the ratio is
$\sigma_{\min}/\sigma_{\max}$.

`dwlsmv_gammadiag(stage1, Y)` reports Γ̂'s rows in a fixed order: asymmetry, not PSD, too few
patterns, threshold aliasing, centring, raw ratio, scale-free ratio. **A flagged asymmetry stops
it** (decided 2026-09-26): rows 2–7 are then not evaluated, because `symeigenvalues()` reads one
triangle and a spectrum of an asymmetric matrix would look valid. Y must be the matrix Stage 1
received; a mismatch is an input-contract violation, rc 3498.

Two further Γ̂ rows carry budgets derived in `@sec-gamma-budgets`:

- **Not PSD:** flagged when $\lambda_{\min} < -4D\varepsilon\,\lambda_{\max}$.
- **Centring:** flagged when $z = \max_k|\bar\psi_k|/\sqrt{\hat\Gamma_{kk}} > 10^{-6}$, i.e. when
  $\bar\psi\bar\psi'$ moves scale-free Γ̂ by more than §7.2's matched-ρ̂ tier. It detects a score
  equation not solved to Module 2's rule (statuses 2, 3, 4 and 8, and status 7 unless its cap hit
  while probing a root already found), never a missing propagation term. The rule alone does not
  place a solved pair below the flag a priori: $z \le 2\tau_S/\sqrt J$ per pair, measured at least
  five decades below $10^{-6}$ on the verification's tables (`@sec-gamma-budgets`, 2026-10-02).
  **Repaired one-empty-cell pairs (decided 2026-10-06):**
  - their $\rho$ entries are excluded from $z$;
  - membership is derived from Y by `dwlsmv_emptycell_pairs()`, the same function Stage 1
    repairs by, so the two cannot disagree;
  - no row is added and the 7 × 4 shape is unchanged (decided by the user 2026-10-06): a repair
    is an estimator action, not a diagnostic condition. Stage 1 and Module 11 report the repaired
    pairs (Module 2, *Reporting*).

  Under (a), a repaired pair's own $z_\rho$ is $O(1)$ by construction: 0.45–0.96 in the Monte
  Carlo (`@sec-o6`). That is why it is excluded.

**Small counts are a separate family**, `dwlsmv_countdiag(Y, K)`, detected from category and cell
counts. It can fire while Γ̂ is well conditioned. It is never reported as a Γ-conditioning failure,
and it routes to the sparse-cell policy, not the Γ policy. It reports values only until its
threshold is derived. Its rows are the smallest **boundary**-category count (category 0 or
$K_j-1$) and the smallest **interior**-category count ($1,\dots,K_j-2$), each with its variable
and category position $1,\dots,K_j$, kept apart because they fail differently: a small boundary
category has an exact binomial reference (`@sec-gamma-smallcount`), a small interior one does not
(decided 2026-09-26). Whether pair-cell rows exist at all waits on the S partition: an empty cell
is not by itself damage, and $(Y, K)$ alone cannot tell a structurally empty corner at high
$\rho$ from a damaging one.

**The flag column has three states** (decided 2026-09-26), in the $[\text{flag}, \text{value},
i, j]$ layout of `dwlsmv_admissible()`: 0/1 is a flag evaluated against its criterion; `.` means
**not evaluated**; `.n` means **evaluated, with no flag defined in v1** (the ratio rows and every
`countdiag` row). `.n` is still missing, so a caller testing `flag == 1` never fires on it, and tests
assert `.n` exactly, since `missing()` cannot tell it from `.`.

**Admissibility is a separate channel** (decided 2026-09-26). `dwlsmv_diagadvice()` takes only
Module 10's rows. `dwlsmv_admissible()` and `dwlsmv_admissadvice()` remain Module 7's, because
they answer a different question: whether $\hat\theta$ is a valid point of the parameter space,
not whether the reported SEs and $T_3$ can be trusted. The two disagree in both directions in the
Monte Carlo reference. Module 11 may present both channels together.

**Module 11 — command wrapper.** Owns the surface syntax, its validation, and its normalization
into Section 4.8's canonical model specification. **The surface syntax is sem-style**: a user's
existing `sem` model statement should port across, which is the migration this command exists to
serve.

**Release labels are internal (decided 2026-10-05).** v1a, v1b and v2 are this specification's
phases. Messages and the help file do not use them: they say what dwlsmv does not support, and
"planned for a later version" where a phase is planned.

**The concrete grammar is deliberately not settled here**, and Modules 3/5/6 do not wait on it.
That is the whole purpose of the canonical form — Module 3 never sees the user's text, so the
surface can keep moving while Stage 2 is built. What *is* settled is the family (sem-style, not a
purpose-built option syntax) and the direction of travel (surface → canonical, never the reverse).
An earlier `factors()`/`reg()`/`metric()` option syntax is recorded as superseded in
the development record.

**D1, the parser (decided 2026-10-03): sem parses; Module 11 translates.** Module 11 does not
parse sem's grammar itself. It runs `sem <paths>, <options> noestimate` and serializes the
result — `e(b)`'s column stripes, `e(Cns)`, `e(lxvars)`, `e(lyvars)`, `e(oyvars)` — into the
canonical form above. sem is therefore the authority on its own grammar and defaults (the
first-listed indicator as anchor, exogenous factor covariances free, varlist shorthand, either
arrow direction), and a model ports exactly as sem reads it.
- **sem never sees the user's data.** It runs on a dummy dataset holding only the model's
  observed variables: fixed $N = 200$, a one-factor structure, each variable integer-coded 0–3
  by cutting its latent at $-0.6, 0, 0.6$ (from 2026-10-05; a dummy-construction fix, so that a
  `group()` variable is valid for sem and four groups of about 50 reach Module 11's refusal;
  41 levels gave sem's rc 2001, too few observations per group), and — asserted before sem runs,
  not assumed — no missing values, no zero variance, every $|r| < 0.9$, a well-conditioned
  correlation matrix. The canonical form therefore cannot depend on the data; checking the data
  is Stage 1's job. The user's data, sort order, RNG state and active estimation results are
  restored. (Deliberately bad real data — a duplicated item, a constant item, two complete
  rows, $N = 3$ — did not make `sem, noestimate` fail; the dummy's protection against
  start-value failures is argued, not demonstrated. The data-independence holds regardless.)
- **Which tokens are observed variables (decided 2026-10-05)** is the one thing Module 11 reads
  from the user's text, and it reads it only to build the dummy. The path and option text is split
  on `-> <- ( ) @ * , = .`. Tokens beginning with a capital letter are skipped (sem's latent
  convention) unless `nocapslatent` is given. Every remaining token that `unab` resolves in the
  user's data is an observed variable. The scan classifies tokens and reads no structure: paths,
  constraints and defaults come from sem's `e()` alone. A token that is not a model variable but
  resolves (an option keyword abbreviating a variable name) only adds an unused dummy column.
  If no token resolves, the command stops at rc 198 ("the model names no observed variable in the
  data").
- **Translation.** A single-coefficient `e(Cns)` row fixing a loading at 1, on a factor whose
  variance is not fixed, becomes `metric marker F y`; `var(F)` or `var(e.F)` fixed at 1 becomes
  `metric unitvariance F`; other fixed values become `fixed <v>` lines; a covariance fixed at 0
  is absent from `e(b)`, as it is absent from the canonical form. sem's own constraint and a
  user's `@1` produce the same row, so translation does not record which the user wrote; the
  canonical form is the same either way. Indicator error variances are Θ and are dropped.
- **What sem refuses** reaches the user as sem's message and rc 498, and it already covers two
  of this specification's rules: `var(F@1)` on an endogenous factor (a total variance) and any
  covariance with an error term, including `cov(F1*e.F2)` (D4). **What sem accepts and v1a
  does not** — residual covariances among indicators (v1b), multi-coefficient constraints,
  more than one group, anything else outside the canonical grammar — Module 11 detects in `e()`
  and rejects (rc 198, D8a). **More than one group** is refused (decided 2026-10-05). The
  contract is that `group()` is unsupported, and how the refusal is reached is an
  implementation detail. sem's `group()` runs on the dummy, so the refusal does not depend on
  the user's grouping. The dummy is integer-coded (below) so that sem accepts the group
  variable and Module 11's own message is what the user sees. With the earlier continuous
  dummy, sem refused first ("variable '…' is not nonnegative integer valued", rc 198), which
  described the dummy and not the user's data. `test_translate.do` 4.9 asserts Module 11's
  message. No text-level `group()` rejection is added,
  because a second, syntax-specific path could diverge from what sem recognizes.
- **Carried to other D-items.** A fixed non-unit loading (`y2@.5`) is sem's normalisation for
  its factor, so the factor gets no `metric` line and `dwls_identify()` aborts (rc 3498, "a
  factor has no metric constraint"); Module 11's message for it belongs to D3. `ind` order is
  sem's `e(oyvars)` order (D9).
- Verified on StataMP 18 and StataMP 16.1 (`tools/probe_d1_sem.do`): 15 models, dummy and real
  data agree, and the two versions' output is identical.

**D2, the metric (decided 2026-10-04).** The option is `metric(marker | unitvariance)`, default
`metric(marker)`.
- **`metric(marker)` is sem's own normalisation**, the first-listed indicator (lavaan's and
  Mplus's default too). Module 11 adds nothing: it is what D1's translation emits.
- **`metric(unitvariance)` applies only to factors the user did not normalise.** It never
  overrides a constraint the user wrote. Which factors those are comes from sem: under
  `forcenoanchor`, `e(Cns)` holds the user's constraints and no others. A factor is
  user-normalised if that pass gives it a `metric` line, a fixed loading, or a fixed variance; a
  fixed covariance is not a normalisation. For every other factor Module 11 adds `var(F@1)`,
  or `var(e.F@1)` when F is endogenous (it has a `beta` line), and runs sem again. The
  translation then emits `metric unitvariance F`, which Module 3 expands; endogeneity is read
  from the translated text, not from `e()`, which has been restored to the user's results by
  then. Injecting `var(F@1)` for an endogenous F — the total variance — is refused by sem with
  rc 498, so that error cannot pass silently.
- **The name is deliberate.** It is not sem's `standardized`, which rescales a fitted solution
  for display; it changes which parameters are estimated. A `std…` name was rejected for that
  reason. The output header states the metric ("marker (first-listed indicator)" or "unit
  variance (ψ_FF = 1; the disturbance variance for an endogenous factor)"), and the help file
  says what the option does and does not change.
- **Tested:** the two metrics are reparameterizations of one model, so T3, its df and the
  standardized solution are identical under both (the marker loading nonzero).
  `tests/unit/test_dwlsmv.do` check 3 fits both metrics and asserts this at 1e-6, a tolerance
  chosen after a smoke run.
- Verified on StataMP 18 and StataMP 16.1 (`tools/probe_d2_metric.do`, eight models): the
  `forcenoanchor` constraint counts, the injections and the canonical text are identical.

**D3, normalisations sem accepts (decided 2026-10-04).** Probed through translation and
`dwls_identify()`:
- **A loading fixed at 1 together with `var(F@1)`** translates to `lambda F y fixed 1` plus
  `metric unitvariance F` and is accepted (a second scale constraint the user wrote; sem accepts
  it too). Nothing to do.
- **A non-unit scale fix is rejected by Module 11, in v1a.** sem takes a loading fixed at a value
  other than 1 (`y2@.5`, `y1@2`), or `var(F@v)` / `var(e.F@v)` with $v \ne 1$, as that factor's
  normalisation, so it adds no anchor and the factor reaches Module 3 with no `metric` line.
  Module 11 detects this before Module 3 and says, naming the factor and the constraint, that a
  factor's scale is set by a loading fixed at 1 (the marker) or by its variance fixed at 1 (the
  disturbance variance for an endogenous factor), and points to `metric(unitvariance)`.
  `dwls_identify()`'s rc 3498 ("a factor has no metric constraint") stays the backstop.
  Extending the canonical form to non-unit scale fixes (spec and Module 3 in one commit) is
  possible later work, not v1a. Rescaling in Module 11 is rejected: it would change what is
  estimated.
- **A factor with no observed indicator is rejected by Module 11:** second-order factors are
  outside v1. Under sem's default the higher-order factor is anchored on a latent child
  (`beta F1 G fixed 1`) and has no `metric` line; with `var(G@1)` Module 3 aborts with "a factor
  has no indicators". Supporting them is a scope extension, not a normalisation, and a separate
  later decision.
- **Module 11's rc for these messages** is 198 (D8a).
- **Help file: a factor that sem leaves unanchored.** In
  `(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (G -> y1@1 … y6@1)` (a random intercept G), sem anchors
  neither F1 nor F2: `e(Cns)` holds only G's six constraints. Translation preserves that and
  Module 3 rejects F1. The user anchors F1 and F2, or uses `metric(unitvariance)`, which leaves
  G alone (user-normalised) and fixes F1's and F2's variances. G orthogonal to F1 and F2 is
  written by the user (`cov(G*F1@0) cov(G*F2@0)`); sem's default frees those covariances.
  With unit variance and ψ_F2F1 free this model loses rank at the default starts only; §4.8's
  second point (decided 2026-10-04) accepts it and starts the fit there.

**D4, `cov(F1*e.F2)` (decided 2026-10-04): sem's refusal stands.** sem refuses a covariance
between an exogenous variable and an error term at parse time, rc 498 ("errors (e.varname) may
not be correlated with other variables"), in either spelling (`cov(F1*e.F2)`, `cov(e.F2*F1)`),
`noestimate` or not. Under D1 the statement never reaches Module 11, and Module 11 adds nothing:
it neither pre-parses the text to pass the statement through (which would be parsing sem's
grammar beyond D1's token scan, against D1, and would accept a model sem refuses) nor appends a note to sem's message
(Module 11 cannot tell this rc 498 from sem's others without matching sem's wording).
- **The help file's wording.** The corresponding factor covariance alongside the structural
  relation between the same two factors gives the rank-deficient specification §4.8's rank
  check exists for; sem itself enforces the restriction at its syntax and model-validation
  stage. The help does not present rank deficiency as the reason for sem's refusal: sem's rule
  is general (no error term may covary with an exogenous variable), and its reasons are sem's.
- **Unaffected.** The deficient specification stays expressible in the canonical form
  (`psi F2 F1 free` with `beta F2 F1`), so Module 3's rank-check tests do not depend on Module 11.
  A disturbance covariance between two endogenous factors, `cov(e.F2*e.F3)`, is accepted by sem
  and translates to `psi F3 F2 free` (probed: `dwls_identify()` accepts it, 30 free parameters
  against 29 without).

**D5, `e(b)`'s layout (decided 2026-10-04): every parameter-table row, sem's convention.**
- **`e(b)`** holds every row of the parameter table, free and fixed, in Module 3's θ order
  (thresholds, then λ, β, ψ; Module 11 does not reorder). A fixed row carries its fixed value,
  so `_b[y1:F1]` exists and equals 1 under a marker, as after sem. Elements that are not table
  rows — covariances fixed at 0, which the canonical form omits — are not posted, as sem does
  not post them. Θ is derived and is not in `e(b)`; whether and where it is reported is D7/D8.
- **`e(Cns)`** has one row per fixed parameter, selecting that column alone at its value.
- **`e(V)`** is `dwlsmv_ev(ACOV, N)` embedded with **exact zeros** in the fixed rows and columns.
- **`e(rank)`** is the number of free parameters, posted by Module 11 (`ereturn post` does not
  set it). Its relation to `e(b)` is fixed by the contract: **columns of `e(b)` = free + fixed;
  `e(rank)` = free only.**
- **Probed on StataMP 18 and 16.1, identical** (`tools/probe_d5_post.do`): the constrained post
  is accepted, the display shows `1 (constrained)`, `_b[]` reaches fixed and free rows, `test`
  and `lincom` work, and `test` on a fixed coefficient reports the constraint dropped; posting
  free rows only makes `_b[y1:F1]` fail (rc 111). **A V that is not positive definite is not
  refused:** `ereturn post` prints "variance matrix is nonsymmetric or highly singular", returns
  rc 0, and the display then shows every SE, z and interval as missing. It does not say why, so
  dwlsmv's own checks remain the only diagnosis. A 46-character stripe posts and displays.
- **The membership/rank contract is tested** (`tests/unit/test_dwlsmv.do` check 2; the user's), on a model whose fixed
  and free counts differ (a marker, at least one other fixed row): `colsof(e(b))` = free + fixed;
  `e(rank)` = free; `e(Cns)` has exactly one row per fixed parameter, each selecting one fixed
  column at its value and no free one; `e(V)` is exactly zero in every fixed row and column and
  only there. It fails under a mutation that sets `e(rank)` to `colsof(e(b))`.
- Whether a fit reaches posting depends on the open rank-deficient-Γ question; the layout does
  not.
- **Multi-group posting convention, amending D5 for G > 1 (decided 2026-10-12, user, option
  (a); built the same day, `dwlsmv_m11_post_g`, `test_invariance_post.do`).** For $G > 1$, `e(b)` contains the full group-expanded parameter table at
  every invariance level. Parameters not free at the selected level are posted at their fixed
  values, including the absent-default values $\nu = 0$, $\omega = 1$ and $\alpha = 0$ (G4).
  Parameters fixed by the group-1 identification convention are likewise posted. The canonical
  table (`e(canon)`, Module 3's parameter map) continues to distinguish actual table rows from
  these additional posted fixed values: a posted coefficient is not thereby a row of the table,
  and `e(rank)` counts free parameters only. Consequence: `e(b)` has the same columns at every
  level, so a configural-to-scalar `esttab` lines up row for row. **Not yet settled (user):** the
  sem stripe probe establishes D6's compatibility with sem under `ginvariant(none)` only, not
  with sem under its own invariance classes; before D6 is treated as settled at G > 1, that
  distinction is to be probed and the added fixed values' names and ordering confirmed.
  Declined: (b) table rows plus group 1's fixed rows wherever a later group has them (columns
  vary by level); (c) table rows only (`_b[]` fails for group 1's rows). **Order, S4-2 (decided
  2026-10-12, user, option (a)):** D5's θ order extended by complete blocks: τ, λ, β, ψ as the
  parameter map orders them, then ν (every group, indicator order), ω (every group, indicator
  order), α (every group, factor order), group 1's fixed values at the head of each block; one
  layout at every level. Declined: (b) interleaved; (c) sem's indicator-equation order. **Probe
  (2026-10-12, `tools/probe_sem_ginvariant_stripe.do`, StataMP 18):** sem's classes `mcoef`,
  `mcons`, `mcoef mcons`, `mcoef mcons merrvar` post the same 34 names in the same order as
  `ginvariant(none)`, expressing equality as `e(Cns)` rows; `all` omits `/mean(F)#g` (32); sem
  posts group 1's mean as fixed 0, our α convention. D6's membership check is unaffected (we run
  sem with `ginvariant(none)` only); Stata 16 not yet run.
- **Notes and Module 10 at G > 1, S4-3 (decided 2026-10-12, user, option (c), narrowed; built the same day, `dwlsmv_m11_gnotes`, `test_invariance_post.do` checks 7, 8): groupwise
  display only for notes whose groupwise meaning is already established; anything whose meaning
  or trigger depends on the pooled multi-group calculation is deferred.** Groupwise, printed per
  group under "Group g (grp = v):", nothing stored: the three Stage 1 notes, each a property of one
  group's pairwise estimation exactly as at G = 1: pairs included at the bound under `keepbound`,
  the one-empty-cell repair, polychoric correlations at or beyond 0.99. Deferred (not printed,
  not stored at G > 1 in 0.3.0, each its own later point): the admissibility verdict and advice
  (read from the pooled θ̂; O12's ω ≤ 0 note belongs here), Module 10's diagnostic notes (their
  readings involve the pooled Δ and the fit's SEs), `e(bound)`'s stored form, `e(gamma_ratio)`,
  `e(gamma_ratio_raw)` and `e(gamma_rank)` (one Γ̂ at G = 1; per group or pooled undecided).
  `e(bread_ratio)` is posted (O-T10). `standardized` with `invariance()` was refused at rc 198 until
  2026-11-01, when `dwlsmv_m11_std_g` was built (O17). Declined: (a) no notes at G > 1;
  (b) deciding every shape before step 4.

**D6, parameter names (decided 2026-10-04): built from Module 3's parameter map, checked
against sem.** Module 11 builds every `e(b)` name from `dwls_identify()`'s parameter map (`ptype`,
row, column), not from the user's text and not from sem's stripe:

| Parameter | `e(b)` name | Source in the map |
|---|---|---|
| λ | `y:F` | `ptype` 1 |
| β | `To:From` | `ptype` 3 |
| ψ, variance | `/var(F)` exogenous; `/var(e.F)` endogenous | `ptype` 2, diagonal |
| ψ, covariance | `/cov(Fc,Fr)`, column factor first, `e.` on each endogenous factor | `ptype` 2, lower triangle |
| τ | `/y:cut1` … `/y:cut(K−1)` (gsem's convention for an ordinal outcome) | `ptype` 4 |

- **Endogenous** means the factor has a β row in the parameter map — Module 3's reading, as
  §4.8 requires — so a ψ name always says ψ: an endogenous factor's variance, free or fixed, is
  `/var(e.F)`, never `/var(F)`. **One addition (decided 2026-10-23):** a factor whose only paths
  are fixed at 0 has no β row (a fixed zero is no path, and sem posts no column for it), yet sem
  still reads it as endogenous and posts `/var(e.F)`. Module 11 therefore also treats as endogenous
  any factor for which sem's stripe holds `/var(e.F)`, and D2's injection treats as endogenous every
  factor in sem's `e(lyvars)`. The two readings differ only in naming: with no live path Φ = Ψ for
  that factor, so the model, Δ and T3 equal the `cov(F1*F2@0)` fit's.
- **Runtime check.** sem's `e(b)` stripe is captured during translation (D1). Module 11 compares
  every structural name it built (λ, β, ψ) with that stripe, less `_cons` and the indicator error
  variances, in both directions; any mismatch aborts as an internal error, since it means Module 3
  and sem read the model differently. Thresholds have no sem counterpart and are not checked.
- The internal `par.name` strings (`psi[F2,F1]`, `y1|t1`) are not posted.
- **Probed on StataMP 18** (`tools/probe_d6_names.do`): on seven models — marker and unit
  variance with an endogenous factor, exogenous covariances, a cross-loading, two endogenous
  factors with correlated disturbances, a factor with two predictors, one factor — the built and
  sem names match exactly, 0 unmatched either way; gsem's cutpoints for an `ologit` outcome are
  `/y1:cut1 /y1:cut2 …`, and `_b[/y1:cut1]` resolves. Identical on StataMP 16.1.

**D7, `standardized` (decided 2026-10-04; built 2026-10-05): part of v1.** It was first deferred
from v1a as a release requirement, then built (step (3), below); the v1a deferral is kept as
history.
- **The v1a deferral.** `standardized` was rejected with a clear "not available in this version" message that
  points to the `nlcom` stopgap below. No standardized value, and no Θ table, was displayed in v1a.
  **Superseded 2026-10-05 by step (3), below.**
- **Step (3), the implementation (decided by the user 2026-10-05; built the same day).**
  `standardized` is accepted on the fit and on replay (`dwlsmv, standardized`). It **replaces**
  the coefficient table: $g(\hat\theta)$ in `e(b)`'s order and names, with delta-method SEs
  $\sqrt{\operatorname{diag}(J_g\,\texttt{e(V)}\,J_g^\top)}$ and the usual z, p and interval,
  followed by one Θ row `/var(e.y)` per indicator. **Θ appears in the standardized table only**;
  the default table is unchanged. A row whose standardized value is a constant (SE exactly 0:
  an exogenous factor's `var()`, a covariance fixed at 0) shows as `(constrained)`, as fixed
  rows do in the default table. Header and notes block are unchanged; one line above the table
  says the table is standardized (`std.all`) and that `/var(e.y)` rows are the derived residual
  variances. **Nothing is stored** (D8d unchanged): `dwlsmv_m11_std()` rebuilds the parameter
  table from `e(b)`'s stripe (D6) and `e(indicators)`, so replay needs nothing more; the display
  holds and restores the posted results. Verified against the D7 fixture on Case E's four fits
  (`tests/unit/test_standardized.do`). **Returned in `r()` (decided 2026-10-10, user, option (c), estout
  review; not built):** the standardized display returns `r(b_std)` (1 × k) and `r(V_std)`, named as
  the table rows (`e(b)`'s names plus the Θ rows `/var(e.y)`), still storing nothing in `e()`; a
  user adds them with `estadd matrix b_std = r(b_std)`. An eclass command leaves `r()` set by an
  rclass helper it calls last (`tools/probe_eclass_rreturn.do`, StataMP 18). Declined: (a) display
  only, (b) `e(b_std)` stored (a derivable quantity that can go stale). Not established: estout's
  layout of rows `e(b)` lacks (the Θ rows); the optional `tools/` diagnostic shows it.
- **Release.** Standardized estimates with delta-method standard errors were required before
  release and were built in this order, all three landing 2026-10-05: (1) the standardization map $g(\theta)$ and its Jacobian $J_g$
  in the derivation, with executable verification (**landed 2026-10-05**, derivation
  @sec-standardize; $J_g$ because @sec-gate5 already uses $G=\Gamma W\Delta$); (2) a tier-2 fixture from lavaan's
  `standardizedSolution()` through `fixtures/generate/`; (3) implementation. The derived residual
  variances Θ are reported in the same step, on the same machinery: Θ is derivable from `e(b)`,
  so nothing new is stored for it.
- **The map.** Under the delta parameterization every $y^*$ has unit variance (§4.3), so
  standardizing rescales by the factor standard deviations only: $\lambda_{jf}\sqrt{\Phi_{ff}}$,
  $\beta_{tf}\sqrt{\Phi_{ff}/\Phi_{tt}}$, covariances to correlations; thresholds are unchanged;
  $\theta_{jj} = 1 - (\Lambda\Phi\Lambda')_{jj}$ is already standardized.
  **$\Phi = (I-B)^{-1}\Psi(I-B)^{-\top}$ is derived, and `metric(unitvariance)` affects
  identification; it does not imply that all factors have unit total variance.** An endogenous
  factor's disturbance variance is 1 under it, its total variance is not (Case E, lavaan:
  $\Phi_{F2F2} = 1.280$), so its loadings change on standardization. A shortcut that rescales by
  $\psi_{ff}$, or skips rescaling under unit variance, is wrong for every endogenous factor.
- **The oracle convention, verified** (`tools/probe_d7_stdall.R`, lavaan 0.7.2, delta
  parameterization, Case E, both metrics): `std.all` equals `std.lv` — estimates identically,
  standard errors to $2.5 \times 10^{-12}$ — as unit-variance $y^*$ implies. The fixture may use
  either; the derivation states which (`std.all`).
- **Factor `~~` rows (decided 2026-10-05).** The row is $\psi_{fg}/(s_fs_g)$, lavaan's `std.all`
  convention. For an endogenous factor's own row that is $1-R^2_f$, shown as such. In the
  correlated case (both factors exogenous) it is $r_{fg}$, shown as the factor correlation;
  $r_{fg}$ for a pair involving an endogenous factor is not shown.
- **Fixture writer (decided 2026-10-05).** The D7 fixture is its own generator,
  `fixtures/generate/80_lavaan_standardized.R`, with its own manifest, written with
  `write_csv_lf()` (LF on every platform), columns `est_std` and `se`; it covers Case E's
  structural and equivalent correlated models under both metrics. **Landed 2026-10-05**
  (`fixtures/expected/lavaan_case_E_standardized.csv`, `fixtures/MANIFEST_lavaan_std.json`).
- **Stopgap for beta testers** (help file). Once `e(b)` and `e(V)` are posted (D5), `nlcom`
  applies the delta method itself. For a marker-identified exogenous factor,
  `nlcom _b[y2:F1]*sqrt(_b[/var(F1)])`. For an endogenous factor the total variance must be
  written out, e.g. for `F2 <- F1`, $\Phi_{F2F2} = \beta^2\psi_{F1F1} + \psi_{F2F2}$:
  `nlcom _b[y5:F2]*sqrt(_b[F2:F1]^2*_b[/var(F1)] + _b[/var(e.F2)])`. The help says the
  expressions are the user's to write and that `nlcom` is not checked against lavaan's
  standardized SEs; **owed** once Module 11 posts: one such check on Case E. **Done 2026-10-05**
  (`test_standardized.do` check 5: y2 and y5 against the D7 fixture, inside spec 7.2's tiers).

**D8a, return codes (decided 2026-10-04).**
- **A Stage 1 abort is rethrown as rc 3498.** Module 11 runs Stage 1 under `capture noisily`,
  never plain `capture` (which hides the `errprintf()` guidance), prints the mapping from
  positional labels to variable names, and exits with the caught rc. One code from Mata to ado
  level, so a test's `_rc == 3498` holds at both. The same applies to any other Mata rc 3498
  abort Module 11 calls through.
- **Module 11's own rejections exit rc 198** (invalid syntax): D3's non-unit scale fixes and
  second-order factors, the other D1 rejections, a non-numeric indicator (decided 2026-10-05; the message names the variables, and the
  wrapper converts `confirm`'s own rc to 198), sem's estimation options (D10), and, until D7 step (3)
  on 2026-10-05, `standardized`. The message is the
  user's guidance; it is printed with `display as error` before `exit 198`. Module 3's rc 3498
  stays the backstop for the same specifications. sem's own refusals keep sem's rc (498, D4).
- **Probed** (`tools/probe_d8a_catch.do`, StataMP 18 and 16.1, identical): `capture noisily`
  keeps the `errprintf()` and `_error()` text visible, also inside `quietly`; plain `capture`
  hides both; a legend printed after the catch follows the message; `exit _rc` gives the caller
  3498.

**D8b, positional labels to variable names (decided 2026-10-04).** Stage 1's label `yj` is the
j-th variable of sem's `e(oyvars)` (D1). Stage 1 is still not given names.
- **After a Stage 1 abort** Module 11 prints the full legend, one line per indicator
  (`y3 = q7_anx`), after the caught message and before the rethrow (D8a). No per-pair data
  returns on an abort, so relabelling only the offending pairs would need a Stage 1 change; that
  is rejected.
- **The near-bound note** is computed by Module 11 from `stage1.s` and `stage1.bound`, so it
  prints variable names directly (`q7_anx~~q9_dep`) and carries no legend.

**D8c, the opt-in's name (decided 2026-10-04): `keepbound`.** It sets Stage 1's `keep_bound`
to 1 and says what happens: pairs at the bound are kept rather than aborting. It implies no
repair and no admission of empty-cell 2×2 pairs or statuses 3–8. Since 2026-10-06 it has no interaction with
empty-cell pairs at all: a one-empty-cell pair is repaired by default whether or not `keepbound` is
set. The pairs that still abort (two empty cells on a diagonal, $|\hat\rho^+| \ge B$) are not
admitted by it (Module 2, *One-empty-cell repair*). O6 adds no option. Rejected:
`includebound` (close to repair language), `boundpairs` (names no action), `lavaanbound` (ties
the name to another program), a sub-option `stage1(keepbound)` (would set a pattern before O6).

**D8d, what `e()` carries beyond D5 (decided 2026-10-04).**
- **`e(bound)`**, posted always: `stage1.bound` as a 1 × P matrix, one column per pair in the ρ
  block's order, 1 for a pair included at $\pm B$ under `keepbound`, else 0. Its stripe is
  `varA:varB` (the pair's variable names as equation and column), so
  `colnumb(e(bound), "q7_anx:q9_dep")` finds a pair. It is the one stored copy: a count or a list
  of included pairs is derived from it, not stored. Near-bound pairs are not stored (decided
  2026-10-03). **Why not `varA~~varB`:** a stripe element is one name of at most 32 characters
  (`tools/probe_d8d_stripe.do`: 32 posts, 33 is "invalid name", rc 198, with or without a leading
  slash; identical on StataMP 18 and 16.1), and two variable names may each be 32. The
  near-bound note still displays `varA~~varB`.
- **Fit statistics** under project names: `e(T3)`, `e(df_T3)`, `e(p_T3)`; also `e(N)`,
  `e(metric)` (`marker` or `unitvariance`) and `e(cmd)` = `dwlsmv`. `e(p_T3)` is missing when
  `e(df_T3)` is 0 (recorded 2026-10-05 from the code, `dwlsmv.ado`). **Names kept, 2026-10-10 (user, option (a),
  estout review):** no `e(chi2)`, `e(df_m)` or `e(p)` (in Stata's convention a test against the
  constant-only model, not the saturated one); the help gives an `esttab … stats(N T3 df_T3 p_T3)`
  line. **`e(ll)` is never posted:** `lrtest` then refuses with "does not contain scalar e(ll)",
  rc 498 (`tools/probe_lrtest.do`, StataMP 18), which is right, since a likelihood-ratio test is
  not valid for WLSMV; nested comparisons use DIFFTEST.
- **Fit indices** (Module 12 *Fit indices*, option A, F5 frozen 2026-10-21): `e(rmsea_s)`,
  `e(cfi_s)`, `e(tli_s)` (family S, from $T_3$) and the c1 baseline's `e(T3_b)`, `e(df_T3_b)`, at
  every fit; missing by F3/F4's rules. Printed in the header under T3 as "RMSEA (scaled)", "CFI (scaled)", "TLI (scaled)", %9.4f, "." when missing (user, 2026-10-21, display "yes"; layout chosen at the build); shown only when posted.
- **For `difftest`'s refit** (single-group fits, 2026-10-23; Module 12 *Architecture*):
  `e(canon)`, `e(keepbound)`, `e(datasignature)`, `e(datasignaturevars)`. No Stata version,
  edition or OS is posted (decided 2026-10-23, user). **At $G>1$ (2026-10-31):** `e(canon)` and
  `e(canon_configural)` as before, plus `e(groupvar)`, `e(keepbound)` and the data signature over
  the indicators and the group variable.
- **Also posted** (recorded 2026-10-05 from the code; the help file documents them):
  `e(indicators)`, the indicators in D9's order; `e(cmdline)`; `e(version)`, the release
  (`0.2.0`), equal to `dwlsmv.ado`'s version line (added 2026-10-07, the user's; checked by
  `test_dwlsmv_gaps.do` G14 and the release-staging script); `e(properties)` = `b V` and
  `e(sample)`, set by `ereturn post`. **`e(Cns)` is posted only when some parameter is fixed**
  (D5's rows exist only then); a model with no fixed parameter has no `e(Cns)`. sem's `e(chi2_ms)` and
  `e(df_ms)` are not used: they suggest an ML chi-square, and sem's post-estimation commands would
  read them as one.
- **Module 10's singular-value ratios, stored and not printed (2026-10-07, the user's; this
  reverses their earlier "not stored" for these four values):** `e(gamma_ratio)` (Γ̂'s
  smallest/largest singular value, scale-free, within its rank), `e(gamma_ratio_raw)` (the same,
  raw), `e(gamma_rank)` (that rank) and `e(bread_ratio)` (the smallest/largest singular value of
  $W^{1/2}\Delta$). **At G > 1 (O-T10, decided 2026-10-12, user, option (a); built 2026-10-12, `test_invariance_post.do` check 6, no-π mutant red):** Δ is
  the stacked $(\sum_g D_g) \times$ npar matrix and $W = \mathrm{blockdiag}(\pi_g/\mathrm{diag}(\hat\Gamma_g))$,
  weighting (a), the W of the pooled minimizer and of `e(V)`'s bread; still one scalar. Test
  (user): a G = 2 fit with unequal $N_g$, asserted against an SVD computed directly; a mutant
  dropping $\pi_g$ fails it. Declined: (b) unweighted; (c) one ratio per group (a group's Δ is
  rank-deficient by design under ties). Each is missing when its row was not evaluated or has no supported rank.
  They are diagnostic values with no established cutoff and no practitioner reading, so they
  are kept for specialists rather than printed, and no option exposes them. Checked by
  `test_dwlsmv_gaps.do` G15 against Module 10 computed directly.
- **Not stored:** Module 10's flags and notes and the admissibility verdict (derived, displayed
  only); Θ and standardized estimates (D7).

**D8e, display (decided 2026-10-04).** The header's metric line is D2's; the near-bound note's
names are D8b's; the wording rule is Module 10's (fact first, a reading second, no advice as a
diagnosis).
- **After the coefficient table, one notes block,** in this order: the `keepbound` footnote, the
  near-bound note, the admissibility notes (Module 7), Module 10's diagnostic notes. Text that
  bears on whether to trust the fit sits in one place.
- **What prints on every fit, and the wording (decided 2026-10-07, the user's).**
  - Module 10's notes print only when a condition fires. The singular-value ratios no longer
    print on any fit (stored instead, D8d), and no user-facing text uses an internal release
    label such as "v1".
  - **Category counts print on every fit, and are descriptive diagnostics only.** Under
    "Category counts (the fewest observations in any one category):" they give the smallest
    lowest-or-highest-category count and the smallest middle-category count, with their
    items. In 0.2.0 they are not a warning and not an admissibility criterion: no threshold
    is applied and nothing is flagged from them. What counts as an inadequate category is the
    open `minok` decision, not settled here.
  - Wording rule for every note: keep a technical term where it helps and explain it at once
    in ordinary words; no matrix notation. Notes on an exactly singular Γ̂ state the condition
    and its consequence and stop: the remedy is an open policy.
  - Γ̂ asymmetry and a Γ̂ that is not positive semidefinite cannot come from Stage 1, so their
    notes are internal-check reports: "An internal numerical check failed … Please report this
    result to the package developers."
- **The `keepbound` footnote** names, by variable, the pairs included at $\pm B$
  (`q7_anx~~q9_dep`), built from `e(bound)`. Stage 1's own block, with positional labels, is
  unchanged and still prints; no Stage 1 change.
- **The D7 deferral leaves no trace in the output.** No line says that Θ or standardized
  estimates are not shown; the help file says it, and the refusal of `standardized` (rc 198, D8a)
  points to the `nlcom` stopgap. **Superseded 2026-10-05** by D7 step (3): `standardized` shows
  the standardized table and Θ; the default output still carries no line about them.
- **Replay (decided 2026-10-05).** `dwlsmv [, standardized level(#)]` with no model redisplays
  the header and the coefficient table, standardized if asked, from `e()`. The notes block is
  not redisplayed: its near-bound and diagnostic content comes from Stage 1 and Module 10
  results that `e()` does not carry (D8d). Tested by `test_dwlsmv_gaps.do` G7 and
  `test_dwlsmv.do` check 9.

**D8f, the lavaan/Mplus difference (decided 2026-10-04).**
- **A help section, "Comparison with lavaan and Mplus",** holds every cross-software fact: that
  lavaan 0.7.2 and Mplus 8.10 fit data with bound pairs with a warning where dwlsmv aborts by
  default, and why; that `keepbound` admits only the remaining status-2 pairs, at $\pm B$ as
  lavaan does, with an effect on SEs and T3 that is not established; D2's shared marker default;
  D7's `std.all` = `std.lv`. It describes the one-empty-cell repair (O6, decided 2026-10-06)
  and does not call the opt-in's effect small.
- **After the abort legend** (D8b), Module 11 prints one fixed sentence (shown on two lines) pointing to
  that section; no figures. Stage 1's text is unchanged.
- **Figures** about lavaan's `zero.add` or Mplus's repair value are quoted in help only after
  being re-checked against the most recent decisions row on the point when the help is written
  (the user, 2026-10-05: always the most up-to-date decision; written as the 2026-10-03 row).

**Lines in the command's terms after unchanged messages (decided 2026-10-04).** After a Stage 1
abort, "keep_bound = 1 above is the dwlsmv option keepbound.", only when the user did not give
`keepbound` and Stage 1 rerun with `keep_bound` = 1 passes. With `keepbound` given, Stage 1
has already run with it, so the line would be circular (made explicit 2026-10-05). After a Module 3 abort, for each factor with no `metric` line in the
canonical text: "Factor F has no fixed scale: fix one of its loadings at 1, e.g.
(F -> y@1 ...), or specify metric(unitvariance)." Stage 1's and Module 3's texts are unchanged
(D8a, D8f).

**A Stage 2 fit that does not converge is posted (decided 2026-10-04), as sem does:** `e(b)`,
`e(V)` and T3 at the last iterate (`dwls_inference()` does not check convergence, decided
2026-09-24), `e(converged)` = 0, and "convergence not achieved" in the header. Refusing to post
was rejected. `e(converged)` and `e(iterations)` join D8d's scalars. **Tested only by its header**
(`test_dwlsmv_gaps.do` G10, `e(converged)` set by hand): a fit that actually fails to converge
and is posted is not tested, and no iteration-limit option is added to force one (decided by
the user 2026-10-05).

**D9, indicator order and category recoding (decided 2026-10-04).**
- **Order.** The indicators are sem's `e(oyvars)` (D1), which lists observed indicators in order
  of first appearance in the path specification as written, not dataset order: a cross-loaded
  indicator keeps its first position, and a factor whose paths are split across pieces is not
  regrouped (`tools/probe_d9_oyvars.do`, five models, identical on StataMP 18 and 16.1). That
  order fixes Stage 1's $(Y, K)$ columns, the §4.2 index map, the positional labels `yj` (D8b),
  the threshold block of θ and the `/y:cutk` block of `e(b)`. Module 11 writes the `ind` line
  from `e(oyvars)` and builds $(Y, K)$ in the same order; Module 3 only reads the line. **Tested**
  (`tests/unit/test_dwlsmv_gaps.do` G1; `tests/unit/test_dwlsmv.do` check 8): a model whose path
  order differs from dataset order, asserting that Stage 1's columns, the `ind` line and the
  `/y:cutk` names agree.
- **Category recoding**, by Module 11, per indicator: (1) obtain the distinct observed values;
  (2) sort them; (3) assign consecutive internal category indices $0, \dots, K-1$; (4) reject an
  unobserved interior category under the current policy. $K$ is the number of distinct observed
  values. Recoding the observed support (`1,2,4 → 0,1,2`) is not the same operation as
  accepting a declared four-category variable and collapsing category 3: the first is ordinary
  recoding, the second is the O6 policy question.
- **Nothing declares a category in v1a.** Value labels are not read, and a gap in the codes is
  observed support, so step (4) has no trigger in v1a: it applies once something declares
  categories, and what happens then is O6's to decide.

**D10, sem's estimation options, `level()` and the example data (decided by the user
2026-10-05).**
- **sem's estimation options are refused, rc 198** (D8a): `vce()`, `method()`, `iterate()`,
  `technique()`, `tolerance()`, `ltolerance()`, `nrtolerance()`, `nonrtolerance`, `from()`,
  `startvalues()`, `startgrouping()`, `noivstart`, `difficult`. Before this they reached
  `sem, noestimate` on the dummy data (D1) and changed nothing: probed 2026-10-05, each gave rc 0
  and the plain fit's `e(b)` and `e(V)`, so `vce(cluster id)` returned unclustered standard
  errors without a word. `vce()`'s message says clustering and weights are planned for a later version (v2 internally). Options
  that change how sem reads the model (`covariance()`, `variance()`, `latent()`,
  `nocapslatent`, ...) still pass through; `nolog` passes and has no effect. Weights are already
  refused by Stata's `syntax` (rc 101).
- **`level()` is honoured** in the coefficient table, on the fit and on replay, with or without
  `standardized`; nothing is stored for it.
- **The example data.** The help file's examples use `dwlsmv_example.dta`, distributed with the
  command, not a test fixture: `fixtures/data/case_E.csv` unchanged, with labels, built by
  `tools/build_example.do` into `build/` (not committed). Users do not see the fixture names;
  the help calls it "the example data".
- Tests: `test_dwlsmv_gaps.do` G11–G13.

**D11, marker-metric recovery through unit variance (decided by the user 2026-10-06;
built 2026-10-07; Stata 16.1 milestone green at `b059f29`).** Derivation `@sec-marker-recovery`; evidence and the agreed direction in
E-016, *Marker-metric start trap and recovery through unit variance*. Under the
marker metric every free loading starts at 0.7, Module 3's default (corrected 2026-10-07; was "+1"). With reverse-keyed indicators, the fit can stop at
a degenerate non-converged point: on 200 of 200 replicates of the trap design, where unit variance
and a reverse-coded marker fit both reach the same T3.

- **Scope.** Applies only to `metric(marker)`, the default. A `metric(unitvariance)` fit is never
  recovered.
- **Trigger** (decided with the user 2026-10-06). The first marker fit is recovered when any of
  these holds:
  - `converged` = 0;
  - `dwlsmv_admissible()` flags any condition;
  - Module 10's bread row reports numerical rank deficiency.

  Non-convergence alone is present on every trapped fit in the study. Inadmissibility alone
  would miss 47 of 200.
- **Eligibility** (all-or-nothing, the user's rule; derivation Result 2), read from Module 3's
  marker parameter table. All three must hold:
  - every factor has exactly one fixed nonzero loading, equal to 1;
  - no ψ_ff is fixed;
  - every other fixed λ, β or ψ is exactly 0.

  If the model is ineligible, recovery is not attempted: the first fit is reported, with a note
  naming the reasons.
- **Procedure.** Every step reuses existing machinery, with no change to the locked API:
  1. Module 11 builds the unit-variance canonical text from the marker text. For each
     `metric marker F y`, three edits (corrected 2026-10-07, the user's; the swap alone fails):
     - replace that line with `metric unitvariance F`;
     - add `lambda F y free`, because the marker loading exists only through the metric line;
     - drop `psi F F free`, which Module 11 always writes. Kept, it collides with the metric
       line, and `dwls_identify()` refuses two constraints on one element.

     This is `mr_uv_table()` in text form. Module 3 expands the result, as it expands any
     `metric` line (§4.8).
  2. `dwls_identify()`, then `dwls_fit()` on the same Stage 1 output. Stage 1 is not re-run.
  3. The transform (derivation Result 1). With c_f the marker's unit-variance loading:
     - λ ↦ λ/c_f;
     - β_fg ↦ β_fg·c_f/c_g;
     - ψ_fg ↦ ψ_fg·c_f·c_g.

     Thresholds are unchanged.
  4. The transformed values become the starts of the free λ, β and ψ rows in the marker
     canonical text (`free <start>`, §4.8). Each start must read back as the identical
     double; the written format is decided, below. Thresholds cannot be carried in the text:
     the grammar has no threshold directive (2026-10-07, the user's). Then `dwls_identify()`
     runs through the normal marker path.

     **The final fit skips ULS** (decided 2026-10-07, the user's, after measuring). It is
     DWLS only, `dwlsmv_fit_from()` (Module 7), from the table's λ/β/ψ starts and with every
     threshold started at Stage 1's sample threshold s_τ. Δ's threshold block is the
     identity and σ_τ(θ) = τ, so with a diagonal weight the optimum's thresholds are s_τ
     exactly (derivation @sec-marker-recovery, Result 5); the start puts them there. Then
     `dwls_inference()`. Measured on the trap and
     case E (E-016): 2 DWLS iterations against 7 after ULS on the trap, θ̂
     within 2.1e-7 and F 1.3e-14 lower. The derivation's one-iteration prediction is not
     asserted.

     If `dwls_identify()` finds Δ rank deficient at the transformed starts and replaces them
     by its second point θ₂ (§4.8), the recovery goes on from θ₂. That is not a failure. The
     note says so, without the parameter values (2026-10-07, the user's).

  The reported fit is that final marker fit. The unit-variance fit is a recovery device only:
  nothing is chosen by comparing F across starts.
- **Success and failure.**
  - Recovery **succeeds** when the final marker fit converges. Its admissibility and Module 10
    notes are reported as for any fit: a genuine Heywood case stays inadmissible, correctly.
  - Recovery **fails** when the unit-variance fit does not converge, when some c_f is 0 or
    non-finite (the map does not exist), or when the final fit does not converge. Then the first
    marker fit is reported, unchanged, with its own warnings and a note that recovery was tried
    and why it failed.
  - **An abort during recovery is a failed recovery, not a failed command** (2026-10-07, the
    user's). This covers any `_error()` in `dwls_identify()`, `dwls_fit()`,
    `dwls_inference()` or the transform, including θ₂ also being rank deficient. Mata cannot
    catch an abort, so the recovery fits run in a nested `_stata("mata: …")` call. Its rc
    reaches the caller, and the note names the step that stopped. The abort's own message is
    not shown.
  - **A Break during recovery is a failed recovery too** (2026-10-07, the user's). The
    nested call returns rc 1 like any other nonzero rc; the first fit is reported with
    the failed note, and the command goes on.
  - **An abort before the trigger is read does not fire it** (2026-10-07, the user's).
    The first fit's `dwls_inference()` and `dwlsmv_breaddiag()` run outside the nested
    call; an abort there ends the command, as before D11, and recovery is not attempted.
  - A unit-variance fit that converges but is inadmissible is still transformed
    (2026-10-07, the user's).
- **The near-zero c_f cutoff: closed 2026-10-07 (user).** Only an exactly zero or non-finite
  c_f stops recovery. There is no positive floor; no evidence sets one (derivation, *c_f → 0*).
- **Version: 0.2.0** (decided 2026-10-07, the user's). `v0.1.0` stays at `aca3ac8`, unchanged. D11 adds a
  stored result and changes the output, so it ships as a minor release.
- **Iterations.** Each of the three fits has the existing budget (300 per `dwlsmv_minimize()`
  call). There is no new option (D10 refuses `iterate()`).
- **Reporting.**
  - `e(recovery)`, a string, always posted: `none`, `recovered`, `ineligible` or `failed`.
    - It is `none` when the trigger does not fire, whatever the eligibility.
    - It is `none` under `metric(unitvariance)`.
    - Eligibility is judged only once the trigger has fired.

    (2026-10-07, the user's.)
  - `e(converged)` and `e(iterations)` are the reported fit's. The first fit's and the
    unit-variance fit's iteration counts are not stored (2026-10-07, the user's).
  - The trigger is read as follows (2026-10-07, the user's):
    - an admissibility flag of exactly 1 fires it; a missing flag (a condition that could not
      be evaluated) does not;
    - only bread row 3, numerical rank, fires it. Row 2, a zero column, does not.

    A converged but inadmissible first fit, such as a genuine Heywood case, is recovered too.
    Its note says what fired.
  - The notes block, first in D8e's order (wording 2026-10-07, the user's). `<what fired>` lists
    the causes: "did not converge", "was not admissible", "could not determine all parameters
    reliably (the information in the fit was insufficient)" (the bread's numerical rank;
    reworded 2026-10-07 from "had a numerically rank-deficient bread", the user's).
    - Recovered: "The marker-metric fit `<what fired>`; it was restarted from the
      unit-variance solution, transformed to the marker metric." When θ₂ was used, a second
      line: "Recovery used the second identification start because the transformed start
      was rank deficient."
    - Ineligible: "The marker-metric fit `<what fired>`; recovery through the unit-variance
      solution was not attempted, because the model is not eligible for it:", then one
      line per reason, naming factors and parameters.
    - Failed: "The marker-metric fit `<what fired>`; a restart from the unit-variance solution,
      transformed to the marker metric, was attempted. Recovery attempt failed; first fit
      retained:", then the reason.
  - Replay does not redisplay it (D8e); `e(recovery)` carries it.
- **Open for implementation, not decided here.**
  - `dwls_fit()` always runs ULS from the starts, then DWLS from the ULS point (Module 7).
    The final fit therefore does not begin DWLS exactly at the transformed point. The
    derivation's prediction of about one DWLS iteration is not asserted. **Decided
    2026-10-06 (user): keep ULS for the first build** and measure the final fit on the trap
    fixture. Skipping ULS when starts are supplied would change Module 7's contract, and is
    revisited only with that evidence. **2026-10-07 (user): investigate skipping ULS for the
    recovery's final fit.** First measure threshold initialization on the trap fixture:
    Stage 1's sample thresholds as starts, or skipping ULS only when every free λ/β/ψ row has
    a start. Then amend Module 7 and D11 together. **Decided 2026-10-07 (user), after
    measuring: the recovery's final fit skips ULS**, with thresholds at s_τ (Procedure step
    4). "Starts supplied" is not a general notion at the `dwls_fit()` boundary: only the
    recovery's final fit skips ULS, by calling `dwlsmv_fit_from()`. The locked API is
    unchanged. Every other fit, the first marker fit and the unit-variance fit included,
    still runs ULS, and is unchanged bit for bit (measured on the trap and case E).
  - **Decided 2026-10-07 (user): the precision used to write the starts.** Each start is written
    in the first of three forms that reads back as the identical double: `strofreal(x)`,
    `%21.0g`, then `%21x`. The read-back is asserted. The same rule now writes fixed values
    in `dwlsmv_tr_num()` (Module 11's translate step), which until then stopped at
    `%21.0g`. Measured 2026-10-07 on StataMP 18, with random doubles, 40,000 per decade:
    - `%21.0g` reads back for every value with |x| ≥ 10⁻². It fails for 22% of values near
      10⁻³ and 92% at 10⁻⁴ and below, where it switches to e-notation with about 15
      significant digits.
    - `%21x` reads back for every value tested.
- **Tests (built 2026-10-07):** `tests/unit/test_marker_recovery.do` and
  `test_translate.do` check 7. They cover:
  - the tier-1/2 marker-trap fixture, on which the trigger fires and the recovered fit matches
    lavaan's marker solution and the transformed reference;
  - case E unaffected (its marker fit converges, so `e(recovery)` = `none`);
  - an ineligible model with the trap pattern (`ineligible`, no recovery);
  - a failure path. A recovery-time abort is forced by giving `dwlsmv_m11_recover()` a text
    that `dwls_identify()` refuses. The c_f = 0 and non-finite refusals are tested on the
    transform. A non-converging unit-variance or final fit, and a θ₂ substitution, have no
    known input and are not exercised.
  - the starts' precision.
  - a Stata 16.1 milestone, green at `b059f29` (2026-10-07).

Two obligations the surface carries whatever grammar it lands on, both of which are surface
conventions and must not leak downward:

- **Marker selection.** If the syntax implies a marker by position — the first-listed indicator,
  say — that convention lives here. The canonical form names the marker explicitly
  (`metric marker F y`), so Module 3 is never asked to reconstruct it.
- **Metric is not expanded here.** Module 11 emits `metric unitvariance F`; it does not decide
  which parameter that fixes. Section 4.8 assigns that expansion to Module 3, because the
  disturbance-versus-total reading depends on the model graph.

**Residual covariances among indicators (`residcov(b e)`) are phased to v1b, not v1a.** The
original call was to support them from the start; on reflection, deferring them is the better
call, for a specific implementation-risk reason: residual covariances add a real complication
(keeping $\Theta$'s off-diagonal block admissible, per Module 7's admissibility check) at exactly
the point where the base backend — Modules 2–10 — is least proven. Splitting it out:

- **v1a:** no residual covariances — measurement paths, factor-to-factor regressions and the
  metric choice only. This is what Modules 0–11 should target first, validated end to end
  (Section 7) before anything else is added. Note that the canonical form of Section 4.8 has no
  residual-covariance directive either, so v1a's scope is enforced at the seam rather than only at
  the surface.
- **v1b:** add residual covariances once v1a passes its full validation suite. The syntax and the
  underlying $\Theta$ admissibility handling were already designed above and in Module 7 — v1b is
  "turn it on and validate it," not a redesign — but it only starts once v1a is a proven backend,
  so a residual-covariance bug is never competing for attention with a base-estimator bug.

**Observed exogenous covariates (an observed `x` predicting a factor) are out of scope for v1.** The
tempting argument for allowing them — recorded because it is wrong in an instructive way, and
someone will reach for it again — is that a right-hand-side-only covariate `x` could enter
a structural regression without needing threshold/polychoric treatment, on the reasoning that it's "MIMIC-style"
and never appears in `s`. That reasoning doesn't hold up: estimating a structural coefficient on
`x` (even one nothing in the model predicts) requires knowing how `x` covaries with the ordinal
items whose factor it predicts — e.g. to estimate `F2`'s regression on `x`, the model needs
$\text{Cov}(y_j^*, x)$ for `F2`'s indicators `d`, `e`, `f`, and that covariance *is* a polyserial
correlation between a continuous variable and an ordinal item's latent response. There is no way to
estimate that path from thresholds and polychoric correlations alone.

lavaan's own `conditional.x` option makes this precise and confirms it: with `conditional.x =
TRUE`, the model is conditioned on the x-covariates and the model-implied sample statistics exclude
them entirely (a genuinely different, more involved estimation approach with its own
sample-statistics machinery); with `conditional.x = FALSE` (lavaan's default), covariates are
modeled jointly with everything else, which is exactly the mixed continuous/ordinal sample-statistics
problem Section 3 already deferred to Module 12. There's no third option where a covariate
participates in a structural regression for free. **v1 syntax therefore supports regressions among
factors only; a future covariate extension needs its own explicit choice between the
conditional-on-x and joint-modeling approaches, each with its own extension to `s`, Γ, and
$\sigma(\theta)$, and shouldn't be added as an afterthought to either.**

Investigate — don't assume — how to wire Modules 3–10 into Stata's estimation-command machinery,
but go in with a strong prior against the `ml` evaluator path, for a specific, verified reason
(checked directly against Stata's `ml` manual, not assumed): it's *not* that `ml` categorically
requires the objective to decompose as a sum of per-observation log-likelihood contributions —
that restriction only applies to the `lf` evaluator family (`lf`/`lf0`/`lf1`/`lf2`); the `d0`/`d1`/
`d2` family explicitly does *not* require $\ln L = \sum_j \ln \ell_j$ and can evaluate an arbitrary
dataset-level scalar objective, which would in principle fit $F_{DWLS}$'s structure fine. The actual
blocker is narrower and more decisive: **Stata's own documentation states that `vce(robust)`,
`vce(cluster ...)`, and `pweight` are not supported by `d0`/`d1`/`d2` evaluators and will produce an
error if requested.** Since Module 8's entire inference layer *is* a custom sandwich VCE built from
Γ — not `ml`'s default OIM/OPG information-matrix-based covariance, which assumes actual
likelihood-based asymptotics that don't apply here — routing through `ml` would mean fighting past
a VCE mechanism that explicitly won't accept the kind of variance estimate this estimator needs,
for no compensating benefit. That's a narrower and firmer objection than "it's not technically a
likelihood," and it's the one to cite when deciding.

Given that, compare (1) `moptimize()` + explicit `ereturn post`/`ereturn scalar`/`ereturn matrix`
(entering the pipeline as the leading candidate on the strength of the argument above — it keeps
the Mata code directly in control of inserting Module 8/9's results into `e(V)`, `e(b)`, and scalar
fit statistics without an intermediate abstraction layer fighting it), (2) an `ml` evaluator anyway
(mainly to confirm the objection above actually bites in practice, not to default to it), and (3)
`optimize()` + a fully custom wrapper, on gradient/Hessian support, parameter transformations
(Module 3's identification constraints), convergence reporting, and
`e(b)`/`e(V)`/`predict`/`estat`/`test`/`lincom` compatibility. Frame the decision empirically
rather than as a foregone conclusion: **use `moptimize()` unless testing demonstrates a compelling
postestimation advantage to `ml` that outweighs the need to work around its VCE machinery** — the
argument above is a strong prior, not a substitute for actually checking. Ensure `e(b)` covers
everything Module 3 identified as free and `e(V)` is Module 8's $\text{ACOV}(\hat\theta)/(N-1)$, via `dwlsmv_ev()`.

**Module 12 — extensions (future, not v1).** All of them build on Modules 2–10's backend. By
release (decided 2026-09-26):

- **v2:** DIFFTEST (highest priority), continuous and mixed indicators, modification indices,
  pairwise-present sample statistics (Section 2), and multiple-group models. **Release 0.3.0
  (decided 2026-10-08) carries DIFFTEST with multiple-group invariance**, the first part of v2;
  the rest of v2 is not scheduled.

  **Invariance in 0.3.0 (decided 2026-10-08).** Levels are produced by a generator that emits
  a literal canonical model plus `tie` lines; Module 3 sees only ties and literal directives,
  so "`metric` is the only exception" and check 8's exemption stand as decided. The
  unit-variance metric only; structural paths ($B \ne 0$) are refused with invariance levels.
  A level step with $d = 0$ (an *equivalent* model, @sec-invariance-id Result I2: for example
  the loadings step when a factor has exactly one polytomous item) **[single requested level,
  `invariance(level)`, decided 2026-10-12, user, option (a); not built: the level is fitted and
  posted as generated, and the notes block (D8e) says it is equivalent to the previous level and
  why, detected before the fit from the previous level's generated text through `dwls_identify`
  (integer df difference); nothing new in `e()` (derivable). Test: one polytomous item per factor
  under `invariance(metric)`, note printed, T3/df equal to the thresholds fit (tolerance set at the
  build from the architecture); a no-note mutant fails. Declined: (b) refuse rc 198; (c) no note.
  What "skipped" means in a chain stays with the exemption gate. **Settled 2026-10-13 (user,
  option (c); not built): see check 8's *Zero-df generated pairs*.** The single-level fit and
  its note below are unchanged. **The note's text (decided 2026-10-12, user: (c) where verified, (i) elsewhere, refined to separate equal df from established equivalence):** equal df is all the detector measures, so the equivalence wording is given only where a derivation establishes it. *Verified* = the metric step with every factor having at most one indicator of three or more categories and all its loadings free (Result I2: one polytomous item, the step is empty; none, φ fixed, equivalent to configural): "Note: invariance(metric) is equivalent to invariance(thresholds) for this model, so this step adds no testable restriction: both have df = #, and every factor has at most one indicator with three or more categories, the case in which the equivalence is established." *Any other d = 0* (e.g. the thresholds step with only 3-category or binary items, the scalar step with a single-indicator factor; counted by the main session from G8, not verified): "Note: invariance(L) has the same degrees of freedom as invariance(L−1) for this model (df = #). Equal degrees of freedom do not by themselves make two models equivalent, and equivalence is not established for this pattern." A pattern moves to the verified wording only after its own derivation lands. Built 2026-10-12 (`dwlsmv_m11_fit_g`, `_i2`; `test_invariance_flow.do` check 6; no-note and always-verified mutants red); the note prints before the step-3 stop until step 4 moves it into D8e. Test tolerance: |T_raw − T_raw,prev| ≤ ε = (N − G) nrtol (check 10's band, its κ ≤ 2 allowance chosen after measuring); df exact; T3 difference printed, not asserted. Declined: (ii) no note for unverified patterns; (iii) refuse them.]** is reported as equivalent
  to the previous level and skipped; the output names the skipped level and why. The
  singular all-binary case (@sec-invariance-id Result I5) gets **one guard**, before the
  intercepts-level fit, from the configural estimates (the spread of $\hat\tau_j/\hat\lambda_j$
  per all-binary factor; the cutoff a placeholder until measurement and derivation justify a value, bounded by the confirmatory interval), and, after it,
  **post-fit conditioning diagnostics** (decided 2026-10-09, user): P1 and the Schur-complement conditioning P3 of the $(\phi_g,\kappa_g)$ block (@sec-singular-guards Results P1–P4) are reported as numbers with Module 10's notes, with no cutoff and no warning; they are not a guard and carry no decision authority, and the implementation and help must not imply otherwise. **In 0.3.0 neither runs (decided 2026-10-12, user, T0, T1, T2; the narrow-route cutoff selection was inconclusive):** `invariance(scalar)` with an all-binary factor is refused at rc 198 by Module 11 before any fit (built, `dwlsmv_m11_gcheck`, `test_invariance_flow.do`), and neither P3 nor task T's P1 lands. The "one guard" above is the design for a later version, with the cutoff still a placeholder. **Configural refit, to be revisited:** no configural fit is run inside `invariance(level)` in 0.3.0 (decided 2026-10-12, user, option (a)); its one recorded purpose was feeding this guard, so the 2026-10-12 caller decision (always refit configural, then the guard) applies only once a cutoff exists, and the refit is revisited with the guard cutoff work (E-017). The pooled multiple-group scaling is derived in
  `@sec-multigroup-scaling` and decided (2026-10-09): $T = (N-G)F$, `e(V)` $= \mathrm{ACOV}/(N-G)$.

  **Fit indices: scope open, 0.3.0 or 0.4.0 (recorded 2026-10-09, user).** The spec defines no
  fit index today; the only absolute test is $T_3$. For the invariance workflow the coherent unit is
  **RMSEA, CFI and TLI, all from the same $T_3$ convention, with ΔCFI and ΔRMSEA computed across the
  generated levels against one common baseline model** (user). Nothing is built or derived. Before
  any code, a derivation must settle at least: which RMSEA for a scaled-and-shifted statistic (the
  naive form from $T_3$, or a corrected "robust" form), its multiple-group convention, the baseline
  (independence) model's canonical form under groups and why one baseline serves every level, and
  how CFI/TLI use the baseline's $T_3$ and df. SRMR and WRMR are not part of this unit.

  **Scope decided 2026-10-10 (user): option A, in 0.3.0.** In: the *scaled* RMSEA, CFI and TLI
  computed from $T_3$ and its df (family S of the cloud task P memo; lavaan's `*.scaled` form),
  with the **c1 baseline** (thresholds free per group, every ρ fixed at 0), computed in closed
  form from Stage 1 with no additional fit; ΔCFI and ΔRMSEA across generated levels against that
  one baseline. Validation: single-group lavaan 0.7.2 tier-2 fixtures, and multiple-group tests
  that validate the project's own $T = (N-G)F$ convention (which differs from lavaan's
  $\sum_g (N_g-1)F_g$ by about $10^{-4}$ relative at $G>1$, so those comparisons are
  comparison-specific, never a loosened tolerance). Documentation must distinguish a fit-index
  comparison from DIFFTEST (ΔCFI is a difference of $T_3$ excesses, not a test; $T_3$ need not
  rise along a chain, so ΔCFI can be positive) and must state that no published ΔCFI/ΔRMSEA cutoff
  is validated for DWLS-$T_3$. **The known R2 failure is preserved:** the draft's preregistered
  assertion that its catML RMSEA reproduces lavaan 0.7.2's `rmsea.robust` fails (lavaan 0.7.2
  divides by $N$, the draft by $N-G$), and whatever lands (derivation chunk, fixture or test)
  must record that failure as a failure, so no scratch workaround can read as a passing
  preregistered check. **Deferred to 0.4.0 or later:** the robust families (R1, T3/a; R2,
  catML), any claim about their population targets or finite-sample behaviour with $T_3$, any
  recommended DWLSMV ΔCFI/ΔRMSEA threshold, and any claim that robust indices are validated
  against Mplus WLSMV. **Still open within A** (derivation first): the multiple-group factor
  ($\sqrt G$ or none), TLI clipping at 1, the $D = 0$ and negative-excess conventions, the
  `e()` names, and the treatment of skipped $d = 0$ levels in the Δ table. Until built, nothing
  may imply these indices exist; Module 11 does not print them.

  **Derived 2026-10-21, not built (`@sec-fit-indices`).** Family S, exactly lavaan 0.7.2's
  `*.scaled` (Result FI1): $\mathrm{RMSEA}=\kappa_G\sqrt{\max(T_3-D,0)/(nD)}$ with $n=N-G$ (the
  divisor of the decided $T=(N-G)F$, and lavaan's), $\mathrm{CFI}=1-\max(E,0)/\max(E,E_B,0)$,
  $\mathrm{TLI}=1-E\,D_B/(E_B\,D)$, $E=T_3-D$. **c1 baseline (Result FI2), from Stage 1 only:**
  $T_B=(N-G)\sum_g\pi_g\sum_k\hat\rho_{gk}^2/\hat\Gamma_{g,\rho_k\rho_k}$, $D_B=G\,p(p-1)/2$,
  $\operatorname{Tr}(M_B)=D_B$ and $\operatorname{Tr}(M_B^2)=\sum_g\lVert R_g\rVert_F^2$ ($R_g$ the
  correlation matrix of $\hat\Gamma_{\rho\rho,g}$), both free of $\pi_g$, so
  $T_{3B}=a_BT_B+D_B(1-a_B)$. One baseline serves every generated level because it depends on the
  shared Stage 1 alone. **Finding:** c1 is nested in the configural model but not in levels that
  constrain thresholds (with $K\ge4$), so it is a common reference point, not a model nested in
  every level (lavaan's default baseline has the same property). **ΔCFI/ΔRMSEA (Result FI3)** are
  differences of consecutive levels' indices; ΔCFI $=(E_{L-1}^+-E_L^+)/E_B$ when $E_B$ dominates,
  and can be positive. Verified: the closed form against the generic $U$/$M$ construction and a
  numerical minimum (≤ 3e-16 relative), and against lavaan 0.7.2 `baseline.chisq`,
  `baseline.chisq.scaled`, `rmsea/cfi/tli.scaled` at $G=1,2$ on lavaan's own $T$ convention
  (≤ 2e-16); $n=N$ and $\kappa_G=1$ each miss lavaan. The R2 failure is recorded there as a failure.
  **F1–F6 decided and frozen 2026-10-21 (user; F4 amended):** F1 $\kappa_G=\sqrt G$, $n=N-G$;
  F2 TLI unclipped (raw $E_L$ in TLI); F3 $D_L=0$ or $T_{3L}$ missing: all three missing;
  **F4, two independent conditions:** the *model side* $E_L=T_{3L}-D_L$ decides only the model's
  clamp (RMSEA $=0$ when $E_L\le0$; $\max(E_L,0)$ in CFI) and never whether CFI/TLI are defined;
  the *baseline side* $E_B=T_{3B}-D_B$ decides only whether CFI and TLI are defined (both missing,
  with a note, when $E_B\le0$, $T_{3B}$ missing or $D_B=0$), from the baseline's own $T_3$ and df,
  and never affects RMSEA; for $E_B>0$, CFI $=1-\max(E_L,0)/\max(E_L,E_B)$; exact comparisons;
  F5 `e(rmsea_s)`, `e(cfi_s)`, `e(tli_s)`, `e(T3_b)`, `e(df_T3_b)`; F6 a skipped $d=0$ level's row
  shown, marked "skipped (d = 0)", Δ blank, the next Δ against the last fitted level.
  **Built 2026-10-21** (StataMP 16.1): `dwlsmv_fitidx_base`, `dwlsmv_fitidx`,
  `dwlsmv_fitidx_delta`, `dwlsmv_m11_fitidx`; every fit, $G=1$ and $G>1$, posts the five `e()`
  scalars (Module 11 *Fit statistics*); printed in the header (2026-10-21); the help's Fit indices section, Output and Stored results entries approved by the user 2026-10-22.
  `dwlsmv_fitidx_delta`'s caller is `difftest` (`dwlsmv_dtfront`, G = 1 from 2026-10-23, G > 1
  from 2026-10-31) through `dwlsmv_dtfront_delta`: the restricted model's value minus the less
  restricted one's (H0 minus H1, FI3's sign; user 2026-10-31). Tests: `test_fitidx.do` against
  `fixtures/expected/lavaan_fitindices.csv` (`99f_lavaan_fitindices.R`, lavaan 0.7.2); mutants
  M1–M8 red, M1 the merged-condition F4 mutant. **G > 1 tier 2 (2026-10-31, PASS):**
  `test_mg_fitidx.do` against `99j_lavaan_mg_fitidx.R` (lavaan 0.7.2 converted to $(N-G)F$;
  designs E, D, P, N, four levels; indices, then ΔCFI/ΔRMSEA, both signs;
  E-018).
- **0.3.0 addition (user, 2026-10-26): single-group EFA with oblique geomin rotation, point
  estimates only.** The unrotated fit is v1 unchanged (`test_efa_unrotated.do`); new are the
  rotation (derivation first: geomin criterion and gradient, oblique GPA, convergence, random
  starts), the rotated $\Lambda$ and $\Phi$, and Module 11 surface. **Rotated-loading SEs are not
  in 0.3.0** and are suppressed, not shown unrotated. Derivation: `@sec-geomin`.
  **Settled conventions (user):**
  - *geomin ε (2026-10-29):* Mplus's rule (Browne 2001), $\varepsilon = 10^{-4}$ at $m=2$,
    $10^{-3}$ at $m=3$, $10^{-2}$ at $m\ge4$, overridable by the user. lavaan's default (0.001 for
    every $m$) differs at $m\ne3$; tier-2 fixtures set lavaan's `geomin.epsilon` to this rule.
  - *Row normalization (2026-10-29):* none, and no Kaiser option in 0.3.0 (lavaan's and
    Mplus's default).
  - *Starts, local optima, convergence (2026-10-29):* 30 deterministic starts (identity plus a
    fixed grid of rotations; no RNG state used or changed), count set by `rstarts(#)`. Each start
    iterates oblique GPA until step halving fails (the stall, `@sec-geomin`); a start is converged
    iff its final $\|G_p\|_F < 10^{-5}$. The best converged start is selected; no converged start
    is an error. Output notes the number of distinct minima and of starts reaching the best.
    Stored in `e()`: starts attempted, converged starts, distinct minima among converged starts,
    converged starts reaching the selected minimum, the best criterion $Q$, its final
    $\|G_p\|_F$, the ε used, and the gap in $Q$ between the best and second-best distinct
    minimum (missing when only one is found). Two converged starts are the same minimum iff their
    $Q$ differ by at most $10^{-8}+10^{-6}Q_{\text{best}}$ (the distinct-minimum tolerance, set
    2026-10-29 in Unit 5; separate from the coverage check's $\delta_{\text{cov}}$). Its basis: a
    converged start ($\|G_p\|_F<10^{-5}$) sits within about $\|G_p\|^2/(2h)\approx5\times10^{-9}$
    of its minimum for curvature $h\ge10^{-2}$, so the tolerance exceeds that error while staying
    well below any $Q$ difference that changes a loading at the θ̂ tier.
  - *Start grid (2026-10-29):* pinned as `@sec-geomin` Result GM4 (identity, then products of
    Givens rotations at angles $\pi\{k(q+1)(\sqrt5-1)/2\}$). `set seed` has no effect and the
    help file says so; the grids are nested in `rstarts()`. Changing the grid is a change of
    estimator default. **Random starts: contingent extension, not in 0.3.0.** Triggered only by
    specific evidence from Unit 4: on a fixture case, random starts (lavaan's, or a seeded
    random-start run of the same GPA) reach a criterion lower than the grid's best by more than
    the distinct-minima tolerance. The grid's unproven coverage alone is not a trigger.
  - *Sign and column order (2026-10-29):* lavaan's rule. Reflect each column of Λ whose loading
    sum is negative; then order columns ascending by the mean row index of the loadings with
    $|\lambda|\ge0.8\max|\lambda|$ in that column (lavaan's `"index"`). Φ is reflected and
    permuted with Λ. Replaces `dwlsmv_signcanon()`'s first-free-loading anchor for EFA. Mplus
    documents no rule; on Case E it gave the same order and signs.
  - *Surface (2026-10-29):* an EFA block in the existing path syntax, after Mplus's
    `F1-F2 BY y1-y6 (*1)`: `dwlsmv (F1 F2 -> y1-y6), efa(F1 F2) [epsilon(#) rstarts(#)]`.
    `efa(namelist)` marks the block's factors; Module 11 writes the echelon canonical form
    (anchor zeros, `unitvariance` on the block, no factor covariance) and the rotation runs after
    the fit. `epsilon()` overrides the Mplus-rule ε; `rstarts()` sets the start count (default
    30). In 0.3.0 the model must be exactly one complete block: every latent variable in `efa()`,
    every block factor loading on every block item. Refused at rc 198: other paths, a second or
    partial block, `group()`, `invariance()`, `partial()`, `ginvariant()`, `metric()`,
    `standardized` (the rotated Λ is already std.all), negative df; `difftest` refuses an EFA
    result. df = 0 and $m=1$ are allowed. v3 ESEM lifts the one-block refusal without new syntax.
  - *Posted results (2026-10-29, user):* `e(b)` holds the rotated loadings, the factor
    correlations (under sem's `/cov(F,G)` names) and the thresholds; **no `e(V)`**, so `test` and
    `lincom` cannot run on an arbitrary basis. `e(Lambda)` and `e(Phi)` as matrices; the eight
    `rot_*` scalars; `e(efa)`, `e(rotation)`. T3, df and the scaled fit indices are reported as
    for any fit (rotation-invariant, GM1); at df = 0 the indices are missing (F3) and a note says
    the model fits exactly. `difftest` refuses an EFA result at rc 3498, its own convention.
  - *No converged rotation start (2026-10-29, user):* rc 3498, nothing posted, and the reason in
    the command's terms: whether the unrotated fit converged and is admissible, and that this
    usually means more factors than the data support (seen on Case E at $m=3$, where lavaan's own
    fit is degenerate).

  No `group()`/`invariance()`/`difftest` with EFA in 0.3.0.
- **v3:** ESEM (rotated SEs, EFA blocks inside SEM, target rotation, multiple groups).

  **Architecture principles for the full multiple-group module (direction only; recorded
  2026-10-10, user).** These set direction. They reopen nothing decided for 0.3.0: O2, O3 and O4
  stay as decided until their own derivations land (O2's 0.3.0 scope (B) is recorded at O2).
  1. *Single model system.* Generated and user-written models share one canonical representation
     and one estimation and inference pipeline. Nothing downstream of the canonical text knows
     whether a model was generated (§4.8 G8). The one current exception is check 8's
     generated-chain exemption; it stays until principle 5 provides a derived replacement.
     Check 11's configural preflight, the singular-case guard, and reporting or skipping a
     $d = 0$ step are actions of the generator's caller, not properties of a model and not
     exceptions to canonical-representation parity.
  2. *Generator parity.* The generator emits only models a user could write in the canonical
     form, subject to the same admissibility, identification and validation rules (e.g. O4's
     whitelist applies to every `tie` line, whoever wrote it). Once Module 11's surface can express
     a construct, the generator's use of it must be expressible there too.
  3. *Derivation-gated capabilities.* The refusals in O2 (β with groups), O4 (tie types; β at
     least, for DIFFTEST on slopes) and O3 (cross-row ties) are lifted independently. Each needs,
     before it is lifted: an identification or rank result; a reference-tool check; DIFFTEST
     validity for the capability, which includes re-establishing ĥ's weighting for it (the
     unweighted-ĥ result, E-019, holds for generated, $B = 0$
     chains at $G = 2$ only). Each lifted capability states where it falls outside the
     singular-case guard's calibrated range, which covers generated, B = 0 chains only.
     Flexibility comes from relaxing refusals, never from a second system.
  4. *Extensible Module 11 surface.* Specify the surface's shape now: group declarations,
     parameter references, where constraints are written, and carry-forward of a generated
     level's canonical text. Reserve syntax for later capabilities but define no semantics for a
     capability before its derivation lands; until then it is refused clearly at rc 198. Each
     release exposes only validated capabilities. Designing the group surface for user ties is
     the concrete 0.3.0 action (cheap early, expensive to retrofit). **Long-term interface versus
     what 0.3.0 exposes (user, 2026-10-10):** Module 11 is designed to support user-written
     constraints through the public surface, with semantics enabled only after the relevant
     derivations. Under scope (B), 0.3.0 delivers carry-forward (a generated metric- or
     scalar-level model into a user-written β model) only if Module 11 exposes the required
     group and tie syntax **and** the β derivation and rank checks pass. If that surface is not
     ready, B is not claimed merely because canonical text can be edited: either a documented
     canonical-text workflow is explicitly supported, or 0.3.0 ships (A) and B is deferred.
     **Shape of the group surface, decided 2026-10-10 (user, option (a); built 2026-10-12):** sem's own
     group grammar through D1 (`group(grp)`, group-specific paths `(1: …)`, `@name`), plus one
     dwlsmv option naming the invariance level, **`invariance(configural|thresholds|metric|scalar)`**,
     minimum abbreviation `inv()` (decided 2026-10-10, user, option (a); declined `ginvariance()`,
     one letter from the refused sem option, and flags such as `metric`, which the canonical form
     already uses). **Partial invariance, reserved (decided 2026-10-10, user, option (a)):**
     `partial(namelist)` frees the listed parameters from the level's ties in every group (lavaan's
     `group.partial`), named in D6's group-free form: `y3:F`, `/y3:cut2`, `y3:_cons` (ν), `/scale(y3)`
     (ω). In 0.3.0 any use is refused at rc 198 ("partial invariance is not available in this
     version"); the semantics need a derivation (freeing τ or ν changes identification, I5). Declined:
     (b) group-specific members, (c) no reservation. Combining `invariance()` with user `@` ties needs
     no rule while both are refused; it is part of the tie derivation. Neither `invariance()` nor
     `partial()` is a sem option (rc 198, `tools/probe_partial_option.do`), so Module 11 removes both
     before calling sem and before D1's observed-token scan. **`group()` without `invariance()`
     (decided 2026-10-12, user, option (a); not built)** is refused at rc 198 as a version limit
     ("group() without invariance() is not available in this version. Specify
     invariance(configural), invariance(thresholds), invariance(metric) or invariance(scalar). No
     model was fitted."), reserving no meaning: unconstrained multigroup models (sem's
     `ginvariant(none)` use, e.g. structural models among factors) may be admitted later once
     their surface mapping onto G8 is specified and validated. **Fewer than two groups (decided
     2026-10-12, user, option (a) twice; not built):** `invariance()` without `group()` is
     refused at rc 198 ("invariance() needs group(); it compares groups. No model was
     fitted."), and so is `group()` with one value in the estimation sample ("group() has only one
     value in the estimation sample (grp = 0). No model was fitted."). **Group-specific paths under
     `invariance()` (decided 2026-10-12, user, option (a); not built):** Module 11 builds each
     group's section from sem's per-group `e()`, checks every group mirrors group 1 (L1), and
     refuses a difference at rc 198 before the generator, naming the path and the groups ("Path
     F -> y4 is in group 2 (grp = 1) but not in group 1 (grp = 0). With invariance(), every group
     must have the same model. No model was fitted."); the generator's rc 3498 is the backstop.
     Module 11 always runs sem with
     `ginvariant(none)`: sem's default otherwise ties loadings and intercepts silently (8 `e(Cns)`
     rows against 3). A user's `ginvariant()` is refused at rc 198, pointing to the level option,
     since sem's classes do not map onto Wu–Estabrook levels (no thresholds, no ω). A cross-group
     tie reaches `e(Cns)` as a row with more than one nonzero coefficient (`y2@a` and `(1: …@b)
     (2: …@b)` alike); it is refused at rc 198 until its derivation lands (principle 4, R-Q2).
     Parameter references are the decided names (`_b[y1:2.grp#c.F]`) in `test`, `lincom` and `e()`.
     **Help (user):** state that dwlsmv's invariance levels are not sem's `ginvariant()` classes,
     and specify how: for each level, which parameters are tied, freed or fixed (τ, λ, ν, ω, α,
     ψ), set against the nearest `ginvariant()` class and what that class lacks. Evidence:
     `tools/probe_sem_group_surface.do` (StataMP 18; Stata 16 not yet run).
  5. *Restriction-based nesting, within what the tables can verify.* Check 8 generalizes from
     row-by-row free/fixed subsets to restrictions it can verify from the tables: a parameter
     fixed, or tied, in $H_0$ that is free in $H_1$ ("tied in $H_0$, untied in $H_1$" is a nested
     pair). Nesting in general (containment of the implied σ sets, or pairs related by
     reparameterization) is a derivation gate, not a mechanical check. Wu & Estabrook's level
     chains keep the generated-chain exemption until a derived equivalent replaces it.
  6. *Separate extension gates.* Covariates and longitudinal invariance each need their own
     derivation and scope decision, and each touches a locked invariant. **Covariates** need a
     conditional Stage 1, which touches "Stage 1 takes $(Y, K)$"; use groups for the comparison of
     interest and covariates for nuisance variables, and do not cross grouping variables.
     **Longitudinal invariance is a release of its own** (version not fixed; likely 0.4.0): it
     needs residual covariances, a derivation gate against "Θ is derived, never free"; O3's
     within-group cross-row ties; and occasion-specific scale and location (ν, ω per occasion
     within one group).

The order follows feasibility, not preference: ESEM is the more wanted of the two, but it needs a
rotation module and an SE layer through the rotation, while multiple groups reuse Stage 1
unchanged.
- **Continuous/mixed indicators (v2):** extend Section 4.2's index map with variance/covariance
  `stat_type` entries for continuous variables, extend Module 2 to polyserial/Pearson estimation,
  extend Module 4's sub-phase 4c (mixed continuous/ordinal), and give Section 4.3's scale
  convention its own full treatment for the mixed case (this is the piece v1 deliberately deferred
  — see Section 3). This is what's needed the moment a continuous variable stops being a pure
  right-hand-side covariate (Module 11's note) and becomes endogenous — e.g. a mediation model
  F1 → X → F2 with X continuous: X is then a modeled outcome, not an unmodeled covariate, so it
  needs its own variance and covariances with the ordinal items in `s`, exactly like a continuous
  factor indicator would. Not supported in v1; the same MIMIC-style covariate treatment that
  handles a purely exogenous continuous `x` does not extend to a continuous mediator.
- **Robust chi-square difference testing, DIFFTEST (v2, highest priority of these).** Nested-model
  comparison under a scaled statistic. **Prioritise this above modification indices**: two scaled
  chi-squares cannot simply be differenced, and nested comparison is the backbone of measurement
  invariance testing, which is the single most common use of ordinal CFA in psychometrics. A
  package that cannot do it stops being useful to a large fraction of its users at exactly the
  point they need it most. Source: Asparouhov & Muthén (2006), Web Note 10 (Section 2).

  With $H_0$ nested in $H_1$, $d_i$ parameters in model $i$, $T = T_0 - T_1$ and $d = d_1 - d_0$:

$$P_1 = \Delta_1' W_D \Delta_1, \qquad
  V_1 = P_1^{-1}\Delta_1' W_D \Gamma W_D \Delta_1 P_1^{-1}, \qquad
  H = \partial\theta_1/\partial\theta_0$$

$$M_1 = \big(P_1 - P_1 H (H'P_1H)^{-1} H' P_1\big)\, V_1$$

  then apply Section 4.9's $T_3$ with $\operatorname{Tr}(M_1)$, $\operatorname{Tr}(M_1^2)$ and
  $D = d$. Three things make this cheaper than it looks. $M_1$ is $d_1 \times d_1$, not
  $p^* \times p^*$. It never touches $\Gamma$ directly — only through $V_1$, whose displayed
  expression is the same sandwich `dwlsmv_acov()` computes (Section 4.10's ACOV, on the
  $\sqrt N$ scale, before any sample-size divisor), so the backend already produces the
  ingredient. Whether DIFFTEST wants that quantity as is or rescaled was left open on
  2026-09-23 (the earlier wording, "ACOV scaled by $N$", was ambiguous and was withdrawn);
  **settled 2026-10-08: as is** (*The statistic, frozen*, below). When Module 12 is
  built, verify the explicit formula, scaling included, against the v2 derivation and a
  reference implementation. Mplus's DIFFTEST file stores $\Delta_1$, $P_1$ and $V_1$ from the
  less restricted model between fits; dwlsmv stores none of them (*Architecture*, below).

  **Evaluation point (decided 2026-10-08).** $\Delta_1$, $W_D$, $\Gamma$, $P_1$ and $V_1$ are
  evaluated at the unrestricted estimate $\hat\theta_1$. The restricted fit enters only through
  $H$, which is *estimated* as $\hat h = \Delta_1(\hat\theta_1)^{+}\,\Delta_0(\hat\theta_0)$
  (Moore–Penrose inverse), not taken as the parameter-selection matrix. Under $H_0$ the two
  agree asymptotically; under alternatives they do not, and only $\hat h$ reproduces lavaan
  (`lavTestLRT`, `method = "satorra.2000"`, `A.method = "delta"`) and Mplus's DIFFTEST. The
  selection form (lavaan's `A.method = "exact"`) is not used. Evidence:
  E-003, 2026-10-08; derivation @sec-difftest.

  **The statistic, frozen (decided 2026-10-08).** This settles the scaling question above:
  $V_1$ enters **as is**, the $\sqrt N$-scale ACOV `dwlsmv_acov()` computes, with no
  sample-size factor (@sec-difftest Step 2; Mplus's saved $V_1$ is the same quantity). The
  difference is $T_0 - T_1$ with each $T_j = (N-1)F_j(\hat\theta_j)$, §4.9's convention, so
  the statistic equals lavaan's and differs from Mplus's only by $N$ against $N-1$. Then
  $a_d$, $b_d$ and $T_{d3} = a_d(T_0 - T_1) + b_d$ on $d = q_1 - q_0$ degrees of freedom,
  with $\operatorname{Tr}(M_1)$, $\operatorname{Tr}(M_1^2)$ from $M_1$ built on $\hat h$.

  **Architecture (decided 2026-10-08).** One design serves both ordinary nested comparison and
  multiple-group invariance testing (Wu & Estabrook, 2016, is the reference for the latter's
  identification; the invariance levels and their expansion are a later decision). Three
  layers:

  1. *Statistic.* One internal routine: from $\Delta_1$, $\Delta_0$, $\Gamma$ and
     $T_d = T_0 - T_1$ it returns $T_{d3}$, $d$ and its p-value by $\hat h$. It knows nothing
     of labels, groups' meaning or syntax. $W_D = \operatorname{diag}(\Gamma)^{-1}$ and $V_1$
     are derived inside it, and $T_0$, $T_1$ separately are not inputs: the sign and
     $\varepsilon$ check on $T_d$ (check 10) is validation. $N$ is not an input: $\Gamma$'s
     scale cancels ($\Gamma\to c\Gamma$ sends $W_D\to W_D/c$, $P_1\to P_1/c$, $V_1\to cV_1$,
     so $M_1$ is unchanged), and $T_d$ is on the $(N-1)F$ scale (revised 2026-10-08 from $(\Delta_1, \Delta_0, \Gamma, T_0, T_1, N)$;
     E-003).
  2. *Validation.* Before the routine runs: the models share the data and the estimation
     sample, one dwlsmv version produced them, and $H_0$ is nested in $H_1$. Nesting is
     checked on the in-memory parameter tables; it never enters the statistic. The exact
     checks are the next decision.
  3. *Front end.* A post-estimation command over stored estimates, `lrtest`-style (name not
     decided): the user fits each model, `estimates store`s it, and names two or more stored
     results. The command **refits** every named model from its stored specification over
     **one shared Stage 1**, so a shared Stage 1 holds by construction and no Stage-1
     fingerprint is stored or compared. The refit reproduces the stored fit exactly
     (estimation is deterministic), and that equality is asserted. An in-call form (for
     example a chain of invariance levels in one `dwlsmv` call) may come later and feeds the
     same routine.
     **Built 2026-10-23, single-group comparisons only (user, 2026-10-23; design
     E-020, U1–U8).** `difftest name1 [name2 …]`, its own
     `difftest.ado`, `lrtest`-style: `.` the active results, one name compared with the active
     results, three or more a chain ordered by check 9 and tested as adjacent steps; the active
     results are restored afterwards. Only `e(cmd) == "dwlsmv"` results. **Multiple-group fits
     admitted 2026-10-31 (user), under the narrow claim:** check 6 at $G>1$ is certified for
     refits in the same Stata version, edition and processor count only
     (E-021: StataMP 18, counts 1, 2, 4, 6, every separate-process
     refit bitwise identical); its fault detection is supported on StataMP 18 at 6 processors on
     the 0.3.0 estimator (E-022, 2026-11-03: 0 of 3,639 material
     deciding rows accepted), and at 16.1 rests on half L (E-023). Refits
     across processor counts are not certified: that study's cross-count part was INCONCLUSIVE
     under its frozen rule, and the claim stays narrow (user). The help asks for the same
     processor count at $G>1$. Check 4 also requires the
     same `e(groupvar)` (both empty or equal); check 2's data signature covers the group
     variable. At $G>1$ the shared Stage 1 is one per group (`dwlsmv_m11_gsupport`,
     `dwlsmv_m11_stage1_g`), each model is refitted on `dwls_identify(e(canon))` with
     `dwls_fit` and `dwls_inference` (the invariance() flow's path; recovery is "none"), and
     the refit's `e(b)`, `e(V)` come from `dwlsmv_m11_gbv`, the code that posts them; check 8
     receives each model's `e(canon_configural)`, `e(canon)`, `e(invariance)` and `e(K)` for
     the generated-chain exemption; $c=N-G$, $\Omega$ and the fit indices over all groups.
     `difftest.ado` runs checks 1–5 on the stored results;
     `dwlsmv_dtfront` runs one shared Stage 1 on the common estimation sample, refits each model
     through `dwlsmv_m11_estimate` (the code `dwlsmv` itself runs, recovery included) with full
     inference, then check 7, check 6 **only through `dwlsmv_dt_check6`** (today:
     `dwlsmv_dtrefit` on `e(b)`, `e(T3)`, `e(bound)`, `e(df_T3)` at `dwlsmv_dt_reftol()` = 1e-10;
     `e(V)` carried, not compared, so check 6's open contract can change there alone), checks
     8–9 (`dwlsmv_dtchain`), and per adjacent pair checks 10–11 and the statistic
     (`dwlsmv_dtpair`, $c=N-G$, $N-1$ at $G=1$, nrtol $10^{-12}$). Output (U5 (b)): a models table (npar, $T_3$,
     df, RMSEA, CFI, TLI from the refits) and a steps table (d, $T_{d3}$, p, ΔCFI, ΔRMSEA from
     `dwlsmv_fitidx_delta`, each the restricted model's value minus the less restricted one's); a step within check 10's band prints "within eps; not tested" with
     no $T_{d3}$ or p (U6 (a)); a $d=0$ step "skipped (d = 0)"; the routines' notes verbatim; the
     approved cutoff sentence. `r(models)` (rows by model, most restricted first: npar, T3, df,
     rmsea, cfi, tli), `r(steps)` (rows H0:H1: d, T_d3, p, dCFI, dRMSEA, status), `r(N)`
     (U7 (a)). Tests: `test_difftest_front.do`. Help: `difftest.sthlp`, approved by the user 2026-10-23.

  **What `e()` must carry for the refit (a D8d change, made in the commit that builds it).**
  No statistical payload: $\Delta$, $\operatorname{diag}\Gamma$, $\hat s$ and `T_raw` stay in
  Mata. Only what reconstructs the fit: the canonical specification as `dwlsmv_m11_fit()`
  receives it (the single-line, `;`-separated text, as an `e()` local; name not decided),
  beside `e(cmdline)`, which stays as provenance; `keepbound`, which is a Stage-1 option and is
  not posted today; and whatever else `dwlsmv_m11_fit()` takes that `e(indicators)` and the
  canonical text do not already determine (to be established when it is built). The same-data
  check uses Stata's `signestimationsample` / `checkestimationsample`, together with identical
  `e(sample)`. Both run on the Stata 16 floor (verified 2026-10-08 on StataMP 16.1 and 18:
  unchanged data rc 0, one changed value rc 459, after `estimates restore` rc 0, after
  re-sorting rc 0; the same signature on both; E-024).
  **Built 2026-10-23 for single-group fits:** `e(canon)` (the single-line, `;`-separated text
  `dwlsmv_m11_fit()` receives), `e(keepbound)` (`keepbound` or empty) and the data signature
  (`signestimationsample` over `e(indicators)`: `e(datasignature)`, `e(datasignaturevars)`).
  Nothing else was needed: `e(metric)` selects the recovery, `e(indicators)` the Stage 1
  variables and order. Multiple-group fits already post `e(canon)` (2026-10-12) and post no
  data signature, since `difftest` refuses them.

  **Validation contract (decided 2026-10-08).** The command never silently repairs a
  mismatch, substitutes a different model, or reports a statistic whose reference
  distribution it knows to be inappropriate. Every refusal is rc 3498 with a message naming
  the check and the models.

  *A. One Stage 1 can serve every model.*
  1. Every named model was fitted by `dwlsmv`, by the installed version (`e(version)`). A
     version mismatch is refused here, never absorbed by check 6's tolerance.
     **No Stata-version provenance in `e()` (decided 2026-10-23, user):** no Stata version,
     edition or OS is posted; check 1 compares `e(version)` (the dwlsmv release) only.
     `difftest` refits every model in the current session over one shared Stage 1, so the
     statistic is always computed in one Stata version; a stored fit from another version or
     edition can only make check 6 refuse (visible) or pass within its tolerance. StataCorp
     (2026-10-23, E-025): results can differ across Stata
     versions and MP/SE editions and across processor counts. Check 6's refusal adds one line:
     a model fitted in a different version or edition of Stata should be refitted in this one.
     The help asks users to fit the models and run `difftest` in the same Stata version and
     edition, not with the same processor count.
  2. Same data: each fit posts a data signature (`signestimationsample`, a D8d addition);
     `checkestimationsample` passes for every model.
  3. Identical `e(sample)` across the models.
  4. Same indicator set, same `keepbound`, same `metric`.
  5. Same indicator **order**. A difference is refused, with a message telling the user to
     write the paths in the same order; nothing is reordered.

  *B. The refit reproduces the stored fit.*
  6. `e(b)`, `e(T3)` and `e(bound)` of the refit agree with the stored values within a
     relative tolerance set from the measured processor-count residue (class (k):
     θ̂ by one ulp, $T_3$ by about 1,500 ulps, $2.6\times10^{-13}$ relative;
     E-011 2026-10-02). Not exact, and not a structural comparison only. The
     value is fixed when the check is built, citing that measurement, and labelled as chosen
     after measuring. The measure (user, 2026-10-08) is $\max_i |x_i - y_i|/(1 + |y_i|)$, $y$
     the stored value, $x$ the refit's, for `e(b)` and `e(T3)`. Also decided (user,
     2026-10-08): `e(bound)`, a discrete 0/1 status, must be **identical**, never loosened
     by the tolerance; `e(df_T3)` must be equal (the refit is the same model); a missing
     `e(T3)` in either result is refused, since `dwlsmv` posts one for every fit, a
     saturated one included ($T_3 = T_{\text{raw}}$, Section 4.9), which the measure
     compares near 0 as an absolute difference. Whether a DIFFTEST against a saturated
     model is meaningful is the front end's question, not this check's. Internal:
     `dwlsmv_dtrefit`; the tolerance is `dwlsmv_dt_reftol` = **1e-10, chosen after measuring**
     (2026-10-09): 100 × the observed maximum 8.05e-13 across processor counts 1, 2, 4, 6, rounded up to a power of ten, a rule frozen before the run; same Stata version only, G = 1 only (E-026). **At G > 1 (decided 2026-10-14, user, option (c)):** a frozen threshold gate, not a recalibration. The same rule is applied to G > 1 fits (E-028). If its value is ≤ 1e-10, `dwlsmv_dt_reftol` stays 1e-10 for every G and its signature is unchanged. If its value exceeds 1e-10, or a precondition fails (a repeat at one processor count not bitwise identical, a fit not converged, `e(df_T3)` not identical across counts), the outcome is **stop and return**: the gate stays not complete, `dwlsmv_dt_reftol` is unchanged, and the result goes back to the user with a residue diagnosis before any change to the implementation or the preregistration. Declined: (a) one value, the larger of the two; (b) a value per G. **Executed 2026-10-14: STOP AND RETURN** (V = 1e-9; E-028). The tolerance stays 1e-10 and the gate is not complete. Next (user, 2026-10-14): a sensitivity derivation, then a preregistered diagnostic whose thresholds come from it. Neither outcome changes the tolerance; that takes a separate decision. The derivation (`@sec-mg-sensitivity`, Results S1–S4) landed 2026-10-14. **Amended the same day (user):**
- **(D) is added.** It is derivative-computation error: production Δ's correlation block is a central difference (`dwlsmv_delta_cols`), so the computed gradient is $\hat g(\theta)=g(\theta)+e(\theta)$ and the optimizer follows a perturbed Gauss–Newton iteration. (D) is a third mechanism beside (A) and (C). It may interact with (C), and the attribution does not force that interaction into either category without evidence. A deterministic $e$ cannot alone explain a cross-count difference under identical arithmetic, so how $e$ couples to changes in inputs or iteration path must be established.
- **A controlled diagnostic runs before the 4 × 4 grid,** on fixture E at the metric level (consistency data):
  1. the gradient noise, at a saved θ, against a higher-accuracy reference derivative;
  2. the stopping floor separated from cross-count variation, with Stage 1 and the starts held fixed, production Δ against Δ frozen;
  3. an intervention: a diagnostic-only, independently checked, more accurate derivative;
  4. only then, the grid.

  The accurate derivative is a diagnostic control, not a production change. Steps 2–3 run only under a
  separate protocol frozen by the user (2026-10-14).
  7. Both fits converged, before and after the refit. D11-recovered fits are accepted; the
     refit reproduces the recovery.

  *C. Nesting, on the in-memory parameter tables.*
  8. Rows are identified by (type, row, column); a row absent from a table (D5 omits
     covariances fixed at 0) counts as fixed at its type's absent-row default (§4.8 G4, decided
     2026-10-09: ω 1, every other type 0; at $G = 1$ every built type reads 0, as before). Every parameter free in $H_0$ is free in
     $H_1$; every parameter fixed in both is fixed at the same value; $d = q_1 - q_0 \ge 1$,
     and $d$ equals the difference in `e(df_T3)`. Pairs nested only through a
     reparameterization (a different marker, for example) are refused; the help says so.
     **The models declare the same factors and indicators in the same order** (user,
     2026-10-08): a row present in both tables carries the same name in each, names holding
     the factor and indicator names; otherwise the pair is refused, since rows matched by
     index would be different parameters. Nothing is reordered. Internal: `dwlsmv_dtnest`.
     **At $G > 1$ (decided 2026-10-12, user; built 2026-10-12, `test_dtnest_groups.do`):** rows are keyed by (type, group,
     row, column) (§4.8 G6) and absent rows read by G4. **A pair in which either table has a tie
     (two rows sharing one `fidx`) is refused, rc 3498, fail-closed:** the free/fixed rules
     above do not read ties, so without this a pair tying loadings in $H_0$ and thresholds in
     $H_1$ would be accepted though not nested. Counts ($d = q_1 - q_0$) are necessary
     evidence, not proof. A tie rule (every tie of $H_1$ holds in $H_0$), its rank condition
     and an executed validation are a deferred gate, after the help file and Module 11; the
     gate's proposed focus (user) is the level generator and the generated-chain exemption below. The β guard
     is unchanged and runs first.
     **Generated invariance chains are exempt (decided 2026-10-08).** Wu & Estabrook's levels
     are each identified on their own: a level adds equality constraints *and* frees
     identification constraints, so its parameter table is not a subset of the previous
     level's, and this check would refuse every valid step. A pair is exempt only when both
     fits were produced by dwlsmv's level generator and the comparison command, regenerating
     each level's canonical specification from the posted configural specification and level
     name, matches the stored canonical text exactly. Nesting then holds by construction,
     established in the derivation and tests. Anything else, including a hand-edited
     specification that names a level, gets this check unchanged.
     **Limit (decided 2026-11-02, user):** the derivation does not cover a configural
     specification with a loading or factor covariance fixed at a nonzero value. A fixed
     row cannot absorb the rescaling that nests the levels, so such a pair is not exempt and
     reaches the tie refusal, with a note. This over-refuses pairs that are in fact nested
     until a derivation covers them. The Result I2 equivalence note is withheld likewise.
     The surface already refuses a fixed nonzero loading under `invariance()`, so from
     the surface only the fixed covariance reaches this limit.
     **Non-adjacent generated pairs (decided 2026-10-13, user, option (a); built 2026-10-13, `test_dtgenerated.do` g2 and c1):** the
     exemption may admit any pair of generated levels from the same canonical base (e.g.
     configural against scalar), subject to its four conditions, check 9's ordering and the
     established domain requirements (`@sec-inv-nesting`: nesting holds where $\phi_g>0$; check
     10's diagnostic otherwise). The comparison is a joint test of the endpoint models; it is
     not decomposed into adjacent steps. **Preconditions before any implementation relies on
     it:** `inv-nesting` gains a multi-factor assertion; the $G>1$ statistic's derivation
     establishes how the group proportions $\pi_g$ enter the routine; the tier-2 DIFFTEST
     fixture covers design D (a zero-df step). Not settled by this decision: the reporting
     policy for a comparison whose total df difference is zero (settled separately, below).
     Condition (4) (check 6's refit) is check 6's, run by the front end before this check
     (§4.8 G8, decided 2026-10-13, option (ii)); `dwlsmv_dtgenerated` checks (1)–(3).
     **A failed generated claim (decided 2026-10-13, user, option (a); built 2026-10-13, `test_dtgenerated.do`):** when both
     models post a non-empty `e(invariance)` and one of (1)–(3) fails, the exemption prints one
     line saying the pair is not treated as a generated chain and naming the first failed
     condition; check 8 then runs unchanged (at $G>1$, the tie refusal, rc 3498). No line when
     only one model posts a level. A configural text the generator rejects keeps the generator's
     own rc 3498. Tests: each failed condition asserts its line and the tie message; a
     single-level pair prints none. Mutants: line removed; wrong condition named; line printed
     for single-level pairs. Declined: (b) refuse inside the exemption; (c) no line.
     Declined: (b) adjacent levels only; (c) non-adjacent only across $d=0$ steps. All three
     preconditions met 2026-10-13 (`inv-nesting-mf`, `mg-difftest-omega`,
     `99e_lavaan_mg_difftest.R`).
     **Zero-df generated pairs (decided 2026-10-13, user, option (c); built 2026-10-13, `dwlsmv_dtstep`, `test_dtstep.do`).** A pair
     exempt under the conditions above with $d = q_1 - q_0 = 0$ (adjacent or not) has no test:
     the statistic routine is not called, and the pair yields status "no test" with $d = 0$
     and $T_d$ (reported as 0 when $|T_d|\le\varepsilon$, check 10's band), no $T_{d3}$ and no
     p-value, carrying the decided wording: the established-equivalence text where Result I2
     establishes it (the metric step with every factor having at most one indicator of three
     or more categories, all its loadings free; the `_i2` detector), the equal-df text
     elsewhere. In a chain the step is reported as skipped and the remaining adjacent pairs
     are tested as usual. **Consistency check where equivalence is established:** for an I2
     pair, $|T_d|>\varepsilon$ aborts, rc 3498, naming both fits (equivalent models must fit
     identically to the fits' precision, so one has missed its minimum). For any other $d=0$
     pair, $|T_d|>\varepsilon$ is reported as information (same df, the fits differ, so
     equivalence could not be established; amended 2026-10-13, user, from "the models are not
     equivalent"), never refused.
     **Texts (approved 2026-10-13, user; written for an applied reader, ε stated as the precision
     the fits are computed to; T0 − T1 is called the difference in fit, never a chi-squared
     difference):** prefix of every $d=0$ note "No test for *H0* against *H1*: the two models have
     the same degrees of freedom (d = 0).", then the decided note text (the equal-df text naming
     H1's level for L−1). Information line: "The two models have the same degrees of freedom but
     do not fit the data equally well (difference in fit #, larger than the # to which the fits
     are computed), so equivalence of the two models could not be established. With no
     difference in degrees of freedom there is nothing to test." I2 abort: "dwlsmv check 8 (*H0*,
     *H1*): for this model, invariance(metric) is mathematically the same model as
     invariance(thresholds) (every factor has at most one indicator with three or more
     categories), so both fits should agree to the precision they are computed to. They differ by
     # (precision #). The likeliest cause is that one of the two fits stopped short of its best
     solution, though another numerical problem cannot be ruled out. Refit both models; no
     comparison is reported." **Who prints (decided 2026-10-13, user, option (ii)):** aborts are
     printed by the routine before rc 3498; the notes are returned by `dwlsmv_dtstep` and the
     front end prints them verbatim, with a front-end test asserting each appears. Declined: (i)
     the routine prints the notes (they would precede the front end's table). No new tolerance: $\varepsilon$ is check 10's.
     Fixture evidence (one dataset each, a measurement, not a result): B's thresholds step
     $T_d=6.0\times10^{-11}$, D's metric step $-3.6\times10^{-10}$, both within
     $\varepsilon\approx4\times10^{-9}$. **Tests (at the build):** an I2 pair within the band
     (no test, verified wording); the same pair handed $|T_d|>\varepsilon$ (rc 3498); an
     unverified $d=0$ pair inside and outside the band (no test, equal-df wording; the outside
     case carries the information line, rc 0); a chain with a $d=0$ middle step (skipped, the
     other pairs tested). Mutants: the routine called at $d=0$; the I2 abort removed; the
     abort applied to unverified pairs; the band replaced by exact zero. Declined: (a) refuse
     every $d=0$ pair; (b) report without the consistency check; (d) refuse unverified $d=0$
     pairs. **Non-adjacent $d=0$ pairs (decided 2026-10-13, user, option (a); built 2026-10-13, `test_dtstep.do` z6):** the
     verified wording and the I2 abort apply only to the thresholds–metric pair with `_i2`
     holding; every other $d=0$ pair, adjacent or not (configural–metric included), gets the
     equal-df wording and the information line, never the abort. Composition of verified steps
     is not stated in `@sec-invariance-id`; extending the verified set waits for the derivation
     that first verifies a second pattern. Test: a configural–metric $d=0$ pair, equal-df
     wording, outside the band rc 0 with the information line. Mutant: verified wording for any
     pair containing the metric step. Declined: (b) any pair whose every step is verified; (c)
     any pair containing the metric step.
  9. The models are ordered by their free-parameter count, so the user may name them in any
     order. More than two must form a chain under check 8; adjacent pairs are tested.
     Internal: `dwlsmv_dtchain`. **Equal free-parameter counts (decided 2026-10-13, user, option
     (a); built 2026-10-13, `test_dtgenerated.do` z1):** when both tied models are generated (non-empty `e(invariance)`), the tie is
     broken by level rank, most restricted first (scalar, metric, thresholds, configural), so a
     $d=0$ generated step leaves the chain, and which pairs are tested, independent of the order
     the models are named in; any other tie keeps input order and check 8 refuses it as before.
     Tests: a chain with a $d=0$ step named in every order gives one ordering; a mixed tie keeps
     input order and is refused; mutant: tie-break dropped. Declined: (b) input order, documented;
     (c) refuse ties not given in level order. **The exemption is built (2026-10-13, gate step 2b,
     `test_dtgenerated.do`):** `dwlsmv_dtnest` and `dwlsmv_dtchain` take each model's
     `e(canon_configural)`, `e(canon)`, `e(invariance)` and `e(K)` as optional inputs, call
     `dwlsmv_dtgenerated` after the β guard and before the tie refusal, and `dwlsmv_dtchain`
     returns each adjacent pair's $d$ and exempt flag; without those inputs every pair gets check 8
     as written. The per-pair reporting (zero-df pairs, check 10's diagnostic, check 11's callers)
     is `dwlsmv_dtstep` (built 2026-10-13, step 2c).

  *D. At the statistic.*
  10. $T_0 - T_1 < -\varepsilon$ aborts: a premise (a fit at its global minimum) has
      failed, and the message names the suspect fit. $|T_0 - T_1| \le \varepsilon$ is a
      restriction that does not bind to the fits' precision: reported as 0 with a note
      saying so. That is a statement of precision, not truncation. **$\varepsilon =
      c\,\mathrm{nrtol}$, $c$ the multiplier of $F$ in $T$** (decided 2026-10-09; it was
      $(N-1)\,\mathrm{nrtol}$ from 2026-10-08): $T = (N-G)F$, so $c = N - G$, which is $N - 1$
      for one group. The pooled band requires the minimizer to be handed weighting (a),
      $w = \pi_g w_g^0$ (@sec-multigroup-scaling MG7); `dwls_fit()` minimizes under it (decided
      2026-10-10, user). The band comes from the fits' convergence
      criterion, never from results: @sec-difftest derives it as
      $\tfrac12\kappa\,c\,\mathrm{nrtol}$, $\kappa = \lambda_{\max}(H_t^{-1}H_{GN})$, from
      `dwlsmv_minimize()`'s Newton-decrement stop (`nrtol` = $10^{-12}$ in $F$), and the
      allowance $\kappa \le 2$ is **chosen after measuring** $\kappa$ = 1.000–1.037 on four
      Case E fits. The test that builds this check asserts the decrement at the returned
      $\hat\theta$ is at most `nrtol`. Internal: `dwlsmv_dtband`. What "reported as 0" means
      for $T_{d3}$ and its p-value is open (the front end's).
      **Generated pairs with a negative $H_0$ factor variance (decided 2026-10-13, user,
      option (c) with five refinements; built 2026-10-13 in `dwlsmv_dtband`, `test_dtstep.do`).** A generated pair is nested only where
      $\phi_g>0$ (`@sec-inv-nesting`; @sec-invariance-id Result I6): at $\hat\phi_g<0$ $H_0$'s
      fitted $\sigma$ lies outside $H_1$'s image and $T_0\ge T_1$ is not guaranteed, though the
      asymptotic reference is unaffected. Estimation, Heywood treatment, the exemption's
      conditions and check 11 are unchanged. Added: when $T_0-T_1<-\varepsilon$ **and** the pair
      passed check 8 through the generated-chain exemption **and** some free later-group
      $\hat\psi_{gff}$ of $H_0$ is $<0$ ($B=0$, so $\psi=\phi$; only $H_0$ at the metric or
      scalar level can trigger it, adjacent or not), the abort (rc 3498) carries a distinct
      message: the nesting premise failed at $H_0$'s estimate, naming group, factor and value,
      and saying it is not a convergence finding; it does not name $H_1$ as suspect (wording
      decided at the build). A user-written pair keeps the existing message (its nesting is
      global). $\hat\phi_g=0$ is left to the existing level rules (every 0.3.0 map is defined
      there; the all-binary scalar map, which divides by $\sqrt\phi$, is refused). Whether the
      branch lives in `dwlsmv_dtband` or its caller is decided at the build. **Tests:**
      constructed inputs (an exempt pair with $\hat\phi<0$: new message; the same values
      user-written: old message; $\hat\phi=0$: old message; $\hat\phi<0$ with $T_0\ge T_1$: no
      abort), each rc 3498 where refused; mutants: branch removed, branch unscoped, `<=` for
      `<`. The trigger condition is derived; its frequency is not measured. **Message (approved
      2026-10-13, user):** "dwlsmv check 10: *H0* fits the data better than *H1* (by #, beyond the
      # to which the fits are computed), which a more restricted model should not do. *H0*'s
      estimates include a factor variance below 0: psi[F,F] in group *g* is estimated at #. The
      invariance levels are nested only when every factor variance is positive, so at this
      estimate the comparison is not valid, whatever else may contribute to the difference. This
      is a property of *H0*'s solution, not a sign that the fit failed to converge. No statistic is
      reported." (one line per negative variance). No next-step advice here (user): practical
      Heywood-case guidance is a later, separately reviewed help-file note, not part of check 10.
  11. Restrictions on the boundary of the parameter space are refused, since the $\chi^2$
      reference does not hold there (restated and decided 2026-10-08). The test is whether
      $H_0$'s **fixed values force** the boundary whatever the free parameters are, never
      where $\hat\theta_0$ happens to land (an estimate near the boundary is not a boundary
      restriction). Two conditions, exact, with no cutoff:
      (a) some principal submatrix of $\Psi$ whose entries are all fixed in $H_0$ is
      singular (a variance fixed at 0; two variances and their covariance fixed so that the
      correlation is $\pm1$) **or not positive semidefinite** (a covariance fixed at 1.2
      between unit variances: outside the parameter space, refused too; user, 2026-10-08).
      Equivalently, some *maximal* all-fixed block is not positive definite: a PD block has
      only PD principal submatrices, and every all-fixed block lies in a maximal one. Absent
      rows count as fixed at their type's absent-row default (§4.8 G4: ω 1, every other type
      0; decided 2026-10-09), so an absent ω enters rule (b)'s product as 1. $\Phi = (I-B)^{-1}\Psi(I-B)^{-\top}$ with $I-B$ invertible, so
      $\Phi$ is singular exactly when $\Psi$ is and needs no separate check.
      (b) some $\Theta_{jj} = 1 - (\Lambda\Phi\Lambda')_{jj}$ is a function of fixed values
      only and is $\le 0$ (for example a unit-variance factor with a loading fixed at 1).
      "A function of fixed values only" is structural: every loading in row $j$, every
      $\beta$ into the factors that row reaches (with their ancestors), and every $\psi$
      among those factors is fixed. A boundary forced *through* free parameters is not
      refused. "Exact" is literal (user, 2026-10-08): a strict pivot $> 0$ and
      $\Theta_{jj} > 0$ in floating point, with no band; values not exact in binary can
      decide an exactly singular case by one ulp. The table checked is $H_0$'s, and a boundary
      $H_1$ fixes too is refused as well, as the rule reads. Internal: `dwlsmv_dtboundary`.
      The "two factors are one" test (correlation fixed at 1) is expressible only in the
      unit-variance metric, where (a) refuses it; v1a's grammar cannot fix a marker-metric
      correlation, which needs an equality constraint. The help states the refusal and why.
      **Generated chains: a configural preflight (0.3.0; user, 2026-10-08).** When the
      level generator is invoked, the configural table it is given is checked under this
      rule, group by group, before the configural fit and before any level is generated, and
      a failure is refused there, naming the configural specification. This is a workflow
      extension of the contract, not a change to the rule: the same routine, the same
      structural test, no estimate read. Each generated level's $H_0$ is still checked as
      above, as a backstop; under derivation `@sec-inv-check11` Result C11-4 it fires only if
      the generator itself introduces a forcing value, and its message says so. A
      user-written pair is unchanged: only its $H_0$ is checked, and a user-written $H_1$
      that forces the boundary is not refused by this preflight. That asymmetry is a 0.3.0
      scope choice, not a judgment of admissibility; such an $H_1$ goes to the fit, where
      Modules 7 and 10 act as for any model. Extending the preflight needs its own decision.
      **Call site (decided 2026-10-13, user, P1 (i); built 2026-10-13 with the $G>1$ rule
      below, `test_dtboundary_groups.do`, mutants B1–B8 and P red):** Module 11's `invariance()`
      flow, at the top of `dwlsmv_m11_fit_g`, on `dwls_identify` of the generator's configural
      output, before the level is generated and before any fit; the refusal is check 11's rc
      3498 naming the configural specification. The generated-level backstop's own wording is
      built with the exemption (step 2; P2); until then a level's refusal carries the plain
      message, naming the group at $G>1$. **At $G>1$** (`@sec-inv-check11` C11-1–C11-3): rule (a)
      runs on each $\Psi_g$; rule (b) on $\Theta_{gjj}=1-\omega_{gj}^2(\Lambda_g\Phi_g\Lambda_g')_{jj}$,
      a free $\omega_{gj}$ leaving the row undetermined, an absent one reading 1 (G4).
      **Messages (decided 2026-10-13, user; not built):** the routine states the condition, the
      caller the consequence (mechanism (a): the caller passes its sentence, the routine prints it
      and aborts rc 3498 where it does now; declined (b), a returned status, since a forgetful
      caller fails open). Condition text approved: rule (a) "dwlsmv check 11 (*model*, group *g*):
      the fixed values in Psi force a boundary. Every entry of this block is fixed and the block is
      singular or not positive semidefinite:" then the block's fixed entries, then "(a covariance
      with no directive is fixed at 0; a factor correlation fixed at 1, two factors made one, is
      such a block)."; rule (b) "dwlsmv check 11 (*model*, group *g*): the fixed values force the
      residual variance of *item* to 1 − ω²(ΛΦΛ′)ⱼⱼ = *value*, which is zero or negative (ω =
      *value*). Its loadings:" then the loading rows, then "(Φ from the fixed β and ψ values; ω
      reads 1 when the table has no ω row for the item.)". The item is named from its loading rows'
      names, never an index; "group *g*" only at $G>1$. **Caller sentences (approved 2026-10-13,
      user):** preflight "The configural specification is refused before fitting: no invariance
      level can be tested from it."; user-written comparison "This restriction lies on the
      boundary of the parameter space, where the chi-squared reference of the difference test
      does not hold, so the comparison is refused."; generated-level backstop (P2) "The level
      generator introduced this fixed value: the configural specification passed this check, so
      the forcing value comes from the generated level, not from the model you wrote. This is an
      internal error; the comparison is refused." An empty consequence is refused, rc 3498.

  **Caveat on the source, verified numerically.** Web Note 10's eq. (9) gives the equivalent
  large form $M = W^{-1}\Delta_1'(P_1^{-1} - H(H'PH)^{-1}H')\Delta_1W^{-1}\Gamma$, which **does
  not conform dimensionally as printed** — $W^{-1}$ is $p^*\times p^*$ and $\Delta_1'$ is
  $d_1 \times p^*$. Conformability forces $\Delta_1 \cdots \Delta_1'$. Read that way, the note's
  claimed identity $\operatorname{Tr}(M_1)=\operatorname{Tr}(M)$ and
  $\operatorname{Tr}(M_1^2)=\operatorname{Tr}(M^2)$ holds to machine precision on synthetic
  matrices, and $\operatorname{rank}(M_1) = d_1 - d_0$ as it must. Use $M_1$ and the problem
  does not arise. Remember also that the paper's $W$ is the inverse of this spec's $W_D$
  (Section 4.9).

  Note Web Note 10 predates the 2010 note and applies the *Satterthwaite* adjustment (its eq. 10).
  For consistency with a $T_3$ v1, apply the scale-and-shift instead, as the 2010 note's Section 4
  directs.

- **Modification indices (v2).** For a scaled statistic the modification index is the ordinary
  Sörbom (1989) index **multiplied by $a$** — the same $a$ Module 9 already computes. The
  justification, from the 2010 note's Section 5, is that $a$ and $b$ are functions of averages of
  eigenvalues and freeing one parameter essentially removes one eigenvalue, which barely moves an
  average.

  The cost is not in the scaling but in the ordinary index: Sörbom's method needs first and second
  derivatives of the fit function with respect to parameters that are **fixed** in the current
  model. That is the requirement the Module 3 and Module 6 design notes exist to keep reachable.

  Two caveats to carry into the command's documentation, both from the note itself. The constant
  $a$/$b$ assumption "may not be appropriate, such as when the degrees of freedom is small." And
  more bluntly: modification indices, *even unscaled*, do not provide a precise value nor bound
  the approximation error, so their use "should be limited to an exploratory phase in the model
  building rather than used as evidence for model fit." If this ships, that belongs in the help
  text, not buried in a methods appendix.

- **Multiple-group models (v2; 0.3.0 grammar specified in Section 4.8, not built).** What
  follows is the original framing; the design now lives in Section 4.8 (**Multiple groups**). Stage 1 runs unchanged on each group's own $(Y_g, K)$, so $s$ stacks the groups'
  statistics and $\hat\Gamma$ is block diagonal across independent groups. How each block is
  scaled by its group's share of $N$ in the pooled objective, ACOV and $T_3$ is derived in
  `@sec-multigroup-scaling` (2026-10-09). Decided the same day after the check 10 review:
  $T = (N-G)F$, `e(V)` $= \mathrm{ACOV}/(N-G)$, and check 10's $\varepsilon = c\,\mathrm{nrtol}$
  with $c = N - G$. Cross-group equality constraints are the canonical form's `tie` directive, a Module
  3/Module 11 contract change specified in Section 4.8 (**Multiple groups**) and not yet built. Identification is the hard part. Which thresholds, loadings and scale
  parameters may differ across groups under the delta parameterization, at each level of ordinal
  measurement invariance, must be derived and checked against the literature and against lavaan's
  `group.equal` behaviour before any code. This extension and DIFFTEST are what make invariance
  testing usable, so they should ship together.
- **ESEM (v3):** an unrotated starting-value CFA (all indicators loading on all factors, reusing
  Module 3's free/fixed table with cross-loadings freed), a rotation module (geomin default,
  target rotation per Asparouhov & Muthén 2009), and a delta-method SE propagation layer (Jennrich's
  approach) reusing the same Δ/Γ machinery from Modules 4 and 8.

### 5.1 API shape — **locked (Gate 6)**

The signatures below are binding. Types live in `src/mata/dwlsmv_types.mata`;
each function's body lives in its own file. Mata passes arguments by reference,
so routines that fill a struct take it as an argument and return `void`.

`dwls_polychoric()`'s public signature is unchanged by Module 2's root-finder (2026-09-30).
`struct dwlsmv_pairfit` carries `real scalar status`, the pair's outcome coded as in Module 2's
status table; `converged` is kept for compatibility and set from it (`converged = (status == 1)`),
never maintained independently. `iterations` counts evaluations, the endpoint check's included.

`struct dwlsmv_stage1` carries `real colvector bound` (decided 2026-10-03): $P\times1$, entry $k$
for pair $k$ in the map's order (the statistic at `s[T + k]`), 1 if the pair is a status-2 pair
included at $\pm B$ under `keep_bound`, else 0. It is all zeros by default, because without
`keep_bound` such a pair aborts Stage 1. It is stored rather than derived because $|\hat\rho| = B$
does not identify it: a status-1 root may be accepted at $\pm B$ itself. Module 11 maps the marked
positions to variable names for `e()` and the output.

```text
// ---- the index map (Gate 1) --------------------------------------------
struct dwlsmv_statmap scalar
    dwlsmv_statmap_build(real colvector K)

real scalar    dwlsmv_pos_tau(struct dwlsmv_statmap scalar map,
                              real scalar j, real scalar c)
real scalar    dwlsmv_pos_rho(struct dwlsmv_statmap scalar map,
                              real scalar i, real scalar j)
real rowvector dwlsmv_pos_inv(struct dwlsmv_statmap scalar map, real scalar s)
                              // returns (stat_type, var_i, var_j, thr_k)

// ---- Module 2: one pair's polychoric and its information (Gate 3) ------
struct dwlsmv_pairfit scalar
    dwls_polychoric(real colvector yi,    real colvector yj,
                    real colvector tau_i, real colvector tau_j)

// ---- Modules 1-4: Stage 1 as a whole -----------------------------------
struct dwlsmv_stage1 scalar
    dwls_sample_stats(real matrix Y, real colvector K,
                      | real scalar keep_bound)  // optional, default 0

real matrix
    dwls_gamma(real matrix Y,
               struct dwlsmv_statmap  scalar map,
               struct dwlsmv_pairfit  vector pairs,
               real matrix            psi_tau)   // the STACKED per-variable
                                                 // lookup tables, sum(K) x T
                                                 // -- not the gathered N x T

// ---- Module 3: identification and the parameter map --------------------
struct dwlsmv_parmap scalar
    dwls_identify(string scalar modelspec,
                  struct dwlsmv_statmap vector maps,   // one per group (O14)
                  | struct dwlsmv_stage1 vector st)    // required with 'groups':
                                                       // threshold starts (O5)

// ---- Module 5: model-implied moments (Gate 2) --------------------------
real colvector
    dwls_model(real colvector theta,
               struct dwlsmv_parmap  scalar par,
               struct dwlsmv_statmap vector maps)   // stacked sigma at G > 1

// ---- Module 6: the Jacobian (Gate 2) -----------------------------------
real matrix
    dwls_delta(real colvector theta,
               struct dwlsmv_parmap  scalar par,
               struct dwlsmv_statmap vector maps)   // stacked Delta at G > 1

// ---- Module 7: fitting --------------------------------------------------
struct dwlsmv_fit scalar
    dwls_fit(struct dwlsmv_stage1 vector sts,    // one per group (O14)
             struct dwlsmv_parmap scalar par)

// ---- Modules 8-9: inference; fills acov, T3, shift_a, shift_b ----------
void
    dwls_inference(struct dwlsmv_fit scalar f)
```

**What the lock actually constrains**, beyond the names:

*Stage 1 takes `(Y, K)` and nothing else* — except the optional scalar
`keep_bound` (decided 2026-10-03), which chooses how Stage 1 treats status-2
pairs other than empty-cell 2×2 pairs: 0 (default) aborts, 1 includes them at $\pm B$ and marks
them in `stage1.bound`. It is an estimator setting, not model input, and with it
set Stage 1 reproduces lavaan's own treatment of those pairs, so the paired
validation design is unaffected. No model specification, no starting
values, no `theta`. Any signature change that lets a structural model reach
`dwls_sample_stats()` or `dwls_gamma()` breaks the paired validation design of
Section 7.1, which depends on Stage 1 being checkable against lavaan with no
structural model specified at all.

*`dwls_gamma()` receives the pair fits, it does not compute them.* Module 2
estimates; Module 4 assembles. Keeping the assembly in one place is what makes
the two-stage propagation term auditable — it appears once.

*`dwls_delta()` and `dwls_model()` take `theta` explicitly* rather than reading
it from a fitted object, because Module 7 calls them at trial values and
Module 3 calls `dwls_delta()` at *starting* values for the rank check of
Section 4.8, and at its second point when the first is deficient. A signature that presumed a converged fit would not support
either caller.

*`dwls_identify()`'s `modelspec` is the canonical model specification of
Section 4.8, not the user's text.* The signature says `string scalar`, and the
grammar of that string is as much a part of the lock as the type is — Module 11
parses the sem-style surface and normalizes into the canonical form; Module 3
parses only the canonical form. A `dwls_identify()` that accepted surface syntax
would couple Stage 2 to a grammar deliberately left open, and would move
malformed-user-input diagnostics into Mata, where `_error()` truncates
destructively at 100 characters.

*`dwls_identify()` returns fixed parameters as well as free ones.* Required by
Sörbom modification indices (Module 12, v2), which differentiate with respect
to currently-fixed parameters. This is a shape constraint accepted now because
it costs one column now and a refactor of Modules 5–7 later.

*There is no Θ anywhere in the API.* Under delta it is derived, not estimated
(Section 4.6). Its absence is the enforcement.

#### Internals below the locked boundary

These are **recorded, not locked**. They exist here so a later module does not
reinvent one under a second name — which needs their names and purposes, not
their types. **Their signatures live in their own file headers and nowhere
else**, deliberately: a signature carried in two places can disagree, where a
signature three lines above its function cannot.

| Routine | Owner | Purpose |
|---|---|---|
| `dwlsmv_listwise` | 1 | Section 3's missing-data rule, across all variables |
| `dwlsmv_catcounts` | 1 | Category counts, rejecting missing, non-integer and out-of-range codes |
| `dwlsmv_catdegen` | 1 | Degeneracy class: 0 none, 1 boundary, 2 interior, 3 both |
| `dwlsmv_cumprop` | 1 | Cumulative sample proportions |
| `dwlsmv_thresholds` | 1 | τ̂ = `invnormal(q̂)`, and nothing more |
| `dwlsmv_psi_tau_table` | 1 | The §4.5.1 staircase, K × (K−1), one variable |
| `dwlsmv_cellprob` | 2 | Cell probabilities by rectangle differencing |
| `dwlsmv_cellprob_raw` | 2 | The same arithmetic without the abort on a non-positive cell; the root-finder's evaluator (`@sec-emptycells`, `@sec-rootfind`) |
| `dwlsmv_dcellprob_rho` | 2 | Their ρ-derivative, via Plackett |
| `dwlsmv_rootfind` | 2 | The internal root finder. Performs the walk, frontier bracketing, safeguarded Fisher iteration and endpoint verification; accepts the production score evaluator and, for internal testing only, an optional supplied bracket. Returns ρ̂, `status` and the evaluation count. `dwls_polychoric()` calls it with the production evaluator and no bracket. Its steps are `dwlsmv_rf_side` (walk and frontier), `dwlsmv_rf_end` and `dwlsmv_rf_check` (endpoint check), `dwlsmv_rf_taur` and `dwlsmv_rf_ulp` |
| `dwlsmv_rf_eval`, `dwlsmv_rf_evalp` | 2 | One evaluation at ρ: usability, S̄, the masked expected information c, and τ_S. `_evalp` takes the derivative routine as an argument (production passes `dwlsmv_dcellprob_rho`; only tests pass another). `dwlsmv_rf_dP` returns δ_P |
| `dwlsmv_wmat` | 2 | The w-weights |
| `dwlsmv_crosstab` | 2 | One pair's observed contingency table |
| `dwlsmv_emptycell1`, `dwlsmv_emptycell_repair`, `dwlsmv_emptycell_pairs` | 2 | The one-empty-cell pair (O6): whether a table is one (binary 2×2, exactly one empty cell); its keep-margins repair, c = 0.5; the P × 1 membership from `(Y, map)`, used by Stage 1, Module 10 and Module 11 |
| `dwlsmv_pairinfo` | 2 | Fills u, I_rr, I_ra, I_rb at an explicit weight matrix — n/N ships, π is the oracle (§4.5.2) |
| `dwlsmv_koff` | 4 | Row offsets into the stacked `psi_tau` block |
| `dwlsmv_psi_tau_stack` | 4 | Builds that block from (Y, map) |
| `dwlsmv_psi_rho_table` | 4 | One pair's ψ_ρ lookup — **the only place the propagation term lives** |
| `dwlsmv_psi_matrix` | 4 | Gathers the N × D influence matrix by lookup |
| `dwlsmv_psibar` | 4 | The §4.5.5 centring diagnostic, from cell counts |
| `dwlsmv_canon_lines` | 3 | Splits a canonical specification into directive lines; strips `#` comments and `char(13)` |
| `dwlsmv_canon_find` | 3 | Position of a declared name in the factor or indicator list — where a user-facing name becomes a matrix subscript |
| `dwlsmv_canon_status` | 3 | Reads `free [<start>] \| fixed <value>` off a directive, with the per-kind start-value defaults |
| `dwlsmv_nfac` | 3/5/6 | The factor count, **derived** from the table rather than stored — see below |
| `dwlsmv_parvals` | 5/6 | nall × 1 current values, free from θ and fixed from the table |
| `dwlsmv_matrices` | 5/6 | Fills Λ, Ψ, B from those values. Owns the index conventions: λ is (indicator, factor), β is (to, from) |
| `dwlsmv_rystar` | 5 | Λ, Ψ, B → R_y*, with the (I−B) singularity check |
| `dwlsmv_corrvec` | 5/6 | The correlation block of σ(θ) in §4.2 order — the thing Module 6 differentiates |
| `dwlsmv_delta_cols` | 6 | Δ columns for a chosen **set** of table rows, free or fixed (the Module 12 shape) |
| `dwlsmv_fdwls` | 7 | The `optimize()` evaluator: $(s-\sigma)'\operatorname{diag}(w)(s-\sigma)$, its analytic gradient (exact given Δ, which is itself a central difference), and the Gauss–Newton Hessian. $w$ stays a vector |
| `dwlsmv_minimize` | 7 | One weighted least-squares minimization from a start; returns converged 0/1 and never aborts on non-convergence. `dwls_fit()` calls it for ULS, then DWLS. Not polished itself: `dwlsmv_fit_from()` polishes a converged DWLS run |
| `dwlsmv_polish` | 7 | The post-convergence polish (adopted 2026-10-23, `@sec-polish`): `(s, w, par, maps, theta)` → 1 × 3 (status, steps, last relative step); theta by reference. Pure Gauss–Newton steps, no extension and no F comparison; stops at a step ≤ 1e-11 or 20 steps; abandoned (theta unchanged) on a missing step, a first step > 1e-6, or growth above 1e-11. Called only by `dwlsmv_fit_from()`, after a converged run. In `dwls_fit.mata` |
| `dwlsmv_fit_from` | 7 | The DWLS run from a given start, sign canonicalization and the fit struct: `dwls_fit()` after its ULS run. D11's recovery calls it directly for its final fit, skipping ULS (2026-10-07) |
| `dwlsmv_fit_w` | 7 | W_D = diag(Γ)⁻¹ as a vector, after checking Stage 1's output; formed only here |
| `dwlsmv_m11_recstart` | 11 | The recovery's final-fit start: the table's free values, each threshold at s_τ |
| `dwlsmv_signcanon` | 7 | Canonical factor signs after convergence: reverses each reflectable factor whose first free loading is negative, exactly as $\Lambda D$, $D\Psi D$, $DBD$ |
| `dwlsmv_admissible` | 7 | The §5 admissibility conditions at any θ: flag, worst value and location per condition. Derived, so never stored on the fit |
| `dwlsmv_admissadvice` | 7 | Those flags as notes a user can act on — the same flags/prose split as `dwlsmv_diagadvice`; practitioner wording 2026-10-07 (D8e) |
| `dwlsmv_acov` | 8 | Section 4.10's ACOV from (Δ, Γ): expected-information bread, $W_D$ a vector, `cholinv()` on the bread. Divides by no sample size. `dwls_inference()` calls it |
| `dwlsmv_ev` | 8 | $e(V) = \operatorname{ACOV}/(N-1)$ and nothing else; called where `e(V)` is posted (Module 11) |
| `dwlsmv_traces` | 9 | $(\operatorname{Tr}M, \operatorname{Tr}M^2)$ from (Δ, Γ) by Section 4.9's expansion; $U$ and $M$ never formed. Takes matrices so it can be fed lavaan's own Δ and Γ |
| `dwlsmv_scaleshift` | 9 | $(a, b)$ from the traces and $D$; $(.,.)$ when $D = 0$ |
| `dwlsmv_dtboundary` | 12 | DIFFTEST check 11 on one parameter table ($H_0$'s), no data and no free value read: (a) a maximal all-fixed block of Ψ not positive definite (LDL′, strict pivot, no cutoff); (b) a $\Theta_{jj}$ determined by fixed λ, β, ψ alone that is $\le 0$. Absent rows at the G4 default; rc 3498. Group by group at $G>1$ (2026-10-13): rule (b) is $1-\omega_{gj}^2(\Lambda_g\Phi_g\Lambda_g')_{jj}$, a free ω leaving the row undetermined, an absent one reading 1; messages name the group. Also the configural preflight in `dwlsmv_m11_fit_g`. Signature `(par, model, consequence)` (2026-10-13, mechanism (a)): the routine prints the condition, then the caller's consequence sentence, then aborts; an empty consequence is refused rc 3498; rule (b) names the item from its loading rows' names. Helpers `dwlsmv_dtb_cliques` (maximal cliques), `dwlsmv_dtb_pd`, `dwlsmv_dtb_grp`, `dwlsmv_dtb_group`, `dwlsmv_dtb_item` |
| `dwlsmv_dtband` | 12 | DIFFTEST check 10 from $(T_0, T_1, c, \text{nrtol})$, $c$ the multiplier of $F$ in $T$ ($N - G$; computed by the caller): $\varepsilon = c\,\text{nrtol}$; returns (status, $T_d$ or 0, $\varepsilon$), status 1 proceed, 0 within the band; aborts rc 3498 when $T_0 - T_1 < -\varepsilon$, naming $H_1$ as the suspect fit, and on a missing input, $c \le 0$, or $\text{nrtol} \le 0$. `nrtol` is passed in, never restated. Signature changed 2026-10-09 (took $N$). Optional `(exempt, par0, theta0)` (2026-10-13, step 2c): for an exempt pair whose $H_0$ has a free later-group $\hat\psi_{gff} < 0$ ($B = 0$), the abort carries check 10's approved negative-variance message instead (helper `dwlsmv_dtb_negphi`) |
| `dwlsmv_dtstep` | 12 | One adjacent pair after checks 8–9 (2026-10-13, step 2c): `(par0, theta0, d, exempt, T0, T1, c, nrtol, df0, cc, K, inv0, inv1, name0, name1 | note)` → (status, $d$, $T_d$, $\varepsilon$). Check 11 on $H_0$ with the backstop sentence (exempt) or the comparison sentence; at $d \ge 1$ `dwlsmv_dtband` (status 1 test, 2 within the band); at $d = 0$ (exempt only) status 0 "no test", the I2 case (thresholds–metric with `_i2`) verified wording within the band and rc 3498 outside, any other pair the equal-df wording plus the information line outside the band. Prints only aborts; `note` returns the approved text, printed verbatim by the front end (option (ii)). Does not compute the statistic. `test_dtstep.do` |
| `dwlsmv_dtpair` | 12 | One adjacent pair, gate step 3 (decided 2026-10-14, user, option (a); built 2026-10-14, `test_dtpair.do`): `dwlsmv_dtstep`'s inputs plus $\Delta_1$, $\Delta_0$ (stacked, at $\hat\theta_1$ and $\hat\theta_0$) and $\Omega$ from `dwlsmv_omega` → dtstep's row with $(T_{d3}, p, a_d, b_d)$ appended. Calls `dwlsmv_difftest` only when `dwlsmv_dtstep` returns status 1; at status 0 or 2 the four are missing and the statistic is not called. Aborts rc 3498 when the Jacobians' $q_1 - q_0$ is not dtstep's $d$ (another pair's Jacobians). Declined: (b) the call inside `dwlsmv_dtstep`; (c) left to the front end |
| `dwlsmv_dtnest` | 12 | DIFFTEST check 8 on two parameter tables ($G=1$): rows identified by (ptype, row, col), a row absent counts as fixed at its type's absent-row default (§4.8 G4; 0 for every type built at $G=1$); rows present in both carry the same name (factor and indicator order); every parameter free in $H_0$ free in $H_1$; every parameter fixed in both fixed at the same value (exact); $d = q_1 - q_0 \ge 1$ and $d$ = the difference in `e(df_T3)`. Returns $d$; rc 3498. Reads no free row's value. Helper `dwlsmv_dtn_find`. Optional `(cc0, can0, inv0, K0, cc1, can1, inv1, K1 | exempt)` (2026-10-13): after the β guard, `dwlsmv_dtgenerated`; an exempt pair skips rules 0–2 and the tie refusal, needs $d \ge 0$ and $d$ = the `e(df_T3)` difference, and sets `exempt` = 1 |
| `dwlsmv_dtchain` | 12 | DIFFTEST check 9: orders $n \ge 2$ tables by `par.npar` (ties by level rank when every tied model posts a level, scalar first; otherwise by input position, then refused by check 8's $d \ge 1$), tests every adjacent pair with `dwlsmv_dtnest`, returns the permutation, least restricted last. Optional `(cc, can, inv, Kp | info)`: the generated inputs per model (`Kp` pointers to each `e(K)`), and `info`, $(n-1)\times 2$, each pair's $d$ and exempt flag. rc 3498 |
| `dwlsmv_dtgenerated` | 12 | Check 8's exemption (built 2026-10-13): `(cc0, can0, inv0, K0, cc1, can1, inv1, K1, name0, name1)` → 1 iff both post a level and (1) the same `e(canon_configural)`, (2) the same `e(K)`, (3) `dwlsmv_levelgen` reproduces each `e(canon)` exactly; else 0, with one line naming the first failed condition when both post a level. Condition (4) is check 6's (option (ii)). Helper `dwlsmv_dtg_note` |
| `dwlsmv_levelgen` | 12 | The invariance level generator (§4.8 G8; L1 (b), L2 (a), 2026-10-12): `(configural text, K, level)` → that level's canonical text with its ties, in the normative order. Validates that later groups mirror group 1 and that K (G × p) is equal across groups past configural; refuses `beta`, `nu`, `omega`, `alpha` and `tie` in its input. rc 3498 |
| `dwlsmv_m11_estimate` | 11 | One single-group fit with D11's recovery, factored out of `dwlsmv_m11_fit` (2026-10-23) so `difftest`'s refit runs the same code: `(canon, metric, st)` → `par`, `fit` (inference done), `A`, `Bd`, `rec`, `rnote` of the reported fit |
| `dwlsmv_dt_check6` | 12 | `difftest`'s only entry to check 6 (2026-10-23): stored and refit `(b, V, T3, bound, df)`; today `dwlsmv_dtrefit` at `dwlsmv_dt_reftol()`, `V` carried for the open contract and not compared; rc 3498 |
| `dwlsmv_dtfront` | 12 | `difftest`'s Mata side (2026-10-23, G = 1): `_init`, `_add` (reads the active `e()`), then `(touse, ind, kb)` → shared Stage 1, refits, checks 6–11, statistic, fit indices and Δ, both tables, `r(models)`, `r(steps)`, `r(N)`; check 6 under `_stata()` so its refusal gains the Stata version/edition line; helpers `_post`, `_c6`, `_print`, `_lab`, `_delta` (ΔCFI/ΔRMSEA per step, H0 minus H1) |
| `dwlsmv_fitidx_base` | 12 | The c1 baseline from Stage 1 alone (`@sec-fit-indices` FI2): `st` (one per group) → $(T_{3B}, D_B, T_B, a_B)$; $D_B=0$ gives $(.,0,.,.)$; a missing polychoric or Γ entry, or a non-positive polychoric variance, rc 3498 |
| `dwlsmv_fitidx` | 12 | Family S with F1–F4 as frozen 2026-10-21: $(T_3, D, n, G, T_{3B}, D_B)$ → (RMSEA, CFI, TLI); model side and baseline side independent (F4); $n\le0$ or missing, or $G$ not a positive integer, rc 3498 |
| `dwlsmv_geomin` | 12 | Oblique geomin of an unrotated EFA $A$ (`@sec-geomin`): GM4 grid starts, stall-exit GPA, converged iff $\|G_p\|_F<10^{-5}$; returns rotated Λ, sets Φ and the 8 `e()` diagnostics; lavaan's sign and index order; invalid $A$, ε, start count, or no converged start, rc 3498 |
| `dwlsmv_fitidx_delta` | 12 | ΔCFI/ΔRMSEA across a chain (FI3, F6): $X$ (CFI, RMSEA per level), skip flags → differences against the last non-skipped earlier level, skipped rows missing; called by `dwlsmv_dtfront`; malformed input or a skipped first level, rc 3498 |
| `dwlsmv_m11_fitidx` | 11 | Sets the locals posted as `e(rmsea_s)`, `e(cfi_s)`, `e(tli_s)`, `e(T3_b)`, `e(df_T3_b)`, with $n=N-G$ |
| `dwlsmv_dtrefit` | 12 | DIFFTEST check 6's comparison of a refit's (`e(b)`, `e(T3)`, `e(bound)`, `e(df_T3)`) with the stored values: $\max_i |x_i-y_i|/(1+|y_i|) \le$ tol for `e(b)` and `e(T3)`, `e(bound)` identical, `e(df_T3)` equal, a missing `e(T3)` refused; rc 3498. tol is an argument; a missing or non-positive tol is refused |
| `dwlsmv_sens` | 12 | **Diagnostic only** for the check 6 G > 1 investigation (built 2026-10-14, `test_sens.do`). No production path calls it, and it changes no tolerance. `(f, step, epsV, A, c, nu \| ds, eps)` → npar × 6 in check 6's units: κˢ and κʷ (S1 per unit perturbation); $S^{(C)}$ (S3(4), needs the measured last step); the S4(b) stopping envelope; $S^{(A)}$ from measured $(ds, \varepsilon)$, else missing; and the $S^{(D)}$ envelope (S5(c), with epsV the evaluation accuracy of $V$, an argument with no default). Fills $A=J^{-1}\Delta'W$, $c$ and $\nu$. $J=P-R$, with $R$ a central difference of $\Delta'(w\circ r)$. rc 3498 on an unconverged fit, bad step/ds/eps/epsV, iterations < 1, or $c\ge1$. `@sec-mg-sensitivity` |
| `dwlsmv_sens_t3` | 12 | **Diagnostic only**, the `e(T3)` side of the check 6 G > 1 investigation (built 2026-10-16 for E-029, approved item ⚑5; `test_sens_t3.do`). No production path calls it, and it changes no tolerance. `(f, step, epsV, rhoT \| ds, eps, dOm)` → 1 × 8: $T_3$; then, in check 6's units, $2\beta_1$, $\beta_2$ (with $\Theta$ the a-priori S1–S5 envelope from `dwlsmv_sens`), $\beta_3$, $\beta_{\rm eval}$ (TS6, with the **assumed** constant $c_{\rm ch}=6(q+1)$ for `cholinv()`, and since 2026-10-20 SE6's Q4 form $2\gamma_{n-1}\sum\lvert\text{terms}\rvert$ per `quadcross()`, $n=D$, or $q$ for the $E^{-1}$ products; Q1's $2u\sum\lvert\text{terms}\rvert$ was refuted per run by the exact-reference reading, E-030, user's decision), $\beta_{\rm elem}$, $B_{T3}$ (TS7) and $B_{\rm eval}=\beta_{\rm eval}/(1+\lvert T_3\rvert)+r_T$. (ds, eps, dOm) are all three or none; none is the Stage-1-fixed case. Needs `dwls_inference()` to have run. rc 3498 on an unconverged fit, inference not run, $D=0$, $\lvert T_3\rvert<10^{-2}$, negative or missing epsV or rhoT, a partial or malformed (ds, eps, dOm); `dwlsmv_sens`'s own refusals apply. Helpers, each diagnostic only and pinned by the test: `dwlsmv_sens_ebar` (S5's $\bar E$, `dwlsmv_sens`'s loop copied), `dwlsmv_sens_jmat` ($J=P-R$, copied), `dwlsmv_sens_ddelta` ($\partial\Delta/\partial\theta_l$ by central differences of production $\Delta$, $h=10^{-4}\max(1,\lvert\theta_l\rvert)$). `@sec-mg-sensitivity-t3`. Since 2026-10-19 the routine is `dwlsmv_sens_t3_core` with `exact` = 0; its arithmetic is unchanged |
| `dwlsmv_sens_hess` | 12 | **Diagnostic only** (plan E-031 step 3, built 2026-10-19; `test_sens_hess.do`). `(theta, par, maps)` → $D\times q^2$: Result MS6's closed-form $\partial^2\sigma/\partial\theta_i\partial\theta_l$ in column $(l-1)q+i$ (`dwlsmv_sens_ddelta`'s layout). **Domain:** `dwlsmv_sens_hdom(par, G, caller)` certifies the table: (c1) no ptype 3 row, (c2) one factor, (c3) at most one row per θ element per group. Outside the domain it refuses at rc 3498, naming each failing clause, and has **no finite-difference fallback**; models with β are not covered. `@sec-mg-hessian` |
| `dwlsmv_sens_ddelta_x`, `dwlsmv_sens_jmat_x`, `dwlsmv_sens_x`, `dwlsmv_sens_t3_x` | 12 | **Diagnostic only, opt-in exact variants** (2026-10-19; `test_sens_hess.do`). They have the same arguments, outputs and refusals as `dwlsmv_sens_ddelta`, `dwlsmv_sens_jmat`, `dwlsmv_sens` and `dwlsmv_sens_t3`, with two differences: $\partial\Delta/\partial\theta$ and $R$ in $J=P-R$ come from `dwlsmv_sens_hess`, and every variant refuses outside Result MS6's domain. $P$ keeps production $\Delta$. The finite-difference routines are unchanged and stay the default, because the frozen grid and T3 studies ran them. `test_sens_t3` does not use the exact path (plan steps 4–5). `@sec-mg-hessian` |
| `dwlsmv_dt_reftol` | 12 | Check 6's tolerance, **1e-10, chosen after measuring** (2026-10-09): 100 × the observed processor-count maximum 8.05e-13, rounded up to a power of ten (rule frozen before the run); re-checked at G > 1 |
| `dwlsmv_difftest` | 12 | The DIFFTEST statistic for one nested pair, any $G$: from $(\Delta_1, \Delta_0, \Gamma, T_d)$; at $G > 1$ (2026-10-14, MG7(3), `mg-difftest-omega`, signature unchanged) the caller passes the stacked $\Delta_j$, $\Omega = \operatorname{blockdiag}(\Gamma_g/\pi_g)$ from `dwlsmv_omega` in the $\Gamma$ slot (so $W = \operatorname{diag}(\Omega)^{-1}$ is weighting (a)) and $T_d$ with $T = (N-G)F$; $\hat h = \Delta_1^{+}\Delta_0$, $M_1$ on $\hat h$ with $V_1$ from `dwlsmv_acov`, then $(T_{d3}, d, p, a_d, b_d, \operatorname{Tr}M_1, \operatorname{Tr}M_1^2)$. No validation (nesting, sign, ε) and no $N$; $d = 0$ is an error, the caller's to skip |
| `dwlsmv_asym` | 10 | $(a, 4D\varepsilon, i, j)$, $(i, j)$ where $a$ is attained ($i<j$; $(0,0)$ when $a=0$; added 2026-09-26 for `gammadiag`'s row 1), with $a = \max|A_{ij}-A_{ji}|/\sqrt{A_{ii}A_{jj}}$: the **one** definition of Γ's symmetry measure and budget. `dwlsmv_acov` and `dwlsmv_traces` abort with rc 3498 when $a$ exceeds the budget, and Module 10's diagnostic reports the same comparison; none carries its own tolerance |
| `dwlsmv_npatterns` | 10 | $n_{\text{pat}}$ from $(Y, K)$ as §4.5.6 defines it, with Stage 1's per-column validation (rc 3498 on a missing, non-integer or out-of-range code) |
| `dwlsmv_thralias` | 10 | Pairwise relations among the columns of $[\mathbf 1\mid Z]$ (identical, complementary, constant), found from the exact integer counts of $Z'Z$ with no cutoff; rows $(a, b, \text{type})$ in $s$'s positions. Relations among three or more columns are a recorded v1 limit |
| `dwlsmv_countdiag` | 10 | The small-count family from $(Y, K)$: 2 × 4, the smallest boundary- and interior-category counts with (variable, category position), every flag `.n`. Same validation as `dwls_sample_stats()`; an empty declared category reads 0 rather than aborting. No pair-cell rows until the S partition decides their role |
| `dwlsmv_gammadiag` | 10 | Γ̂'s seven rows from `(stage1, Y)`, 7 × 4: asymmetry, not PSD, too few patterns, threshold aliasing, centring, raw ratio, scale-free ratio (§5 Module 10). A flagged asymmetry stops it: rows 2–7 are `.`. rc 3498 when Y is not the matrix Stage 1 was built on, or Γ/ψ̄ do not match the map. Folds in the `dwlsmv_symdiag()` of the 2026-09-25 shape row |
| `dwlsmv_dalias` | 10 | $d_{\text{alias}}$ from `dwlsmv_thralias()`'s relations, by connected components: each constant column 1, a component of $m$ columns $m-1$. Its own function so the count is testable on relations no Stage 1 could be run on |
| `dwlsmv_numrank` | 3, 10 | The one numerical-rank criterion: the number of singular values above $10^{-8}$ of the largest, with the singular values and right vectors (as rows) returned for the caller. `dwls_identify()` aborts on it at the starting values; `dwlsmv_breaddiag()` reports it at $\hat\theta$ |
| `dwlsmv_breaddiag` | 10 | The bread's four rows from `(Delta, Gamma)`, 4 × 4: npar > D; an identically zero Δ column; numerical rank of $W_D^{1/2}\Delta$ (`dwlsmv_numrank`); its $\sigma_{\min}/\sigma_{\max}$, flag `.n`. npar is `cols(Delta)`. Rows 3–4 judge only the non-zero columns and are `.` when row 1 flags or no column is non-zero. rc 3498 on an empty or missing Δ, or a Γ that is not D × D with a positive diagonal |
| `dwlsmv_diagadvice` | 10 | The three row matrices as prose, `(G, C, B, ind, slab, pname)`. Flag notes by family, only when a row fires: Γ̂ (the sampling covariance of the thresholds and correlations) and the identification of the parameters at the estimate (the bread), in the 2026-10-07 practitioner wording (D8e). Then **always** the category counts, as descriptive diagnostics with no threshold (D8e), never worded as Γ conditioning. The singular-value ratios are not printed (stored by Module 11, D8d). Unevaluated rows read "not evaluated". It prints nothing and aborts only on a contract violation (rc 3498): wrong shapes, or a location outside its labels. Admissibility is not an input (separate channel) |

Modules 1 and 2 keep their routines taking a single column or a single pair
rather than `Y`, which is what keeps `dwls_sample_stats()` the only place the
loops live.

**Three of these constrain the caller, not just the implementer:**

- **`dwlsmv_listwise()` is deliberately not a step inside `dwls_sample_stats()`.**
  Deletion changes $N$, and $N$ is what every influence function is divided by;
  buried inside Stage 1 it would make the $N$ that Γ was built with invisible to
  the caller reporting it. The caller deletes, sees `rows()` change, and hands
  Stage 1 a matrix whose row count is the $N$ it will report. Module 11 still
  performs the authoritative deletion at the dataset level, since `e(sample)`
  cannot be marked from a matrix whose rows have already been dropped.
- **`dwls_polychoric()` does not set the `i` and `j` fields** of the pairfit it
  returns. The locked signature is not told which variables it was handed, which
  is deliberate — a pair estimator has no business knowing its position in the
  model. The assembler fills them, and `dwlsmv_psi_matrix()` checks them against
  the map.
- **`dwlsmv_catdegen()` classifies and does not act.** Section 4.5.6 requires the
  interior/boundary distinction to survive as far as Module 10; a routine that
  repaired, collapsed or warned would destroy it here.
- **`dwlsmv_nfac()` is a function and not a field on `dwlsmv_parmap`**, for the
  reason `dwlsmv_koff()` is: a derived quantity stored in a struct is a second
  copy that can disagree with the first. It is safe to derive *because*
  `dwls_identify()` requires every declared factor to carry exactly one metric
  constraint, and both forms of that constraint emit a row indexed by the
  factor. Relaxing that invariant means storing the count, not patching the
  maximum.

Both future extensions (Module 12) should only need to touch
`dwls_sample_stats()`/`dwls_gamma()` (continuous indicators) or
`dwls_model()`/`dwls_delta()`/`dwls_inference()` (ESEM's rotation) — neither
should need to know how the other parts of the backend work.

**One efficiency claim, stated precisely because it is easy to overstate.** The
lookup tables mean no *bivariate normal evaluation* is repeated per observation
— those are per cell, $K_iK_j$ per pair. Forming Γ still touches every
observation, because the cross-product of two statistics resting on different
pairs needs the joint over as many as four variables. Do not read "lookup
table" as "no $O(N)$ work".

## 6. Risk matrix

| Module | Risk | Severity | Mitigation |
|---|---|---|---|
| Module 4 (Γ) | Two-stage influence function implemented without the threshold-propagation term, silently producing a different (wrong) estimator that may still look plausible | Critical | Explicit derivation and unit tests for $\partial\rho_{jk}/\partial\tau_j$ propagation terms (Section 4.5); validate against Monte Carlo $N\cdot\text{Var}(\hat s)$, not just lavaan. |
| Module 9 (test statistic) | Wrong WLSMV variant, or an under-derived scale-and-shift correction | Critical | Hard gate: derivation appendix from Muthén, du Toit & Spisic (1997) checked algebraically before coding; validate against `test="scaled.shifted"` specifically. |
| Module 0/4 | Validating Γ against the wrong lavaan object (`wls.v` instead of `gamma`) | High | Use `lavInspect(fit, "gamma")` for Γ; keep it explicitly distinct from `wls.v` in the harness code (Section 2, Module 0). |
| Module 4/4.4 | A stray factor of $N$ migrating between Γ, the objective, and `e(V)` | High | Fix the $\sqrt{N}(\hat s - s) \to N(0,\Gamma)$ convention explicitly and unit-test the final $/(N-1)$ scaling step in isolation, requiring lavaan's SEs to match $/(N-1)$ and reject $/N$. |
| Modules 2/4 | Sparse data / empty contingency cells driving correlations to ±1 or failing to converge | High | Fail/warn by default; continuity correction as an explicit opt-in mode added after validation. |
| Modules 1/8/9 | An empty **interior** category collides two thresholds at a finite value, making Γ silently rank-deficient. DWLS point estimates are unaffected (W_D uses only the diagonal), so the command converges and returns plausible $\hat\theta$ with wrong robust SEs and a wrong test statistic | High | Distinguish interior from boundary empty categories in Module 1's category validation; promote Module 10's Γ rank/condition check from diagnostic to a precondition of Modules 8–9 (Section 4.5). Point-estimate agreement with lavaan will *not* surface this. |
| Module 4 | Γ built by plugging estimates into closed forms instead of averaging influence functions — mathematically defensible, but disagrees with lavaan at $O(N^{-1/2})$ and will be reported as an SE bug | Medium–High | Empirical construction is a compatibility requirement (Section 4.5, Section 3 matrix); closed forms are unit-test oracles only. |
| Module 3 | Factor-metric and ordinal-response-scale identification conflated, producing a $\theta$/Δ mismatch with validation targets | High | Keep them as separate rows in the parameter table (Section 4.8); record the choice in the parameterization compatibility matrix. |
| Module 7 | Optimizer converges to an inadmissible solution | Medium–High | Explicit admissibility check on every fitted $\hat\theta$. |
| Module 7 | $F_{DWLS}$ local minima from poor starting values | Medium | Use the ULS solution as the starting value. |
| Module 2 | Root-finder cost across many pairs, dominated by the endpoint check (54–105 evaluations per pair against 7–17 for the root) | Medium: the mitigation is concrete and extensively exercised, but the module stays numerically delicate and its errors propagate downstream (user, 2026-09-30) | Superseded 2026-09-29: `optimize()` is replaced by a bracketed root-finder on the analytic score (`@sec-rootfind`), which on the failing high-ρ pairs takes 14–17 evaluations against `optimize()`'s 6,457–12,351. The endpoint check's four-ulp width was kept knowing its cost (2026-09-30). Section 7.5's timing re-measures it. |
| Module 11 | Reimplementing `sem`'s path parser by hand, or routing through `ml`'s `d0`/`d1`/`d2` evaluators despite their documented refusal to support `vce(robust)`/custom sandwich VCEs — exactly what Module 8 needs | Medium (schedule/architecture risk) | Simplified syntax for v1; strong prior toward `moptimize()` + `ereturn post` (verified reason in Module 11), confirmed by investigation rather than assumed. |

## 7. Validation plan

Three components, in increasing scope: a paired replication design that runs both implementations
on identical data (7.1), the numerical tolerances those comparisons are judged against (7.2), and
an independent Monte Carlo check on Γ that never touches lavaan (7.3). The tolerance table lives in
7.2 only. It has been carried in two verbatim copies before, which is precisely the kind of thing
that drifts out of sync once one copy is edited.

**How the lavaan reference is fitted — a requirement, not a suggestion.** Every relevant option
is stated explicitly rather than left to package defaults, so that a future change to a lavaan
default surfaces as a spec violation instead of as reference numbers that quietly moved:

```r
lavaan(..., estimator = "WLSMV", parameterization = "delta", test = "scaled.shifted",
       zero.add = c(0.5, 0.0), zero.keep.margins = TRUE)
```

`parameterization = "delta"` currently matches lavaan's default; stating it makes the
compatibility target immune to that changing. `zero.add = c(0.5, 0.0)` is lavaan's documented
empty-cell handling — 0.5 for 2×2 tables, nothing for larger ones. `zero.keep.margins = TRUE`
is lavaan's effective default (set by `lav_options_mimic()`) and is pinned for the same reason.
**dwlsmv applies the same repair to a one-empty-cell pair** (decided 2026-10-06; Module 2,
*One-empty-cell repair*), so lavaan's default fit is the reference there: s, Γ̂ (construction
(a)), θ̂, SEs and T3 (`90_lavaan_emptycell.R`, `test_lavaan_emptycell.do`). A larger table is not
repaired by either. A perfect-association pair is repaired by neither, and dwlsmv aborts where
lavaan fits at its 0.999 cap. Comparisons that rerun lavaan at `zero.add = c(0, 0)`, such as 70_'s
`ref` arm, match production only where no pair has exactly one empty cell (70_'s header).
Unpinned, a real discrepancy and a zero-cell mismatch are indistinguishable.

**The reference objects, named precisely.** Two of these have been swapped or mislabeled in a
draft of this spec before, which is the mistake this table exists to prevent — check the names
against `lavInspect`'s documentation rather than from memory:

```text
lavInspect(fit, "wls.obs")  → observed sample-statistics vector (s)
lavInspect(fit, "wls.est")  → model-implied sample statistics (σ(θ̂)) — NOT the observed vector
lavInspect(fit, "gamma")    → Γ (N × asymptotic covariance matrix of the sample statistics)
lavInspect(fit, "wls.v")    → lavaan's stored WLS weight object; never use it as Γ
lavInspect(fit, "th")       → thresholds
```

**`names(lavInspect(fit, "wls.obs"))` is the authoritative ordering of Γ**, not
`dimnames(gamma)`. The reason is structural: Γ is by definition the asymptotic covariance of the
WLS statistics vector, so that vector's order is Γ's order by construction, where `dimnames` is a
convenience a future release could populate differently. If the two disagree, `wls.obs` wins.
Aligning comparisons to §4.2's index map by name gives a permutation-free comparison that still
fails loudly if lavaan reorders.

### 7.1 Primary validation design: paired Monte Carlo replication

A single shared dataset run through both Stata and lavaan only establishes that the two agree (or
don't) on that one draw — and a one-off comparison has two specific blind spots. A small
*systematic* discrepancy (e.g. a slightly wrong scaling term in Γ) can be masked by ordinary
sampling variability on any single dataset, so a real bug looks like "close enough." And a single
dataset that happens to hit a sparse-cell or near-boundary case can produce a large discrepancy
that's really just that one draw being unusual, not a bug at all. Neither failure mode is
distinguishable from a one-off check. The fix is to generate many replicate datasets from the same
population model and run *both* pipelines on every one of them, so discrepancies can be examined as
a distribution rather than a single point.

**Design:**

1. **Generate once, in one place.** For each population model / sample-size combination in
   Section 7.3's coverage list, generate a few hundred to a thousand replicate raw item-level
   datasets in R (e.g. `lavaan::simulateData()` or a custom generator), under a fixed, recorded
   seed for reproducibility.
2. **Write each replicate to a format both tools read identically.** `.dta` is the natural choice
   — Stata reads it natively, R reads it via `haven`. Explicitly verify the round-trip is
   lossless (no silent precision loss in the write/read cycle) before trusting anything
   downstream — this is exactly the kind of silent-mismatch source the whole exercise exists to
   catch, so it's worth a dedicated check rather than an assumption.
3. **Run both estimators on every identical replicate.** Fit `dwlsmv` (Stata) and
   `estimator = "WLSMV"` (lavaan) on each of the same replicate files; save parameter estimates,
   SEs, the test statistic, convergence status, and admissibility flags (Module 7) from both.
4. **Compute per-replicate paired discrepancies**, not just a pooled average: $\hat\theta_{\text{Stata}}
   - \hat\theta_{\text{lavaan}}$, $\widehat{SE}_{\text{Stata}} - \widehat{SE}_{\text{lavaan}}$, and the
   test-statistic difference, for every replicate, every parameter.
5. **Inspect the discrepancy distribution, not just its mean.** A systematic bug shows up as a
   discrepancy distribution centered away from zero even when no individual replicate looks
   alarming. A conditional bug (only manifesting under sparse cells, small $N$, or specific
   parameter regions) shows up as a pattern when discrepancy is plotted or regressed against
   per-replicate dataset characteristics (cell sparsity, $N$, true parameter values) — this is
   real diagnostic power a single-dataset comparison structurally cannot provide, since there's
   nothing to correlate against with a sample size of one.

This single design subsumes what used to be two separate tiers: because the population parameters
used to generate each replicate are known, the same replicate set also gives, for each estimator
independently, parameter bias, empirical SD of $\hat\theta$ vs. mean reported SE (should match
under correct inference), coverage of nominal 95% CIs, Type I error rate of the model test
statistic under a correctly-specified population model, convergence rate, and the frequency of
inadmissible solutions (Module 7) — all computed on the same replicate set as the paired Stata-vs-
lavaan comparison, rather than as a separate simulation run.

Apply the Section 7.2 tolerances to the paired discrepancies at the per-module numerical-validation
stage (Modules 0–10, on a handful of hand-picked or small-$N$ replicates where every intermediate
object can be checked); apply them more loosely, and mainly as an "is the discrepancy distribution
centered near zero with no systematic pattern" check rather than a per-replicate pass/fail, at the
full end-to-end stage across hundreds of replicates.

### 7.2 Tiered numerical tolerances, absolute + relative

A pure relative-error criterion breaks down near zero (some Γ/Δ entries legitimately are near
zero), so specify both: accept when $|a - b| \le \epsilon_{\text{abs}} + \epsilon_{\text{rel}}
|b|$.

| Component | Suggested tolerance |
|---|---|
| Thresholds | $\epsilon_{\text{abs}} \approx 10^{-8}$, $\epsilon_{\text{rel}} \approx 10^{-8}$ |
| Bivariate correlations | $\epsilon_{\text{abs}} \approx 10^{-8}$, $\epsilon_{\text{rel}} \approx 10^{-8}$ |
| Γ — see the split below | (one row is not enough; a single ^{-5} cannot cover them) |
| DWLS objective value | $\epsilon_{\text{rel}} \approx 10^{-8}$ |
| Parameter estimates ($\hat\theta$) | $\epsilon_{\text{abs}} \approx 10^{-5}$, $\epsilon_{\text{rel}} \approx 10^{-4}$ |
| Standard errors | $\epsilon_{\text{abs}} = 10^{-8}$, $\epsilon_{\text{rel}} = 10^{-5}$ |
| $T_3$ | $\epsilon_{\text{abs}} = 10^{-6}$, $\epsilon_{\text{rel}} = 10^{-5}$ |
| $\operatorname{Tr}(M)$, $\operatorname{Tr}(M^2)$, $a$ | $\epsilon_{\text{rel}} = 10^{-5}$ |
| $b$ | $\epsilon_{\text{abs}} = 10^{-6}$, no relative term |
| df | exact |
| Tier 3 (Mplus), SEs and $T_3$ | $\epsilon_{\text{abs}} \approx 10^{-4}$, $\epsilon_{\text{rel}} \approx 10^{-3}$, except where a convention check sets its own (Section 4.9, 4.10) |

**The standard-error and test-statistic rows were set on 2026-09-24 from the architecture, not
from what passes** (E-032). They replace $10^{-4}$ abs + $10^{-3}$ rel for SEs and
$10^{-3}$ rel for the statistic. Those looser rows were only placeholders. At SE sizes near 0.02 the
absolute term alone made the SE row about 0.5% relative. As a result the row could not see the
$N$-versus-$(N-1)$ conventions of Sections 4.9 and 4.10 ($\approx 2\times10^{-4}$ and
$1.0\times10^{-4}$), nor a unit weight or observed information on case D.

Each row now follows from two things: the source of the end-to-end gap, and what can be zero by
construction.

- **Relative terms.** Against lavaan the gap is lavaan's own convergence. Its ρ̂ carries a score
  of up to $8.8\times10^{-5}$ on case E. This reaches the SEs at $\le 1.3\times10^{-7}$ and $T_3$
  at $\le 1.0\times10^{-7}$ relative, which `test_acov.do` check 3 and `test_inference.do`
  check 3 print. $10^{-5}$ gives about 75–100× headroom over that, while still catching both
  $1/N$ conventions: $T_3$'s by 14–18× (asserted, `test_inference.do` check 4). The SE divisor's
  shift is 9.2–9.7× the tier on the current SEs, measured once in R by the Module 9 review; this
  was first written as "10× or more", which it is not. At the tier's edge that is enough, and
  `test_acov.do` checks 2 and 4 pin the divisor exactly anyway.
- **Absolute terms: which quantities have one, and why.**
  - $a$ and both traces need no absolute term: they are bounded away from 0 whenever $D > 0$.
  - $b$ is a cancellation and can be arbitrarily small (case D: 0.001). It therefore takes an
    absolute tolerance only, set at about 70× the $1.4\times10^{-8}$ measured.
  - $T_3 = aT + b$ tends to $b$ as the fit becomes perfect, so it inherits $b$'s absolute term.
  - SEs are positive but not bounded away from zero by construction. Their absolute term never
    binds on the current fixtures (smallest SE 0.012). It keeps the criterion meaningful on a
    fixture where it would.
- **The θ̂ row is unchanged.** Estimates can legitimately be zero, and its absolute term is what
  handles that.
- **Mplus keeps the old row.** Its θ̂ differs from lavaan's by about $10^{-5}$, a cross-tool
  difference that is not ours to close.

Loosen a given tier only when the competing implementations are known to use genuinely different
numerical approximations at that step, documented explicitly rather than loosened silently to make
a test pass.

#### 7.2.1 Attribution to the reference's stopping point (validation policy, user, 2026-10-11)

**Status: defined and committed before its first general use (910c1db); first executed 2026-10-11 in
`test_fit_groups_betag1.do` check 6 (5 elements attributed, worst 0.46 of the bound; condition 2 at
most 8.6e-4 of τ). Negative control: a 3e-5 relative error planted in ACOV failed every level
through condition 2; condition 3 alone did not separate it, because its bound carries the
condition-2 residual. The policy's protection against a formula error therefore rests on
condition 2, which is why condition 2 covers every element and not only the missing ones.** It
generalises the rule `test_fit_groups_beta.do` check 4 introduced for $a$ and $b$ (2026-10-10, chosen
after a failure there) and was adopted as policy after `test_fit_groups_betag1.do` check 6 failed on
$e(V)$ (2026-10-11). It is a controlled policy, not an exception for one test. A tier is never
loosened by it: every element is first judged by its own tier, and the rule can only attribute a miss
whose size and cause it predicts.

**Why.** An end-to-end comparison of an inference quantity mixes two things: our formula, and the
point where each optimizer stopped. When the reference stops short of the minimum, a quantity that
depends on $\hat\theta$ moves by roughly its derivative times the gap in $\hat\theta$. The θ̂ tier
already admits that gap, so an inference tier set from the formula's precision can miss on it
without any error in the formula. The rule asks the comparison to show, at the reference's own
$\hat\theta$, that the formula agrees, and that the miss at ours is the predicted consequence of
our better minimum.

**Scope.**
- *Applies to:* end-to-end comparisons against a **tier-2 (lavaan)** reference fit, of inference
  quantities that are functions of $\hat\theta$ through `dwls_inference()`: $a$, $b$, $T_3$, standard
  errors and the covariances of $e(V)$. It applies in a test that declares it, per check and per
  level.
- *Does not apply to:* $\hat\theta$ itself (its own tier covers the convergence gap); Stage 1
  statistics (their own Newton-step rule); $T_{\text{raw}}$ (its propagated bound); tier-3 (Mplus)
  comparisons; population, Monte Carlo or sampling comparisons; anything compared at matched inputs,
  where agreement is expected to machine precision.
- *Existing instances stand as written.* `test_fit_groups_beta.do` check 4 (for $a$, $b$) keeps its
  own constants, which are at least as strict as the general ones below.

**The rule.** Let $\theta_o$ be our $\hat\theta$, and $\theta_l$ the reference's, sign-matched and
mapped to our θ order (a tie is one element). Let $d = \theta_o - \theta_l$. For an element $e$ of
the compared quantity:
- $x_e(\theta)$ is our value from `dwls_inference()` at $\theta$, with Δ rebuilt there and Γ, π
  and $D$ the fit's;
- $x_{L,e}$ is the reference's value;
- $\tau_e$ is the element's tier evaluated at the reference: $\epsilon_{\text{abs}} +
  \epsilon_{\text{rel}}|x_{L,e}|$, and for a covariance of $e(V)$, $10^{-8} + 10^{-5}
  \sqrt{V_{L,ii} V_{L,jj}}$.

An element that misses its plain tier at $\theta_o$ passes only if all three hold:
1. **Better minimum (per level):** $F(\theta_o) \le F(\theta_l)(1 + 10^{-12})$, with $F$ our DWLS
   objective.
2. **Formula agrees at the reference's point (per level, every element of the quantity, not only the
   missing ones):** $|x_e(\theta_l) - x_{L,e}| \le 10^{-2}\,\tau_e$.
3. **Predicted size (per missing element):** $|x_e(\theta_o) - x_{L,e}| \le 2P_e + |x_e(\theta_l) -
   x_{L,e}|$. Here $P_e = \sum_k |\partial x_e/\partial\theta_k|\,|d_k|$ is a triangle bound on the
   first-order change, so terms cannot cancel. The factor 2 is the Newton-step rule's mean-value
   factor. The second term is the residual already bounded by condition 2.

**Numerical method.** Each $\partial x_e/\partial\theta_k$ is a central difference at $\theta_o$,
with step $h_k = 10^{-6}\max(1, |\theta_{o,k}|)$. For $e(V)$ one call of `dwls_inference()` per
perturbed θ gives every element at once: $2q$ calls per level. A standard error's derivative is
taken on $\sqrt{V_{ii}}$ directly, not propagated from $V_{ii}$.

**Edge cases.**
- $d_k = 0$ contributes nothing, so $d = 0$ gives $P = 0$ and condition 3 demands agreement within
  condition 2's residual.
- Anything missing or non-finite ($F$, $x$, a derivative, a reference value) makes the rule fail. It
  never passes by default.
- The rule is evaluated and printed at every level the test declares, including levels where nothing
  missed. It decides only elements that missed their plain tier.
- The test reports how many elements missed, and how many the rule attributed.

**Consequences.**
- A pass under the rule supports attribution to the reference's stopping point. It is not
  agreement within the plain tier, and the test says which it was.
- A failure under the rule is a failure and is kept.
- The rule is never used to compare $\hat\theta$, and no constant in it is changed after a run.

**Γ needs four tolerances, not one, because "compare Γ" means four different comparisons whose
achievable agreement differs by ten orders of magnitude.** A single ^{-5} is
simultaneously far too loose for two of these comparisons and unachievable for a third. Chasing an
unachievable tolerance, or passing a test that should have been tight, are both expensive.

| Comparison | What is being compared | Achievable | Tolerance |
|---|---|---|---|
| Within-variable threshold block, empirical vs closed form | The same arithmetic (Section 4.5 identity) | Machine precision | $\epsilon_{\text{abs}} \approx 10^{-12}$, no relative term. A failure is a coding error, never noise. |
| `dwlsmv` Γ̂ vs `lavInspect(fit,"gamma")`, **at a matched $\hat\rho$** | Two *empirical* constructions from the same $\psi$ | Floating-point | $\epsilon_{\text{abs}} \approx 10^{-12}$, $\epsilon_{\text{rel}} \approx 10^{-12}$. This is the real test of whether the two $\psi_\rho$ agree; $10^{-5}$ would let a genuinely different influence function pass. Measured 1.4e-13 on `case_D`. |
| `dwlsmv` Γ̂ vs `lavInspect(fit,"gamma")`, **end to end** | The above, *plus* two independently converged $\hat\rho$ | Capped by $\hat\rho$ | $\epsilon_{\text{abs}} \approx 10^{-8}$, $\epsilon_{\text{rel}} \approx 10^{-8}$ — the **bivariate-correlation** tier above, not a loosened Γ tier. Measured 9.1e-10. |
| Cross-variable entries, empirical Γ̂ vs closed form at the same estimates | A sample proportion estimating $\Phi_2$ | $O(N^{-1/2})$ | Not a fixed number — calibrate per design (see below). At $n=5000$ on a four-indicator model the max gap runs ≈ 0.06. |
| Analytic Γ vs Monte Carlo $N\cdot\text{Var}(\hat s)$ (Section 7.3) | Sampling variability across $R$ replications | $O(R^{-1/2})$ | Not a fixed number — use $z$ statistics (see 7.3). **$10^{-5}$ here would need $R \approx 2\times10^{10}$ replications.** |

**Why the lavaan comparison is two rows and not one, since the split looks like a way of
excusing a miss.** A single ^{-10} would fail: the end-to-end comparison measures
9.1e-10. The cause was established by substitution rather than assumed: rebuilding
Γ̂ at *lavaan's own* $\hat\rho$ collapses the gap to 1.4e-13, the same agreement the threshold
block reaches. So the arithmetic matches and what is left is that the two implementations
stop at different points — and the direction settles who should move. Our score vanishes at
$\hat\rho$ to 4.8e-16; lavaan's to 2.6e-10. **lavaan is the less converged of the two, so
tightening our optimizer cannot close the gap.** A Γ comparison that inherits $\hat\rho$'s
uncertainty therefore inherits $\hat\rho$'s tolerance tier, which is $10^{-8}$. Both rows are
asserted in `tests/unit/test_gamma.do`; keeping the tight one is what stops the loose one
from hiding an arithmetic difference.

For the two sampling comparisons, do not invent a percentage. Raw percentage deviation is an
actively misleading criterion on a covariance matrix: the small off-diagonal entries carry large
*relative* Monte Carlo error, so the worst percentage in the matrix is dominated by noise on
entries nobody cares about, and a threshold set to accommodate it would be far too permissive for
the entries that matter. Report a $z$ statistic against each entry's own Monte Carlo standard
error, $\operatorname{SE}(\hat\sigma_{ab}) \approx
\sqrt{(\sigma_{aa}\sigma_{bb} + \sigma_{ab}^2)/R}$, and judge the maximum $|z|$ against the number
of entries being compared.

### 7.3 Independent Γ validation (Monte Carlo, no lavaan dependency)

For a small binary/ordinal system where repeated resampling is cheap, compute a Monte Carlo
estimate of $N \cdot \text{Var}(\hat s)$ directly from many replicated samples of a known
population, and compare it to Module 4's analytic Γ. This is a distinct check from Section 7.1 and
worth keeping separate, specifically because *matching lavaan is not sufficient if the goal is
validating the mathematics rather than reproducing another piece of software* — two
implementations can agree numerically while sharing the same bug, or while one has silently
dropped the threshold-propagation term from Section 4.5 in a way that happens not to move the
diagonal much. Because this check never touches lavaan at all, it doesn't have that failure mode.

**The Monte Carlo estimator must reproduce Module 2's exact finite-sample algorithm, not just draw
from the same population.** Same threshold estimator, same pairwise correlation estimator, same
zero-cell/sparse-cell treatment (Module 2's fail/warn default, or whatever correction mode is being
tested), same convergence and failure handling, same sample size, same category thresholds as the
production code. If the Monte Carlo harness quietly uses a simpler or different estimator (e.g. a
textbook tetrachoric formula instead of Module 2's actual bivariate-ML routine) it's validating a
different estimator's Γ, not the one being shipped — this check only has teeth if it's checking the
literal code path Module 4 depends on.

**Judge this check against a simulated reference distribution, not a tolerance.** The quantity
being compared is itself a random variable, so "is the gap small" is not a well-posed question
until you know how large the gap is *supposed* to be. The procedure: simulate the same comparison
$R'$ times under the known population, record the maximum absolute discrepancy each time, and ask
what percentile the observed value falls at. A value near the median is a pass; a value in the
extreme upper tail is a signal. This is strictly more informative than any fixed threshold, it
needs no judgement call about what number to pick, and it automatically adapts to $N$, $R$, the
number of categories and the number of variables — all of which move the expected gap.

Worked example of the procedure, from the lavaan cross-check (Module 0): four indicators at
$n = 5000$, observed maximum gap 0.0580; reference distribution from 300 replications had median
0.0565, IQR [0.0492, 0.0659] and a 90% interval of [0.0394, 0.0873]; observed value at the 55th
percentile, so: pass. As a cross-check on the reference distribution itself, the analytic standard
error of the single noisiest entry was 0.0650 — larger on its own than the observed maximum across
all 100 entries, which is the sanity check that the reference distribution is not itself wrong.

**Two things this check should test that a scalar tolerance cannot.** First, *where* the
discrepancy sits: by the Section 4.5 identity, the within-variable entries must be exact regardless
of $N$ or $R$, so any discrepancy concentrated there is structural rather than stochastic no matter
how small it is. Second, *what shape* it has: a constant factor across all entries, a transposed
block, or a discrepancy that fails to shrink as $R$ grows are all real mismatches that a
well-centred maximum could otherwise hide.

### 7.4 End-to-end model coverage

At least three population models (single-factor CFA, correlated two-factor CFA, a small SEM with a
structural path between factors), each simulated at a couple of sample sizes — this is the
population-model list Section 7.1's paired replication design runs against.

**Edge cases:** sparse contingency tables/empty cells for pairs with few categories, and binary
indicators (tetrachoric special case) specifically, since that's sub-phase 4a's validation target —
include these as their own conditions in the paired replication design (Section 7.1) rather than
as one-off checks, so their discrepancy patterns, not just pass/fail, are visible.

**The second and third population models are one data-generating process and two fits.** A path
$F_1 \to F_2$ with a disturbance on $F_2$ induces exactly a correlation between the factors, so the
observed data is distributionally identical to the correlated two-factor model. They differ in what
is *fitted*, not in what is generated, and generating a third replicate set would produce one
indistinguishable from the second while inviting the belief that a third process had been covered.

That also supplies a check that needs no external reference at all: fitted to the same data the two
models must agree on $npar$, $df$, $\sigma(\hat\theta)$ and $T_3$ (Section 4.8). Module 7 should
assert it directly, on its own output. It is asserted on lavaan's output in
`tests/unit/test_structural_metric.do` against `fixtures/expected/lavaan_case_E_equivalence.csv`,
and — since Module 7, 2026-09-23 — on our own fits in the same file's checks 12–15, for $npar$,
$df$, $\sigma(\hat\theta)$ and $T_{\text{raw}}$, and since Module 9 (2026-09-24) for $T_3$ in check 15.

### 7.5 Timing against lavaan, Mplus and `gsem` (reported, not a validation gate)

*Added 2026-09-26.* Timing belongs in Section 7 because it reuses Section 7's machinery: the same
designs, the same pinned lavaan options (`70_lavaan_mc_reference.R`, `20_lavaan_reference.R`) and the
same Mplus setup (`30_mplus_reference.R`). **It is not a gate.** A timing difference cannot make an
estimate wrong, and mixing performance into pass/fail criteria would blur what the gates mean.

- **What is timed.** Stage 1 (sample statistics and Γ, `dwls_sample_stats()`) and Stage 2 (fit and
  inference, `dwls_identify()` through `dwls_inference()`) are timed **separately**, because they
  scale differently. Stage 1's cost is known to be concentrated in Module 2's per-pair
  root-finding. Timing is wall clock over a stated number of repetitions, reported as the median,
  and excludes one-time loading and compilation (Mata definition files, lavaan's package load,
  Mplus start-up where it can be separated). *Settled 2026-09-28:* the separate stages are
  **dwlsmv's own**. lavaan, Mplus and `gsem` are timed as **totals**, through their public
  interfaces only, with no internal stage boundary; dwlsmv reports the same total beside its
  Stage 1 / Module 2 / Stage 2 breakdown. Splitting lavaan is a follow-up, built only if the
  totals show a gap on which §6's Module 2 decision turns.
- **Designs.** A grid over $p$, $K$ and $N$ that includes the size already quoted for
  `quadcross()`'s cost (E-072; $p = 40$, $K = 5$, $N = 10{,}000$), plus a few Monte Carlo conditions
  of Section 7.1's reference (for example `mc_d`, a P2 point and a b1 point), so that timings sit
  next to known statistical behaviour.
- **Fairness conditions.** Record each tool's processor use: StataMP's thread count, R/lavaan
  (largely single-threaded), and Mplus's own `PROCESSORS` setting. Record every version. Use
  identical data files, as the paired run already does. Report the development machine (6 cores,
  Stata 18) and the milestone machine (4 cores, Stata 16) separately, never pooled.
- **Where things live.** The driver goes in `tools/`, like `tools/time_fit.do`: a diagnostic, not
  a test. Results go in a dated performance note, and any decision they prompt goes in
  the development record. The first candidate was the Module 2 row of Section 6's risk matrix; it was
  settled on correctness grounds instead (`@sec-rootfind`, 2026-09-29), and the timing re-measures
  the new algorithm. **Nothing goes in the suite.**
- **A third comparator: `gsem` with ordinal probit links (`oprobit`), estimator versus
  estimator.** *Added 2026-09-27.* This is **not** a validation comparison, as lavaan's is. It
  answers the question Stata users will ask first, "why not `gsem` with ordinal probit?", and its
  numbers belong in the install page's motivation.
  - *Expected shape.* `gsem` is full-information ML and integrates over the factors numerically, so
    its cost grows with $N$ and with the number of factors (quadrature points raised to that
    power). dwlsmv's cost is dominated by Stage 1's per-pair work and hardly depends on the
    number of factors. The expectation is that `gsem` is competitive with one factor and falls
    far behind with several. This is a prediction to measure, not a result.
  - *Timing protocol.* The same as above: median wall clock over stated repetitions, loading
    excluded. For `gsem` only the total can be timed; dwlsmv keeps its two stages separate.
  - *Conditions to pin.* Complete data only, because `gsem`'s default for missing data differs
    from dwlsmv's listwise deletion. Record `intmethod()`, `intpoints()` and convergence status
    for every fit: they drive both `gsem`'s accuracy and its speed.
  - *What an estimate comparison can show.* **Never per-replicate agreement.** ML and DWLS are
    different estimators, so on a given dataset they differ by more than sampling noise, far
    outside Section 7.2's tiers. The fair check is statistical and optional. Over Section 7.1's
    Monte Carlo files, both should be close to unbiased for the true parameters. Their SDs can
    be compared as a rough efficiency comparison, with ML usually a little more efficient.
    Label any such result "estimator versus estimator".
  - *Parameterization.* `gsem`'s ordinal probit fixes the latent response's residual variance
    at 1, as Mplus's theta parameterization does. dwlsmv uses delta, with unit total variance.
    Loadings, thresholds and $\psi$ must be put on one scale before any comparison. The
    conversion is standard, but it is exactly where a sign or scale slip would hide. It
    therefore needs a derivation note in the `.qmd` with executable verification, or a check
    against lavaan's theta and delta output, before any number is read.
  - *No shared test statistic.* `gsem` gives no WLSMV-style $\chi^2$ for ordinal data, so $T_3$
    is not compared.

### 7.6 What the tests in this release establish

**What runs.** `tests/run_all.do` runs `tests/unit/test_harness.do` first, which checks the
comparator itself, then every `tests/unit/test_*.do`. Each test exits 0 (pass), 8 (a deliberate
skip) or 9 (failure). The result is the summary in `run_all.log`, ending `suite green` when nothing
failed; on Windows, Stata's process exit code is not the result.

**What the tests compare against.** Committed reference files in `fixtures/expected/`, read with
the tolerances of Section 7.2:

- **tier 1:** base-R mathematical references, such as closed forms and population quantities;
- **tier 2:** lavaan 0.7.2 on committed data: the sample statistics s, Γ̂, estimates, standard
  errors, T3, fit indices, standardized solutions, DIFFTEST, multiple-group fits and EFA;
- **tier 3:** Mplus, on one reference model.

Fixtures are written only by `fixtures/generate/*.R`, and each family has a manifest
(`fixtures/MANIFEST*.json`) recording its package versions and checksums; `fixtures/verify.R`
checks the files against the manifests. The tests also check internal contracts: the index map,
identification and its refusals, error paths (each refusal is asserted by its return code), and
the agreement of analytic derivatives with finite differences.

**What they do not establish.**

- Agreement with lavaan or Mplus is agreement on the committed models and data, not correctness in
  general: two implementations can agree while sharing a mistake (Section 7.3).
- The Monte Carlo designs of Sections 7.1 and 7.3 are not rerun by the suite.
  `fixtures/generate/50_replicates.R` regenerates their data (`fixtures/replication/`, not
  committed).
- Some references were recorded on one Stata version at one processor count, and those tests skip
  elsewhere (for example `test_sens_t3.do`). A skip is not a pass.
- `test_version_floor.do` checks that the Mata routines used exist in Stata 16, not that they behave
  the same. The release was run on Stata 16.1 and 18 on Windows only.
- The R chunks of `spec/derivations/dwlsmv-derivations.qmd` are a separate check. Rendering the
  file executes them (R and lavaan); the Stata suite does not.
- The preregistered studies and decisions cited as E-### (see *Evidence status* at the top) are
  not reproduced by any of the above.

## 8. Pre-implementation gates

Seven statistical design decisions that had to be settled on paper before substantial Mata
code existed, on the reasoning that letting an implementation infer them is how a subtle
estimator mismatch gets embedded in hundreds of lines rather than caught in a derivation.

**All seven are discharged**, as of 2026-09-17. The status table and
the record of when and why each closed are in the development record. What each gate
covered, and where its result now lives:

| Gate | Question | Result lives in |
|---|---|---|
| 0 | Data-type scope | §3 — ordinal/binary only for v1 |
| 1 | The sample-statistic vector and its index map | §4.2; derivation `@sec-gate1`; `statmap` locked in §5.1 |
| 2 | The model-implied vector σ(θ) | §4.6, §4.8; derivation `@sec-gate2` |
| 3 | The two-stage influence function for Γ, including threshold propagation | §4.5.2; derivation `@sec-general` |
| 4 | Γ's normalization | §4.4 — `lavInspect(fit, "gamma")`, never `wls.v` |
| 5 | The exact WLSMV variant | §4.9; derivation `@sec-gate5` |
| 6 | The locked Mata API | §5.1; `src/mata/dwlsmv_types.mata` |

Two results from that work are carried forward as standing cautions rather than history,
because both describe a way of being wrong that looks right:

- **The naive influence function is *also* exactly centred**, so §4.5.5's centring
  diagnostic cannot detect a missing propagation term. It catches derivative errors, not
  omissions. The check that does separate them is the pointwise expansion test — with the
  term present the residual's standard deviation times √N halves as N quadruples; with it
  dropped, that quantity is flat.
- **Gate 5's statistic comes from two sources, not one.** Muthén, du Toit & Spisic (1997)
  supplies U; the shift is Asparouhov & Muthén (2010). An instruction to derive the whole
  statistic from the first paper is impossible to follow — it does not contain one.
## 9. Starting a session

Onboarding lives in the development repository's instructions, which opens with a source-of-truth table naming
which document answers which question, and closes with a numbered workflow. Read that
rather than a prompt reproduced here; a second copy of the reading order is a second
copy that can go stale.

Two points specific to this document:

- **The §8 gates are closed and are not to be reopened.** Their results are recorded in
  the sections above and in `spec/derivations/dwlsmv-derivations.qmd`. If a gate's result
  looks wrong, say so and stop — do not quietly work around it or re-derive it.
- **Where this spec and the derivation file overlap, the derivation file wins.** This
  document summarizes it.

The history of how the prompt itself evolved — including why Gates 1, 3 and 5 were
sequenced first — is in the development record.
