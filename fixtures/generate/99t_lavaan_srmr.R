# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 2: lavaan references for the residual
# diagnostics (SRMR and the largest normalized residuals)
# (derivation @sec-resid-diag, Results SR1-SR5, O1-O7;
#  E-117, sections 3b, 3c)
#
# Drafted by a read-only agent (2026-11-22); edited by the main session for the
# owner's 2026-11-22 decisions 1(b) and 2(a) and the ranking-uncertainty rule,
# then run by the main session from the repository root on 2026-11-22.
#
# These files are frozen (committed) BEFORE any Mata code for the residual
# diagnostics exists; the Mata calculation is tested against them.
#
# WHAT THIS CERTIFIES. Agreement of dwlsmv's calculation with lavaan's
# categorical SRMR, sample/implied correlations and normalized residuals
# (tier 2). The comparator differs from dwlsmv in code and optimizer and SHARES
# the formulas; the formulas themselves are checked in the derivation's tier-1
# chunks. It certifies nothing about the ranking rule (that is dwlsmv's own,
# derived from SR2 and the owner's gate; lavaan has no ranking), about
# repaired or bound pairs being valid scales (SR4: they are not; the flags
# here only say which pairs are excluded), or about multiple groups (G = 1).
#
# VERSION. lavaan 0.7-2, the default library, asserted at run time. Why not
# 0.7-3 (the separate library 99p_ uses): this generator touches no rotation,
# no rotated SE and no information option. What it reads (srmr family,
# sample and implied correlations, Gamma, the normalized residual) is
# identical in the two versions (package section 1), and 0.7-2 is the version
# of every fit-index, unrotated-EFA and empty-cell fixture dwlsmv is already
# tested against (99f_, 99g_, 90_), so a residual compared here sits beside
# the same fit's T3 and sigma elsewhere. 0.7-3's differences (default bread
# "observed", the "bb" GPA algorithm, the std_ov expression) enter none of
# the quantities exported.
#
# FITS (options as 99f_/99g_: WLSMV, delta, scaled.shifted, std.lv = TRUE,
# zero.add = c(0.5, 0), zero.keep.margins = TRUE; the last two are the
# one-empty-cell repair dwlsmv adopted, spec Module 2):
#   E_cfa1   case_E  y1-y6   one factor
#   E_cfa2   case_E  y1-y6   two correlated factors (20_'s M_COR: f1 =~ y1-y3,
#                            f2 =~ y4-y6, f1 ~~ f2)
#   E_efa1   case_E  y1-y6   EFA m = 1 (the same model as E_cfa1; kept as its
#                            own tag because dwlsmv reaches it by another command)
#   E_efa2   case_E  y1-y6   EFA m = 2, df = 4
#   F_efa1-3 case_F  y1-y9   EFA m = 1, 2, 3 (df 27, 19, 12)
#   G_efa1,2 case_G  y1-y6   EFA m = 1, 2 (binary; the m = 2 case is the one
#                            99h_/99p_ use for case G), df 9 and 4
#   R_cfa1   emptycell y1-y4 one factor; y1~y2 has exactly one empty cell and
#                            is repaired by lavaan (the repaired-pair fixture)
#   E3_df0   case_E  y1-y3   one factor, df = 0 (O4)
#   G3_df0   case_G  y1-y3   one factor, binary, df = 0 (O4)
# ATTEMPTED, recorded in the fits file only, nothing exported from them
# (dwlsmv fails at these m on case_F, package 3a; whether lavaan does is the
# question): F_efa4, F_efa5.
#
# EFA IS FIT AS AN EXPLICIT ECHELON CFA (as 99g_): m(m-1)/2 anchor zeros
# (factor k has loadings 0 on items 1..k-1), factors uncorrelated, std.lv.
# Why: it is the model dwlsmv's identify layer is given, so the anchors are
# the same; and residuals are rotation-invariant (SR6: Lambda Lambda' is
# unchanged by any rotation), so any identification gives the same
# implied correlations. As a stop-the-line check each EFA fit is also run
# as lavaan's own efa("b") block with rotation = "none" and must agree on the
# implied correlations and on SRMR (tolerance below).
#
# EXPORTS
#   lavaan_srmr.csv  long: fit, quantity, i, j, item_i, item_j, value
#     Scalar quantities (i, j, item_i, item_j empty):
#       srmr, srmr_bentler_nomean, crmr_nomean   fitMeasures(), lavaan's own
#       srmr_hand, crmr_hand                     sqrt(sum e^2/(P+p)), sqrt(sum e^2/P)
#       N (ntotal), p, P, df, npar, n_thresholds
#       max_th_resid     largest |tau_sample - tau_implied| (lavaan's own convergence)
#       n_repaired, n_at_bound, n_small_e (|e| < 1e-5), n_eligible
#       top3_margin      (|z|(3) - |z|(4)) / |z|(3) among eligible pairs, df > 0 only:
#                        how far the third-ranked pair is from the fourth; a ranking
#                        compared in Stata is only as stable as this is large
#     Pair quantities, one row per pair i < j in item order (the spec's
#     column-major lower-triangle order is i outer, j inner, the same order):
#       r          sample polychoric (lavInspect "sampstat", the repaired value
#                  for a repaired pair)
#       rho        model-implied correlation (lavInspect "implied")
#       e          r - rho, observed minus implied, signed
#       gamma_kk   lavaan's Gamma at that polychoric entry (label "yi~~yj")
#       se_r       sqrt(gamma_kk / N), N = ntotal (not N - 1)
#       z_hand     e / se_r: dwlsmv's definition (SR2)
#       z_lavaan   residuals(type = "normalized"); equals e itself, NOT z, where
#                  |e| < 1e-5 (O7: lavaan divides by 1 there)
#       lav_div1   1 where |e| < 1e-5, i.e. where z_lavaan is not comparable to z_hand
#       nzero_cells  empty cells of the pair's contingency table
#       repaired   1 if both items binary and exactly one empty cell (lavaan's repair)
#       at_bound   1 if |r| >= 0.999 (spec Module 2's B). NONE in any committed
#                  dataset: no fixture exercises keepbound (see below)
#       eligible   1 - max(repaired, at_bound)
#       rank       rank of |z_hand| among eligible pairs (1 = largest), ties by
#                  smaller i then smaller j; empty where df = 0 (not ranked, O4)
#       rank_incl  the same ranking with ineligible pairs left in, so a test can
#                  show what the exclusion changed; empty where df = 0
#       top3       1 where rank <= 3 (empty where df = 0)
#   lavaan_srmr_fits.csv  fit, data, items, model, m, attempted, converged, df, N,
#     npar, admissible, lavaan_version, note. Every attempted fit, converged or not.
#     admissible is lavaan's post.check (no negative variance, no |r| > 1), for O3.
#
# Read with stringcols(_all) and strtoreal(), delimiters(","); empty = missing.
#
# PAIRS HELD AT THE BOUND UNDER keepbound: no committed dataset produces one.
# The only keepbound test (tests/unit/test_dwlsmv.do, section 4) builds its
# duplicate variable (q7 = y1) inside the test. lavaan has no keepbound: a
# perfectly associated pair is not a fit it reports an SRMR for. So no
# at-bound fixture exists, none is invented, and the at_bound column is
# a data-derived zero. The Mata test for that exclusion is hand-built (the
# package's step 3), not lavaan-referenced.
#
# STOP-THE-LINE CHECKS, ALL SET BEFORE THE FIRST RUN (tolerances taken from the
# derivation's sr-lavaan chunk, which set them before ITS first run, except where
# a line says otherwise):
#   [1 srmr]   hand srmr, crmr_nomean equal fitMeasures() to 1e-12 abs; srmr equals
#              srmr_bentler_nomean to 1e-12; raw e equals residuals(type = "raw")
#              to 1e-12; srmr_bentler^2 (T + P + p) = srmr^2 (P + p) to 1e-10 and
#              max_th_resid <= 1e-5 (thresholds drop out at G = 1, SR1)
#   [2 z]      z_hand equals z_lavaan to 1e-8 (relative to max(1, |z|)) at every
#              pair with |e| >= 1e-5 (O7); z computed with N - 1 must MISS by more
#              than 1e-6 on every fit with such a pair and df > 0 (N is ntotal)
#   [3 flags]  repaired pairs: exactly the pairs with one empty cell in a 2x2
#              table; no non-binary pair has an empty cell and no binary pair has
#              two (either would be a case the flags do not cover); R_cfa1 has
#              exactly one repaired pair, y1~y2, and every other fit none
#   [4 efa]    DIAGNOSTIC, PRINTED, NOT ASSERTED (owner, 2026-11-22). Each EFA fit against lavaan's
#              efa("b") rotation = "none": the max implied-rho gap and the srmr gap are printed. History:
#              the drafting agent first asserted rho within 1e-6 abs (99g_'s sigma gap); in its %TEMP%
#              run that FAILED at F_efa1 (rho gap 1.28e-6, srmr gap 1.1e-8: two lavaan fits of one
#              model stopping at lavaan's own floor at different points, same fmin), and it then
#              loosened the criterion. The owner chose neither the old nor the loosened tolerance:
#              this compares lavaan with lavaan, not dwlsmv with lavaan, so it is a diagnostic. The
#              deciding check is [1] (hand formula against lavaan's output).
#   [5 df0]    df = 0 fits: |srmr| < 1e-6 and max |e| < 1e-6 (a real solution);
#              set before the first run, a plausibility floor, not a tier
#   [6 df]     the df of every exported fit equals the one this header lists
#   [7 data]   every dataset complete (dwlsmv listwise-deletes; none applies)
# A failing check stops the script; no fixture value is moved to make one pass.
#
# TOLERANCES for the Stata test (recorded in the manifest; every one an
# existing tier, none invented):
#   rho (implied)          theta-hat tier, 1e-5 abs + 1e-4 rel on lavaan's value
#                          (spec 7.2; test_fit.do applies it to sigma)
#   r                      bivariate-correlation tier, 1e-8 abs + 1e-8 rel
#                          (spec 7.2; test_lavaan_emptycell.do check 1)
#   gamma_kk               end to end, 1e-8 abs + 1e-8 rel (test_gamma.do 5b)
#   PRINTED, NOT ASSERTED end to end where lavaan's own rho-hat misses the tier (the data
#   column of lavaan_srmr_fits.csv says which):
#     case_E    r and gamma_kk  test_gamma.do 5d (2026-09-25): lavaan's case E rho-hat is
#                               1.67e-8 from the root, r ratio 1.12
#     case_F    r and gamma_kk  owner 2026-11-22, after the cause check in
#                               E-117: dwlsmv's r within
#                               8.2e-14 of each pair's score root, lavaan's up to 3.3e-8 off;
#                               r ratio 2.83, gamma_kk 4.09 (measured in test_resid_diag.do)
#     emptycell gamma_kk only   test_lavaan_emptycell.do check 3; r passes and stays asserted
#   Where gamma_kk is printed, gamma_kk at lavaan's own rho-hat (these fixtures' r) is
#   asserted at 1e-12 abs + rel, spec 7.2's matched-rho row (as test_gamma.do 5c).
#   e = r - rho            not a tier of its own: its allowed error is the SUM of the two
#                          above, DE below. A tier on e relative to |e| would be tighter
#                          than the rho tier allows, since |e| << |rho|.
#   srmr, crmr             the bound the e tolerances imply (srmr_tol, crmr_tol, exported):
#                          by the triangle inequality on the Euclidean norm,
#                          |srmr_a - srmr_b| <= ||e_a - e_b|| / sqrt(P + p)
#                                            <= sqrt(sum DE^2 / (P + p)), and with P for crmr
#   N, p, P, df, nzero_cells, flags, rank exact
#   z                      NOT COMPARED DIRECTLY (owner, 2026-11-22, option 2(a)). No
#                          existing tier applies to a ratio with a denominator near
#                          0.01. The Stata test compares its parts: rho, r and gamma_kk.
#   rank / top3            asserted ONLY where the interval rule below says the order
#                          is determined; otherwise printed, not asserted.
#
# RANKING UNCERTAINTY (owner, 2026-11-22). Every tolerance below is an existing tier
# (owner's choice (a), 2026-11-22); none is new. A correct implementation may differ from
# lavaan in both parts of z = e / se_r:
#   numerator   |e_dwlsmv - e| <= DE = (TOL_RHO_ABS + TOL_RHO_REL |rho|)   theta-hat tier, on rho
#                                     + (TOL_G    + TOL_G |r|)            correlation tier, on r
#               the sum, not the larger: the two errors can add, since e = r - rho
#   denominator |gamma_dwlsmv - gamma_kk| <= TOL_G + TOL_G gamma_kk, so gamma lies within the
#               relative factors (1 - EPSG) and (1 + EPSG), EPSG = (TOL_G + TOL_G gamma_kk) / gamma_kk
# Each pair's |z| then lies in the worst-case interval, numerator and denominator perturbed
# together in the direction that widens it:
#   zlo = max(0, |e| - DE) / (se_r sqrt(1 + EPSG))
#   zhi =       (|e| + DE) / (se_r sqrt(1 - EPSG))
# PRINTED PARTS: where r or gamma_kk is not asserted (above), the tier is not a guaranteed
# bound. The Stata test therefore rebuilds those fits' intervals with the OBSERVED |r| and/or
# |gamma_kk| gaps of its own run in place of the printed parts' tier terms (a measured
# quantity, not a tolerance; the rho term stays at its asserted tier) and asserts the ranking
# only where those intervals are determined. The z_lo, z_hi exported here are the tier
# versions, for every fit.
# Because the denominator moves each pair's |z| by a different amount, it can reorder two
# pairs with close |z| even when their e agree exactly; the intervals carry that. Two pairs
# are in a determined order only if their intervals do not overlap: zlo(a) > zhi(b).
#   top3_set_determined    min zlo over the top 3 eligible > max zhi over the other eligible
#   top3_order_determined  additionally zlo(rank 1) > zhi(rank 2) and zlo(rank 2) > zhi(rank 3)
# Where either is 0 the test does not assert that part of the ranking. The rule is worst-case
# and per pair; it is not a probability, and it says nothing about sampling variation, only
# about two implementations computing the same quantity.
# History: the first version of this rule (commit bf1fc0e, never used by a Stata test) took a
# relative 1e-6 on gamma_kk, chosen after measuring test_gamma.do 5d. That was a new tolerance,
# not an existing tier; the owner chose the existing tier instead (2026-11-22). A sensitivity
# run on bf1fc0e's fixtures (main session, R, one-off) found every fit's set and order flags
# unchanged under the tier, under 1e-6, and under 1e-4 on gamma_kk with 1e-6 on r.
#
# Own manifest (MANIFEST_lavaan_srmr.json); fixtures/verify.R needs an entry
# (the main session adds it; this draft does not edit verify.R).
# Run from the repository root:  Rscript fixtures/generate/99t_lavaan_srmr.R
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))

if (!requireNamespace("lavaan", quietly = TRUE))
  stop("lavaan is required for tier-2 fixtures.")
library(lavaan)
# the default library's version, the one the fit-index and unrotated-EFA
# fixtures were generated with (see VERSION above); a different version stops
# the line instead of moving a fixture
stopifnot(as.character(packageVersion("lavaan")) == "0.7.2")
LAVAAN_VER <- as.character(packageVersion("lavaan"))

OUT <- file.path("fixtures", "expected")
MAN <- file.path("fixtures", "MANIFEST_lavaan_srmr.json")

B_BOUND   <- 0.999   # spec Module 2, the declared hard limit on |rho|
SMALL_E   <- 1e-5    # lavaan's normalized residual divides by 1 below this (O7)
TOL_HAND  <- 1e-12   # [1] hand vs fitMeasures / raw residuals
TOL_DIL   <- 1e-10   # [1] thresholds-drop-out identity
TOL_TH    <- 1e-5    # [1] lavaan's own threshold convergence
TOL_Z     <- 1e-8    # [2] z_hand vs z_lavaan, relative to max(1, |z|)
TOL_ZN1   <- 1e-6    # [2] the N - 1 control must miss by more than this
TOL_DF0   <- 1e-6    # [5]
TOL_RHO_ABS <- 1e-5  # ranking rule: theta-hat tier on rho (spec 7.2; test_efa_unrotated.do)
TOL_RHO_REL <- 1e-4
TOL_G       <- 1e-8  # ranking rule: correlation tier on r and end-to-end tier on gamma_kk,
                     # abs and rel (spec 7.2; test_lavaan_emptycell.do 1; test_gamma.do 5b)

rd <- function(f) {
  d <- read.csv(file.path("fixtures", "data", f))
  stopifnot(all(complete.cases(d)))                          # [7]
  d[] <- lapply(d, ordered)
  d
}
DATAFILES <- c("case_E.csv", "case_F.csv", "case_G.csv", "emptycell.csv")
dat <- list(E = rd("case_E.csv"), F = rd("case_F.csv"), G = rd("case_G.csv"),
            R = rd("emptycell.csv"))

yv <- function(p) sprintf("y%d", seq_len(p))

# --- model strings -----------------------------------------------------------
M_ONE <- function(p) sprintf("f =~ %s", paste(yv(p), collapse = " + "))
M_COR <- "f1 =~ y1 + y2 + y3
          f2 =~ y4 + y5 + y6
          f1 ~~ f2"                                          # 20_'s M_COR
M_ECH <- function(p, m) {                                     # echelon EFA, as 99g_
  L <- vapply(seq_len(m), function(k) {
    it <- yv(p); it[seq_len(k - 1)] <- paste0("0*", it[seq_len(k - 1)])
    sprintf("f%d =~ %s", k, paste(it, collapse = " + "))
  }, "")
  cv <- if (m > 1) unlist(lapply(2:m, function(a) sprintf("f%d ~~ 0*f%d", seq_len(a - 1), a))) else NULL
  paste(c(L, cv), collapse = "\n")
}
M_EFAB <- function(p, m)                                      # lavaan's own efa block
  sprintf('%s =~ %s', paste(sprintf('efa("b")*f%d', seq_len(m)), collapse = " + "),
          paste(yv(p), collapse = " + "))

fit_srmr <- function(model, d, ...)   # options as 99f_/99g_
  lavaan::sem(model, data = d, ordered = names(d),
              estimator = "WLSMV", parameterization = "delta",
              test = "scaled.shifted", zero.add = c(0.5, 0.0),
              zero.keep.margins = TRUE, std.lv = TRUE, ...)

#            tag       data  items model_type  m  model                       efa-check  required
specs <- list(
  E_cfa1 = list("E", 1:6, "cfa", 1, M_ONE(6),     FALSE, TRUE),
  E_cfa2 = list("E", 1:6, "cfa", 2, M_COR,        FALSE, TRUE),
  E_efa1 = list("E", 1:6, "efa", 1, M_ECH(6, 1),  TRUE,  TRUE),
  E_efa2 = list("E", 1:6, "efa", 2, M_ECH(6, 2),  TRUE,  TRUE),
  F_efa1 = list("F", 1:9, "efa", 1, M_ECH(9, 1),  TRUE,  TRUE),
  F_efa2 = list("F", 1:9, "efa", 2, M_ECH(9, 2),  TRUE,  TRUE),
  F_efa3 = list("F", 1:9, "efa", 3, M_ECH(9, 3),  TRUE,  TRUE),
  G_efa1 = list("G", 1:6, "efa", 1, M_ECH(6, 1),  TRUE,  TRUE),
  G_efa2 = list("G", 1:6, "efa", 2, M_ECH(6, 2),  TRUE,  TRUE),
  R_cfa1 = list("R", 1:4, "cfa", 1, M_ONE(4),     FALSE, TRUE),
  E3_df0 = list("E", 1:3, "cfa", 1, M_ONE(3),     FALSE, TRUE),
  G3_df0 = list("G", 1:3, "cfa", 1, M_ONE(3),     FALSE, TRUE),
  F_efa4 = list("F", 1:9, "efa", 4, M_ECH(9, 4),  FALSE, FALSE),   # attempted only
  F_efa5 = list("F", 1:9, "efa", 5, M_ECH(9, 5),  FALSE, FALSE))
DF_EXPECT <- c(E_cfa1 = 9, E_cfa2 = 8, E_efa1 = 9, E_efa2 = 4, F_efa1 = 27, F_efa2 = 19,
               F_efa3 = 12, G_efa1 = 9, G_efa2 = 4, R_cfa1 = 2, E3_df0 = 0, G3_df0 = 0)

quiet_fit <- function(model, d, ...) {
  warn <- character(0)
  fit <- tryCatch(withCallingHandlers(fit_srmr(model, d, ...),
           warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") }),
           error = function(e) { warn <<- c(warn, conditionMessage(e)); NULL })
  list(fit = fit, warn = unique(warn))
}

# --- hand quantities (derivation sr-lavaan's sr_hand, with the labels checked) ----
pair_stats <- function(fit, yn, d) {
  ss <- lavInspect(fit, "sampstat"); im <- lavInspect(fit, "implied")
  Ro <- ss$cov[yn, yn]; Rm <- im$cov[yn, yn]
  lt <- which(lower.tri(Ro), arr.ind = TRUE)                  # column-major, row > col
  i <- lt[, 2]; j <- lt[, 1]                                  # i < j
  ord <- order(i, j); i <- i[ord]; j <- j[ord]
  lab <- paste0(yn[i], "~~", yn[j])
  G <- unclass(lavInspect(fit, "gamma"))
  stopifnot(all(lab %in% rownames(G)))
  N <- lavInspect(fit, "ntotal")
  r <- Ro[cbind(j, i)]; rho <- Rm[cbind(j, i)]
  zl <- residuals(fit, type = "normalized")$cov[yn, yn][cbind(j, i)]
  rl <- residuals(fit, type = "raw")$cov[yn, yn][cbind(j, i)]
  nz <- vapply(seq_along(i), function(k) {
    t <- table(d[[yn[i[k]]]], d[[yn[j[k]]]]); c(zero = sum(t == 0),
      bin = as.integer(nrow(t) == 2 && ncol(t) == 2),
      binlev = as.integer(nlevels(d[[yn[i[k]]]]) == 2 && nlevels(d[[yn[j[k]]]]) == 2))
  }, c(zero = 0, bin = 0, binlev = 0))
  list(i = i, j = j, lab = lab, r = r, rho = rho, e = r - rho, gk = unname(diag(G)[match(lab, rownames(G))]),
       N = N, zl = zl, rl = rl, nzero = nz["zero", ], bin = nz["binlev", ],
       th = ss$th - im$th, p = length(yn))
}

rel_gap <- function(a, b) abs(a - b) / pmax(1, abs(b))

# --- fit everything -----------------------------------------------------------
fits <- list(); fitrows <- list()
for (tag in names(specs)) {
  s <- specs[[tag]]
  d <- dat[[s[[1]]]][, paste0("y", s[[2]]), drop = FALSE]
  r <- quiet_fit(s[[5]], d)
  fit <- r$fit
  conv <- !is.null(fit) && isTRUE(lavInspect(fit, "converged"))
  if (length(r$warn)) cat(sprintf("%-7s lavaan said: %s\n", tag, paste(r$warn, collapse = " | ")))
  fm <- if (conv) fitMeasures(fit, c("df", "npar", "ntotal")) else c(df = NA_real_, npar = NA_real_, ntotal = NA_real_)
  adm <- if (conv) as.integer(isTRUE(tryCatch(lavInspect(fit, "post.check"), error = function(e) NA))) else NA_integer_
  fitrows[[tag]] <- data.frame(fit = tag, data = if (s[[1]] == "R") "emptycell" else sprintf("case_%s", s[[1]]),
    items = sprintf("y%d-y%d", min(s[[2]]), max(s[[2]])), model = s[[3]], m = s[[4]],
    attempted = 1L, converged = as.integer(conv), df = fm[["df"]], N = fm[["ntotal"]],
    npar = fm[["npar"]], admissible = adm, lavaan_version = LAVAAN_VER,
    note = if (s[[7]]) "" else "attempted only; nothing exported")
  cat(sprintf("%-7s converged %d  df %s  N %s  npar %s\n", tag, conv, format(fm[["df"]]),
              format(fm[["ntotal"]]), format(fm[["npar"]])))
  if (conv && s[[7]]) fits[[tag]] <- list(fit = fit, d = d, spec = s)
}
req <- names(specs)[vapply(specs, function(s) s[[7]], TRUE)]
if (!all(req %in% names(fits))) stop("did not converge: ", paste(setdiff(req, names(fits)), collapse = ", "))

# --- rows ---------------------------------------------------------------------
rows <- list()
add_s <- function(tag, q, v) rows[[length(rows) + 1L]] <<- data.frame(fit = tag, quantity = q,
  i = NA_integer_, j = NA_integer_, item_i = NA_character_, item_j = NA_character_, value = v)
add_p <- function(tag, q, h, yn, v) rows[[length(rows) + 1L]] <<- data.frame(fit = tag, quantity = q,
  i = h$i, j = h$j, item_i = yn[h$i], item_j = yn[h$j], value = as.numeric(v))

summ <- list()
for (tag in names(fits)) {
  fit <- fits[[tag]]$fit; d <- fits[[tag]]$d; s <- fits[[tag]]$spec
  yn <- names(d); h <- pair_stats(fit, yn, d)
  P <- length(h$e); p <- h$p; N <- h$N
  fm <- fitMeasures(fit, c("srmr", "srmr_bentler_nomean", "srmr_bentler", "crmr_nomean", "df", "npar", "ntotal"))
  df <- fm[["df"]]
  stopifnot(P == p * (p - 1) / 2, !anyNA(h$e), !anyNA(h$gk), all(h$gk > 0),
            df == DF_EXPECT[[tag]])                                           # [6]
  # hand quantities
  srmr_h <- sqrt(sum(h$e^2) / (P + p)); crmr_h <- sqrt(sum(h$e^2) / P)
  g1 <- c(srmr = abs(srmr_h - fm[["srmr"]]), crmr = abs(crmr_h - fm[["crmr_nomean"]]),
          nomean = abs(fm[["srmr"]] - fm[["srmr_bentler_nomean"]]),
          raw = max(abs(h$e - h$rl)),
          dilute = abs(fm[["srmr_bentler"]]^2 * (length(h$th) + P + p) - fm[["srmr"]]^2 * (P + p)))
  mth <- max(abs(h$th))
  stopifnot(g1[c("srmr", "crmr", "nomean", "raw")] < TOL_HAND, g1["dilute"] < TOL_DIL, mth <= TOL_TH)  # [1]
  # z
  se_r <- sqrt(h$gk / N); z <- h$e / se_r; zN1 <- h$e / sqrt(h$gk / (N - 1))
  big <- abs(h$e) >= SMALL_E
  gz <- if (any(big)) max(rel_gap(z[big], h$zl[big])) else NA_real_
  gzN1 <- if (any(big)) max(rel_gap(zN1[big], h$zl[big])) else NA_real_
  if (any(big)) stopifnot(gz < TOL_Z)                                          # [2]
  if (any(big) && df > 0) stopifnot(gzN1 > TOL_ZN1)
  if (any(!big)) stopifnot(max(abs(h$zl[!big] - h$e[!big])) < TOL_HAND)         # O7 as found
  # flags
  nzero_nonbin <- sum(h$nzero[h$bin == 0]); nzero2 <- sum(h$nzero[h$bin == 1] >= 2)
  stopifnot(nzero_nonbin == 0, nzero2 == 0)                                    # [3]
  repaired <- as.integer(h$bin == 1 & h$nzero == 1)
  at_bound <- as.integer(abs(h$r) >= B_BOUND)
  elig <- 1L - pmax(repaired, at_bound)
  if (tag == "R_cfa1") stopifnot(sum(repaired) == 1, h$lab[repaired == 1] == "y1~~y2")
  else stopifnot(sum(repaired) == 0)
  # ranking, df > 0 only; ties by smaller i, then smaller j
  rank_e <- rank_a <- rep(NA_real_, P); top3 <- rep(NA_real_, P); margin <- NA_real_
  # ranking uncertainty: worst-case |z| interval per pair (header, RANKING UNCERTAINTY)
  DE   <- (TOL_RHO_ABS + TOL_RHO_REL * abs(h$rho)) + (TOL_G + TOL_G * abs(h$r))
  EPSG <- (TOL_G + TOL_G * h$gk) / h$gk
  zlo  <- pmax(0, abs(h$e) - DE) / (se_r * sqrt(1 + EPSG))
  zhi  <- (abs(h$e) + DE) / (se_r * sqrt(1 - EPSG))
  srmr_tol <- sqrt(sum(DE^2) / (P + p)); crmr_tol <- sqrt(sum(DE^2) / P)
  set_det <- ord_det <- NA_integer_
  if (df > 0) {
    o <- function(idx) idx[order(-abs(z[idx]), h$i[idx], h$j[idx])]
    ie <- which(elig == 1); oe <- o(ie); rank_e[oe] <- seq_along(oe)
    oa <- o(seq_len(P)); rank_a[oa] <- seq_len(P)
    top3 <- as.numeric(!is.na(rank_e) & rank_e <= 3); top3[is.na(rank_e)] <- 0
    za <- abs(z[oe]); margin <- if (length(za) >= 4) (za[3] - za[4]) / za[3] else NA_real_
    k3 <- oe[seq_len(min(3, length(oe)))]; rest <- setdiff(oe, k3)
    set_det <- as.integer(length(rest) == 0 || min(zlo[k3]) > max(zhi[rest]))
    ord_det <- as.integer(set_det == 1 &&
      all(vapply(seq_len(length(k3) - 1), function(k) zlo[k3[k]] > zhi[k3[k + 1]], TRUE)))
  } else {                                                                    # [5]
    stopifnot(abs(fm[["srmr"]]) < TOL_DF0, max(abs(h$e)) < TOL_DF0)
  }
  # lavaan's own EFA block, rotation none: DIAGNOSTIC, printed, not asserted     # [4]
  gefa <- c(rho = NA_real_, srmr = NA_real_)
  if (s[[6]]) {
    re <- quiet_fit(M_EFAB(p, s[[4]]), d, rotation = "none")
    if (is.null(re$fit) || !isTRUE(lavInspect(re$fit, "converged"))) {
      cat(sprintf("%-7s diagnostic [4]: lavaan's efa() block did not converge; no gap printed\n", tag))
    } else {
      he <- pair_stats(re$fit, yn, d)
      gefa <- c(rho = max(abs(h$rho - he$rho)), srmr = abs(fm[["srmr"]] - fitMeasures(re$fit, "srmr")[[1]]))
    }
  }
  # rows
  for (q in c("srmr", "srmr_bentler_nomean", "crmr_nomean")) add_s(tag, q, fm[[q]])
  add_s(tag, "srmr_hand", srmr_h); add_s(tag, "crmr_hand", crmr_h)
  add_s(tag, "N", N); add_s(tag, "p", p); add_s(tag, "P", P); add_s(tag, "df", df)
  add_s(tag, "npar", fm[["npar"]]); add_s(tag, "n_thresholds", length(h$th))
  add_s(tag, "max_th_resid", mth)
  add_s(tag, "n_repaired", sum(repaired)); add_s(tag, "n_at_bound", sum(at_bound))
  add_s(tag, "n_small_e", sum(!big)); add_s(tag, "n_eligible", sum(elig))
  if (df > 0) { add_s(tag, "top3_margin", margin)
                add_s(tag, "top3_set_determined", set_det); add_s(tag, "top3_order_determined", ord_det) }
  add_p(tag, "z_lo", h, yn, zlo); add_p(tag, "z_hi", h, yn, zhi); add_p(tag, "e_tol", h, yn, DE)
  add_s(tag, "srmr_tol", srmr_tol); add_s(tag, "crmr_tol", crmr_tol)
  add_p(tag, "r", h, yn, h$r); add_p(tag, "rho", h, yn, h$rho); add_p(tag, "e", h, yn, h$e)
  add_p(tag, "gamma_kk", h, yn, h$gk); add_p(tag, "se_r", h, yn, se_r)
  add_p(tag, "z_hand", h, yn, z); add_p(tag, "z_lavaan", h, yn, h$zl)
  add_p(tag, "lav_div1", h, yn, as.integer(!big))
  add_p(tag, "nzero_cells", h, yn, h$nzero); add_p(tag, "repaired", h, yn, repaired)
  add_p(tag, "at_bound", h, yn, at_bound); add_p(tag, "eligible", h, yn, elig)
  if (df > 0) { add_p(tag, "rank", h, yn, rank_e); add_p(tag, "rank_incl", h, yn, rank_a)
                add_p(tag, "top3", h, yn, top3) }
  topk <- if (df > 0) which(top3 == 1)[order(rank_e[top3 == 1])] else integer(0)
  incl_hit <- if (df > 0) sum(repaired == 1 & rank_a <= 3) else 0
  cat(sprintf("%-7s srmr %.7f (gaps %.1e)  P %d  repaired %d  small|e| %d  z_gap %s  N-1 control %s  efa(): %s\n  top3: %s  margin %s%s\n",
      tag, fm[["srmr"]], max(g1[c("srmr", "crmr", "raw")]), P, sum(repaired), sum(!big),
      format(gz, digits = 2), format(gzN1, digits = 2),
      if (s[[6]]) sprintf("rho gap %.1e srmr gap %.1e", gefa[1], gefa[2]) else "n/a",
      paste(sprintf("%s (e %+.4f z %+.2f)", h$lab[topk], h$e[topk], z[topk]), collapse = "; "),
      format(margin, digits = 3), if (incl_hit) "  [a repaired pair would be in the top 3 if kept]" else ""))
  if (df > 0) cat(sprintf("        ranking determined: set %d  order %d\n", set_det, ord_det))
}
long <- do.call(rbind, rows); rownames(long) <- NULL
stopifnot(!anyNA(long$value[!long$quantity %in% c("rank", "rank_incl", "top3", "top3_margin")]))

# --- write --------------------------------------------------------------------
files <- file.path(OUT, c("lavaan_srmr.csv", "lavaan_srmr_fits.csv"))
write_csv_lf(long, files[1], na = "")
fr <- do.call(rbind, fitrows); rownames(fr) <- NULL
write_csv_lf(fr, files[2], na = "")

dmd5 <- tools::md5sum(file.path("fixtures", "data", DATAFILES))
write_manifest(MAN, files,
  seeds = "none (deterministic given case_E/F/G.csv and emptycell.csv)",
  extra = list(packages = c(
    paste("lavaan", LAVAAN_VER, "(the default library, not 99p_'s lib-lavaan-0.7-3; asserted == 0.7.2 at run time)"),
    "generator: fixtures/generate/99t_lavaan_srmr.R",
    "WLSMV, delta, scaled.shifted, std.lv TRUE, zero.add c(0.5, 0), zero.keep.margins TRUE (as 99f_, 99g_)",
    "EFA fit as explicit echelon CFA; lavaan's efa('b') rotation 'none' gaps printed as a diagnostic, not asserted (owner 2026-11-22)",
    paste0("data md5: ", paste(sprintf("%s %s", basename(names(dmd5)), dmd5), collapse = "; ")),
    "hand-vs-fitMeasures stop-the-line: srmr, crmr, raw e 1e-12; z vs normalized 1e-8 rel where |e| >= 1e-5 (derivation sr-lavaan, set before running)",
    "Stata tolerances, all existing tiers (owner 2026-11-22, choice (a)): rho theta-hat tier 1e-5 abs + 1e-4 rel (spec 7.2; test_fit.do sigma); r 1e-8 abs + rel (test_lavaan_emptycell.do 1); gamma_kk 1e-8 abs + rel end to end (test_gamma.do 5b)",
    "printed, not asserted end to end (lavaan's rho-hat misses the tier): case_E r and gamma_kk (test_gamma.do 5d); case_F r and gamma_kk (owner 2026-11-22, cause check docs/notes/residual-diagnostics/caseF-cause); emptycell gamma_kk only (test_lavaan_emptycell.do 3); where gamma_kk is printed it is asserted at lavaan's rho-hat at 1e-12 (as test_gamma.do 5c)",
    "e: allowed error e_tol = rho tier + r tier (the sum); srmr, crmr: srmr_tol, crmr_tol = sqrt(sum e_tol^2 / (P + p)), / P (triangle inequality); z not compared directly (owner 2026-11-22, 2(a)); df, N, flags exact",
    "ranking: asserted only where top3_set_determined / top3_order_determined = 1; worst-case |z| interval per pair, numerator (e_tol) and denominator (gamma_kk within its tier) perturbed together; z_lo, z_hi exported; on case E the Stata test rebuilds the intervals from its observed r and gamma_kk gaps (generator header)",
    "keepbound: no committed dataset has a pair at the bound; at_bound is all 0 and no at-bound fixture exists",
    R.version.string)))

cat(sprintf("%d lavaan SRMR reference files written (%d fits attempted, %d exported); %s written.\n",
            length(files), length(specs), length(fits), MAN))
