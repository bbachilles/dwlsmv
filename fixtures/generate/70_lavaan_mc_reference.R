# Executed end to end under R 4.6.0 and lavaan 0.7.2: the 14 P1-era conditions
# on 2026-09-25, the other 28 on 2026-09-26 (in two runs; the second after
# bad.marker.crit was pinned, see FITTING ARMS), then the committed-fixture
# path and MANIFEST_mc.json. fixtures/verify.R passes on all three files. The
# pin was shown inert on a fit it does not switch: mc_a1n_40, refitted under
# it, gave a byte-identical slim table.
# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 2: lavaan's Monte Carlo SE reference
#
# E-054, Rank-deficient Gamma, step 3: "Step 3
# needs a reference, not a trend." For every mc_* condition of
# 50_replicates.R, fit lavaan to each replicate and compare the Monte Carlo SD
# of theta-hat with the mean fitted SE, per parameter, by a per-entry z
# (spec 7.3's approach; never a percentage -- a project rule).
#
# THE REFERENCE IS A PROPERTY OF THE ESTIMATOR AND IS ESTABLISHED BY LAVAAN
# ALONE. Nothing here depends on production, on Module 2's behaviour near
# rho = 1, or on its eventual fix. Production is compared later, on the same
# replicate files, as it stands then; its failures are reported, never
# excluded.
#
# WHAT IS COMMITTED (fixtures/expected/, tier 2, CRLF like the rest of tier 2)
#   lavaan_mc_conditions.csv  the generating specification and seed scheme of
#                             every condition and fitting arm
#   lavaan_mc_summary.csv     per condition x arm x parameter: population
#                             value, MC mean, MC SD, mean fitted SE, the paired
#                             variance-scale z, n successful, failure and
#                             warning counts. MC SD, mean SE and z are given
#                             over ALL converged fits and over the UNCAPPED
#                             ones (no rho-hat at the 0.999 bound), with the
#                             capped fraction, because at the boundary points
#                             (mc_c_999, mc_B_9975) the two differ in kind
#   lavaan_mc_condition_stats.csv  per condition x arm: the distribution
#                             (min, 1%, 5%, 25%, 50%) over converged fits of
#                             each candidate diagnostic -- min category count,
#                             min cell count, raw and scale-free Gamma ratios,
#                             and the centring z (its upper tail, and the
#                             count over the 1e-6 flag)
#                             -- so a count-based flag and a ratio-based flag
#                             can be set against the SE failure at ANY cutoff
#                             without refitting (E-060, 2026-09-25)
# WHAT IS NOT (gitignored, regenerable)
#   fixtures/replication/<condition>/lavaan_fits.csv   one slim row per
#       replicate and arm: condition, arm, replicate, seed, status, warning
#       flags, theta-hat, fitted SE, Gamma ratios, n_pat, and the OBSERVED
#       min category count, min cell count and number of empty cells, and
#       the centring measure z = max|psibar_k| / sqrt(Gamma_kk) (spec 5)
# The replicate .dta files are TRANSIENT (50_replicates.R's header): generate
# a condition, run this file on it, keep it for the production paired run,
# then delete it. The per-condition slim tables are what let this file build
# the summary after the data are gone.
#
# FITTING ARMS
#   ref         the reference fit, on EVERY condition: zero.add = c(0, 0) and
#               zero.keep.margins = TRUE. That is production's estimator,
#               which adds nothing to any cell. zero.add[1] acts only on 2x2
#               tables (lav_bvord.R) and zero.add[2] is 0 by default, so on a
#               condition with no pair of binary items this is identical to
#               lavaan's default.
#   lavdefault  every condition WITH a pair of binary items (a 2x2 table
#               exists), decided from the condition's thresholds, never from
#               its name: lavaan's default zero.add = c(0.5, 0), lavaan's own
#               result. Secondary tier (S); never a reference for our SEs.
#
# SINCE 2026-10-06 (O6: production repairs a one-empty-cell pair as lavaan's
# default does; spec Module 2, *One-empty-cell repair*) the ref arm's premise,
# "production's estimator, which adds nothing to any cell", holds only on
# replicates with NO binary pair whose 2x2 table has exactly one empty cell.
# On a replicate that has one, production equals lavdefault, not ref, and ref
# is no longer production's estimator. Nothing committed changes: the summary
# files were built before O6, and replicates without such a pair are
# unaffected. A rerun of the paired comparison on 2x2 conditions
# (mc_a2n_40, mc_a2n_60, ...) must take, per replicate, lavdefault where the
# replicate has a one-empty-cell pair, or the comparison measures the repair
# rather than the estimator. Identify those replicates from the data
# (dwlsmv_emptycell_pairs()'s rule), not from lavaan_fits.csv's empty-cell
# count, which counts every empty cell of every table. Not regenerated here;
# this header records the premise change.
# Identification is marker-variable, first indicator, as each condition's
# `ident` records; asserted on every fit's parameter table.
# bad.marker.crit = 0: lavaan 0.7.2 otherwise moves the marker off a "poor"
# first indicator (warning "switching to another marker item"), which fits a
# different parameter vector. Found 2026-09-26 on mc_a2n_40 replicate 86
# (lambda-hat_1 = 0.075, y1 -> y6) by that assertion. Pinned, as
# zero.keep.margins is: a reference must not depend on a default.
#
# REPLICATES THAT ARE NOT FITTED, AND WHY
#   status "empty_category": a declared category has no observations.
#   Production's Module 1 aborts; lavaan's ordered() would silently drop the
#   level and fit a different parameter vector. Counted, never fitted.
# Everything that IS fitted stays in the summary whatever happens to it:
# non-convergence is its own status and is counted; an inadmissible solution
# is counted and INCLUDED (Module 7 reports, never enforces). The all-converged
# columns are the PRIMARY result. *_admissible columns (inadmissible fits
# removed) sit alongside as a SENSITIVITY result, with n_admissible, as do
# *_uncapped (rho-hat at the 0.999 bound removed). Neither replaces the primary.
#
# THE PER-ENTRY STATISTIC. For parameter k over the R_k successful fits,
#   d_r = SE_r^2 - (theta_r - mean(theta))^2 * R_k / (R_k - 1)
#   z_var = mean(d) / (sd(d) / sqrt(R_k))
# This carries the Monte Carlo error of both the SD and the mean SE, and does
# not assume theta-hat is normal (b1's thresholds are Poisson-driven and
# heavily kurtotic, where SD / sqrt(2(R - 1)) is wrong). z_sd, the
# normal-theory (mean SE - SD) / (SD / sqrt(2 (R - 1))), is written alongside.
# Judge max |z| against the number of entries compared, not entry by entry.
#
# RUN ORDER. This file owns fixtures/MANIFEST_mc.json and nothing else, so it
# is independent of 20_ and 60_ and of MANIFEST_lavaan.json's merge hazard:
#   DWLSMV_CONDITIONS=mc_c_99 Rscript fixtures/generate/50_replicates.R
#   DWLSMV_CONDITIONS=mc_c_99 Rscript fixtures/generate/70_lavaan_mc_reference.R
#   ... per condition; delete fixtures/replication/<condition>/*.dta after the
#   production paired run ...
#   Rscript fixtures/generate/70_lavaan_mc_reference.R     (summarize only)
#
#   DWLSMV_CONDITIONS  fit only these conditions (default: none; summarize)
#   DWLSMV_WORKERS     PSOCK workers for fitting (default 1)
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))
for (pkg in c("lavaan", "haven", "pbivnorm"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop(pkg, " is required.")

REP  <- file.path("fixtures", "replication")
OUT  <- file.path("fixtures", "expected")
FIT  <- strsplit(Sys.getenv("DWLSMV_CONDITIONS", ""), ",")[[1]]
NWRK <- as.integer(Sys.getenv("DWLSMV_WORKERS", "1"))
CAP  <- 0.999                      # lavaan's (and production's) polychoric bound

popf <- file.path(REP, "population.rds")
if (!file.exists(popf)) stop("run fixtures/generate/50_replicates.R first: ", popf)
pop <- readRDS(popf)
pop <- pop[grepl("^mc_", pop$condition), ]
if (!nrow(pop)) stop("population.rds holds no mc_* condition; 50_replicates.R predates them.")
if (length(FIT) && !all(FIT %in% pop$condition))
  stop("unknown condition(s): ", paste(setdiff(FIT, pop$condition), collapse = ", "))

# A 2x2 table exists iff at least two items are binary (one threshold each).
# Keyed on the data-generating thresholds, not the condition's name, so a
# renamed or added binary condition cannot silently lose its reference arm.
has_2x2 <- function(nm) {
  thr <- pop$thresholds[[match(nm, pop$condition)]]
  sum(lengths(thr) == 1L) >= 2L
}
arms_of <- function(nm) if (has_2x2(nm)) c("ref", "lavdefault") else "ref"
zero_of <- function(nm, arm) if (arm == "ref") c(0, 0) else c(0.5, 0)
tier_of <- function(part, arm) if (arm == "lavdefault") "S" else part

# population theta under marker identification (lambda_1 = 1):
# thresholds as generated, lambda_j / lambda_1, psi = lambda_1^2
pop_theta <- function(row) {
  lam <- as.vector(row$loading[[1]]); thr <- row$thresholds[[1]]; p <- length(lam)
  th <- unlist(lapply(seq_len(p), function(j)
    setNames(thr[[j]], sprintf("y%d|t%d", j, seq_along(thr[[j]])))))
  c(setNames(lam[-1] / lam[1], sprintf("f=~y%d", 2:p)), th, "f~~f" = lam[1]^2)
}

# --- the centring row's axis (spec 5, Module 10) ------------------------------
# z = max_k |psibar_k| / sqrt(Gamma_kk). psibar_tau = 0 exactly; for a pair,
# psibar_rho = S(rho-hat) / (N I_rhorho) (derivation @sec-gamma-rank, Lemma 2),
# i.e. sum(W u) / sum(W u^2) with W the observed proportions and u = P_rho / P,
# at lavaan's thresholds (the sample ones) and lavaan's rho-hat. Gamma_kk is
# lavaan's own diagonal. pbivnorm (a lavaan import) gives the corners, which is
# ample for a flag at 1e-6; it is not the 1e-12 arithmetic of the Gamma tests.
centring_z <- function(Y, K, rho, gdiag) {
  N <- nrow(Y); p <- ncol(Y); T <- sum(K - 1L); k <- 0L; zmax <- 0
  th <- lapply(seq_len(p), function(j) qnorm(cumsum(tabulate(Y[, j] + 1L, K[j]))[-K[j]] / N))
  for (a in seq_len(p - 1L)) for (b in seq(a + 1L, p)) {
    k <- k + 1L; r <- rho[b, a]
    A <- pmin(pmax(c(-Inf, th[[a]], Inf), -40), 40); B <- pmin(pmax(c(-Inf, th[[b]], Inf), -40), 40)
    g <- expand.grid(x = A, y = B)
    F <- matrix(pbivnorm::pbivnorm(g$x, g$y, r), length(A))
    f <- matrix(exp(-(g$x^2 - 2 * r * g$x * g$y + g$y^2) / (2 * (1 - r^2))) /
                  (2 * pi * sqrt(1 - r^2)), length(A))
    cell <- function(M) M[-1, -1] - M[-nrow(M), -1] - M[-1, -ncol(M)] + M[-nrow(M), -ncol(M)]
    P <- cell(F); Pr <- cell(f)
    W <- table(factor(Y[, a], levels = 0:(K[a] - 1)), factor(Y[, b], levels = 0:(K[b] - 1))) / N
    use <- W > 0 & P > 0; u <- Pr[use] / P[use]
    zmax <- max(zmax, abs(sum(W[use] * u) / sum(W[use] * u^2)) / sqrt(gdiag[T + k]))
  }
  zmax
}

# --- one replicate ------------------------------------------------------------
fit_one <- function(file, K, zadd) {
  Y <- as.matrix(haven::read_dta(file))
  storage.mode(Y) <- "integer"
  cnt <- lapply(seq_len(ncol(Y)), function(j) tabulate(Y[, j] + 1L, K[j]))
  # The candidate count-based diagnostic, OBSERVED, over every declared level:
  # the smallest category count, and the smallest cell and number of empty
  # cells over every pair's full K_j x K_k table.
  cells <- unlist(lapply(seq_len(ncol(Y) - 1L), function(a)
    lapply(seq(a + 1L, ncol(Y)), function(b)
      as.vector(table(factor(Y[, a], levels = 0:(K[a] - 1)),
                      factor(Y[, b], levels = 0:(K[b] - 1)))))))
  # Threshold-block deficiency d_alias = T + 1 - rank[1 | Z] (derivation
  # @sec-gamma-rank, Theorem 3), counted as that section's Level 3 counts it.
  # qr() on a 0/1 matrix is a numerical rank, but it is the derivation's own
  # check, and on integer columns it has never been ambiguous here.
  Z <- do.call(cbind, lapply(seq_len(ncol(Y)), function(j)
    outer(Y[, j], seq_len(K[j] - 1L), "<") * 1))
  base <- list(n_pat = nrow(unique(Y)), d_alias = ncol(Z) + 1L - qr(cbind(1, Z))$rank,
               min_cat = min(unlist(cnt)),
               min_cell = min(cells), n_empty_cells = sum(cells == 0L))
  if (any(unlist(cnt) == 0L)) return(c(base, status = "empty_category"))
  d <- as.data.frame(lapply(seq_len(ncol(Y)), function(j) ordered(Y[, j], levels = 0:(K[j] - 1))))
  names(d) <- sprintf("y%d", seq_len(ncol(Y)))
  w <- character(0)
  fit <- withCallingHandlers(
    tryCatch(lavaan::sem(paste("f =~", paste(names(d), collapse = " + ")),
                         data = d, ordered = names(d), estimator = "WLSMV",
                         parameterization = "delta", test = "scaled.shifted",
                         zero.add = zadd, zero.keep.margins = TRUE,
                         bad.marker.crit = 0),
             error = function(e) e),
    warning = function(m) { w <<- c(w, conditionMessage(m)); invokeRestart("muffleWarning") })
  flags <- list(warn_nearly1 = any(grepl("nearly\\) 1", w)),
                warn_negvar  = any(grepl("variances are negative", w)),
                warn_other   = sum(!grepl("nearly\\) 1|variances are negative", w)))
  if (inherits(fit, "error")) return(c(base, flags, status = "error"))
  pt <- lavaan::parTable(fit)
  mk <- pt$op == "=~" & pt$rhs == "y1"
  if (!(sum(mk) == 1L && pt$free[mk] == 0L && pt$est[mk] == 1))
    stop("marker identification did not hold on ", file)
  if (!lavaan::lavInspect(fit, "converged")) return(c(base, flags, status = "nonconverged"))
  fr  <- pt$free > 0L                              # coef() is an S4 method; read the table
  est <- setNames(pt$est[fr], paste0(pt$lhs, pt$op, pt$rhs)[fr]); se <- setNames(pt$se[fr], names(est))
  G   <- lavaan::lavInspect(fit, "gamma")
  # The ratio rows (spec 5, Module 10) report sigma_{r*} / sigma_1 with
  # r* = min(D - d_alias, n_pat - 1): the smallest singular value that is NOT
  # structurally zero. Past r*, singular values are zero by Theorem 1 (too few
  # patterns) or Theorem 3 (a null vector of Gamma_tautau is one of Gamma,
  # since Gamma is PSD), and that exact deficiency is the structural rows' to
  # report. Folding it into the ratio would mix the two kinds of row.
  ks  <- min(ncol(G) - base$d_alias, base$n_pat - 1L)
  sv  <- function(M) { s <- svd(M, nu = 0, nv = 0)$d; s[ks] / s[1] }
  rho <- lavaan::lavInspect(fit, "sampstat")$cov
  c(base, flags, status = "ok", r_star = ks,
    ratio_raw  = sv(G),
    ratio_corr = sv(G / sqrt(outer(diag(G), diag(G)))),
    centring_z = centring_z(Y, K, rho, diag(G)),
    rho_cap    = any(abs(rho[lower.tri(rho)]) >= CAP - 1e-9),
    inadmissible = min(diag(lavaan::lavInspect(fit, "est")$theta)) < 0,
    as.list(setNames(est, paste0("est:", names(est)))),
    as.list(setNames(se,  paste0("se:",  names(se)))))
}

# --- fit the selected conditions, one slim table each --------------------------
for (nm in FIT) {
  row <- pop[pop$condition == nm, ]
  K   <- vapply(row$thresholds[[1]], length, integer(1)) + 1L
  idx <- read.csv(file.path(REP, nm, "index.csv"), stringsAsFactors = FALSE)
  if (nrow(idx) != row$nrep) stop(nm, ": index holds ", nrow(idx), " replicates, expected ", row$nrep)
  if (!all(file.exists(idx$file))) stop(nm, ": replicate files missing; regenerate with 50_replicates.R")
  rows <- list()
  for (arm in arms_of(nm)) {
    za <- zero_of(nm, arm)
    one <- function(i) c(list(condition = nm, arm = arm, replicate = idx$replicate[i],
                              seed = idx$seed[i]), fit_one(idx$file[i], K, za))
    if (NWRK > 1L) {
      cl <- parallel::makeCluster(NWRK)
      parallel::clusterExport(cl, c("fit_one", "centring_z", "idx", "K", "za", "nm", "arm", "CAP"))
      res <- tryCatch(parallel::parLapply(cl, seq_len(nrow(idx)), one),
                      finally = parallel::stopCluster(cl))
    } else res <- lapply(seq_len(nrow(idx)), one)
    rows <- c(rows, res)
  }
  cols <- unique(unlist(lapply(rows, names)))
  tab <- do.call(rbind, lapply(rows, function(r) {
    r <- r[cols]; names(r) <- cols
    as.data.frame(lapply(r, function(v) if (is.null(v)) NA else v),
                  check.names = FALSE, stringsAsFactors = FALSE) }))
  write.csv(tab, file.path(REP, nm, "lavaan_fits.csv"), row.names = FALSE)
  cat(sprintf("  %-14s %s\n", nm, paste(names(table(tab$status)), table(tab$status),
                                         sep = " ", collapse = ", ")))
}

# --- summarize every fitted condition ------------------------------------------
# The COMMITTED fixtures are written only once every condition has been fitted.
# Before that, the same summary is written to fixtures/replication/ (gitignored)
# as screening_*.csv, so a partition (P1 first) can be looked at before the rest
# is run. Screening files are never committed and never a reference.
have <- file.exists(file.path(REP, pop$condition, "lavaan_fits.csv"))
if (!any(have)) {
  cat("\nno condition fitted yet; nothing to summarize.\n")
  quit(save = "no", status = 0)
}

paired_z <- function(th, se) {
  R <- length(th); if (R < 3L) return(c(z_var = NA, z_sd = NA))
  d <- se^2 - (th - mean(th))^2 * R / (R - 1)
  s <- sd(th)
  c(z_var = mean(d) / (sd(d) / sqrt(R)),
    z_sd  = (mean(se) - s) / (s / sqrt(2 * (R - 1))))
}


qs <- function(x) {                                  # the lower tail is what a flag sees
  x <- x[is.finite(x)]
  if (!length(x)) return(setNames(rep(NA_real_, 5), c("min", "q01", "q05", "q25", "q50")))
  setNames(c(min(x), quantile(x, c(0.01, 0.05, 0.25, 0.50), names = FALSE, type = 7)),
           c("min", "q01", "q05", "q25", "q50"))
}

summ <- list(); spec <- list(); cstat <- list()
for (i in which(have)) {
  row <- pop[i, ]; nm <- row$condition; th0 <- pop_theta(row)
  tab <- read.csv(file.path(REP, nm, "lavaan_fits.csv"), check.names = FALSE,
                  stringsAsFactors = FALSE)
  for (arm in arms_of(nm)) {
    a   <- tab[tab$arm == arm, ]; ok <- a$status == "ok"
    cap <- ok & (a$rho_cap %in% TRUE); unc <- ok & !cap
    # Admissible-only is a SENSITIVITY result. The primary is every converged
    # fit, inadmissible ones included: Module 7 reports inadmissibility and
    # never enforces it, so production's SEs are produced for those fits too.
    adm <- ok & !(a$inadmissible %in% TRUE)
    za  <- zero_of(nm, arm); tr <- tier_of(row$partition, arm)
    spec[[length(spec) + 1L]] <- data.frame(
      condition = nm, arm = arm, tier = tr, n = row$n,
      nrep = row$nrep, seed_first = row$seed_first, seed_last = row$seed_first + row$nrep - 1L,
      seed_rule = "20260919 + position * 1000000 + replicate (50_replicates.R)",
      ident = row$ident, zero_add = paste(za, collapse = ";"), zero_keep_margins = TRUE,
      lambda = paste(format(as.vector(row$loading[[1]]), digits = 17), collapse = ";"),
      thresholds = paste(vapply(row$thresholds[[1]], function(t)
        paste(format(t, digits = 17), collapse = ";"), character(1)), collapse = "/"),
      note = gsub(",", ";", row$note), stringsAsFactors = FALSE)
    # the candidate diagnostics' distributions, for flag-vs-failure comparison
    st <- lapply(c("min_cat", "min_cell", "ratio_raw", "ratio_corr"), function(v) {
      # (centring_z's UPPER tail is the one a flag sees; it is summarized below)
      src <- if (v %in% c("min_cat", "min_cell")) a[[v]][a$status != "empty_category"] else a[[v]][ok]
      q <- qs(src); names(q) <- paste(v, names(q), sep = "_"); q })
    cstat[[length(cstat) + 1L]] <- data.frame(
      condition = nm, arm = arm, tier = tr, n_ok = sum(ok), n_capped = sum(cap),
      frac_capped = if (sum(ok)) sum(cap) / sum(ok) else NA_real_,
      n_inadmissible = sum(a$inadmissible[ok] %in% TRUE),
      as.list(unlist(st)),
      centring_z_q50 = median(a$centring_z[ok]), centring_z_max = max(a$centring_z[ok]),
      n_centring_z_gt_1e6 = sum(a$centring_z[ok] > 1e-6),
      med_n_pat = median(a$n_pat),
      med_n_empty_cells = median(a$n_empty_cells), stringsAsFactors = FALSE)
    for (k in names(th0)) {
      if (is.null(a[[paste0("est:", k)]])) stop(nm, ": no estimate column for ", k)
      one <- function(use) { e <- a[[paste0("est:", k)]][use]; s <- a[[paste0("se:", k)]][use]
        z <- paired_z(e, s)
        c(mc_mean = if (length(e)) mean(e) else NA, mc_sd = if (length(e) > 1) sd(e) else NA,
          mean_se = if (length(s)) mean(s) else NA, z_var = unname(z["z_var"]), z_sd = unname(z["z_sd"])) }
      al <- one(ok); un <- one(unc); ad <- one(adm)
      summ[[length(summ) + 1L]] <- data.frame(
        condition = nm, arm = arm, tier = tr, parameter = k, pop_value = unname(th0[k]),
        as.list(al), as.list(setNames(un, paste0(names(un), "_uncapped"))),
        as.list(setNames(ad, paste0(names(ad), "_admissible"))),
        n_ok = sum(ok), n_uncapped = sum(unc), n_admissible = sum(adm),
        frac_capped = if (sum(ok)) sum(cap) / sum(ok) else NA_real_,
        n_empty_category = sum(a$status == "empty_category"),
        n_nonconverged = sum(a$status == "nonconverged"), n_error = sum(a$status == "error"),
        n_inadmissible = sum(a$inadmissible[ok] %in% TRUE),
        n_warn_nearly1 = sum(a$warn_nearly1 %in% TRUE),
        n_warn_negvar = sum(a$warn_negvar %in% TRUE),
        n_warn_other = sum(a$warn_other > 0, na.rm = TRUE), stringsAsFactors = FALSE)
    }
  }
}

# --- screening only, while the set is incomplete ------------------------------
# Same tables, written under fixtures/replication/ (gitignored). The committed
# fixtures and MANIFEST_mc.json are not touched until every condition exists.
if (!all(have)) {
  for (x in list(list(spec, "screening_conditions.csv"), list(summ, "screening_summary.csv"),
                 list(cstat, "screening_condition_stats.csv")))
    write.csv(do.call(rbind, x[[1]]), file.path(REP, x[[2]]), row.names = FALSE)
  cat(sprintf(paste0("\n%d of %d conditions fitted. SCREENING files written to %s/",
                     "screening_*.csv (not fixtures, never committed).\n",
                     "The committed summary waits for: %s\n"),
              sum(have), nrow(pop), REP, paste(pop$condition[!have], collapse = " ")))
  quit(save = "no", status = 0)
}

# --- write (tier 2 convention: plain write.csv, CRLF, quote = FALSE) -----------
written <- character(0)
emit <- function(df, name) {                        # 60_'s guard, verbatim in intent
  chr <- vapply(df, function(v) is.character(v) || is.factor(v), logical(1))
  bad <- vapply(df[chr], function(v) any(grepl("[,\r\n\"]", as.character(v))), logical(1))
  if (any(bad)) stop(name, ": unquoted field contains a comma, quote or newline in ",
                     paste(names(df)[chr][bad], collapse = " "))
  f <- file.path(OUT, name)
  write.csv(df, f, row.names = FALSE, quote = FALSE)
  written <<- c(written, f)
}
num17 <- function(df) { for (v in names(df)) if (is.double(df[[v]]))
  df[[v]] <- format(df[[v]], digits = 17, trim = TRUE); df }   # read with stringcols(_all)
emit(do.call(rbind, spec), "lavaan_mc_conditions.csv")
emit(num17(do.call(rbind, summ)), "lavaan_mc_summary.csv")
emit(num17(do.call(rbind, cstat)), "lavaan_mc_condition_stats.csv")

# --- manifest: this file's own -------------------------------------------------
# Not merged into MANIFEST_lavaan.json: 20_ rewrites that file from its own
# list, and its single seeds field belongs to 20_ and 60_'s deterministic
# fixtures. fixtures/verify.R checks this manifest as its own optional entry.
write_manifest(file.path("fixtures", "MANIFEST_mc.json"), written,
  seeds = paste("mc replicates seeded 20260919 + position * 1000000 + replicate by",
                "50_replicates.R; each condition's range is in lavaan_mc_conditions.csv"),
  extra = list(packages = c(
    paste("lavaan", as.character(packageVersion("lavaan"))),
    paste("haven", as.character(packageVersion("haven"))),
    "se type: robust.sem", "zero.keep.margins = TRUE (pinned, as 20_)",
    R.version.string)))
cat(sprintf("\n%d Monte Carlo reference files written to fixtures/expected/; MANIFEST_mc.json written.\n",
            length(written)))
