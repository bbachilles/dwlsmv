# Executed 2026-09-26 under R 4.6.0 on the complete paired run (42 conditions).
# Results: E-071, 2026-09-26.
# ---------------------------------------------------------------------------
# dwlsmv -- summary of the production paired run (tools/paired_mc.do)
#
# A diagnostic tool, not a fixture generator: it writes nothing into the
# repository's tracked files. It reads
#   fixtures/replication/<condition>/production_fits.csv   (tools/paired_mc.do)
#   fixtures/replication/<condition>/lavaan_fits.csv       (70_, ref arm)
#   fixtures/expected/lavaan_mc_summary.csv                (committed, tier 2)
# and prints, per condition (decided with the user 2026-09-26):
#
#   A. STATUS: every production outcome counted -- ok, nonconverged,
#      empty_category, error by stage and rc. Nothing is excluded from the
#      count; the statistics below say which set they are over.
#   B. PRIMARY: production's own per-entry z_var, 70_'s statistic verbatim,
#      over production's converged fits (status ok), with SE/SD. Judged
#      against a two-sided Bonferroni bound over the entries of each tier.
#   C. SECONDARY: production's mean SE against lavaan's committed MC SD
#      (descriptive ratio), and the per-replicate paired differences against
#      lavaan on replicates both fitted, as the share inside spec 7.2's tiers
#      (estimates 1e-5 abs + 1e-4 rel, SEs 1e-8 abs + 1e-5 rel; architectural,
#      set 2026-09-24). Descriptive: lavaan's own convergence enters them.
#   D. MODULE 10 IN PRODUCTION: per condition, how often each flag fires, and
#      the centring flag against the replicates where production and lavaan
#      disagree beyond the SE tier -- does the row catch what it exists for?
#
# Also written, gitignored: fixtures/replication/production_summary.csv.
#
#   Rscript tools/paired_mc_summary.R
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

REP <- file.path("fixtures", "replication")
EXP <- file.path("fixtures", "expected")
rd  <- function(f) read.csv(f, colClasses = "character", sep = ",", check.names = FALSE)
num <- function(x) suppressWarnings(as.numeric(x))         # "." and ".n" become NA
bonf <- function(n) qnorm(1 - 0.05 / (2 * n))

paired_z <- function(th, se) {                              # 70_'s statistic, verbatim
  R <- length(th); if (R < 3L) return(c(z_var = NA, z_sd = NA))
  d <- se^2 - (th - mean(th))^2 * R / (R - 1)
  s <- sd(th)
  c(z_var = mean(d) / (sd(d) / sqrt(R)),
    z_sd  = (mean(se) - s) / (s / sqrt(2 * (R - 1))))
}

conds <- rd(file.path(EXP, "lavaan_mc_conditions.csv"))
conds <- conds[conds$arm == "ref", ]
lsum  <- rd(file.path(EXP, "lavaan_mc_summary.csv"))
lsum  <- lsum[lsum$arm == "ref", ]

have <- file.exists(file.path(REP, conds$condition, "production_fits.csv"))
if (!all(have)) cat("NOT YET RUN:", paste(conds$condition[!have], collapse = " "), "\n")
conds <- conds[have, ]

A <- list(); B <- list(); C <- list(); D <- list()
for (i in seq_len(nrow(conds))) {
  nm <- conds$condition[i]; tier <- conds$tier[i]
  p  <- rd(file.path(REP, nm, "production_fits.csv"))
  l  <- rd(file.path(REP, nm, "lavaan_fits.csv")); l <- l[l$arm == "ref", ]
  ok <- p$status == "ok"

  # ---- A. status ----------------------------------------------------------
  err <- p[p$status == "error", ]
  A[[nm]] <- data.frame(condition = nm, tier = tier, n = nrow(p), ok = sum(ok),
    nonconverged = sum(p$status == "nonconverged"), empty_category = sum(p$status == "empty_category"),
    error = nrow(err),
    error_detail = if (nrow(err)) paste(names(table(paste0(err$stage, ":", err$rc))),
                                        table(paste0(err$stage, ":", err$rc)), sep = "x", collapse = " ") else "",
    lavaan_ok = sum(l$status == "ok"), stringsAsFactors = FALSE)

  # ---- B. primary: production's own z_var -------------------------------------
  pars <- sub("^est:", "", grep("^est:", names(p), value = TRUE))
  for (k in pars) {
    e <- num(p[[paste0("est:", k)]][ok]); s <- num(p[[paste0("se:", k)]][ok])
    use <- is.finite(e) & is.finite(s); e <- e[use]; s <- s[use]
    z <- paired_z(e, s)
    lr <- lsum[lsum$condition == nm & lsum$parameter == k, ]
    B[[length(B) + 1L]] <- data.frame(condition = nm, tier = tier, parameter = k, n = length(e),
      mc_sd = sd(e), mean_se = mean(s), se_sd = mean(s) / sd(e),
      z_var = unname(z["z_var"]), z_sd = unname(z["z_sd"]),
      lavaan_mc_sd = num(lr$mc_sd), prod_se_over_lavaan_sd = mean(s) / num(lr$mc_sd),
      stringsAsFactors = FALSE)
  }

  # ---- C. paired, per replicate, on replicates both fitted --------------------
  m <- merge(p[ok, ], l[l$status == "ok", ], by = "replicate", suffixes = c(".p", ".l"))
  est_in <- se_in <- 0; worst_est <- worst_se <- 0
  bad_rep <- rep(FALSE, nrow(m))
  for (k in pars) {
    ep <- num(m[[paste0("est:", k, ".p")]]); el <- num(m[[paste0("est:", k, ".l")]])
    sp <- num(m[[paste0("se:", k, ".p")]]);  sl <- num(m[[paste0("se:", k, ".l")]])
    te <- abs(ep - el) <= 1e-5 + 1e-4 * abs(el)
    ts <- abs(sp - sl) <= 1e-8 + 1e-5 * abs(sl)
    bad_rep <- bad_rep | !ts
    worst_est <- max(worst_est, abs(ep - el) / (1e-5 + 1e-4 * abs(el)), na.rm = TRUE)
    worst_se  <- max(worst_se,  abs(sp - sl) / (1e-8 + 1e-5 * abs(sl)), na.rm = TRUE)
    est_in <- est_in + sum(te, na.rm = TRUE); se_in <- se_in + sum(ts, na.rm = TRUE)
  }
  ne <- nrow(m) * length(pars)
  C[[nm]] <- data.frame(condition = nm, tier = tier, n_paired = nrow(m),
    est_share_in_tier = est_in / ne, se_share_in_tier = se_in / ne,
    worst_est_over_tier = worst_est, worst_se_over_tier = worst_se,
    reps_se_outside_tier = sum(bad_rep), stringsAsFactors = FALSE)

  # ---- D. Module 10 in production ------------------------------------------
  f <- function(col) sum(num(p[[col]][ok]) == 1, na.rm = TRUE)
  cz <- num(m$g5f) == 1
  D[[nm]] <- data.frame(condition = nm, tier = tier, n_ok = sum(ok),
    asym = f("g1f"), notpsd = f("g2f"), fewpat = f("g3f"), alias = f("g4f"),
    centring = f("g5f"), bread_npar = f("b1f"), bread_zero = f("b2f"), bread_rank = f("b3f"),
    ratio_sf_med = median(num(p$g7v[ok]), na.rm = TRUE),
    bnd_count_min = min(num(p$c1v), na.rm = TRUE),
    se_mismatch_reps = sum(bad_rep), mismatch_with_centring_flag = sum(bad_rep & cz, na.rm = TRUE),
    centring_flag_without_mismatch = sum(!bad_rep & cz, na.rm = TRUE), stringsAsFactors = FALSE)
}

A <- do.call(rbind, A); B <- do.call(rbind, B); C <- do.call(rbind, C); D <- do.call(rbind, D)
options(width = 250)
cat("\n== A. Status (every replicate counted) ==\n"); print(A, row.names = FALSE)

cat("\n== B. PRIMARY: production's own z_var over its converged fits ==\n")
for (t in unique(B$tier)) {
  bt <- B[B$tier == t, ]; bd <- bonf(nrow(bt))
  cat(sprintf("\ntier %s: %d entries, Bonferroni bound |z| < %.3f\n", t, nrow(bt), bd))
  agg <- do.call(rbind, lapply(split(bt, bt$condition), function(x) {
    w <- which.max(abs(x$z_var))
    data.frame(condition = x$condition[1], n = x$n[1], max_abs_z = abs(x$z_var[w]), at = x$parameter[w],
               se_sd_min = min(x$se_sd), se_sd_max = max(x$se_sd), exceed = sum(abs(x$z_var) > bd),
               prod_se_over_lav_sd_min = min(x$prod_se_over_lavaan_sd),
               prod_se_over_lav_sd_max = max(x$prod_se_over_lavaan_sd))
  }))
  print(format(agg, digits = 4), row.names = FALSE)
}

cat("\n== C. SECONDARY: paired against lavaan, per replicate, spec 7.2 tiers (descriptive) ==\n")
print(format(C, digits = 4), row.names = FALSE)

cat("\n== D. Module 10's flags in production (counts over production's converged fits) ==\n")
print(format(D, digits = 4), row.names = FALSE)

write.csv(B, file.path(REP, "production_summary.csv"), row.names = FALSE)
cat("\nwritten (gitignored):", file.path(REP, "production_summary.csv"), "\n")
