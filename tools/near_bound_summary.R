# tools/near_bound_summary.R -- summary of tools/near_bound_pairs.do output.
# Run after the tool, from the repository root: Rscript tools/near_bound_summary.R
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek
#
# Summaries fixed with the user
# 2026-10-03 BEFORE the run: per condition, status 1 and status 6 separately,
# count and min / 0.1% / 1% / 5% / 25% / 50% quantiles of
#   dist = B - |rho|                       (forms a, b)
#   z    = (B - |rho|) / sqrt(gamma_rr / N) (form c; Gamma is N x ACOV, so the
#                                           variance is gamma_rr / N)
#   |S/N| at sign(rho) * B, usable pairs only (form d)
# The usable share is reported BEFORE the score quantiles; unusable is
# structurally undefined, rows with k missing (Module 1 aborted) are counted
# separately as not computed. No thresholds.

dir <- file.path(Sys.getenv("TEMP"), "gammaC", "nb", "pairs")
files <- sort(list.files(dir, pattern = "\\.csv$", full.names = TRUE))
probs <- c(0, 0.001, 0.01, 0.05, 0.25, 0.5)
qn <- c("min", "q001", "q01", "q05", "q25", "q50")

qs <- function(x) {
  if (length(x) == 0) return(setNames(rep(NA_real_, 6), qn))
  setNames(quantile(x, probs, names = FALSE, type = 7), qn)
}

out <- list()
for (f in files) {
  d <- read.csv(f, colClasses = "character")
  num <- c("replicate", "k", "B", "rho", "status", "N", "gamma_rr",
           "score_usable", "score_SN")
  for (v in num) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
  cond <- d$condition[1]
  nreps <- length(unique(d$replicate))
  notcomp <- sum(is.na(d$k))
  p <- d[!is.na(d$k), ]
  for (st in c(1, 6)) {
    s <- p[p$status == st, ]
    if (nrow(s) == 0) next
    dist <- s$B - abs(s$rho)
    hasg <- !is.na(s$gamma_rr)
    z <- dist[hasg] / sqrt(s$gamma_rr[hasg] / s$N[hasg])
    us <- s$score_usable == 1
    sc <- abs(s$score_SN[us])
    out[[length(out) + 1]] <- data.frame(
      condition = cond, status = st, reps = nreps, reps_not_computed = notcomp,
      pairs = nrow(s), at_B = sum(dist == 0),
      t(setNames(qs(dist), paste0("dist_", qn))),
      z_n = sum(hasg), t(setNames(qs(z), paste0("z_", qn))),
      usable_share = mean(us),
      t(setNames(qs(sc), paste0("absSN_", qn))))
  }
  others <- table(p$status[!(p$status %in% c(1, 6))])
  if (length(others)) cat(cond, "other statuses:", paste(names(others), others, sep = "=", collapse = " "), "\n")
}
res <- do.call(rbind, out)
write.csv(res, file.path(Sys.getenv("TEMP"), "gammaC", "nb", "nb_summary.csv"), row.names = FALSE)
options(width = 250)
print(res[, c("condition", "status", "reps", "reps_not_computed", "pairs", "at_B",
              "dist_min", "dist_q01", "dist_q05", "dist_q50", "z_n", "z_min", "z_q01",
              "z_q05", "z_q50", "usable_share", "absSN_min", "absSN_q05", "absSN_q50")],
      digits = 3, row.names = FALSE)

# ---- the decided rule's counts (appended 2026-10-03) ----
# Pairs the near-bound note would flag under the rule decided 2026-10-03
# (chosen after measuring): |rho-hat| >= 0.99, pairs entering s with bound = 0,
# i.e. status 1 or 6. Per condition: pairs flagged, replicates with >= 1 flag.
dir <- file.path(Sys.getenv("TEMP"), "gammaC", "nb", "pairs")
for (f in sort(list.files(dir, pattern = "\\.csv$", full.names = TRUE))) {
  d <- read.csv(f, colClasses = "character")
  d <- d[d$k != ".", ]
  rho <- as.numeric(d$rho)
  st <- d$status
  fl <- st %in% c("1", "6") & abs(rho) >= 0.99
  if (sum(fl) > 0) {
    cat(sprintf("%-11s flagged pairs %6d of %6d; replicates with a flag %4d of %4d\n",
                d$condition[1], sum(fl), sum(st %in% c("1", "6")),
                length(unique(d$replicate[fl])), length(unique(d$replicate))))
  }
}
cat("conditions not listed flag nothing\n")
