# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 2: lavaan references for second-order
# factors (release B, unit B4; design E-114)
#
# The metrics are the ones release B offers (B3.1, spec Module 12 *Release B*):
# V (G by its variance, F's by item markers; dwlsmv's metric(marker)) and U
# (every disturbance variance 1; metric(unitvariance), lavaan std.lv). Metric
# M (G by a latent marker) is not offered and has no fixture (B3.14). A model
# whose G has an item (MD) is item-anchored on G, as any factor with an item.
#
#   fit tag       data    items    model
#   H4_V, H4_U    case_H  y1-y12   G over F1-F4
#   H3_V, H3_U    case_H  y1-y9    G over F1-F3, inside the region
#   H3_cor        case_H  y1-y9    correlated F1-F3, marker (SO2.3 partner)
#   HD_M, HD_U    case_H  y1-y13   G over F1-F4 plus y13 (B0's MD); in HD_M
#                                  y13 is G's marker
#   I_V           case_I  y1-y9    G over F1-F3, IMPROPER: psi_F1 < 0
#   I_cor         case_I  y1-y9    correlated F1-F3, marker (SO2.3 partner)
#   I_U           case_I  y1-y9    attempted; expected NOT to converge (SO4.4);
#                                  recorded in the fits file only
#   J_V, J_U      case_J  y1-y18   G1 over F1-F3, G2 over F4-F6, G1 ~~ G2
#
#   lavaan_so_parameters.csv    fit, label, lhs, op, rhs, est, se, z, pvalue
#                               (parameterEstimates, as lavaan_case_E_parameters.csv)
#   lavaan_so_standardized.csv  fit, label, lhs, op, rhs, est_std, se
#                               (standardizedSolution(type = "std.all"))
#   lavaan_so_t3.csv            fit, quantity, value: T df trM trM2 a b T3,
#                               as 20_'s lavaan_t3.csv, rebuilt and asserted
#   lavaan_so_fits.csv          fit, data, items, converged, npar, df, chisq,
#                               chisq_scaled
#
# Only converged fits enter the first three files.
#
# TRAP (as 80_): lavaan's standardized SEs come from a forward-difference
# Jacobian, about 1.5e-8 relative against an analytic one. A comparison
# against that column must allow for it.
#
# Stop-the-line checks, all set BEFORE the first run (not after measuring);
# the optimizer-gap tolerances are 80_'s, which were chosen after measuring on
# Case E. They are grouped by the three conclusions the owner keeps apart
# (2026-11-19); passing one establishes nothing about the others, and each
# printed line carries its group's tag.
#   [1 parameterization] V and U are the same model after transformation:
#     - each V (or MD's item-anchored) fit equals its U fit mapped by
#       @sec-second-order SO3.1 (d_f = the U loading of f's marker; d_G = 1
#       for a G scaled by its variance): relative gap < 1e-5 on lambda,
#       beta, psi;
#     - standardized estimates agree across the two metrics, < 1e-5, SEs
#       < 1e-5. The estimate tolerance was 80_'s 1e-6 before the first run;
#       it is 1e-5, the map check's, CHOSEN AFTER MEASURING (owner,
#       2026-11-19, option (a)): Case J stopped at 1.87e-6, lavaan's own
#       stopping floor (fmin equal to 2.6e-12; no tighter rel.tol converges;
#       E-114, the B4 design).
#   [2 equivalence] G over three and the correlated partner, a different
#     model, have the same df and chi-square, < 1e-6 relative (SO2.3), on
#     H3 (proper) and I (improper).
#   [3 admissibility] the references are the cases they are meant to be:
#     H and J fits proper (every free psi > 0); I_V improper (psi_F1 < 0, its
#     standardized loading on G above 1). This certifies the FIXTURE only.
#     Whether dwlsmv classifies and reports them correctly is B5's test,
#     against these values.
#   [orientation] production's sign rule (B3.10) holds in every fit: in U
#     every F's first loading is positive, and for every G the standardized
#     paths sum to a positive value. A fixture-layout check, none of the three.
# A lavaan release that breaks one of these stops the script rather than
# moving a fixture.
#
# Own manifest (MANIFEST_lavaan_so.json); fixtures/verify.R has the entry.
# Run from the repository root:  Rscript fixtures/generate/99s_lavaan_second_order.R
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))

if (!requireNamespace("lavaan", quietly = TRUE))
  stop("lavaan is required for tier-2 fixtures.")
library(lavaan)

OUT <- file.path("fixtures", "expected")
MAN <- file.path("fixtures", "MANIFEST_lavaan_so.json")

read_case <- function(nm) {
  d <- read.csv(file.path("fixtures", "data", sprintf("case_%s.csv", nm)))
  d[] <- lapply(d, ordered)
  d
}
dat <- list(H = read_case("H"), I = read_case("I"), J = read_case("J"))

fit_so <- function(model, d, items, ...)   # options as 20_'s fit_E()
  lavaan::sem(model, data = d[, items], ordered = items,
              estimator = "WLSMV", parameterization = "delta",
              test = "scaled.shifted", zero.add = c(0.5, 0.0),
              zero.keep.margins = TRUE, ...)

y  <- function(i) sprintf("y%d", i)
fl <- function(k, first) sprintf("F%d =~ %s", k, paste(y(first + 0:2), collapse = " + "))
F_lines <- function(K) paste(vapply(seq_len(K), function(k) fl(k, 3 * k - 2), ""), collapse = "\n")

# V: G's first path freed, its variance fixed; the F's keep lavaan's first-item
# marker. U: plain syntax with std.lv = TRUE.
m_V  <- function(K) paste(F_lines(K),
  sprintf("G =~ NA*F1 + %s", paste(sprintf("F%d", 2:K), collapse = " + ")),
  "G ~~ 1*G", sep = "\n")
m_U  <- function(K) paste(F_lines(K),
  sprintf("G =~ %s", paste(sprintf("F%d", 1:K), collapse = " + ")), sep = "\n")
m_cor <- function(K) F_lines(K)
m_DM <- paste(F_lines(4), "G =~ y13 + F1 + F2 + F3 + F4", sep = "\n")  # y13 marks G
m_JV <- paste(F_lines(6), "G1 =~ NA*F1 + F2 + F3", "G2 =~ NA*F4 + F5 + F6",
              "G1 ~~ 1*G1", "G2 ~~ 1*G2", "G1 ~~ G2", sep = "\n")
m_JU <- paste(F_lines(6), "G1 =~ F1 + F2 + F3", "G2 =~ F4 + F5 + F6",
              "G1 ~~ G2", sep = "\n")

specs <- list(
  H4_V   = list("H", 1:12, m_V(4),  FALSE),
  H4_U   = list("H", 1:12, m_U(4),  TRUE),
  H3_V   = list("H", 1:9,  m_V(3),  FALSE),
  H3_U   = list("H", 1:9,  m_U(3),  TRUE),
  H3_cor = list("H", 1:9,  m_cor(3), FALSE),
  HD_M   = list("H", 1:13, m_DM,    FALSE),
  HD_U   = list("H", 1:13, paste(F_lines(4), "G =~ y13 + F1 + F2 + F3 + F4", sep = "\n"), TRUE),
  I_V    = list("I", 1:9,  m_V(3),  FALSE),
  I_cor  = list("I", 1:9,  m_cor(3), FALSE),
  I_U    = list("I", 1:9,  m_U(3),  TRUE),
  J_V    = list("J", 1:18, m_JV,    FALSE),
  J_U    = list("J", 1:18, m_JU,    TRUE))

fits <- list(); fitrows <- list()
for (tag in names(specs)) {
  s <- specs[[tag]]; items <- y(s[[2]])
  warn <- character(0)
  fit <- tryCatch(withCallingHandlers(
           fit_so(s[[3]], dat[[s[[1]]]], items, std.lv = s[[4]]),
           warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") }),
         error = function(e) { warn <<- c(warn, conditionMessage(e)); NULL })
  conv <- !is.null(fit) && isTRUE(lavInspect(fit, "converged"))
  if (length(warn)) cat(sprintf("%-7s lavaan said: %s\n", tag, paste(unique(warn), collapse = " | ")))
  fm <- if (conv) fitMeasures(fit, c("npar", "df", "chisq", "chisq.scaled")) else rep(NA_real_, 4)
  fitrows[[tag]] <- data.frame(fit = tag, data = paste0("case_", s[[1]]),
    items = sprintf("y%d-y%d", min(s[[2]]), max(s[[2]])), converged = as.integer(conv),
    npar = fm[1], df = fm[2], chisq = fm[3], chisq_scaled = fm[4])
  cat(sprintf("%-7s converged %d  df %s  chisq %s\n", tag, conv,
              format(fm[2]), format(fm[3], digits = 10)))
  if (conv) fits[[tag]] <- fit
}
must <- setdiff(names(specs), "I_U")
if (!all(must %in% names(fits)))
  stop("did not converge: ", paste(setdiff(must, names(fits)), collapse = ", "))
if ("I_U" %in% names(fits))
  cat("NOTE: I_U converged, against SO4.4's expectation; its rows are written. Look before committing.\n")

est <- function(tag) lavInspect(fits[[tag]], "est")
rel <- function(a, b) max(abs(a - b)) / max(abs(b))

# --- U -> target map (SO3.1). marker[f] = the item marking factor f in the
# target, or NA for a factor scaled by its variance (d = 1).
map_check <- function(src, tgt, marker) {
  u <- est(src); t <- est(tgt)
  lv <- colnames(u$lambda)
  stopifnot(identical(lv, colnames(t$lambda)))
  d <- vapply(lv, function(f) if (is.na(marker[[f]])) 1 else u$lambda[marker[[f]], f], 0)
  L <- sweep(u$lambda, 2, d, "/")
  B <- if (is.null(u$beta)) NULL else diag(d, length(d)) %*% u$beta %*% diag(1 / d, length(d))
  P <- diag(d, length(d)) %*% u$psi %*% diag(d, length(d))
  g <- c(lambda = rel(L, t$lambda), beta = rel(B, t$beta), psi = rel(P, t$psi))
  cat(sprintf("[1 parameterization] map %s -> %s: %s\n", src, tgt, paste(sprintf("%s %.1e", names(g), g), collapse = "  ")))
  stopifnot(all(g < 1e-5))
}
mk <- function(K, ng = 1, item_on_G = FALSE) {
  m <- setNames(as.list(y(3 * seq_len(K) - 2)), sprintf("F%d", seq_len(K)))
  G <- if (ng == 1) "G" else sprintf("G%d", seq_len(ng))
  for (g in G) m[[g]] <- if (item_on_G) "y13" else NA
  m
}
map_check("H4_U", "H4_V", mk(4))
map_check("H3_U", "H3_V", mk(3))
map_check("HD_U", "HD_M", mk(4, item_on_G = TRUE))
map_check("J_U",  "J_V",  mk(6, ng = 2))

# --- standardized solution
std <- do.call(rbind, Map(function(tag, fit) {
  s <- standardizedSolution(fit, type = "std.all")
  s <- s[s$op %in% c("=~", "~", "~~", "|"), ]
  data.frame(fit = tag, label = paste0(s$lhs, s$op, s$rhs),
             lhs = s$lhs, op = s$op, rhs = s$rhs, est_std = s$est.std, se = s$se)
}, names(fits), fits))
rownames(std) <- NULL
stopifnot(!anyNA(std$est_std), !anyNA(std$se))
for (pr in list(c("H4_U", "H4_V"), c("H3_U", "H3_V"), c("HD_U", "HD_M"), c("J_U", "J_V"))) {
  a <- std[std$fit == pr[1], ]; b <- std[std$fit == pr[2], ]
  stopifnot(setequal(a$label, b$label))
  b <- b[match(a$label, b$label), ]
  d_est <- max(abs(a$est_std - b$est_std)); d_se <- max(abs(a$se - b$se))
  cat(sprintf("[1 parameterization] std %s vs %s: est %.1e  se %.1e\n", pr[1], pr[2], d_est, d_se))
  stopifnot(d_est < 1e-5, d_se < 1e-5)
}

# --- SO2.3: G over three has the correlated model's chi-square
for (pr in list(c("H3_V", "H3_cor"), c("I_V", "I_cor"))) {
  c1 <- fitMeasures(fits[[pr[1]]], c("chisq", "df")); c2 <- fitMeasures(fits[[pr[2]]], c("chisq", "df"))
  cat(sprintf("[2 equivalence] SO2.3 %s vs %s: df %g %g, chisq rel gap %.1e\n", pr[1], pr[2], c1[2], c2[2], abs(c1[1] / c2[1] - 1)))
  stopifnot(c1[2] == c2[2], abs(c1[1] / c2[1] - 1) < 1e-6)
}

# --- orientation (B3.10) and properness
pe <- function(tag) parameterEstimates(fits[[tag]])
for (tag in names(fits)) {
  st <- std[std$fit == tag, ]
  for (g in intersect(c("G", "G1", "G2"), st$lhs)) {
    paths <- st[st$lhs == g & st$op == "=~" & grepl("^F", st$rhs), "est_std"]
    stopifnot(sum(paths) > 0)
  }
  if (grepl("_U$", tag)) {
    p <- pe(tag)
    for (f in unique(p$lhs[p$op == "=~" & grepl("^F", p$lhs)]))
      stopifnot(p$est[p$lhs == f & p$op == "=~"][1] > 0)
  }
  p <- pe(tag)
  vars <- p[p$op == "~~" & p$lhs == p$rhs & p$lhs %in% c("G", "G1", "G2", sprintf("F%d", 1:6)) & is.na(match(p$se, 0)), ]
  if (tag == "I_V") {
    b1 <- st$est_std[st$label == "G=~F1"]
    cat(sprintf("[3 admissibility] I_V: psi_F1 %.4f (se %.4f), standardized loading %.4f\n",
                vars$est[vars$lhs == "F1"], vars$se[vars$lhs == "F1"], b1))
    stopifnot(vars$est[vars$lhs == "F1"] < 0, b1 > 1)
  } else if (!grepl("^I_", tag)) {
    if (nrow(vars) == 0) {
      # U fixes every factor's psi at 1: nothing to check here; U cannot
      # represent psi <= 0 at all (SO3.3), which is not a check of the data
      stopifnot(grepl("_U$", tag))
      cat(sprintf("[3 admissibility] %s: no free factor variance (every psi fixed at 1)\n", tag))
    } else {
      cat(sprintf("[3 admissibility] %s: smallest free psi %.4f\n", tag, min(vars$est)))
      stopifnot(all(vars$est > 0))
    }
  }
}

# --- T3 components, rebuilt as 20_ does
trace_noU <- function(Gm, Dl) {
  w  <- 1 / diag(Gm)
  WD <- w * Dl
  Ei <- solve(crossprod(Dl, WD))
  G  <- Gm %*% WD
  B  <- crossprod(WD, G)
  C  <- crossprod(G, w * G)
  EB <- Ei %*% B
  c(trM  = sum(w * diag(Gm)) - sum(Ei * B),
    trM2 = drop(crossprod(w, (Gm * Gm) %*% w)) - 2 * sum(Ei * C) + sum(EB * t(EB)))
}
rows_t3 <- list()
for (tag in names(fits)) {
  fa <- fits[[tag]]; ss <- fa@test$scaled.shifted
  stopifnot(!is.null(ss), ss$df == fa@test$standard$df)
  tr <- trace_noU(lavTech(fa, "gamma")[[1]], lavTech(fa, "delta")[[1]])
  a  <- 1 / ss$scaling.factor; b <- ss$shift.parameter; Tr <- fa@test$standard$stat
  g  <- c(trM  = abs(tr[["trM"]]  / ss$trace.UGamma  - 1),
          trM2 = abs(tr[["trM2"]] / ss$trace.UGamma2 - 1),
          a    = abs(sqrt(ss$df / tr[["trM2"]]) / a - 1),
          b    = abs(ss$df - sqrt(ss$df / tr[["trM2"]]) * tr[["trM"]] - b),
          T3   = abs((a * Tr + b) / ss$stat - 1))
  if (any(g > 1e-12))
    stop(tag, ": lavaan's scaled.shifted is not the @sec-gate5-traces expansion: ",
         paste(sprintf("%s %.2e", names(g), g), collapse = ", "))
  rows_t3[[tag]] <- data.frame(fit = tag,
    quantity = c("T", "df", "trM", "trM2", "a", "b", "T3"),
    value = c(Tr, ss$df, ss$trace.UGamma, ss$trace.UGamma2, a, b, ss$stat))
}

# --- write
par <- do.call(rbind, lapply(names(fits), function(tag) {
  p <- pe(tag)
  data.frame(fit = tag, label = paste0(p$lhs, p$op, p$rhs), lhs = p$lhs, op = p$op,
             rhs = p$rhs, est = p$est, se = p$se, z = p$z, pvalue = p$pvalue)
}))
files <- file.path(OUT, c("lavaan_so_parameters.csv", "lavaan_so_standardized.csv",
                          "lavaan_so_t3.csv", "lavaan_so_fits.csv"))
write_csv_lf(par, files[1])
write_csv_lf(std, files[2])
write_csv_lf(do.call(rbind, rows_t3), files[3])
fr <- do.call(rbind, fitrows); rownames(fr) <- NULL
write_csv_lf(fr, files[4])

write_manifest(MAN, files,
  seeds = "none (deterministic given case_H/I/J.csv)",
  extra = list(packages = c(
    paste("lavaan", as.character(packageVersion("lavaan"))),
    paste("se type:", lavInspect(fits[[1]], "options")$se),
    "standardizedSolution(type = 'std.all'); se by lavaan's forward-difference Jacobian",
    "zero.keep.margins = TRUE (pinned, as 20_)",
    R.version.string)))

cat(sprintf("%d second-order reference files written; %s written.\n", length(files), MAN))
