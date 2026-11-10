# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 2: lavaan oracle for sigma(theta) and Delta
#
# Executed 2026-09-22 under R 4.6.0 with lavaan 0.7.2, on the development
# machine. One thing broke and it was worth breaking: the emit() guard below
# rejected the first label convention, `lambda[3,1]`, because the comma inside
# an unquoted field silently widens the CSV. See LABEL CONVENTIONS. Regenerated
# twice, byte-identical both times -- nothing here depends on an optimizer.
#
# What this produces, and why it is separate from 20_lavaan_reference.R
# --------------------------------------------------------------------
# 20_ fits models and records what lavaan ESTIMATED. This file never uses an
# estimate. It evaluates lavaan's model-implied correlation matrix R_y* and its
# Jacobian d sigma / d theta' at matrices WE STATE, with do.fit = FALSE, so
# every number here is a deterministic function of the stated input and not of
# any optimizer's stopping rule. That is what makes it usable as a Module 5/6
# oracle at machine precision rather than at the parameter-estimate tolerance
# tier of spec Section 7.2.
#
# WHAT `lavInspect(fit, "delta")` ACTUALLY IS -- measured, not assumed
# --------------------------------------------------------------------
# It is the JACOBIAN d sigma / d theta', not the delta-parameterization scaling.
# A previous session believed the opposite and believed the Jacobian was
# `lavaan:::computeDelta()`. Both halves are wrong under lavaan 0.7.2:
#
#   lavInspect(fit, "delta")     -> lav_inspect_delta() -> lav_model_delta()
#                                   D_s x n_free.  THE JACOBIAN.
#   lavInspect(fit, "mm.delta")  -> the model matrix named `delta`, p x 1,
#                                   the delta-parameterization scaling. Under
#                                   parameterization = "delta" with Theta
#                                   derived it is all ones.
#   lavaan:::computeDelta        -> DOES NOT EXIST in 0.7.2. It was the 0.6.x
#                                   name; 0.7.x calls it lav_model_delta().
#
# Evidence, reproduced by the FD check below on every configuration: on case D
# `lavInspect(fit,"delta")` is 16 x 14 and equals `lav_model_delta()[[1]]` to
# 0.0, and equals a central-difference Jacobian of lavaan's own
# lav_model_wls_est() to 4.0e-13 at h = 1e-4, degrading as h shrinks (6.6e-12 at
# 1e-5, 5.8e-10 at 1e-7) -- the signature of round-off, not of a wrong object.
# `mm.delta` is 4 x 1. Dimension alone would have been suggestive; the values
# are what settle it.
#
# THE INTERFACE IS MATRICES, NOT theta
# ------------------------------------
# The oracle's input is an explicit (Lambda, Psi, B) plus a threshold vector
# tau, written out in lavaan_modelimplied_input.csv. It is deliberately NOT a
# theta vector: theta's ordering is spec Section 4.8's and belongs to Module 3,
# which is a separate seam still being built. A consumer maps its own theta into
# these matrices and looks the answer up by label.
#
# LABEL CONVENTIONS -- key on these, never on position
# ----------------------------------------------------
# sigma rows (and Jacobian rows) carry the Section 4.2 label:
#
#     y<j>|t<c>      a threshold, variable then ascending threshold index
#     y<i>~~y<j>     a polychoric correlation, column-major lower triangle
#
# identical to the labels in lavaan_case_D_gamma.csv and case_*_gamma_tau.csv.
#
# Jacobian COLUMNS carry the matrix element differentiated with respect to:
#
#     lambda_<j>_<f>  row j = position of the indicator in the declared
#                     indicator order (the Section 4.2 variable order);
#                     col f = position of the factor in the declared factor
#                     order
#     beta_<t>_<f>    to-factor t, from-factor f, both in declared factor order
#     psi_<i>_<j>     LOWER TRIANGLE ONLY, i >= j, declared factor order.
#                     This is Psi, the DISTURBANCE (co)variance matrix. Phi_eta
#                     is derived (spec Section 4.8) and never appears.
#     y<j>|t<c>       a threshold parameter -- the same Section 4.2 label the
#                     row carries, because tau enters sigma verbatim and the
#                     threshold block of Delta is a literal identity
#
# Underscores, and not the obvious bracketed spelling. That spelling puts a COMMA inside a
# field of a comma-separated file that tier 2 writes with quote = FALSE, which
# does not fail -- it silently adds a column, and `import delimited` accepts
# the ragged result. That is how this file was first written and the emit()
# guard below is what caught it; the guard is kept rather than the lesson being
# left as a comment. Brackets would have survived quoting, but quoting a field
# whose content is a comma also feeds `import delimited`'s delimiter
# auto-detection, which this repository has already been bitten by once (see
# the pipe trap below). Underscores sidestep both.
#
# lavaan_modelimplied_input.csv keys its rows on the SAME element labels, so
# input and Jacobian join on one string.
#
# COLUMN ORDER IN THE JACOBIAN FILE CARRIES NO MEANING. It is lavaan's free
# parameter order, which is neither Section 4.8's theta order nor anything a
# consumer should reproduce. The file is long-format (fit, row, col, value)
# precisely so that a positional read is not even possible.
#
# TRAP, from the project's coding rules, and it bites this file harder than most: read these
# with `import delimited ..., delimiters(",") stringcols(_all)` and convert
# with strtoreal(). Every label column here contains `|`, so delimiter
# auto-detection can choose `|` and every lookup then misses at once; and the
# Jacobian is the one fixture in the repository that a `float` import would
# ruin quietly, since an exact 1.0 in the threshold block would survive while
# the correlation block lost its low mantissa bits.
#
# THE OVERPARAMETERIZED CONFIGURATION
# -----------------------------------
# `E_over` is case E's F2 <- F1 model with BOTH `beta_2_1` and `psi_2_1`
# free -- spec Section 4.8's worked error, expressible in the canonical form as
# the extra line `psi F1 F2 free`. Case E is the only fixture with an
# endogenous factor, which is the whole reason this can be asked here and
# nowhere else: in a correlated two-factor model both factors are exogenous and
# psi_2_1 is simply the factor covariance, freely estimated and not redundant
# with anything.
#
# Expected rank and deficiency, derived before measuring:
#
#   Phi = (I-B)^-1 Psi (I-B)^-T  with B_2_1 = beta, Psi = [[1, p],[p, 1]]
#   gives Phi_11 = 1, Phi_12 = beta + p, Phi_22 = beta^2 + 2 beta p + 1.
#   At p = 0:  d Phi_12 / d beta = d Phi_12 / d p = 1, and
#              d Phi_22 / d beta = 2 beta + 2p = 2 beta = d Phi_22 / d p.
#   The two columns of Delta coincide exactly, so
#
#     n_free = 20,  rank(Delta) = 19,  deficiency = 1,
#     null direction = +-1/sqrt(2) on (beta_2_1, psi_2_1), zero elsewhere,
#
#   which is exactly what spec Section 4.8 claims. Note the "at p = 0": the
#   determinant of that 2x2 block is -2p, so the collinearity is EXACT only at
#   psi_2_1 = 0. That is not a limitation -- Section 4.8's check is run at the
#   starting values and 0 is where a disturbance covariance starts -- but a
#   consumer that perturbs psi_2_1 off zero before checking will see the rank
#   come back and should not conclude the check is broken. The model is still
#   unidentified there; the redundancy has merely rotated into the F2 loadings.
#   E_over is therefore stated with psi_2_1 = 0 and beta_2_1 = 0.40, so the
#   collinearity is exact while beta is not.
#
# RUN ORDER, and a hazard worth stating
# -------------------------------------
#   Rscript fixtures/generate/20_lavaan_reference.R
#   Rscript fixtures/generate/60_lavaan_modelimplied.R     <- this file
#
# Both write fixtures/MANIFEST_lavaan.json, tier 2's manifest. 20_ rewrites it
# from its own file list alone, so running 20_ AFTER this file silently drops
# these six entries from the manifest and fixtures/verify.R stops checking
# them. This file merges instead: it re-reads the manifest, VERIFIES every
# pre-existing tier-2 hash still holds, and rewrites the union. So the safe
# order is 20_ then 60_, and if you ever run 20_ again, re-run 60_ after it.
#
# Tier 2 is Windows-generated and consistently CRLF (.gitattributes marks
# fixtures/expected/** as -text so git never renormalizes it). These files use
# plain write.csv(), like 20_ does, NOT tier 1's write_csv_lf().
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))

if (!requireNamespace("lavaan", quietly = TRUE))
  stop("lavaan is required for tier-2 fixtures.")
library(lavaan)

OUT <- file.path("fixtures", "expected")
written <- character(0)

# --- the Section 4.6 mapping, in base R, as an independent oracle ------------
# theta -> Sigma_y* -> R_y*. Under the delta parameterization Theta is DERIVED,
# Theta_jj = 1 - (Lambda Phi Lambda')_jj, so R_y* = Sigma_y* identically and
# the "standardization" of Section 4.6 is a no-op. This function exists so that
# lavaan's answer is checked against the specification rather than adopted from
# it, and so the finite-difference Jacobian below is OUR function's, not a
# second reading of lavaan's.
sigma_spec <- function(Lambda, Psi, B, tau) {
  m   <- nrow(B)
  A   <- solve(diag(m) - B)
  Phi <- A %*% Psi %*% t(A)
  Sig <- Lambda %*% Phi %*% t(Lambda)
  theta_diag <- 1 - diag(Sig)
  if (any(theta_diag <= 0))
    stop("stated matrices are inadmissible under delta: derived Theta_jj <= 0")
  R <- Sig
  diag(R) <- 1
  # lower.tri() in column-major order is exactly Section 4.2's correlation
  # ordering -- the same order 00_common.R's correlation_labels() emits.
  c(tau, R[lower.tri(R)])
}

# --- the stated configurations ----------------------------------------------
# Values are chosen distinct from one another on purpose. A design with two
# equal loadings cannot fail when a row and a column index are swapped, and an
# oracle that cannot fail is not an oracle.
CFG <- list(
  D_marker = list(
    case = "D",
    model = "f =~ y1 + y2 + y3 + y4",
    std.lv = FALSE,
    factors = "f",
    inds = paste0("y", 1:4),
    Lambda = matrix(c(1.00, 0.80, 0.70, 0.60), 4, 1),
    Psi    = matrix(0.50, 1, 1),
    B      = matrix(0, 1, 1),
    tau    = c(-0.60, 0.10, 0.80, -0.50, 0.20, 0.90, -0.40, 0.30, -0.30, 0.40),
    note   = "single factor; marker metric (lambda_1_1 fixed to 1); psi free"
  ),
  D_unitvar = list(
    case = "D",
    model = "f =~ y1 + y2 + y3 + y4",
    std.lv = TRUE,
    factors = "f",
    inds = paste0("y", 1:4),
    Lambda = matrix(c(0.75, 0.65, 0.55, 0.45), 4, 1),
    Psi    = matrix(1, 1, 1),
    B      = matrix(0, 1, 1),
    tau    = c(-0.60, 0.10, 0.80, -0.50, 0.20, 0.90, -0.40, 0.30, -0.30, 0.40),
    note   = "single factor; unit-variance metric (psi_1_1 fixed to 1); all loadings free"
  ),
  E_marker = list(
    case = "E",
    model = "f1 =~ y1 + y2 + y3\n f2 =~ y4 + y5 + y6\n f2 ~ f1",
    std.lv = FALSE,
    factors = c("f1", "f2"),
    inds = paste0("y", 1:6),
    Lambda = rbind(cbind(c(1.00, 0.65, 0.55), 0), cbind(0, c(1.00, 0.60, 0.50))),
    Psi    = diag(c(0.55, 0.40)),
    B      = matrix(c(0, 0.40, 0, 0), 2, 2),
    tau    = c(-0.30, 0.60, -0.25, 0.70, -0.20, 0.80,
               -0.35, 0.55, -0.15, 0.65, -0.40, 0.50),
    note   = "two factors; F2 <- F1; marker metric on both; psi_1_1 and psi_2_2 free"
  ),
  E_unitvar = list(
    case = "E",
    model = "f1 =~ y1 + y2 + y3\n f2 =~ y4 + y5 + y6\n f2 ~ f1",
    std.lv = TRUE,
    factors = c("f1", "f2"),
    inds = paste0("y", 1:6),
    Lambda = rbind(cbind(c(0.75, 0.65, 0.55), 0), cbind(0, c(0.70, 0.60, 0.50))),
    Psi    = diag(2),
    B      = matrix(c(0, 0.40, 0, 0), 2, 2),
    tau    = c(-0.30, 0.60, -0.25, 0.70, -0.20, 0.80,
               -0.35, 0.55, -0.15, 0.65, -0.40, 0.50),
    note   = paste("two factors; F2 <- F1; unit-variance metric: psi_2_2 = 1 is",
                   "the DISTURBANCE variance so implied Var(F2) = beta^2 + 1")
  ),
  E_over = list(
    case = "E",
    model = "f1 =~ y1 + y2 + y3\n f2 =~ y4 + y5 + y6\n f2 ~ f1\n f1 ~~ f2",
    std.lv = TRUE,
    factors = c("f1", "f2"),
    inds = paste0("y", 1:6),
    Lambda = rbind(cbind(c(0.75, 0.65, 0.55), 0), cbind(0, c(0.70, 0.60, 0.50))),
    Psi    = diag(2),                      # psi_2_1 = 0, and FREE -- see header
    B      = matrix(c(0, 0.40, 0, 0), 2, 2),
    tau    = c(-0.30, 0.60, -0.25, 0.70, -0.20, 0.80,
               -0.35, 0.55, -0.15, 0.65, -0.40, 0.50),
    note   = "DELIBERATELY OVERPARAMETERIZED: beta_2_1 and psi_2_1 both free -- rank deficiency 1"
  )
)

# Expected rank = number of free parameters, except for the configuration built
# to violate that. Written down here rather than read off the singular values,
# so the assertion can fail.
EXPECTED_RANK <- c(D_marker = 14L, D_unitvar = 14L,
                   E_marker = 19L, E_unitvar = 19L, E_over = 19L)

DATA <- list(
  D = local({ d <- read.csv(file.path("fixtures", "data", "case_D.csv"))
              d[] <- lapply(d, ordered); d }),
  E = local({ d <- read.csv(file.path("fixtures", "data", "case_E.csv"))
              d[] <- lapply(d, ordered); d })
)

# --- per-configuration evaluation -------------------------------------------
rows_input <- list(); rows_sigma <- list(); rows_jac <- list()
rows_sval  <- list(); rows_rank  <- list(); rows_null <- list()

for (tag in names(CFG)) {
  cfg <- CFG[[tag]]
  dat <- DATA[[cfg$case]]
  p   <- length(cfg$inds)
  m   <- length(cfg$factors)

  fit <- lavaan::sem(cfg$model, data = dat, ordered = names(dat),
                     estimator = "WLSMV", parameterization = "delta",
                     test = "scaled.shifted", zero.add = c(0.5, 0.0),
                     zero.keep.margins = TRUE,   # lavaan's effective default, pinned 2026-09-25
                     std.lv = cfg$std.lv, do.fit = FALSE)

  # --- map lavaan's free parameters onto matrix elements --------------------
  pt <- parTable(fit)
  fr <- pt[pt$free > 0, ]
  fr <- fr[order(fr$free), ]

  ix_ind <- function(v) match(v, cfg$inds)
  ix_fac <- function(v) match(v, cfg$factors)

  spec_of <- lapply(seq_len(nrow(fr)), function(k) {
    r <- fr[k, ]
    if (r$op == "=~") {
      i <- ix_ind(r$rhs); j <- ix_fac(r$lhs)
      list(mat = "lambda", i = i, j = j,
           label = sprintf("lambda_%d_%d", i, j),
           row_name = r$rhs, col_name = r$lhs)
    } else if (r$op == "~") {
      i <- ix_fac(r$lhs); j <- ix_fac(r$rhs)
      list(mat = "beta", i = i, j = j,
           label = sprintf("beta_%d_%d", i, j),
           row_name = r$lhs, col_name = r$rhs)
    } else if (r$op == "~~") {
      a <- ix_fac(r$lhs); b <- ix_fac(r$rhs)
      i <- max(a, b); j <- min(a, b)          # lower triangle, always
      list(mat = "psi", i = i, j = j,
           label = sprintf("psi_%d_%d", i, j),
           row_name = cfg$factors[i], col_name = cfg$factors[j])
    } else if (r$op == "|") {
      i <- ix_ind(r$lhs); c_ <- as.integer(sub("^t", "", r$rhs))
      list(mat = "tau", i = i, j = c_,
           label = sprintf("%s|%s", r$lhs, r$rhs),
           row_name = r$lhs, col_name = r$rhs)
    } else stop(tag, ": unexpected free parameter op '", r$op, "'")
  })
  if (any(vapply(spec_of, function(s) anyNA(c(s$i, s$j)), logical(1))))
    stop(tag, ": a free parameter did not map onto a declared indicator/factor")
  par_labels <- vapply(spec_of, `[[`, character(1), "label")
  if (anyDuplicated(par_labels))
    stop(tag, ": duplicated parameter label -- ",
         paste(par_labels[duplicated(par_labels)], collapse = " "))

  # --- gather/scatter between the stated matrices and lavaan's free vector ---
  gather <- function(Lambda, Psi, B, tau) {
    vapply(spec_of, function(s) switch(s$mat,
      lambda = Lambda[s$i, s$j],
      beta   = B[s$i, s$j],
      psi    = Psi[s$i, s$j],
      tau    = tau[(cumsum(c(0, tau_per_var))[s$i]) + s$j]), numeric(1))
  }
  scatter <- function(x) {
    Lambda <- cfg$Lambda; Psi <- cfg$Psi; B <- cfg$B; tau <- cfg$tau
    for (k in seq_along(spec_of)) {
      s <- spec_of[[k]]
      if (s$mat == "lambda") Lambda[s$i, s$j] <- x[k]
      else if (s$mat == "beta") B[s$i, s$j] <- x[k]
      else if (s$mat == "psi") { Psi[s$i, s$j] <- x[k]; Psi[s$j, s$i] <- x[k] }
      else tau[(cumsum(c(0, tau_per_var))[s$i]) + s$j] <- x[k]
    }
    list(Lambda = Lambda, Psi = Psi, B = B, tau = tau)
  }

  tau_per_var <- vapply(cfg$inds, function(v) nlevels(dat[[v]]) - 1L, integer(1))
  stopifnot(sum(tau_per_var) == length(cfg$tau))

  x0 <- gather(cfg$Lambda, cfg$Psi, cfg$B, cfg$tau)
  stopifnot(!anyNA(x0))

  mod <- lavaan:::lav_model_set_parameters(fit@Model, x = x0)

  # --- did lavaan actually receive the matrices we stated? ------------------
  # Including the FIXED entries: a marker lambda or a unit-variance psi is part
  # of the stated input, and if lavaan disagreed about one of them the oracle
  # would describe a model nobody asked for. Read the model matrices back.
  gl <- mod@GLIST
  stopifnot(identical(lavNames(fit, "ov"), cfg$inds),
            identical(lavNames(fit, "lv"), cfg$factors))
  B_lav <- if (is.null(gl$beta)) matrix(0, m, m) else gl$beta
  stopifnot(max(abs(gl$lambda - cfg$Lambda)) == 0,
            max(abs(gl$psi    - cfg$Psi))    == 0,
            max(abs(B_lav     - cfg$B))      == 0)
  # Theta is derived and the delta scaling is fixed at 1 -- spec Section 4.6.
  stopifnot(max(abs(as.numeric(gl$delta) - 1)) == 0)
  theta_implied <- 1 - diag(cfg$Lambda %*%
                     (solve(diag(m) - cfg$B) %*% cfg$Psi %*% t(solve(diag(m) - cfg$B))) %*%
                     t(cfg$Lambda))
  stopifnot(max(abs(diag(gl$theta) - theta_implied)) < 1e-14)

  # --- sigma(theta) and Delta ----------------------------------------------
  sig_lav <- as.numeric(lavaan:::lav_model_wls_est(mod)[[1]])
  Delta   <- lavaan:::lav_model_delta(mod)[[1]]
  s_labels <- rownames(lavInspect(fit, "delta"))
  stopifnot(length(s_labels) == length(sig_lav), nrow(Delta) == length(sig_lav),
            ncol(Delta) == length(par_labels))

  # Section 4.2 ordering, checked against 00_common.R rather than against
  # lavaan's own dimnames, so the two references have to agree.
  want_labels <- c(threshold_labels(unname(tau_per_var)), correlation_labels(p))
  if (!identical(s_labels, want_labels))
    stop(tag, ": lavaan's sigma ordering differs from Section 4.2.\n  lavaan: ",
         paste(s_labels, collapse = " "), "\n  spec:   ",
         paste(want_labels, collapse = " "))

  # --- check 1: lavaan's sigma IS the Section 4.6 mapping -------------------
  sig_spec <- sigma_spec(cfg$Lambda, cfg$Psi, cfg$B, cfg$tau)
  g_sig <- max(abs(sig_lav - sig_spec))
  if (g_sig > 1e-14) stop(tag, ": sigma disagrees with Section 4.6 at ", g_sig)

  # --- check 2: Delta IS the Jacobian, and our column labels are right ------
  # Finite differences of OUR sigma_spec(), scattered through the label map. If
  # a column were mislabelled -- lambda_3_1 where lambda_1_3 was meant, say --
  # this comparison is what fails. It is therefore a test of the labelling and
  # not only of lavaan.
  h <- 1e-5
  J <- vapply(seq_along(x0), function(k) {
    xp <- x0; xp[k] <- xp[k] + h
    xm <- x0; xm[k] <- xm[k] - h
    (do.call(sigma_spec, scatter(xp)) - do.call(sigma_spec, scatter(xm))) / (2 * h)
  }, numeric(length(sig_lav)))
  g_jac <- max(abs(Delta - J))
  if (g_jac > 1e-8) stop(tag, ": Delta is not d sigma/d theta' -- gap ", g_jac)

  # --- check 3: the threshold block of Delta is a LITERAL identity ----------
  # Spec Section 4.8 puts thresholds first in theta so that this block needs no
  # computation. Asserted at 0.0, not at a tolerance: tau enters sigma verbatim
  # and the delta scaling is 1, so anything else is a structural surprise.
  Tn   <- sum(tau_per_var)
  is_t <- vapply(spec_of, function(s) s$mat == "tau", logical(1))
  stopifnot(max(abs(Delta[seq_len(Tn), is_t] - diag(Tn))) == 0,
            max(abs(Delta[seq_len(Tn), !is_t])) == 0,
            max(abs(Delta[-seq_len(Tn), is_t])) == 0)

  # --- rank ----------------------------------------------------------------
  sv  <- svd(Delta)
  d   <- sv$d
  tol <- max(d) * 1e-10
  rk  <- sum(d > tol)
  if (rk != EXPECTED_RANK[[tag]])
    stop(tag, ": rank ", rk, " but expected ", EXPECTED_RANK[[tag]])

  cat(sprintf("%-10s  sigma %2d x free %2d   rank %2d (expected %2d)   ",
              tag, length(sig_lav), length(x0), rk, EXPECTED_RANK[[tag]]))
  cat(sprintf("sigma gap %.2e  Delta gap %.2e  sval %.4e .. %.4e\n",
              g_sig, g_jac, max(d), min(d)))

  # --- long-format emission -------------------------------------------------
  free_lab <- par_labels
  elem_row <- function(mat, i, j, label, row_name, col_name, value)
    data.frame(fit = tag, matrix = mat, element = label, row_name = row_name,
               col_name = col_name, value = value,
               status = if (label %in% free_lab) "free" else "fixed")

  inp <- list()
  for (i in seq_len(p)) for (j in seq_len(m))
    inp[[length(inp) + 1L]] <- elem_row("lambda", i, j,
      sprintf("lambda_%d_%d", i, j), cfg$inds[i], cfg$factors[j], cfg$Lambda[i, j])
  for (i in seq_len(m)) for (j in seq_len(i))
    inp[[length(inp) + 1L]] <- elem_row("psi", i, j,
      sprintf("psi_%d_%d", i, j), cfg$factors[i], cfg$factors[j], cfg$Psi[i, j])
  for (i in seq_len(m)) for (j in seq_len(m)) if (i != j)
    inp[[length(inp) + 1L]] <- elem_row("beta", i, j,
      sprintf("beta_%d_%d", i, j), cfg$factors[i], cfg$factors[j], cfg$B[i, j])
  tl <- threshold_labels(unname(tau_per_var))
  tv <- rep(seq_len(p), tau_per_var)
  tc <- unlist(lapply(tau_per_var, seq_len))
  for (k in seq_along(tl))
    inp[[length(inp) + 1L]] <- elem_row("tau", tv[k], tc[k], tl[k],
      cfg$inds[tv[k]], sprintf("t%d", tc[k]), cfg$tau[k])
  rows_input[[tag]] <- do.call(rbind, inp)

  rows_sigma[[tag]] <- data.frame(fit = tag, label = s_labels, value = sig_lav)

  rows_jac[[tag]] <- data.frame(
    fit   = tag,
    row   = rep(s_labels, times = length(par_labels)),
    col   = rep(par_labels, each = length(s_labels)),
    value = as.numeric(Delta))

  rows_sval[[tag]] <- data.frame(fit = tag, index = seq_along(d), svalue = d)

  rows_rank[[tag]] <- data.frame(
    fit = tag, case = cfg$case, n_sigma = length(sig_lav), n_free = length(x0),
    expected_rank = EXPECTED_RANK[[tag]], rank = rk,
    deficiency = length(x0) - rk,
    sval_max = max(d), sval_min = min(d), sval_ratio = min(d) / max(d),
    note = cfg$note)

  # --- the null space, where there is one ----------------------------------
  if (length(x0) - rk > 0L) {
    v <- sv$v[, which.min(d)]
    v <- v / sqrt(sum(v^2))
    if (v[which.max(abs(v))] < 0) v <- -v      # sign convention: largest entry +
    rows_null[[tag]] <- data.frame(fit = tag, param = par_labels, value = v)
    big <- par_labels[abs(v) > 1e-8]
    cat(sprintf("%-10s  null direction on: %s\n", tag, paste(big, collapse = " ")))
    # Spec Section 4.8's claim, asserted rather than described.
    stopifnot(setequal(big, c("beta_2_1", "psi_2_1")))
    stopifnot(abs(abs(v[par_labels == "beta_2_1"]) - 1 / sqrt(2)) < 1e-8,
              abs(abs(v[par_labels == "psi_2_1"])  - 1 / sqrt(2)) < 1e-8,
              max(abs(v[!(par_labels %in% big)])) < 1e-8)
    cat(sprintf("%-10s  singular-value gap: %.4e vs %.4e  (ratio %.2e)\n",
                tag, max(d), min(d), min(d) / max(d)))
  }
}

# --- write -------------------------------------------------------------------
# Plain write.csv(): tier 2 is Windows-generated CRLF throughout. Tier 1's
# write_csv_lf() is deliberately NOT used here -- mixing a single LF file into
# this directory is the float-import trap described above, running the other way.
#
# quote = FALSE is tier 2's convention, so no emitted FIELD may contain a comma
# or a newline -- it would silently add a column and the file would still
# import. This bit first, on the `note` column of the rank file, where three
# perfectly ordinary English commas turned a 11-column CSV into a ragged one
# that `import delimited` would have accepted without complaint. Guarded rather
# than fixed once.
emit <- function(df, name) {
  chr <- vapply(df, function(v) is.character(v) || is.factor(v), logical(1))
  bad <- vapply(df[chr], function(v) any(grepl("[,\r\n\"]", as.character(v))),
                logical(1))
  if (any(bad))
    stop(name, ": unquoted field contains a comma, quote or newline in column(s) ",
         paste(names(df)[chr][bad], collapse = " "),
         " -- write.csv(quote = FALSE) would emit a ragged file")
  if (any(grepl("[,\r\n\"]", names(df))))
    stop(name, ": a column NAME contains a comma or quote")
  f <- file.path(OUT, name)
  write.csv(df, f, row.names = FALSE, quote = FALSE)
  written <<- c(written, f)
}

emit(do.call(rbind, rows_input), "lavaan_modelimplied_input.csv")
emit(do.call(rbind, rows_sigma), "lavaan_modelimplied_sigma.csv")
emit(do.call(rbind, rows_jac),   "lavaan_modelimplied_jacobian.csv")
emit(do.call(rbind, rows_sval),  "lavaan_modelimplied_svals.csv")
emit(do.call(rbind, rows_rank),  "lavaan_modelimplied_rank.csv")
emit(do.call(rbind, rows_null),  "lavaan_modelimplied_nullvec.csv")

# --- manifest: merge into tier 2's, do not clobber it ------------------------
# 20_lavaan_reference.R rewrites MANIFEST_lavaan.json from its own file list.
# Re-read it, check every entry it already holds still hashes as recorded --
# blessing a drifted fixture here would be worse than failing -- and write the
# union.
man <- file.path("fixtures", "MANIFEST_lavaan.json")
prior <- character(0)
if (file.exists(man)) {
  lines <- readLines(man, warn = FALSE)
  hit <- regmatches(lines, regexec(
    '"(fixtures/[^"]+)"[[:space:]]*:[[:space:]]*"([0-9a-f]{32})"', lines))
  hit <- Filter(function(x) length(x) == 3L, hit)
  prior <- vapply(hit, `[`, character(1), 2L)
  want  <- vapply(hit, `[`, character(1), 3L)
  prior <- setdiff(prior, written)
  want  <- want[match(prior, vapply(hit, `[`, character(1), 2L))]
  if (length(prior)) {
    if (!all(file.exists(prior)))
      stop("MANIFEST_lavaan.json lists files that are missing: ",
           paste(prior[!file.exists(prior)], collapse = " "),
           "\n  Run fixtures/generate/20_lavaan_reference.R first.")
    got <- unname(tools::md5sum(prior))
    if (any(got != want))
      stop("a tier-2 fixture recorded in MANIFEST_lavaan.json has changed:\n  ",
           paste(prior[got != want], collapse = "\n  "),
           "\n  Refusing to rewrite the manifest around it. Regenerate tier 2",
           " (20_lavaan_reference.R) and re-run this file.")
  }
}

write_manifest(man, c(prior, written),
               seeds = "none (deterministic given case_D.csv, case_E.csv and the stated matrices)",
               extra = list(packages = c(
                 paste("lavaan", as.character(packageVersion("lavaan"))),
                 "se type: robust.sem",
                 R.version.string)))

cat(sprintf("\n%d model-implied reference files written; manifest now holds %d.\n",
            length(written), length(prior) + length(written)))
cat("sigma and Delta here are evaluated at STATED matrices with do.fit = FALSE,\n")
cat("so they are optimizer-independent and comparable at machine precision.\n")
