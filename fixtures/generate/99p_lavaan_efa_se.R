# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 2: rotated-EFA standard errors (release A, A4)
#
# Drafted by a read-only agent (2026-11-17). Executed by the main session from
# the repository root 2026-11-17 (R 4.6.0, lavaan 0.7-3 from the separate
# library): all three cases converged with no warning, constraint residual at
# most 2.7e-13; output identical to the agent's scratch run; verify.R clean.
#
# lavaan's efa(rotation = "geomin") rotated standard errors and the FULL rotated
# covariance matrix, on the same cases and fit options as 99h_ (Cases E, F, G;
# oblique geomin, epsilon by Mplus's rule, row_weights none, 30 seeded starts,
# order_lv_by "index"), with these differences, each recorded in the manifest:
#
#   std_ov    FALSE. A COMPARATOR SETTING, not a dwlsmv convention (spec Module 12,
#             Release A, A3.10). In lavaan 0.7-2, std_ov TRUE rescales the rotated
#             loadings in the bordered constraint by diag(L L' + Theta), which for an
#             oblique solution is not diag(L Phi L' + Theta) (= 1 under the delta
#             parameterization); the constraint residual at lavaan's own estimate was
#             then 6e-3..1e-2 instead of ~1e-7 (E-109, comparator).
#             lavaan 0.7-3 (commit 55786b9) changed that expression; FALSE is kept
#             because it is settled, and the grid in the A4 report separates the
#             two std_ov settings per version. Estimates do not depend on the option.
#   version   lavaan 0.7-3 from a separate library (owner, 2026-11-17); 99h_'s
#             fixture was generated with 0.7-2. The rotated POINT can differ between
#             versions (0.7-3 adds the "bb" GPA algorithm as default), so these
#             estimates are NOT assumed equal to lavaan_efa_geomin.csv's.
#   bread     information.bread = "expected". lavaan 0.7-3 changed the default bread
#             of robust.sem to "observed" (option information.bread, new in 0.7-3).
#             Measured on this project's Cases (diagnostic grid): at the same point
#             and options, 0.7-3's default differs from "expected" by up to 1.4e-2
#             relative in the SEs and 3.8e-2 in the covariance (unrotated CFA too:
#             it is not a rotation effect), and "expected" reproduces 0.7-2's SEs
#             and the hand-computed (D'WD)^-1 D'W Gamma W D (D'WD)^-1 exactly.
#             Set because the spec's estimator uses the expected (DWLS) bread; this
#             is a comparator setting recorded here and in the manifest.
#   se        lavaan's default rotation.se = "bordered" (the delta route's step is
#             hard-wired and is not used).
#   tightness gpa_tol, gpa_algorithm, max_iter and the Jacobian; see below.
#
# TIGHT SETTINGS (spec 7.2 "Rotated EFA (release A ...)": "fixed tight, set and
# recorded by A4"). Chosen on CONVERGENCE grounds only (gpnorm reached, no
# warning, the constraint residual of the selected solution), never on
# agreement with any dwlsmv quantity. lavaan 0.7-3 offers:
#   gpa_tol   accepted range ]0, 0.01]. GPA stops when the Frobenius norm of the
#             projected gradient is below gpa_tol; it warns "GP rotation algorithm
#             did not converge" at max_iter. Experiment (Cases E, F, G; this
#             file's options, seed per case): "bb" converges with no warning at
#             1e-6, 1e-8, 1e-10, 1e-12, 1e-14 (max constraint residual 6.8e-8 at
#             1e-6, 4e-11 at 1e-10, 3e-13 at 1e-12, 3e-15 at 1e-14). The "legacy"
#             algorithm (0.7-2's only one) converges at 1e-6 (residual 3e-7..4e-7)
#             but does not finish within minutes at 1e-8 (stall floor ~1e-7).
#             GPA_TOL = 1e-12: converged with margin, not at the arithmetic floor.
#   gpa_algorithm "bb" (0.7-3 default) set explicitly: the only algorithm that
#             can be driven tight. gpa_fwindow left at its default 0.
#   max_iter  10000 (default), not reached at GPA_TOL with "bb".
#   tol       the PAIRWISE algorithm's tolerance; unused with algorithm "gpa".
#   rstarts   30, seeded per case, as 99h_. The selected start is the minimum
#             criterion; the SEs are conditional on it (comparator.md item 5).
#   Jacobian  the bordered Jacobian is numDeriv::jacobian(method = "Richardson")
#             with numDeriv's defaults, hard-wired in lav_step16_rotation and NOT
#             settable through a lavaan option. This generator MEASURES it: the
#             Jacobian is recomputed outside lavaan with a finer Richardson
#             setting (d 1e-5, r 6) and the largest absolute Jacobian difference
#             is recorded (quantity "diag"). It decides nothing and no fixture
#             value is taken from it.
#   rotation.se "delta" (0.7-3: central differences at a tightened tolerance) is
#             not used here; it is compared in the diagnostic grid only.
#
#   Case E  six 3-category items, m = 2      (fixtures/data/case_E.csv)
#   Case F  nine 4-category items, m = 3     (fixtures/data/case_F.csv)
#   Case G  six binary items, m = 2          (fixtures/data/case_G.csv)
#   Further cases (near-zero loading, near-tie in column order, empty cell) are
#   added to CASES only if the A4 rules in the report call for them.
#
#   lavaan_efa_se.csv  long: case, quantity, label, value
#     test    T (= (N-1) F), df, T3
#     rot     eps, Q, N, gpa_tol, constraint residual (lavaan's border function
#             at its own rotated estimate, one row per constraint "r1", "r2", ...)
#     est     thresholds "y1|t1", rotated loadings "fk=~yj", factor correlations
#             "fk~~fl" (k > l)
#     se      standard errors of the same parameters (std_ov FALSE bordered)
#     order   column-ordering statistic per rotated column, "index_fk" (mean row
#             index of loadings >= .8 max|loading| in the column; the order key),
#             and "minabs_over_sqrteps" per case (min |loading| / sqrt(eps))
#     diag    Jacobian-noise diagnostics (see above)
#   lavaan_efa_se_vcov.csv  long: case, row, col, value
#     the lower triangle (incl. diagonal) of vcov(fit), over every rotated free
#     parameter (thresholds, loadings, factor correlations; labels as above),
#     in lavaan's free-parameter order. lavaan's covariance is singular by
#     construction (rank n_u, spec A3.5) and is NOT regularised here.
#
# Read with stringcols(_all) and strtoreal(), delimiters(",") (the float-import trap).
# Run from the repository root:  Rscript fixtures/generate/99p_lavaan_efa_se.R
# Own manifest (MANIFEST_lavaan_efa_se.json), for 80_'s reason.
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))

if (!requireNamespace("lavaan", quietly = TRUE))
  stop("lavaan is required for tier-2 fixtures.")
if (!requireNamespace("numDeriv", quietly = TRUE))
  stop("numDeriv is required (lavaan's bordered Jacobian).")
# lavaan 0.7-3 from a SEPARATE library (the default library keeps 0.7-2 for the
# other generators). Owner's decision, 2026-11-17: lavaan commit 55786b9
# (2026-08-27, in 0.7-3) (a) makes the bordered border function under
# std_ov = TRUE use Lambda Phi Lambda' + Theta, and (b) changes the delta-route
# Jacobian numerics and adds the "bb" GPA algorithm (gpa_algorithm, new default).
# This file does not attribute 0.7-2/0.7-3 differences to one correction; see
# the diagnostic grid (scratch) that separates version, std_ov and rotation.se.
LAVAAN_LIB <- Sys.getenv("DWLSMV_LAVAAN_LIB",
                         file.path(Sys.getenv("LOCALAPPDATA"), "R", "lib-lavaan-0.7-3"))
stopifnot(dir.exists(LAVAAN_LIB),
          as.character(packageVersion("lavaan", lib.loc = LAVAAN_LIB)) == "0.7.3")
library(lavaan, lib.loc = LAVAAN_LIB)
stopifnot(as.character(packageVersion("lavaan")) == "0.7.3",
          normalizePath(find.package("lavaan")) ==
            normalizePath(file.path(LAVAAN_LIB, "lavaan")))

OUT <- file.path("fixtures", "expected")
MAN  <- file.path("fixtures", "MANIFEST_lavaan_efa_se.json")
SEED <- 20261117

GPA_TOL  <- 1e-12    # see TIGHT SETTINGS; chosen on convergence grounds
GPA_ALG  <- "bb"     # 0.7-3 default, set explicitly
MAX_ITER <- 10000L
STD_OV   <- FALSE    # comparator setting (A3.10), not a dwlsmv convention
BREAD    <- "expected" # lavaan 0.7-3 changed the default bread to "observed"; "expected" is 0.7-2 behaviour and the (Delta'W Delta)^-1 bread

eps_rule <- function(m) if (m == 2) 1e-4 else if (m == 3) 1e-3 else 1e-2
geomin_Q <- function(L, eps) sum(exp(rowSums(log(L^2 + eps)) / ncol(L)))

CASES <- list(E = list(file = "case_E.csv", m = 2),
              F = list(file = "case_F.csv", m = 3),
              G = list(file = "case_G.csv", m = 2))

fit_efa <- function(dat, m, rot, ...) {
  yn  <- names(dat)
  lhs <- paste(sprintf('efa("b")*f%d', seq_len(m)), collapse = " + ")
  lavaan::sem(paste(lhs, "=~", paste(yn, collapse = " + ")), data = dat,
              ordered = yn, estimator = "WLSMV", parameterization = "delta",
              test = "scaled.shifted", zero.add = c(0.5, 0.0), information.bread = BREAD,
              zero.keep.margins = TRUE, rotation = rot, ...)
}

# canonical label of each free parameter of the rotated table, in free order
free_labels <- function(fit, fn) {
  pt <- parTable(fit); pt <- pt[pt$free > 0, ]
  lab <- character(nrow(pt))
  for (i in seq_len(nrow(pt))) {
    if (pt$op[i] == "|") lab[i] <- paste0(pt$lhs[i], "|", pt$rhs[i])
    else if (pt$op[i] == "=~") lab[i] <- paste0(pt$lhs[i], "=~", pt$rhs[i])
    else if (pt$op[i] == "~~") {                       # factor covariance, k > l
      k <- match(c(pt$lhs[i], pt$rhs[i]), fn); k <- k[order(-k)]
      lab[i] <- paste0(fn[k[1]], "~~", fn[k[2]])
    } else stop("unexpected free parameter: ", pt$op[i])
  }
  stopifnot(!anyDuplicated(lab))
  lab
}

rows <- list(); vrows <- list()
seeds <- integer(0)
for (nm in names(CASES)) {
  cs  <- CASES[[nm]]
  dat <- read.csv(file.path("fixtures", "data", cs$file))
  dat[] <- lapply(dat, ordered)
  yn  <- names(dat); m <- cs$m; eps <- eps_rule(m); N <- nrow(dat)

  f0 <- fit_efa(dat, m, "none")
  seeds[nm] <- SEED + utf8ToInt(nm)
  set.seed(seeds[[nm]])
  wrn <- character(0)
  f1 <- withCallingHandlers(
    fit_efa(dat, m, "geomin",
            rotation.args = list(geomin_epsilon = eps, rstarts = 30L,
                                 row_weights = "none", std_ov = STD_OV,
                                 gpa_tol = GPA_TOL, gpa_algorithm = GPA_ALG, max_iter = MAX_ITER,
                                 order_lv_by = "index", orthogonal = FALSE)),
    warning = function(w) { wrn <<- c(wrn, conditionMessage(w))
                            invokeRestart("muffleWarning") })
  stopifnot(lavInspect(f0, "converged"), lavInspect(f1, "converged"))
  # convergence grounds: no GPA / half-stepping warning from the rotation
  stopifnot(!any(grepl("GP rotation|half-stepping", wrn)))
  ra <- lavInspect(f1, "options")$rotation.args
  stopifnot(ra$geomin_epsilon == eps, ra$rstarts == 30, ra$row_weights == "none",
            isFALSE(ra$std_ov), ra$gpa_tol == GPA_TOL, ra$gpa_algorithm == GPA_ALG,
            ra$max_iter == MAX_ITER,
            ra$order_lv_by == "index", isFALSE(ra$orthogonal),
            identical(lavInspect(f1, "options")$rotation.se, "bordered"))

  e0 <- lavInspect(f0, "est"); e1 <- lavInspect(f1, "est")
  L  <- e1$lambda[yn, , drop = FALSE]; P <- e1$psi; fn <- colnames(L)
  stopifnot(max(abs(diag(P) - 1)) < 1e-12, all(colSums(L) > 0),
            max(abs(L %*% P %*% t(L) -
                    e0$lambda[yn, ] %*% e0$psi %*% t(e0$lambda[yn, ]))) < 1e-12)

  # --- standard errors and full vcov, in the rotated table's free order ----
  lab <- free_labels(f1, fn)
  pt  <- parTable(f1); free_idx <- which(pt$free > 0)
  est <- pt$est[free_idx]; se <- pt$se[free_idx]
  V   <- vcov(f1); stopifnot(nrow(V) == length(lab), !anyNA(V), !anyNA(se))
  stopifnot(max(abs(sqrt(diag(V)) - se)) < 1e-12)
  dimnames(V) <- list(lab, lab)
  ll <- which(lower.tri(V, diag = TRUE), arr.ind = TRUE)

  # --- convergence record: lavaan's constraint residual at its own estimate -
  mod <- f1@Model
  xr  <- lavaan:::lav_model_get_parameters(mod)
  rr  <- lavaan:::lav_model_efa_rotate_border_x(xr, lavmodel = mod,
           lavoptions = f1@Options, lavpartable = f1@ParTable)

  # --- Jacobian-noise diagnostic (measures, decides nothing) ----------------
  jac0 <- mod@ceq.efa.JAC
  jac1 <- numDeriv::jacobian(lavaan:::lav_model_efa_rotate_border_x, x = xr,
            lavmodel = mod, lavoptions = f1@Options, lavpartable = f1@ParTable,
            method = "Richardson",
            method.args = list(eps = 1e-5, d = 1e-5, r = 6, v = 2))
  dj <- max(abs(jac1 - jac0))

  ss <- f1@test$scaled.shifted
  qq <- geomin_Q(L, eps)
  idx <- sapply(seq_len(m), function(k)
           mean(which(abs(L[, k]) >= 0.8 * max(abs(L[, k])))))
  names(idx) <- fn
  rows[[nm]] <- rbind(
    data.frame(case = nm, quantity = "test", label = c("T", "df", "T3"),
               value = c(f1@test$standard$stat, ss$df, ss$stat)),
    data.frame(case = nm, quantity = "rot",
               label = c("eps", "Q", "N", "gpa_tol", paste0("r", seq_along(rr))),
               value = c(eps, qq, N, GPA_TOL, rr)),
    data.frame(case = nm, quantity = "est", label = lab, value = est),
    data.frame(case = nm, quantity = "se",  label = lab, value = se),
    data.frame(case = nm, quantity = "order",
               label = c(paste0("index_", fn), "minabs_over_sqrteps"),
               value = c(idx, min(abs(L)) / sqrt(eps))),
    data.frame(case = nm, quantity = "diag", label = "jac_richardson_fine_maxabsdiff",
               value = dj))
  vrows[[nm]] <- data.frame(case = nm, row = lab[ll[, 1]], col = lab[ll[, 2]],
                            value = V[ll])
  cat(sprintf(paste0("case %s: m=%d eps=%g N=%d  T3 %.4f df %d  Q %.10f  ",
                     "max|r| %.2e  min|L|/sqrt(eps) %.3f  index %s  dJ %.2e  warnings %d\n"),
              nm, m, eps, N, ss$stat, ss$df, qq, max(abs(rr)),
              min(abs(L)) / sqrt(eps), paste(sprintf("%.3f", idx), collapse = "/"),
              dj, length(wrn)))
}

out <- do.call(rbind, rows);  rownames(out)  <- NULL
vout <- do.call(rbind, vrows); rownames(vout) <- NULL
stopifnot(!anyNA(out$value), !anyNA(vout$value))
f1p <- file.path(OUT, "lavaan_efa_se.csv")
f2p <- file.path(OUT, "lavaan_efa_se_vcov.csv")
write_csv_lf(out,  f1p)
write_csv_lf(vout, f2p)

write_manifest(MAN, c(f1p, f2p),
  seeds = seeds,
  extra = list(packages = c(
    paste("lavaan", as.character(packageVersion("lavaan")), "from separate library", if (nzchar(Sys.getenv("DWLSMV_LAVAAN_LIB"))) "DWLSMV_LAVAAN_LIB" else "%LOCALAPPDATA%/R/lib-lavaan-0.7-3", "(default library holds 0.7-2); owner decision 2026-11-17; lavaan commit 55786b9 in 0.7-3"),
    paste("numDeriv", as.character(packageVersion("numDeriv"))),
    "estimator WLSMV, parameterization delta, test scaled.shifted",
    "zero.add c(0.5, 0), zero.keep.margins TRUE (as 99g_, 99h_)",
    "rotation geomin, oblique; geomin_epsilon by Mplus's rule (E .0001, F .001, G .0001)",
    paste0("rotation.args: rstarts 30 (seeded per case), row_weights none, ",
           "std_ov FALSE (comparator setting, A3.10; not a dwlsmv convention), ",
           "gpa_tol ", format(GPA_TOL), ", gpa_algorithm ", GPA_ALG, " (0.7-3 default; gpa_fwindow 0 default), max_iter ", MAX_ITER,
           ", order_lv_by index, orthogonal FALSE, algorithm gpa, reflect TRUE, ",
           "jac_init_rot TRUE (default), tol 1e-8 (default, pairwise only)"),
    "rotation.se bordered (default); Jacobian numDeriv::jacobian Richardson, numDeriv defaults (eps 1e-4, d 1e-4, r 4, v 2), not settable in lavaan; diagnostic re-evaluation d 1e-5, r 6 recorded as quantity diag",
    "se robust.sem; information expected; information.bread expected (lavaan 0.7-3 default is observed; set explicitly, see generator header); vcov(fit) is singular by construction (rank n_u), not regularised",
    "inputs: case_E.csv (10_fixtures.R), case_F.csv and case_G.csv (91_efa_data.R)",
    R.version.string)))

cat("lavaan_efa_se.csv, lavaan_efa_se_vcov.csv and", MAN, "written.\n")
