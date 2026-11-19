# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 2: lavaan rotated EFA at a column-order tie
#
# Drafted by a read-only agent 2026-11-17; executed by the main session 2026-11-18.
#
# Release A, A4 near-tie question (design and acceptance criterion:
# E-109, A4 near-tie criterion, declared before any data). lavaan
# 0.7-3's efa(rotation = "geomin") posted rotated solution, rotated standard
# errors and FULL rotated covariance matrix on Cases T (exact tie) and B
# (boundary) from 92_efa_tie_data.R, with options IDENTICAL to 99p_ (this file
# reuses 99p_'s code paths and settings; see its header for the rationale of
# every one):
#
#   oblique geomin, epsilon by Mplus's rule (m = 2: .0001), row_weights none,
#   30 seeded rstarts, order_lv_by "index", std_ov FALSE (comparator setting,
#   A3.10), information.bread "expected", gpa_algorithm "bb", gpa_tol 1e-12,
#   max_iter 10000, rotation.se "bordered", se robust.sem, estimator WLSMV,
#   parameterization delta, test scaled.shifted, zero.add c(0.5, 0),
#   zero.keep.margins TRUE; lavaan 0.7-3 from the separate library.
#
# No reseeding (criterion.md): if lavaan does not converge, or emits ANY
# warning, this script stops with an error naming the case. It does not retry.
#
#   Case T  exact tie: all eight primaries .70; f1 items {1,4,5,8}, f2 {2,3,6,7}
#   Case B  boundary:  as T, y8 on f1 = .56 = 0.8 * .70
#
#   lavaan_efa_tie.csv  long: case, quantity, label, value
#     test    T, df, T3
#     rot     eps, Q (geomin criterion of the posted solution), N, gpa_tol,
#             constraint residual "r1", "r2", ...
#     est     thresholds "y1|t1", rotated loadings "fk=~yj", factor correlation
#             "f2~~f1"
#     se      standard errors of the same parameters
#     order   "index_fk" (mean row index of |loading| >= 0.8 max|loading| in the
#             column of lavaan's POSTED solution), "minabs_over_sqrteps"
#     tie     DIAGNOSTICS, never used for any choice:
#               n_near_cutoff        number of loadings with |l|/(0.8 max) in
#                                    [0.9, 1.1] (max of their own column)
#               ratio_fk=~yj         |l|/(0.8 max of column k), one row per such
#                                    loading
#               both_index_4.5       1 if every column's index is 4.5 (to 1e-12),
#                                    else 0
#               index_gap            |index_f1 - index_f2|
#     diag    Jacobian-noise diagnostic, as 99p_
#   lavaan_efa_tie_vcov.csv  long: case, row, col, value; lower triangle of
#     vcov(fit) over every rotated free parameter, as lavaan_efa_se_vcov.csv.
#
# lavaan 0.7-3 does not expose which of the 30 random starts it selected (the
# rotation's per-start results are not stored on the fit, and its verbose
# per-start log is not emitted from the estimation step), so the selected start
# is not recorded; the criterion value Q of the posted solution is.
#
# Read with stringcols(_all) and strtoreal(), delimiters(",") (the float-import trap).
# Run from the repository root:  Rscript fixtures/generate/99r_lavaan_efa_tie.R
# Own manifest (MANIFEST_lavaan_efa_tie.json), for 80_'s reason.
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))

if (!requireNamespace("lavaan", quietly = TRUE))
  stop("lavaan is required for tier-2 fixtures.")
if (!requireNamespace("numDeriv", quietly = TRUE))
  stop("numDeriv is required (lavaan's bordered Jacobian).")
# lavaan 0.7-3 from a SEPARATE library, exactly as 99p_ (owner, 2026-11-17).
LAVAAN_LIB <- Sys.getenv("DWLSMV_LAVAAN_LIB",
                         file.path(Sys.getenv("LOCALAPPDATA"), "R", "lib-lavaan-0.7-3"))
stopifnot(dir.exists(LAVAAN_LIB),
          as.character(packageVersion("lavaan", lib.loc = LAVAAN_LIB)) == "0.7.3")
library(lavaan, lib.loc = LAVAAN_LIB)
stopifnot(as.character(packageVersion("lavaan")) == "0.7.3",
          normalizePath(find.package("lavaan")) ==
            normalizePath(file.path(LAVAAN_LIB, "lavaan")))

OUT <- file.path("fixtures", "expected")
MAN  <- file.path("fixtures", "MANIFEST_lavaan_efa_tie.json")
SEED <- 20261117

GPA_TOL  <- 1e-12    # as 99p_
GPA_ALG  <- "bb"
MAX_ITER <- 10000L
STD_OV   <- FALSE    # comparator setting (A3.10), not a dwlsmv convention
BREAD    <- "expected"

eps_rule <- function(m) if (m == 2) 1e-4 else if (m == 3) 1e-3 else 1e-2
geomin_Q <- function(L, eps) sum(exp(rowSums(log(L^2 + eps)) / ncol(L)))

CASES <- list(T = list(file = "efatie_T.csv", m = 2),
              B = list(file = "efatie_B.csv", m = 2))

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

  # criterion.md "or warns" covers both fits: the unrotated fit's warnings are
  # collected too (main-session review fix; 99p_ left f0 outside the handler).
  wrn <- character(0)
  f0 <- withCallingHandlers(fit_efa(dat, m, "none"),
    warning = function(w) { wrn <<- c(wrn, paste("unrotated:", conditionMessage(w)))
                            invokeRestart("muffleWarning") })
  seeds[nm] <- SEED + utf8ToInt(nm)
  set.seed(seeds[[nm]])
  f1 <- withCallingHandlers(
    fit_efa(dat, m, "geomin",
            rotation.args = list(geomin_epsilon = eps, rstarts = 30L,
                                 row_weights = "none", std_ov = STD_OV,
                                 gpa_tol = GPA_TOL, gpa_algorithm = GPA_ALG, max_iter = MAX_ITER,
                                 order_lv_by = "index", orthogonal = FALSE)),
    warning = function(w) { wrn <<- c(wrn, conditionMessage(w))
                            invokeRestart("muffleWarning") })
  # criterion.md: no reseeding. Failure or ANY warning stops the script, naming the case.
  if (!isTRUE(lavInspect(f0, "converged")) || !isTRUE(lavInspect(f1, "converged")))
    stop(sprintf("case %s: lavaan did not converge; reported as failed, not reseeded.", nm))
  if (length(wrn))
    stop(sprintf("case %s: lavaan warned (%s); reported as failed, not reseeded.",
                 nm, paste(wrn, collapse = " | ")))
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

  # --- Jacobian-noise diagnostic (measures, decides nothing; as 99p_) -------
  jac0 <- mod@ceq.efa.JAC
  jac1 <- numDeriv::jacobian(lavaan:::lav_model_efa_rotate_border_x, x = xr,
            lavmodel = mod, lavoptions = f1@Options, lavpartable = f1@ParTable,
            method = "Richardson",
            method.args = list(eps = 1e-5, d = 1e-5, r = 6, v = 2))
  dj <- max(abs(jac1 - jac0))

  ss <- f1@test$scaled.shifted
  qq <- geomin_Q(L, eps)
  cmax <- apply(abs(L), 2, max)
  idx <- sapply(seq_len(m), function(k)
           mean(which(abs(L[, k]) >= 0.8 * cmax[k])))
  names(idx) <- fn

  # --- tie diagnostics (recorded, never used for any choice) ---------------
  ratio <- abs(L) / matrix(0.8 * cmax, nrow(L), m, byrow = TRUE)
  near  <- which(ratio >= 0.9 & ratio <= 1.1, arr.ind = TRUE)
  near_lab <- if (nrow(near))
    sprintf("ratio_%s=~%s", fn[near[, 2]], yn[near[, 1]]) else character(0)
  near_val <- if (nrow(near)) ratio[near] else numeric(0)
  both45 <- as.numeric(all(abs(idx - 4.5) < 1e-12))
  tie_rows <- data.frame(case = nm, quantity = "tie",
    label = c("n_near_cutoff", near_lab, "both_index_4.5", "index_gap"),
    value = c(nrow(near), near_val, both45, abs(idx[[1]] - idx[[2]])))

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
    tie_rows,
    data.frame(case = nm, quantity = "diag", label = "jac_richardson_fine_maxabsdiff",
               value = dj))
  vrows[[nm]] <- data.frame(case = nm, row = lab[ll[, 1]], col = lab[ll[, 2]],
                            value = V[ll])
  cat(sprintf(paste0("case %s: m=%d eps=%g N=%d  T3 %.4f df %d  Q %.10f  ",
                     "max|r| %.2e  index %s  both 4.5: %d  near-cutoff %d  dJ %.2e  warnings %d\n"),
              nm, m, eps, N, ss$stat, ss$df, qq, max(abs(rr)),
              paste(sprintf("%.3f", idx), collapse = "/"), as.integer(both45),
              nrow(near), dj, length(wrn)))
}

out <- do.call(rbind, rows);  rownames(out)  <- NULL
vout <- do.call(rbind, vrows); rownames(vout) <- NULL
stopifnot(!anyNA(out$value), !anyNA(vout$value))
f1p <- file.path(OUT, "lavaan_efa_tie.csv")
f2p <- file.path(OUT, "lavaan_efa_tie_vcov.csv")
write_csv_lf(out,  f1p)
write_csv_lf(vout, f2p)

write_manifest(MAN, c(f1p, f2p),
  seeds = seeds,
  extra = list(packages = c(
    paste("lavaan", as.character(packageVersion("lavaan")), "from separate library", if (nzchar(Sys.getenv("DWLSMV_LAVAAN_LIB"))) "DWLSMV_LAVAAN_LIB" else "%LOCALAPPDATA%/R/lib-lavaan-0.7-3", "(default library holds 0.7-2); owner decision 2026-11-17; lavaan commit 55786b9 in 0.7-3"),
    paste("numDeriv", as.character(packageVersion("numDeriv"))),
    "design and criterion: docs/notes/efa-se/a4-neartie/criterion.md (declared 2026-11-17 before any data); options identical to 99p_lavaan_efa_se.R",
    "estimator WLSMV, parameterization delta, test scaled.shifted",
    "zero.add c(0.5, 0), zero.keep.margins TRUE (as 99g_, 99h_, 99p_)",
    "rotation geomin, oblique; geomin_epsilon by Mplus's rule (m = 2: .0001)",
    paste0("rotation.args: rstarts 30 (seeded per case), row_weights none, ",
           "std_ov FALSE (comparator setting, A3.10; not a dwlsmv convention), ",
           "gpa_tol ", format(GPA_TOL), ", gpa_algorithm ", GPA_ALG, " (0.7-3 default; gpa_fwindow 0 default), max_iter ", MAX_ITER,
           ", order_lv_by index, orthogonal FALSE, algorithm gpa, reflect TRUE, ",
           "jac_init_rot TRUE (default), tol 1e-8 (default, pairwise only)"),
    "rotation.se bordered (default); Jacobian numDeriv::jacobian Richardson, numDeriv defaults, not settable in lavaan; diagnostic re-evaluation d 1e-5, r 6 recorded as quantity diag",
    "se robust.sem; information expected; information.bread expected (lavaan 0.7-3 default is observed; set explicitly); vcov(fit) is singular by construction (rank n_u), not regularised",
    "quantity tie: diagnostics only (n_near_cutoff, ratio_fk=~yj within 10% of 0.8*column max, both_index_4.5, index_gap); not used for any choice",
    "selected random start is not exposed by lavaan 0.7-3 and is not recorded; criterion value Q of the posted solution is (quantity rot)",
    "no reseeding: non-convergence or any warning stops the generator",
    "inputs: efatie_T.csv, efatie_B.csv (92_efa_tie_data.R)",
    R.version.string)))

cat("lavaan_efa_tie.csv, lavaan_efa_tie_vcov.csv and", MAN, "written.\n")
