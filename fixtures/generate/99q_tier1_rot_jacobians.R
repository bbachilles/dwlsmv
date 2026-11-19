# ---------------------------------------------------------------------------
# Drafted by a read-only agent (2026-11-17). Executed by the main session from
# the repository root 2026-11-17 (R 4.6.0): deciding criteria passed on E, F, G
# (R 1.4e-16 to 4.0e-16, J 2.0e-16 to 9.3e-16 relative); output identical to
# the agent's scratch run (md5 9c864cb9...); fixtures/verify.R clean.
#
# dwlsmv fixture generation -- tier 1 (base R only, no lavaan): the constraint
# Jacobian R and the rotation-map Jacobian J of the geomin-rotated EFA
# (release A, unit A4; E-108, unit A4; spec/derivations/
# dwlsmv-derivations.qmd @sec-rot-se RS3, RS4).
#
# WHAT IT CERTIFIES (docs/prereg-template.md 3a)
#   R  the m(m-1) x (p m + m(m-1)/2) Jacobian of the stationarity function r
#      (RS1) with respect to (vec Lambda, phi), at a stated rotated point.
#      FIXTURE VALUE: the complex step of r (h = 1e-30).
#      COMPARATOR: the closed form of RS4 (rs_R). Differs in representation
#      (derivative formula vs differentiated r) and algorithm (analytic vs
#      complex step); shares r, i.e. GM2's gradient gm_dQ, and the .qmd chunks
#      that define both. The fixture value does not use rs_R.
#   J  the Jacobian of the rotation map rho (RS3 ii), (vec Lambda, phi) with
#      respect to the free echelon loadings a, at the Newton-polished point.
#      FIXTURE VALUE: the complex step of an actual re-rotation (Newton on
#      r = 0 from the perturbed unrotated loadings; R is not used, the
#      iteration matrix is a forward-difference one).
#      COMPARATOR: the closed form (I - G_T (R G_T)^-1 R) G_u (rs_J). Differs
#      in representation (implicit root of r = 0 vs a linear formula) and
#      algorithm; shares r, g and the unrotated input; uses R.
#   Deciding criterion (fixed before running): R closed vs complex step
#   <= 1e-12 relative to max|R|, at both points; J closed vs re-rotation
#   <= 1e-10 relative to max|J|. The generator writes nothing if either fails.
#   EVERYTHING ELSE is a diagnostic: printed and written, never fails:
#   R J = 0, Delta* J = Delta_u (both by complex step, independent of the
#   closed forms), a real central/Richardson re-rotation of J, the stall-vs-
#   root J gap, residuals, comparison with lavaan's committed rotated solution.
#   The bordered-versus-delta covariance is NOT used: it holds for any R
#   (RS3 i) and is no check of R.
#
# EVALUATION POINTS (spec Module 12, Release A, A3.15: the selected GPA stall
# point, no Newton polish). Per case, two points, both through the sign and
# order rule:
#   "stall"  gm_geomin()'s selected start, as production computes it. R is
#            written here (R is exact at any point, so a test feeds this
#            Lambda, Phi to production's R).
#   "root"   the same start Newton-polished to r ~ 1e-15. R and J are written
#            here (J is the derivative of rho, defined at the root).
#   diag:J_stall_vs_root records how far the closed-form J at the stall point
#   is from the root J; a production test of J at the stall point needs a
#   tolerance that covers it (the owner's choice, made after this number).
#
# INPUTS (all committed; neither R packages nor lavaan are called):
#   fixtures/expected/lavaan_efa_geomin.csv  rotated lambda and phi, E/F/G
#       (tier 2, std_ov TRUE). sigma = Lambda Phi Lambda' is invariant to the
#       rotation, so it determines the echelon A.
#   A (p x m, A[i,j] = 0 for i < j) is DERIVED from them: A = Lambda C Q, C the
#   Cholesky factor of Phi, Q from the LQ decomposition of the top m x m block
#   of Lambda C, diagonal made positive. A is a deterministic function of
#   committed numbers and is itself written to the fixture.
#   fixtures/expected/lavaan_efa_case_E.csv  Case E's own unrotated lavaan
#       loadings: a cross-check of the derivation of A (diagnostic).
#   fixtures/data/case_{E,F,G}.csv  only to count thresholds (nth), which the
#       R/J fixture omits (R's threshold columns are 0, J's threshold rows are
#       [I, 0]).
#   p, m, eps: E 6, 2, 1e-4; F 9, 3, 1e-3; G 6, 2, 1e-4 (Mplus's rule).
# The mathematics is the derivation's: chunks gm-impl, rs-impl, and the one
# function rs_rerot_cs of rs-3-rerot are extracted from the CURRENT .qmd (md5
# in the manifest); this file defines none.
#
# Output (OUTROOT/expected/rotse_tier1_jacobians.csv, numbers %.17g), long form
#   case,point,quantity,row,col,value
# quantities: dims (point "both"; rows 1..7 = p, m, eps, nth, n_free_a, rows of
# R, cols of R), A, Lambda, Phi, T, R, J, r (residual vector), and "diag:*"
# scalars. Matrix indices are 1-based. R rows: the off-diagonal (j,k) of an
# m x m matrix, column by column; R cols: vec Lambda column-major, then the
# strict lower triangle of Phi column-major. J rows: R's columns; J cols: the
# free entries of A (A[i,j], i >= j), column-major.
# and OUTROOT/MANIFEST_rotse_tier1.json.
#
# Run from the repository root:
#   Rscript fixtures/generate/99q_tier1_rot_jacobians.R            # OUTROOT = fixtures
#   Rscript fixtures/generate/99q_tier1_rot_jacobians.R <outroot>  # elsewhere (testing)
# ---------------------------------------------------------------------------
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

args <- commandArgs(TRUE)
OUTROOT <- if (length(args) >= 1) args[1] else "fixtures"
QMD    <- file.path("spec", "derivations", "dwlsmv-derivations.qmd")
COMMON <- file.path("fixtures", "generate", "00_common.R")
IN_GEO <- file.path("fixtures", "expected", "lavaan_efa_geomin.csv")
IN_E   <- file.path("fixtures", "expected", "lavaan_efa_case_E.csv")
DATA   <- file.path("fixtures", "data", sprintf("case_%s.csv", c("E", "F", "G")))
CASES  <- list(E = list(m = 2, eps = 1e-4), F = list(m = 3, eps = 1e-3),
               G = list(m = 2, eps = 1e-4))
TOL_R  <- 1e-12   # deciding: R closed form vs complex step (as RS4)
TOL_J  <- 1e-10   # deciding: J closed form vs re-rotation (as RS3 ii)
stopifnot(file.exists(QMD, COMMON, IN_GEO, IN_E, DATA),
          dir.exists(file.path(OUTROOT, "expected")))

source(COMMON)                                  # write_csv_lf, write_manifest
qmd <- readLines(QMD, encoding = "UTF-8")
chunk <- function(lab) {
  s <- grep(paste0("^```\\{r ", lab, "\\}"), qmd); stopifnot(length(s) == 1)
  e <- s + which(qmd[(s + 1):length(qmd)] == "```")[1]
  qmd[(s + 1):(e - 1)]
}
is_fun_def <- function(ex) is.call(ex) && identical(ex[[1]], as.name("<-")) &&
  is.name(ex[[2]]) && is.call(ex[[3]]) && identical(ex[[3]][[1]], as.name("function"))
for (lab in c("gm-impl", "rs-impl"))            # function definitions only
  for (ex in parse(text = chunk(lab), keep.source = FALSE)) if (is_fun_def(ex)) eval(ex, globalenv())
for (ex in parse(text = chunk("rs-3-rerot"), keep.source = FALSE))   # rs_rerot_cs only
  if (is_fun_def(ex) && identical(ex[[2]], as.name("rs_rerot_cs"))) eval(ex, globalenv())
stopifnot(exists("gm_geomin"), exists("rs_R"), exists("rs_J"), exists("rs_rerot_cs"),
          exists("rs_polish"), exists("rs_canon"), exists("rs_csjac"), exists("gm_align"))
f17 <- function(x) sprintf("%.17g", x)

# ---- inputs -----------------------------------------------------------------
# numeric text read as character and converted, as the project's fixture rule
geo <- read.csv(IN_GEO, colClasses = "character"); geo$value <- as.numeric(geo$value)
unrot_E <- read.csv(IN_E, colClasses = "character"); unrot_E$value <- as.numeric(unrot_E$value)

# echelon A from sigma-determining rotated (Lambda, Phi): see the header
echelon_A <- function(L, Phi) {
  C <- t(chol(Phi)); M <- L %*% C; m <- ncol(L)
  Q <- qr.Q(qr(t(M[1:m, , drop = FALSE])), complete = TRUE)
  A <- M %*% Q
  A <- sweep(A, 2, sign(diag(A[1:m, , drop = FALSE])), "*")
  stopifnot(max(abs(A[row(A) < col(A)])) < 1e-13)
  A[row(A) < col(A)] <- 0
  stopifnot(max(abs(tcrossprod(A) - L %*% Phi %*% t(L))) < 1e-13)
  A
}

out <- list(); dg <- list()
add <- function(cs, pt, qn, M) {                # matrix or vector -> long rows
  M <- if (is.matrix(M)) M else matrix(M, ncol = 1)
  ix <- expand.grid(row = seq_len(nrow(M)), col = seq_len(ncol(M)))
  out[[length(out) + 1]] <<- data.frame(case = cs, point = pt, quantity = qn,
    row = ix$row, col = ix$col, value = f17(as.vector(M)))
}
dgs <- function(cs, pt, nm, v)                  # named scalar diagnostic
  dg[[length(dg) + 1]] <<- data.frame(case = cs, point = pt, quantity = paste0("diag:", nm),
    row = 1L, col = 1L, value = f17(v))

# complex-step references, with NO use of the closed forms
cs_R <- function(L, Phi, eps, p, m) {
  rf <- function(z) { u <- rs_unpack(z, p, m, 0); rs_r(u$L, u$Phi, eps) }
  rs_csjac(rf, c(as.vector(L), Phi[lower.tri(Phi)]))
}
cs_Dstar <- function(L, Phi, p, m) {
  sf <- function(z) { u <- rs_unpack(z, p, m, 0); rs_sig(u$L, u$Phi) }
  rs_csjac(sf, c(as.vector(L), Phi[lower.tri(Phi)]))
}
cs_Du <- function(A, p, m) {
  fr <- which(rs_free(p, m))
  rs_csjac(function(z) { B <- matrix(0 * z[1], p, m); B[fr] <- z; rs_sig(B, diag(m)) }, A[fr])
}

fails <- character(0)
for (cs in names(CASES)) {
  m <- CASES[[cs]]$m; eps <- CASES[[cs]]$eps
  g0 <- geo[geo$case == cs & geo$quantity == "lambda", ]
  p <- nrow(g0) / m
  stopifnot(identical(g0$label, sprintf("f%d=~y%d", rep(1:m, each = p), rep(1:p, m))))
  Llav <- matrix(g0$value, p, m)
  ph <- geo[geo$case == cs & geo$quantity == "phi", ]
  lw <- which(lower.tri(diag(m)), arr.ind = TRUE); lw <- lw[order(lw[, 2], lw[, 1]), , drop = FALSE]
  stopifnot(identical(ph$label, sprintf("f%d~~f%d", lw[, 1], lw[, 2])))
  Plav <- diag(m); Plav[lower.tri(Plav)] <- ph$value; Plav <- Plav + t(Plav) - diag(m)
  A <- echelon_A(Llav, Plav)
  dat <- read.csv(DATA[match(cs, names(CASES))])
  nl <- unique(sapply(dat, function(v) length(unique(v)))); stopifnot(length(nl) == 1)
  nth <- (nl - 1) * ncol(dat)
  fr <- which(rs_free(p, m))
  if (cs == "E") {
    Au <- matrix(unrot_E$value[unrot_E$quantity == "lambda"], p, m)
    dgs(cs, "input", "A_derived_vs_committed_unrotated_E_maxabsdiff", max(abs(abs(A) - abs(Au))))
  }

  # ---- the two points ------------------------------------------------------
  g  <- gm_geomin(A, eps = eps, nstart = 30)
  pz <- rs_polish(A, g$T, eps)
  pts <- list(stall = rs_canon(gm_rot(A, g$T), crossprod(g$T), g$T),
              root  = rs_canon(pz$L, pz$Phi, pz$T))
  stopifnot(max(abs(pts$stall$L - pts$root$L)) < 1e-5)       # same minimum, same labelling
  dgs(cs, "both", "n_distinct_start_minima", length(unique(round(g$fall, 6))))
  dgs(cs, "both", "gpa_Gp_norm", g$s)
  dgs(cs, "both", "gpa_stalled_flag", as.numeric(g$stalled))
  dgs(cs, "both", "newton_its", pz$it)

  Rcs <- list()
  for (pt in names(pts)) {
    P <- pts[[pt]]; L <- A %*% solve(t(P$T)); Phi <- crossprod(P$T)
    Rcs[[pt]] <- cs_R(L, Phi, eps, p, m)
    eR <- rs_relmax(rs_R(L, Phi, eps, 0), Rcs[[pt]])
    if (!(eR < TOL_R)) fails <- c(fails, sprintf("R %s %s %.3g", cs, pt, eR))
    add(cs, pt, "A", A); add(cs, pt, "Lambda", L); add(cs, pt, "Phi", Phi); add(cs, pt, "T", P$T)
    add(cs, pt, "R", Rcs[[pt]]); add(cs, pt, "r", rs_r(L, Phi, eps))
    dgs(cs, pt, "R_closed_vs_cs_rel", eR)
    dgs(cs, pt, "r_max_abs", max(abs(rs_r(L, Phi, eps))))
    dgs(cs, pt, "rank_R_over_rows", rs_nrank(Rcs[[pt]]) / nrow(Rcs[[pt]]))
    dgs(cs, pt, "cond_R_GT", rs_cond(Rcs[[pt]] %*% rs_GT(L, Phi, 0)))
  }

  # ---- J at the root: re-rotation (fixture) against the closed form --------
  P <- pts$root; Minv <- solve(t(P$T)); L <- A %*% Minv; Phi <- crossprod(P$T)
  x <- list(cs = cs, m = m, eps = eps, p = p, A = A, T = P$T)
  Jcs <- rs_rerot_cs(x)                          # rows (vec L, phi), cols free a
  Jcf <- rs_J(rs_R(L, Phi, eps, 0), rs_Gu(L, Minv, 0), rs_GT(L, Phi, 0))
  eJ <- rs_relmax(Jcf, Jcs); if (!(eJ < TOL_J)) fails <- c(fails, sprintf("J %s root %.3g", cs, eJ))
  add(cs, "root", "J", Jcs)
  dgs(cs, "root", "J_closed_vs_rerotation_rel", eJ)

  # ---- diagnostics (never fail) --------------------------------------------
  Rr <- Rcs$root
  dgs(cs, "root", "RcsJcs_zero_rel", max(abs(Rr %*% Jcs)) / (max(abs(Rr)) * max(abs(Jcs))))
  dgs(cs, "root", "DstarJ_minus_Du_rel", rs_relmax(cs_Dstar(L, Phi, p, m) %*% Jcs, cs_Du(A, p, m)))
  Jr <- sapply(seq_along(fr), function(c) {      # real central Richardson re-rotation
    cen <- function(h) {
      Ap <- A; Am <- A; Ap[fr[c]] <- Ap[fr[c]] + h; Am[fr[c]] <- Am[fr[c]] - h
      z1 <- rs_polish(Ap, P$T, eps, tol = 1e-15); z2 <- rs_polish(Am, P$T, eps, tol = 1e-15)
      (c(as.vector(z1$L), z1$Phi[lower.tri(z1$Phi)]) -
         c(as.vector(z2$L), z2$Phi[lower.tri(z2$Phi)])) / (2 * h)
    }
    (4 * cen(5e-5) - cen(1e-4)) / 3
  })
  dgs(cs, "root", "J_vs_real_richardson_rel", rs_relmax(Jr, Jcs))
  S <- pts$stall; Ls <- S$L; Ps <- crossprod(S$T)
  Jst <- rs_J(rs_R(Ls, Ps, eps, 0), rs_Gu(Ls, solve(t(S$T)), 0), rs_GT(Ls, Ps, 0))
  dgs(cs, "stall", "J_stall_vs_root", rs_relmax(Jst, Jcs))
  al <- gm_align(L, Llav, Phi, Plav)             # lavaan's committed solution: sign/order
  dgs(cs, "root", "own_vs_lavaan_dLambda", al$dL); dgs(cs, "root", "own_vs_lavaan_dPhi", al$dPhi)
  add(cs, "both", "dims", c(p, m, eps, nth, length(fr), nrow(Rr), ncol(Rr)))
}
if (length(fails))
  stop("99q: deciding criterion failed, nothing written: ", paste(fails, collapse = "; "), call. = FALSE)

dgd <- do.call(rbind, dg)
print(dgd[, c("case", "point", "quantity", "value")], row.names = FALSE)
csv <- file.path(OUTROOT, "expected", "rotse_tier1_jacobians.csv")
write_csv_lf(rbind(do.call(rbind, out), dgd), csv)
write_manifest(file.path(OUTROOT, "MANIFEST_rotse_tier1.json"), csv,
  seeds = "none (no RNG: gm_starts is the deterministic GM4 grid)",
  extra = list(packages = c(
    "base R only (no lavaan, no mvtnorm)", R.version.string,
    paste("qmd md5", QMD, unname(tools::md5sum(QMD)), "chunks gm-impl, rs-impl, rs-3-rerot (rs_rerot_cs only), function definitions only"),
    paste("common md5", COMMON, unname(tools::md5sum(COMMON))),
    paste("input md5", c(IN_GEO, IN_E, DATA), unname(tools::md5sum(c(IN_GEO, IN_E, DATA)))),
    "cases E (p 6, m 2, eps 1e-4), F (p 9, m 3, eps 1e-3), G (p 6, m 2, eps 1e-4); GM4 grid, 30 starts; sign and order rule rs_canon",
    "points: stall = GPA stall (A3.15), root = Newton-polished (r ~ 1e-15)",
    "R: complex step h 1e-30 of rs_r; comparator closed form RS4; deciding 1e-12 relative to max|R|",
    "J: complex step h 1e-20 of re-rotation rs_rerot_cs (Newton on r = 0, forward-difference iteration matrix, R not used); comparator closed form RS3(ii); deciding 1e-10 relative to max|J|",
    "A derived from committed lavaan_efa_geomin.csv (Lambda, Phi) by Cholesky + LQ, positive diagonal",
    "threshold columns of R are 0 and threshold rows of J are [I, 0]: omitted; nth in dims",
    "not certified: calibration, other data, m >= 4, several minima, near ties, empty binary cells, the Stata route, production tolerances; a pass here is not a pass at any other layer")))
cat("99q: wrote", csv, "\n")
