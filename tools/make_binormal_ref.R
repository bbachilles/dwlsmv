# make_binormal_ref.R -- reference values for measuring delta_Phi of Stata's
# binormal() (tools/measure_binormal.do). Tool data, NOT a fixture.
#
# Output (in the directory this script is run from, or DPHI_OUT):
#   binormal_ref_corners.csv  one row per (h, k, rho) corner, plus margin rows
#   binormal_ref_cells.csv    one row per cell of the rectangles that
#                             dwlsmv_cellprob_raw() forms
#
# Reference method (M1). With r = sin(t), dr / sqrt(1 - r^2) = dt, so
#   Phi2(h,k;rho) = Phi(h)Phi(k) + (1/2pi) int_0^{asin rho} g(t) dt,
#   g(t) = exp(-(h^2 - 2hk sin t + k^2) / (2 cos^2 t)),
# a smooth, bounded integrand on a finite interval (no 1/sqrt(1-r^2)
# singularity), integrated by tanh-sinh quadrature in 160-bit MPFR.
# Margins Phi(h) are MPFR pnorm. Cells are differenced in MPFR, so the
# reference cell has no cancellation error at double resolution.
#
# Independent check (M2), different integration range and identity:
#   rho >= 0: Phi2 = Phi(min(h,k)) - (1/2pi) int_{asin rho}^{pi/2} g(t) dt
#   rho <  0: Phi2 = Phi(h) - Phi2(h, -k; -rho), with the inner term by M2
# Double-precision cross-check (M3): mvtnorm::pmvnorm(algorithm = TVPACK).
# The script stops if M1 and M2, or M1 at two quadrature steps, disagree
# beyond 1e-20 (= 4.5e-5 x 2^-52, far below what is being measured) on any
# corner. Observed (2026-09-29): see the REF| lines of the run.
#
# Each reference is written as a double-double pair ref_hi + ref_lo
# (ref_hi = nearest double, ref_lo = remainder), both %.17g. Reading only
# ref_hi would put up to 0.25 x 2^-52 of rounding into the very error being
# measured; the do-file forms (computed - ref_hi) - ref_lo instead.
#
# Grid choices (thresholds, rho values) were CHOSEN AFTER LOOKING at the
# lavaan threshold estimates in fixtures/replication/{mc_c_998, mc_c_999,
# mc_B_995, mc_B_9975}/lavaan_fits.csv (range -1.112 .. 0.919; medians
# -1.037, 0, 0.842 / -0.253, 0.675).

suppressMessages({ library(Rmpfr); library(mvtnorm) })
PREC <- 160
out_dir <- Sys.getenv("DPHI_OUT", ".")
two_pi <- 2 * Const("pi", PREC)
one <- mpfr(1, PREC)

# tanh-sinh (double-exponential) nodes on [-1,1], MPFR, built once. For an
# integrand analytic on the closed interval the error falls like
# exp(-c/step); convergence is CHECKED, not assumed: M1 at step 1/64 is
# compared with the same rule at step 1/32 (DE_CHECK) and with M2.
de_rule <- function(step, tmax = 4.5) {
  t <- mpfr(seq(-tmax, tmax, by = step), PREC)
  u <- Const("pi", PREC) / 2 * sinh(t)
  list(x = tanh(u), w = step * Const("pi", PREC) / 2 * cosh(t) / cosh(u)^2)
}
DE <- de_rule(1/64); DE_CHECK <- de_rule(1/32)
g_int <- function(h, k, lo, hi, rule = DE) {
  if (lo == hi) return(mpfr(0, PREC))
  h <- mpfr(h, PREC); k <- mpfr(k, PREC)
  t <- (lo + hi) / 2 + (hi - lo) / 2 * rule$x
  g <- exp(-(h * h - 2 * h * k * sin(t) + k * k) / (2 * cos(t)^2))
  g[!is.finite(g)] <- 0   # cos(t) == 0 exactly: integrand limit is 0 (h != k)
  (hi - lo) / 2 * sum(rule$w * g)
}
Pn <- function(x) pnorm(mpfr(x, PREC))

phi2_m1 <- function(h, k, rho) {
  Pn(h) * Pn(k) + g_int(h, k, mpfr(0, PREC), asin(mpfr(rho, PREC))) / two_pi
}
phi2_m2 <- function(h, k, rho) {
  if (rho >= 0) {
    Pn(min(h, k)) - g_int(h, k, asin(mpfr(rho, PREC)), Const("pi", PREC) / 2) / two_pi
  } else {
    Pn(h) - phi2_m2(h, -k, -rho)
  }
}

# ---- grid ------------------------------------------------------------------
thr_fixture <- c(-1.12, -1.10, -1.07, -1.05, -1.037, -1.02, -1.00, -0.97,
                 -0.95, -0.31, -0.253, -0.20, -0.06, 0, 0.06, 0.61, 0.675,
                 0.74, 0.78, 0.842, 0.85, 0.88, 0.92)
thr_wide <- c(-2.5, -1.5, 1.5, 2.5)
rho_tail <- c(0.9999, 0.9995, 0.999, 0.998, 0.9975, 0.995, 0.99, 0.985,
              0.98, 0.975, 0.97, 0.95)
rho_grid <- sort(unique(c(-rho_tail, rho_tail, -0.9, -0.8, -0.6, -0.3, 0,
                          0.3, 0.6, 0.8, 0.9)))
thr <- c(thr_fixture, thr_wide)
pairs <- subset(expand.grid(h = thr, k = thr), h <= k)   # Phi2 symmetric in (h,k)

regime_rho <- function(r) ifelse(abs(r) >= 0.97, ifelse(r > 0, "rho_hi", "rho_lo"), "rho_mid")

rows <- list(); n <- 0; maxdiff12 <- 0; maxdiff13 <- 0; maxdiffde <- 0
# A true Phi2 or cell far below the reference's own absolute accuracy (<= 5e-22
# observed) can come out as tiny negative noise (min seen -3.6e-41) from the
# Phi(h)Phi(k) + integral cancellation; it is written as 0. Only absolute
# error is meaningful for such rows; the do-file takes relative error only
# where the reference exceeds 1e-18.
hilo <- function(x) { if (x < 0) x <- mpfr(0, PREC); hi <- asNumeric(x); c(hi, asNumeric(x - hi)) }
for (i in seq_len(nrow(pairs))) for (rho in rho_grid) {
  h <- pairs$h[i]; k <- pairs$k[i]
  m1 <- phi2_m1(h, k, rho); m2 <- phi2_m2(h, k, rho)
  d12 <- abs(asNumeric(m1 - m2)); maxdiff12 <- max(maxdiff12, d12)
  mde <- Pn(h) * Pn(k) + g_int(h, k, mpfr(0, PREC), asin(mpfr(rho, PREC)), DE_CHECK) / two_pi
  maxdiffde <- max(maxdiffde, abs(asNumeric(m1 - mde)))
  m3 <- pmvnorm(upper = c(h, k), corr = matrix(c(1, rho, rho, 1), 2),
                algorithm = TVPACK(abseps = 1e-16))[1]
  maxdiff13 <- max(maxdiff13, abs(asNumeric(m1) - m3))
  n <- n + 1
  rows[[n]] <- data.frame(kind = "corner",
    regime = paste0(regime_rho(rho), ifelse(h %in% thr_wide | k %in% thr_wide, "_wide", "_fix")),
    h = h, k = k, rho = rho, ref_hi = hilo(m1)[1], ref_lo = hilo(m1)[2], m2_absdiff = d12)
}
if (maxdiff12 > 1e-20) stop("M1 and M2 disagree: ", maxdiff12)
if (maxdiffde > 1e-20) stop("quadrature not converged: ", maxdiffde)
for (h in thr) {   # infinite-threshold margins: evaluator uses normal()
  n <- n + 1
  rows[[n]] <- data.frame(kind = "margin", regime = "margin", h = h, k = NA,
    rho = NA, ref_hi = hilo(Pn(h))[1], ref_lo = hilo(Pn(h))[2], m2_absdiff = 0)
}
corners <- do.call(rbind, rows)

# ---- cells: rectangles dwlsmv_cellprob_raw() forms ---------------------------
sets <- list(                       # threshold vectors (a ; b), NA = unused
  med43   = list(a = c(-1.037, 0, 0.842),     b = c(-0.253, 0.675)),
  med44   = list(a = c(-1.037, 0, 0.842),     b = c(-1.037, 0, 0.842)),
  lo44    = list(a = c(-1.10, -0.06, 0.78),   b = c(-1.12, -0.06, 0.78)),
  hi44    = list(a = c(-0.97, 0.06, 0.92),    b = c(-0.95, 0.06, 0.88)),
  mix44   = list(a = c(-1.10, 0.06, 0.92),    b = c(-0.97, -0.06, 0.78)),
  c33     = list(a = c(-0.31, 0.74),          b = c(-0.20, 0.61)),
  wide44  = list(a = c(-2.5, 0, 2.5),         b = c(-1.5, 0, 1.5)))
PHI2 <- function(x, y, rho) {           # MPFR corner, handling +-Inf
  if (x == -Inf || y == -Inf) return(mpfr(0, PREC))
  if (x == Inf && y == Inf) return(one)
  if (x == Inf) return(Pn(y)); if (y == Inf) return(Pn(x))
  phi2_m1(x, y, rho)
}
crow <- list(); m <- 0
for (sn in names(sets)) for (rho in rho_grid) {
  a <- c(-Inf, sets[[sn]]$a, Inf); b <- c(-Inf, sets[[sn]]$b, Inf)
  C <- lapply(a, function(x) lapply(b, function(y) PHI2(x, y, rho)))
  for (r in 1:(length(a) - 1)) for (s in 1:(length(b) - 1)) {
    P <- C[[r + 1]][[s + 1]] - C[[r]][[s + 1]] - C[[r + 1]][[s]] + C[[r]][[s]]
    m <- m + 1
    av <- c(sets[[sn]]$a, NA, NA)[1:3]; bv <- c(sets[[sn]]$b, NA, NA)[1:3]
    crow[[m]] <- data.frame(set = sn, regime = regime_rho(rho), rho = rho,
      a1 = av[1], a2 = av[2], a3 = av[3], b1 = bv[1], b2 = bv[2], b3 = bv[3],
      r = r, s = s, ref_hi = hilo(P)[1], ref_lo = hilo(P)[2])
  }
}
cells <- do.call(rbind, crow)

# ---- write: every number %.17g, LF line endings --------------------------
fmt <- function(x) ifelse(is.na(x), "", sprintf("%.17g", x))
write_lf <- function(df, file) {
  for (j in names(df)) if (is.numeric(df[[j]])) df[[j]] <- fmt(df[[j]])
  con <- file(file, "wb")
  writeLines(paste(names(df), collapse = ","), con, sep = "\n")
  writeLines(do.call(paste, c(df, sep = ",")), con, sep = "\n")
  close(con)
}
write_lf(corners, file.path(out_dir, "binormal_ref_corners.csv"))
write_lf(cells, file.path(out_dir, "binormal_ref_cells.csv"))

cat(sprintf("REF|corners=%d cells=%d\n", nrow(corners), nrow(cells)))
cat(sprintf("REF|max |M1-M2| (MPFR, independent ranges) = %.3g\n", maxdiff12))
cat(sprintf("REF|max |M1-TVPACK| (double) = %.3g = %.3f x 2^-52\n",
            maxdiff13, maxdiff13 / 2^-52))
cat(sprintf("REF|max |M1(step 1/64)-M1(step 1/32)| = %.3g
", maxdiffde))
cat(sprintf("REF|min cell ref = %s\n", min(cells$ref_hi)))
cat(sprintf("REF|cells <= 1.6e-15: %d\n", sum(cells$ref_hi <= 1.6e-15)))
cat(sprintf("REF|R %s, Rmpfr %s, mvtnorm %s\n", getRversion(),
            packageVersion("Rmpfr"), packageVersion("mvtnorm")))
