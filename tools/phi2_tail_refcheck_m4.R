# phi2_tail_refcheck_m4.R -- an independent check of the Phi2 tail study's MPFR
# reference (E-063, open decision (e);
# E-063). Tool, not a fixture, not the frozen run.
#
# Why: the reference's M1 and M2 share their representation (the arcsine-form
# integrand g(t)), their quadrature (tanh-sinh, step 1/64) and their margins.
# Their agreement does not rule out an error common to both. M4 differs in both:
#   M4  Phi2(h, k; rho) = int_{-40}^{h} phi(x) Phi((k - rho x)/s) dx,
#       s = sqrt(1 - rho^2), by adaptive 20-point Gauss-Legendre in MPFR
#       (nodes by Newton on P_20), panels split at integers and at the step
#       x0 = k / rho. Truncation at -40: the integrand is <= phi(x), so the
#       omitted mass is < Phi(-40) ~ 4e-350.
#   M4c M4 at a 1e3-fold tighter adaptive tolerance: M4's own convergence.
#   TV  mvtnorm's TVPACK (double precision only): a third algorithm, ~1e-16.
#
# Scope (user, 2026-10-25): the 569 corners whose half-step check failed in the
# full reference run, plus a fixed-seed sample of the passing corners, 20 per
# rho value. Reads phi2tail_refcheck.csv from the run's directory (PHI2T_OUT).
# Writes phi2tail_m4check.csv there, and prints M4| summary lines.
#
# Run from the repository root:
#   PHI2T_OUT=<dir of the reference run> Rscript tools/phi2_tail_refcheck_m4.R

suppressMessages({ library(Rmpfr); library(mvtnorm) })
PREC <- 256
dir <- Sys.getenv("PHI2T_OUT", ".")

# M1 and M2 exactly as the reference defines them
src <- readLines("tools/make_phi2_tail_ref.R")
a <- grep("^# ---- MPFR primitives", src); b <- grep("^cache <- new.env", src) - 1
eval(parse(text = src[a:b]))

# ---- M4 ---------------------------------------------------------------------
gl_rule <- function(n) {
  x <- mpfr(cos(pi * (seq_len(n) - 0.25) / (n + 0.5)), PREC)
  for (it in 1:100) {
    p0 <- rep(one, n); p1 <- x
    for (j in 2:n) { p2 <- ((2 * j - 1) * x * p1 - (j - 1) * p0) / j; p0 <- p1; p1 <- p2 }
    dp <- n * (x * p1 - p0) / (x * x - 1)
    dx <- p1 / dp; x <- x - dx
    if (max(abs(asNumeric(dx))) < 1e-70) break
  }
  list(x = x, w = 2 / ((1 - x * x) * dp * dp))
}
GL <- gl_rule(20)
stopifnot(abs(asNumeric(sum(GL$w)) - 2) < 1e-60)          # the rule integrates 1 exactly
sqrt2pi <- sqrt(2 * Const("pi", PREC))
gl_panel <- function(a, b, f) { m <- (a + b) / 2; r <- (b - a) / 2; r * sum(GL$w * f(m + r * GL$x)) }
gl_adapt <- function(a, b, f, whole, abs_tol, rel_tol, depth = 0) {
  m <- (a + b) / 2; L <- gl_panel(a, m, f); R <- gl_panel(m, b, f)
  if (abs(asNumeric(L + R - whole)) <= max(abs_tol, rel_tol * abs(asNumeric(L + R))) || depth >= 50)
    return(L + R)
  gl_adapt(a, m, f, L, abs_tol, rel_tol, depth + 1) + gl_adapt(m, b, f, R, abs_tol, rel_tol, depth + 1)
}
phi2_m4 <- function(h, k, rho, abs_tol = 1e-45, rel_tol = 1e-28, lo = -40) {
  r <- mpfr(rho, PREC); s <- sqrt(1 - r * r); km <- mpfr(k, PREC)
  f <- function(x) exp(-x * x / 2) / sqrt2pi * pnorm((km - r * x) / s)
  br <- sort(unique(c(lo, seq(ceiling(lo), floor(h)), h,
                      if (rho != 0 && k / rho > lo && k / rho < h) k / rho)))
  out <- zero
  for (q in seq_len(length(br) - 1)) {
    aa <- mpfr(br[q], PREC); bb <- mpfr(br[q + 1], PREC)
    out <- out + gl_adapt(aa, bb, f, gl_panel(aa, bb, f), abs_tol, rel_tol)
  }
  out
}

# ---- scope ------------------------------------------------------------------
rc <- read.csv(file.path(dir, "phi2tail_refcheck.csv"))
rc$fail <- pmax(rc$d_m2, rc$d_half) > rc$tol
set.seed(20261025)
pass_idx <- unlist(lapply(split(which(!rc$fail), rc$rho[!rc$fail]), function(ix)
  if (length(ix) <= 20) ix else sample(ix, 20)))
grp <- function(d, g) if (nrow(d)) cbind(d, group = g) else NULL   # a group may be empty
scope <- rbind(grp(rc[rc$fail, ], "half_step_fail"), grp(rc[pass_idx, ], "pass_sample"))

# ---- run ----------------------------------------------------------------------
rows <- lapply(seq_len(nrow(scope)), function(q) {
  h <- scope$h[q]; k <- scope$k[q]; rho <- scope$rho[q]
  m1 <- phi2_m1(h, k, rho); m2 <- phi2_m2(h, k, rho)
  m4 <- phi2_m4(h, k, rho); m4c <- phi2_m4(h, k, rho, 1e-48, 1e-31)
  tv <- pmvnorm(upper = c(h, k), corr = matrix(c(1, rho, rho, 1), 2),
                algorithm = TVPACK(abseps = 1e-16))[1]
  tol <- max(1e-30, 1e-20 * abs(asNumeric(m4)))
  data.frame(group = scope$group[q], h = h, k = k, rho = rho, m4 = asNumeric(m4),
             d41 = abs(asNumeric(m4 - m1)), d42 = abs(asNumeric(m4 - m2)),
             d4c = abs(asNumeric(m4 - m4c)), dtv = abs(tv - asNumeric(m4)),
             d_half = scope$d_half[q], tol = tol)
})
out <- do.call(rbind, rows)
out$agree41 <- out$d41 <= out$tol
out$agree42 <- out$d42 <= out$tol
fmt <- function(x) if (is.numeric(x)) sprintf("%.17g", x) else as.character(x)
con <- file(file.path(dir, "phi2tail_m4check.csv"), "wb")
writeLines(paste(names(out), collapse = ","), con, sep = "\n")
writeLines(do.call(paste, c(lapply(out, function(v) vapply(v, fmt, "")), sep = ",")), con, sep = "\n")
close(con)

for (g in unique(out$group)) {
  z <- out[out$group == g, ]
  cat(sprintf("M4|%s|n=%d|M4~M1 agree=%d|M4~M2 agree=%d|max|M4-M1|=%.3e|max|M4-M2|=%.3e|max|M4-M4c|=%.3e|max|TV-M4|=%.3e (%.2f x 2^-52)\n",
              g, nrow(z), sum(z$agree41), sum(z$agree42), max(z$d41), max(z$d42), max(z$d4c),
              max(z$dtv), max(z$dtv) / 2^-52))
}
cat(sprintf("M4|tolerance used: max(1e-30, 1e-20 |M4|) (the reference's own, unchanged)\n"))
cat(sprintf("M4|R %s, Rmpfr %s, MPFR %s, mvtnorm %s, PREC %d\n", getRversion(),
            packageVersion("Rmpfr"), mpfrVersion(), packageVersion("mvtnorm"), PREC))
