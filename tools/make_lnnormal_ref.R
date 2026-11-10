# make_lnnormal_ref.R -- MPFR reference for the lnnormal()/normal() study (A-L)
# (docs/notes/lnnormal-cert/prereg.md).
#
# EXECUTED on the frozen design: StataMP 18 and 16.1 gates (docs/notes/lnnormal-cert/result.md).
# Banner changed 2026-10-30, after execution; frozen hash in prereg §6.
#
# Usage: LNREF_OUT=<dir> [LNREF_SMOKE=1] Rscript tools/make_lnnormal_ref.R
# Writes <dir>/lnnormal_ref.csv and <dir>/normal_ref.csv, and prints REF| lines.
#
# Columns: comp (design component), beta_x (Stata %21x, authoritative), beta_d (%.17g),
# then for lnnormal: L_hi_x, L_lo_x (L = ln Phi(beta) = hi + lo); for normal: P_hi_x and
# P_lorel_x, with P = hi * (1 + lorel), since lo itself can be subnormal there.
# A reference value below 2^-1022 in magnitude is sent as 0: that changes a normalized
# error u(1 + |L|) by at most 2^-1022 / u, about 2e-292 (prereg §3).

suppressPackageStartupMessages(library(Rmpfr))
PREC <- 256
out <- Sys.getenv("LNREF_OUT")
if (out == "") stop("set LNREF_OUT")
smoke <- Sys.getenv("LNREF_SMOKE") == "1"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ---- adjacent doubles, exactly (prereg §3) ----------------------------------
binade <- function(x) {                       # e with 2^e <= x < 2^(e+1), x > 0 normal
  e <- floor(log2(x))
  while (2^e > x) e <- e - 1
  while (2^(e + 1) <= x) e <- e + 1
  e
}
succ <- function(x) {
  if (x < 0) return(-pred(-x))
  if (x == 0) stop("succ(0) is subnormal; excluded")
  x + 2^(binade(x) - 52)
}
pred <- function(x) {
  if (x < 0) return(-succ(-x))
  if (x == 0) stop("pred(0) is subnormal; excluded")
  e <- binade(x)
  if (x == 2^e) x - 2^(e - 53) else x - 2^(e - 52)
}
nbhd <- function(x) {                         # x and its three nearest doubles each side
  if (x == 0) return(0)                       # neighbours of 0 are subnormal (excluded)
  up <- x; dn <- x; v <- x
  for (j in 1:3) { up <- succ(up); dn <- pred(dn); v <- c(v, up, dn) }
  v
}
check_adjacent <- function(x) {               # pred(x) < x < succ(x), nothing in between
  if (x == 0) return(TRUE)
  p <- pred(x); s <- succ(x)
  p < x && x < s && ((p + x) / 2) %in% c(p, x) && ((x + s) / 2) %in% c(x, s)
}

# ---- Stata %21x -------------------------------------------------------------
hex21 <- function(x) {
  if (!is.finite(x)) stop("non-finite value to hex21")
  if (x == 0) return("+0.0000000000000X-3ff")
  a <- abs(x)
  if (a < 2^-1022) stop("subnormal value to hex21")
  e <- binade(a)
  v <- (a / 2^e - 1) * 2^52                   # exact integer < 2^52
  d <- character(13)
  for (k in 13:1) { d[k] <- sprintf("%x", as.integer(v %% 16)); v <- (v - v %% 16) / 16 }
  sprintf("%s1.%sX%s%03x", if (x < 0) "-" else "+", paste(d, collapse = ""),
          if (e < 0) "-" else "+", abs(e))
}
flush <- function(x) if (abs(x) < 2^-1022) 0 else x

# ---- the design (prereg §3) -------------------------------------------------
Lmin <- -200; Lmax <- 200
one <- mpfr(1, PREC)
lnPhi_m <- function(b) {                      # MPFR ln Phi; log1p(-Phi(-b)) above 0
  bm <- mpfr(b, PREC)
  if (b > 0) log1p(-pnorm(-bm)) else log(pnorm(bm))
}
target <- -1022 * log(mpfr(2, PREC))
bmin <- uniroot(function(b) pnorm(b, log.p = TRUE) + 1022 * log(2), c(-40, -35), tol = 1e-14)$root
while (log(pnorm(mpfr(bmin, PREC))) < target) bmin <- succ(bmin)
while (log(pnorm(mpfr(pred(bmin), PREC))) >= target) bmin <- pred(bmin)

design <- function(lo, hi, seed) {
  if (smoke) {
    set.seed(20261030, kind = "Mersenne-Twister", normal.kind = "Inversion")
    pts <- list(grid = seq(lo, hi, by = 0.37), random = runif(200, lo, hi),
                integer = unlist(lapply(c(-3.3, 1.7), nbhd)),     # off-design, to exercise
                binade = unlist(lapply(c(-0.75, 0.5), nbhd)),     # the adjacency code (0.5 is
                handoff = unlist(lapply(c(0, -0.125), nbhd)))     # a power of 2: pred halves)
  } else {
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
    ints <- seq(ceiling(lo), floor(hi))
    bins <- unlist(lapply(-10:7, function(e) c(2^e, -2^e)))
    bins <- bins[bins >= lo & bins <= hi]
    pts <- list(grid = seq(lo, hi, by = 2^-6), random = runif(20000, lo, hi),
                integer = unlist(lapply(ints, nbhd)), binade = unlist(lapply(bins, nbhd)),
                handoff = unlist(lapply(c(0, bmin), nbhd)))
  }
  d <- do.call(rbind, lapply(names(pts), function(nm) {
    v <- pts[[nm]]
    if (length(v)) data.frame(comp = nm, beta = v) else NULL
  }))
  d[d$beta >= lo & d$beta <= hi, ]
}
dL <- design(Lmin, Lmax, 20261029)
dP <- design(bmin, 0, 20261029)
adj <- unique(c(dL$beta[dL$comp %in% c("integer", "binade", "handoff")],
                dP$beta[dP$comp %in% c("integer", "binade", "handoff")]))
stopifnot(all(vapply(adj, check_adjacent, logical(1))))
stopifnot(all(dL$beta == 0 | abs(dL$beta) >= 2^-1022), all(dP$beta == 0 | abs(dP$beta) >= 2^-1022))

# ---- the reference ----------------------------------------------------------
refL <- t(vapply(dL$beta, function(b) {
  L <- lnPhi_m(b); hi <- asNumeric(L); lo <- asNumeric(L - hi)
  c(flush(hi), flush(lo))
}, numeric(2)))
refP <- t(vapply(dP$beta, function(b) {
  P <- pnorm(mpfr(b, PREC)); hi <- asNumeric(P)
  c(hi, asNumeric((P - hi) / hi))
}, numeric(2)))
stopifnot(all(is.finite(refL)), all(is.finite(refP)), all(refP[, 1] >= 2^-1022))

wr <- function(file, d, ref, nm) {
  lines <- c(paste(c("comp", "beta_x", "beta_d", nm), collapse = ","),
             sprintf("%s,%s,%s,%s,%s", d$comp, vapply(d$beta, hex21, ""), sprintf("%.17g", d$beta),
                     vapply(ref[, 1], hex21, ""), vapply(ref[, 2], hex21, "")))
  con <- file(file, "wb"); writeLines(lines, con, sep = "\n"); close(con)
}
wr(file.path(out, "lnnormal_ref.csv"), dL, refL, c("L_hi_x", "L_lo_x"))
wr(file.path(out, "normal_ref.csv"), dP, refP, c("P_hi_x", "P_lorel_x"))

# ---- diagnostics (prereg §3a; never decide) ---------------------------------
u <- 2^-53
# (a) Mills-ratio series for beta <= -40: ln Phi = ln phi - ln|b| + ln sum_k (-1)^k (2k-1)!! / b^(2k)
mills <- function(b, K = 12) {
  bm <- mpfr(b, PREC); s <- one; t <- one
  for (k in 1:K) { t <- -t * (2 * k - 1) / (bm * bm); s <- s + t }
  -bm * bm / 2 - log(sqrt(2 * Const("pi", PREC))) - log(-bm) + log(s)
}
iM <- which(dL$beta <= -40)
dm <- if (length(iM)) max(vapply(iM, function(i) {
  asNumeric(abs(mills(dL$beta[i]) - (mpfr(refL[i, 1], PREC) + refL[i, 2]))) / (u * (1 + abs(refL[i, 1])))
}, 0)) else NA
# (b) R's pnorm(log.p = TRUE)
dr <- max(abs((pnorm(dL$beta, log.p = TRUE) - refL[, 1]) - refL[, 2]) / (u * (1 + abs(refL[, 1]))))

cat(sprintf("REF|smoke=%s|lnnormal_points=%d|normal_points=%d|beta_min=%s (%.17g)\n",
            smoke, nrow(dL), nrow(dP), hex21(bmin), bmin))
for (nm in unique(dL$comp)) cat(sprintf("REF|L|comp=%s|n=%d\n", nm, sum(dL$comp == nm)))
for (nm in unique(dP$comp)) cat(sprintf("REF|P|comp=%s|n=%d\n", nm, sum(dP$comp == nm)))
cat(sprintf("REF|diag|mills_vs_ref_max_norm=%.3e (beta <= -40, n=%d)|R_pnorm_log_max_norm=%.3e\n",
            dm, length(iM), dr))
cat(sprintf("REF|%s, Rmpfr %s, MPFR %s, PREC %d\n", R.version.string,
            packageVersion("Rmpfr"), mpfrVersion(), PREC))
