# make_normal_ref_r2.R -- MPFR reference for the revised normal() study (A-L revision 2)
# (docs/notes/lnnormal-cert-r2/prereg.md). An extension of make_lnnormal_ref.R, whose frozen
# file and hash are untouched; its helpers are copied here, unchanged.
#
# EXECUTED on the frozen design: StataMP 18 and 16.1 gates (docs/notes/lnnormal-cert-r2/result.md).
# Banner changed 2026-10-30, after execution; frozen hash in prereg §6.
#
# Usage: LNREF_OUT=<dir> LNREF_R1=<frozen normal_ref.csv> [LNREF_SMOKE=1] \
#        Rscript tools/make_normal_ref_r2.R
# Writes <dir>/normal_ref_r2.csv (same columns as the frozen normal_ref.csv) and REF| lines.

suppressPackageStartupMessages(library(Rmpfr))
PREC <- 256
out <- Sys.getenv("LNREF_OUT")
r1 <- Sys.getenv("LNREF_R1")
if (out == "" || r1 == "") stop("set LNREF_OUT and LNREF_R1")
smoke <- Sys.getenv("LNREF_SMOKE") == "1"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ---- helpers, copied unchanged from make_lnnormal_ref.R ---------------------
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
check_adjacent <- function(x) {               # pred(x) < x < succ(x), nothing in between
  if (x == 0) return(TRUE)
  p <- pred(x); s <- succ(x)
  p < x && x < s && ((p + x) / 2) %in% c(p, x) && ((x + s) / 2) %in% c(x, s)
}
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

# ---- new in r2 --------------------------------------------------------------
unhex21 <- function(s) {                      # exact inverse of hex21 (normal or zero)
  m <- regmatches(s, regexec("^([+-])([01])\\.([0-9a-f]{13})X([+-])([0-9a-f]{3})$", s))[[1]]
  if (length(m) != 6) stop("bad %21x: ", s)
  if (m[3] == "0") return(0)
  f <- 0
  for (ch in strsplit(m[4], "")[[1]]) f <- f * 16 + strtoi(ch, 16L)
  e <- strtoi(m[6], 16L) * (if (m[5] == "-") -1 else 1)
  (if (m[2] == "-") -1 else 1) * (1 + f / 2^52) * 2^e
}
ring <- function(x, from, to) {               # the from-th..to-th nearest doubles each side
  up <- x; dn <- x; v <- numeric(0)
  for (j in 1:to) {
    up <- succ(up); dn <- pred(dn)
    if (j >= from) v <- c(v, up, dn)
  }
  v
}
above <- function(x, n) {                     # x and its n nearest doubles above it
  v <- x
  for (j in 1:n) { x <- succ(x); v <- c(v, x) }
  v
}

# beta_h: the largest double with Phi(beta_h) <= 2^-55 (@sec-tphi R-h), found in MPFR
t55 <- mpfr(2, PREC)^-55
bh <- qnorm(2^-55)
while (pnorm(mpfr(bh, PREC)) > t55) bh <- pred(bh)
while (pnorm(mpfr(succ(bh), PREC)) <= t55) bh <- succ(bh)
stopifnot(sprintf("%.17g", bh) == "-8.3743889230674569")

lo <- bh; hi <- 0
if (smoke) {
  set.seed(20261032, kind = "Mersenne-Twister", normal.kind = "Inversion")
  pts <- list(grid = seq(lo + 0.011, hi, by = 0.29), random = runif(200, lo, hi),
              integer = unlist(lapply(c(-6.6, -2.2), ring, from = 4, to = 6)),  # off-design
              binade = unlist(lapply(c(-0.375, -2), ring, from = 4, to = 6)),   # -2: pred halves
              handoff = c(above(-5.5, 6), -3))        # -3 is a frozen point: must be removed
} else {
  set.seed(20261031, kind = "Mersenne-Twister", normal.kind = "Inversion")
  pts <- list(grid = seq(lo + 2^-10, hi, by = 2^-8), random = runif(20000, lo, hi),
              integer = unlist(lapply(-8:-1, ring, from = 4, to = 6)),
              binade = unlist(lapply(-2^(-10:3), ring, from = 4, to = 6)),
              handoff = above(bh, 6))
}
d <- do.call(rbind, lapply(names(pts), function(nm) data.frame(comp = nm, beta = pts[[nm]])))
d <- d[d$beta >= lo & d$beta <= hi, ]
adj <- unique(d$beta[d$comp %in% c("integer", "binade", "handoff")])
stopifnot(all(vapply(adj, check_adjacent, logical(1))))

# remove any point evaluated in the frozen design (prereg §2-3), and report the count
r1x <- read.csv(r1, colClasses = "character")$beta_x
r1b <- vapply(r1x, unhex21, 0, USE.NAMES = FALSE)
stopifnot(identical(vapply(r1b, hex21, ""), r1x))       # the parser is exact on every point
seen <- d$beta %in% r1b
removed <- table(factor(d$comp[seen], levels = names(pts)))
d <- d[!seen, ]
stopifnot(nrow(d) > 0)

refP <- t(vapply(d$beta, function(b) {
  P <- pnorm(mpfr(b, PREC)); h <- asNumeric(P)
  c(h, asNumeric((P - h) / h))
}, numeric(2)))
stopifnot(all(is.finite(refP)), all(refP[, 1] >= 2^-1022))
lines <- c("comp,beta_x,beta_d,P_hi_x,P_lorel_x",
           sprintf("%s,%s,%s,%s,%s", d$comp, vapply(d$beta, hex21, ""), sprintf("%.17g", d$beta),
                   vapply(refP[, 1], hex21, ""), vapply(refP[, 2], hex21, "")))
con <- file(file.path(out, "normal_ref_r2.csv"), "wb"); writeLines(lines, con, sep = "\n"); close(con)

u <- 2^-53
dr <- max(abs(((pnorm(d$beta) - refP[, 1]) / refP[, 1]) - refP[, 2]) / u)
cat(sprintf("REF|smoke=%s|normal_points=%d|beta_h=%s (%.17g)|frozen_points_read=%d\n",
            smoke, nrow(d), hex21(bh), bh, length(r1b)))
for (nm in names(pts)) cat(sprintf("REF|P|comp=%s|n=%d|removed_as_seen=%d\n",
                                   nm, sum(d$comp == nm), removed[[nm]]))
cat(sprintf("REF|diag|R_pnorm_max_norm=%.3e\n", dr))
cat(sprintf("REF|%s, Rmpfr %s, MPFR %s, PREC %d\n", R.version.string,
            packageVersion("Rmpfr"), mpfrVersion(), PREC))
