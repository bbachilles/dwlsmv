# phase2_binormal_ref.R -- delta_Phi phase 2 references and gates
# (docs/notes/deltaphi-phase2/prereg.md sections 3, 3a, 5). Tool data, NOT a
# fixture.
#
# STUDY RUN NOT YET EXECUTED. Drafted 2026-10-31; smoke-run 2026-10-31 (docs/notes/deltaphi-phase2/smoke/).
#
# Reads p2_points.csv and p2_counts.csv written by tools/phase2_binormal.do;
# writes p2_ref.csv beside them and prints P2R| lines. Run from the repo root:
#   P2_OUT=<dir> [P2_CORES=<n>] Rscript tools/phase2_binormal_ref.R
#
# Comparators (3a.1): M4 (conditional 1-D integral, adaptive Gauss-Legendre,
# 256-bit MPFR, unchanged from tools/phi2_tail_refcheck_m4.R) DECIDES both
# gates. M1 (arcsine form, tanh-sinh, unchanged from tools/make_phi2_tail_ref.R)
# is a diagnostic: |M1 - M4| is reported and never decides.
#
# Gates (3, 5; c = 1e-15, e = 2.5 x 2^-52, compared in MPFR):
#   G-in  (Claim 3): every G-in designated point has M4 <= c + e.
#   G-out (Claim 2): every G-out point has |(stata - M4_hi) - M4_lo| <= e.
#   A point M4 cannot resolve (non-finite or negative) is a domain failure, never
#   an exclusion; the test is made in MPFR, so a true P below double range
#   (about 1e-2700 at rho = -0.999) counts as resolved. A gate whose own quota fell short (an ordinary G-in band
#   below its quota at its cap, or a G-out stratum below quota) is reported
#   as floor-bitten; the terminal band never is.
# Verdict per gate: FAIL if any point violates or is unresolved; otherwise
# INCONCLUSIVE if its floor bit; otherwise PASS. (Fail before inconclusive:
# prereg section 5, decided by the user 2026-10-31.)
# Rows with stratum "forced" (the do-file's smoke-only stop-rule check) abort
# the run.

suppressMessages({ library(Rmpfr); library(parallel) })
PREC <- 256
dir <- Sys.getenv("P2_OUT", ".")
cores <- as.integer(Sys.getenv("P2_CORES", "1"))

load_refs <- function() {
  src <- readLines("tools/make_phi2_tail_ref.R")
  a <- grep("^# ---- MPFR primitives", src); b <- grep("^cache <- new.env", src) - 1
  eval(parse(text = src[a:b]), envir = globalenv())
  src4 <- readLines("tools/phi2_tail_refcheck_m4.R")
  a <- grep("^# ---- M4 ", src4); b <- grep("^# ---- scope", src4) - 1
  eval(parse(text = src4[a:b]), envir = globalenv())
  invisible(NULL)
}
load_refs()

pts <- read.csv(file.path(dir, "p2_points.csv"), colClasses = c(gate = "character",
                stratum = "character"))
cnt <- read.csv(file.path(dir, "p2_counts.csv"), colClasses = c(gate = "character",
                stratum = "character"))
if (any(pts$stratum == "forced")) stop("forced (smoke-only) rows in p2_points.csv")
if (any(pts$gate == "Gin" & pts$stata != 0)) stop("a G-in designated point is not flushed")

one_pt <- function(i) {
  h <- pts$h[i]; k <- pts$k[i]; r <- pts$rho[i]
  m4 <- phi2_m4(h, k, r, abs_tol = 0)
  m1 <- phi2_m1(h, k, r)
  hi <- asNumeric(m4)
  c(hi, asNumeric(m4 - hi), asNumeric(abs(m1 - m4)), !is.finite(m4) || m4 < 0)
}
if (cores > 1) {
  cl <- makeCluster(cores)
  clusterExport(cl, c("PREC", "pts", "load_refs", "one_pt"))
  clusterEvalQ(cl, { suppressMessages(library(Rmpfr)); load_refs() })
  res <- parLapply(cl, seq_len(nrow(pts)), one_pt)
  stopCluster(cl)
} else res <- lapply(seq_len(nrow(pts)), one_pt)
res <- do.call(rbind, res)
pts$m4_hi <- res[, 1]; pts$m4_lo <- res[, 2]; pts$d14 <- res[, 3]
pts$unresolved <- res[, 4] == 1   # tested in MPFR: a value below double range is resolved

cpe <- mpfr(1e-15, PREC) + mpfr(2.5, PREC) * mpfr(2, PREC)^-52
e <- mpfr(2.5, PREC) * mpfr(2, PREC)^-52
p4 <- function(i) mpfr(pts$m4_hi[i], PREC) + mpfr(pts$m4_lo[i], PREC)
pts$viol <- FALSE
pts$err52 <- NA
for (i in seq_len(nrow(pts))) {
  if (pts$unresolved[i]) next
  if (pts$gate[i] == "Gin") {
    pts$viol[i] <- p4(i) > cpe
  } else {
    er <- mpfr(pts$stata[i], PREC) - p4(i)
    pts$viol[i] <- abs(er) > e
    pts$err52[i] <- asNumeric(er / mpfr(2, PREC)^-52)
  }
}

fmt <- function(x) if (is.numeric(x)) ifelse(is.na(x), "", sprintf("%.17g", x)) else as.character(x)
con <- file(file.path(dir, "p2_ref.csv"), "wb")
writeLines(paste(names(pts), collapse = ","), con, sep = "\n")
writeLines(do.call(paste, c(lapply(pts, fmt), sep = ",")), con, sep = "\n")
close(con)

verdict <- function(z, short) {
  nv <- sum(z$viol); nu <- sum(z$unresolved)
  if (nv + nu > 0) "FAIL" else if (short) "INCONCLUSIVE" else "PASS"
}
gi <- pts[pts$gate == "Gin", ]
ci <- cnt[cnt$gate == "Gin", ]
short_in <- any(ci$terminal == 0 & ci$accepted < ci$quota)
cat(sprintf("P2R|G-in|points=%d violations(M4 > c+e)=%d unresolved=%d max M4=%.6g (c+e=%.6g)\n",
            nrow(gi), sum(gi$viol), sum(gi$unresolved),
            if (nrow(gi)) max(gi$m4_hi) else NA, asNumeric(cpe)))
cat(sprintf("P2R|G-in|ordinary bands short=%d; terminal band accepted=%s (reduced coverage)\n",
            sum(ci$terminal == 0 & ci$accepted < ci$quota),
            paste(ci$accepted[ci$terminal == 1], collapse = ",")))
cat(sprintf("P2R|G-in|VERDICT %s\n", verdict(gi, short_in)))

go <- pts[pts$gate == "Gout", ]
co <- cnt[cnt$gate == "Gout", ]
short_out <- any(co$accepted < co$quota)
for (s in unique(go$stratum)) {
  z <- go[go$stratum == s, ]
  cat(sprintf("P2R|G-out|%s points=%d violations=%d unresolved=%d max |err|=%.3f x 2^-52\n", s,
              nrow(z), sum(z$viol), sum(z$unresolved), max(abs(z$err52), na.rm = TRUE)))
}
cat(sprintf("P2R|G-out|strata short=%d\n", sum(co$accepted < co$quota)))
cat(sprintf("P2R|G-out|VERDICT %s\n", verdict(go, short_out)))

cat(sprintf("P2R|diagnostic|max |M1 - M4| = %.3g over %d points\n", max(pts$d14), nrow(pts)))
cat(sprintf("P2R|R %s, Rmpfr %s, MPFR %s, PREC %d, cores %d\n", getRversion(),
            packageVersion("Rmpfr"), mpfrVersion(), PREC, cores))
