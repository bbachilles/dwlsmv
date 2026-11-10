# probe_binormal_ref.R -- MPFR references and summaries for the delta_Phi
# phase 1 probes (E-085). Tool data, NOT a
# fixture. Phase 1 has no pass/fail: every line is a diagnostic.
#
# EXECUTED 2026-10-30 (R 4.6.0, Rmpfr 1.1.3); results: E-085.
#
# Reads pb_points.csv written by tools/probe_binormal.do; writes pb_ref.csv
# beside it and prints PBR| lines. Run from the repository root:
#   PB_OUT=<dir holding pb_points.csv> Rscript tools/probe_binormal_ref.R
#
# Comparators (design note section 3):
#   M1  arcsine-form integrand, tanh-sinh, 256-bit MPFR, taken unchanged from
#       tools/make_phi2_tail_ref.R. Every row with want_ref = 1.
#   M2  the same integrand over the complementary range: shares M1's
#       representation, so M1 ~ M2 is a DIAGNOSTIC only. Same rows.
#   M4  conditional one-dimensional integral, adaptive Gauss-Legendre, taken
#       unchanged from tools/phi2_tail_refcheck_m4.R: differs from M1 in
#       representation and algorithm. Every P1 row, and every 20th of the
#       others. Shared by all three: Rmpfr's pnorm (MPFR erfc).
# The error is (stata - ref_hi) - ref_lo, ref = M4 where computed, else M1.
# A reference below the absolute floor M1 can resolve (about 1e-60 here,
# cancellation in Phi(h)Phi(k) + integral) is taken from M4 only; rows where
# M4 was not computed and M1 is below 1e-50 are counted, not used.

suppressMessages(library(Rmpfr))
PREC <- 256
dir <- Sys.getenv("PB_OUT", ".")

src <- readLines("tools/make_phi2_tail_ref.R")
a <- grep("^# ---- MPFR primitives", src); b <- grep("^cache <- new.env", src) - 1
eval(parse(text = src[a:b]))
src4 <- readLines("tools/phi2_tail_refcheck_m4.R")
a <- grep("^# ---- M4 ", src4); b <- grep("^# ---- scope", src4) - 1
eval(parse(text = src4[a:b]))

pts <- read.csv(file.path(dir, "pb_points.csv"), colClasses = c(probe = "character",
                tag = "character"), na.strings = ".")
pts$id <- seq_len(nrow(pts))
use <- pts$want_ref == 1
use_m4 <- use & (grepl("^P1", pts$probe) | (pts$id %% 20 == 0))

hilo <- function(x) { hi <- asNumeric(x); c(hi, asNumeric(x - hi)) }
res <- lapply(which(use), function(i) {
  h <- pts$h[i]; k <- pts$k[i]; r <- pts$rho[i]
  m1 <- phi2_m1(h, k, r); m2 <- phi2_m2(h, k, r)
  m4 <- if (use_m4[i]) phi2_m4(h, k, r, abs_tol = 0) else NULL
  ref <- if (!is.null(m4)) m4 else m1
  hl <- hilo(ref)
  data.frame(id = i, ref_hi = hl[1], ref_lo = hl[2], src = if (is.null(m4)) "M1" else "M4",
             d12 = abs(asNumeric(m1 - m2)),
             d14 = if (is.null(m4)) NA else abs(asNumeric(m1 - m4)),
             rel14 = if (is.null(m4) || m4 == 0) NA else abs(asNumeric((m1 - m4) / m4)))
})
ref <- merge(pts, do.call(rbind, res), by = "id")
ref$unresolved <- ref$src == "M1" & ref$ref_hi < 1e-50
ref$err <- (ref$stata - ref$ref_hi) - ref$ref_lo
ref$err[ref$unresolved] <- NA
ref$abs52 <- abs(ref$err) / 2^-52
ref$rel <- ifelse(ref$ref_hi > 0, abs(ref$err) / ref$ref_hi, NA)

fmt <- function(x) if (is.numeric(x)) ifelse(is.na(x), "", sprintf("%.17g", x)) else as.character(x)
con <- file(file.path(dir, "pb_ref.csv"), "wb")
writeLines(paste(names(ref), collapse = ","), con, sep = "\n")
writeLines(do.call(paste, c(lapply(ref, fmt), sep = ",")), con, sep = "\n")
close(con)

mx <- function(x) if (all(is.na(x))) NA else max(x, na.rm = TRUE)
cat(sprintf("PBR|refs=%d M4=%d unresolved=%d | max|M1-M2|=%.3g max|M1-M4|=%.3g max rel|M1-M4|=%.3g (diagnostics)\n",
            nrow(ref), sum(ref$src == "M4"), sum(ref$unresolved), mx(ref$d12), mx(ref$d14), mx(ref$rel14)))

# P1: the true value at each flush boundary (last nonzero point, first zero point)
p1 <- ref[grepl("^P1", ref$probe), ]
for (side in c("P1a", "P1b")) {
  z <- p1[p1$probe == side, ]
  cat(sprintf("PBR|%s|n=%d|stata==0: %d|true P range [%.3g, %.3g]|inside D: %d\n", side, nrow(z),
              sum(z$stata == 0), min(z$ref_hi), max(z$ref_hi), sum(z$inD == 1)))
}
b <- p1[p1$probe == "P1b" & p1$stata == 0, ]          # first flushed point per boundary
b$rho_r <- round(b$rho, 3)
for (r in sort(unique(b$rho_r))) {
  z <- b[b$rho_r == r, ]
  cat(sprintf("PBR|P1|rho=%6.3f|flushed true P: min %.3g max %.3g|inside D: %d\n",
              r, min(z$ref_hi), max(z$ref_hi), sum(z$inD == 1)))
}

# P2: absolute or relative error near the flush
p2 <- ref[ref$probe == "P2" & !is.na(ref$err), ]
cat(sprintf("PBR|P2|n=%d|max abs err=%.3f x 2^-52|max rel err=%.3g|max lowbit exponent=%s\n",
            nrow(p2), mx(p2$abs52), mx(p2$rel), mx(p2$lowbit)))

# P3: error against rho, by 0.05 bin and by band
p3 <- ref[grepl("^P3", ref$probe) & !is.na(ref$err), ]
p3$bin <- floor(p3$rho / 0.05) * 0.05
for (tg in unique(p3$tag)) {
  z <- p3[p3$tag == tg & p3$probe == "P3u", ]
  s <- tapply(z$abs52, z$bin, mx)
  cat(sprintf("PBR|P3u|%s|(h,k)=(%g,%g)|max abs err x 2^-52 by rho bin: %s\n", tg, z$h[1], z$k[1],
              paste(sprintf("%.2f:%.2f", as.numeric(names(s)), s), collapse = " ")))
  for (c0 in c(-0.925, -0.75, -0.3, 0.3, 0.75, 0.925)) {
    w <- p3[p3$tag == tg & p3$probe == "P3b" & abs(p3$rho - c0) <= 0.0020001, ]
    lo <- w[w$rho < c0, ]; hi <- w[w$rho >= c0, ]
    cat(sprintf("PBR|P3b|%s|band %+.3f|max abs err x 2^-52 below %.3f above %.3f\n",
                tg, c0, mx(lo$abs52), mx(hi$abs52)))
  }
}

# P4: symmetry flags, by probe
for (p in c("P1a", "P1b", "P2", "P3u", "P3b", "P5")) {
  z <- pts[pts$probe == p, ]
  if (!nrow(z)) next
  cat(sprintf("PBR|P4|%s|n=%d|swap=%d refl1=%d refl2=%d refl3=%d (bitwise equal)\n", p, nrow(z),
              sum(z$swap_eq), sum(z$refl1_eq), sum(z$refl2_eq), sum(z$refl3_eq)))
}

cat(sprintf("PBR|R %s, Rmpfr %s, MPFR %s, PREC %d\n", getRversion(),
            packageVersion("Rmpfr"), mpfrVersion(), PREC))
