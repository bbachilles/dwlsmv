# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- Investigation D's flush events (cell 180)
#
# The 14 pair rows of the 12 inside-D event runs of Investigation D
# (E-100; design E-099
# design.md). In each run a binormal() flush at a usable point moved the
# pairwise solver's status and the class of the low end.
#
# Inputs (committed evidence, E-099):
#   events_insideD.csv  the study's rows: arm 0 (production) and arm 1 (flushed
#                       in-D corners replaced by the instrument's Q), with the
#                       M4 replay verdict (Investigation D, 7bc7e15)
#   pairs.csv           extract_pairs.do: each pair rebuilt from the frozen
#                       generator (seed 2026200000 + 1000*180 + run) and run
#                       through the production solver; thresholds and the
#                       estimate as %21x hex; match = 1 iff arm 0 reproduced
#                       the study's row exactly
#   corners.csv         the same: in-D corners where binormal() is exactly 0,
#                       at each arm's ends and estimate
#
# Outputs:
#   phi2_flush_events.csv   one row per pair
#     run, k, vi, vj          study keys (cell 180); pair k is (vi, vj)
#     a1, a2, b1, b2          thresholds, %.17g (exact doubles)
#     n11..n33                the crosstab; f11..f33 the solver input (after
#                             the empty-cell repair, if any)
#     s0, cl0, ch0, xl0, xh0, rho0   arm 0, production: status, end classes
#                             (1 positive, -1 negative, 0 ambiguous), end x, rho-hat
#     s1, cl1, ch1, xl1, xh1, rho1   arm 1, study evidence (Q corners). NOT
#                             certified as a decision: it ran in the study's
#                             double-precision path and is supported by the
#                             M4 replay (replay_same), not by two methods
#     replay_same             1 if the M4 replay matched arm 1
#   phi2_flush_corners.csv  one row per (pair, point) flushed corner
#     run, k, at              at: xl0, xh0, r0, xl1, xh1 or r1
#     i, j, h, kc, rho        the corner, %.17g
#     ref_hi, ref_lo          Method C as a double-double (the reference)
#     m4_hi, m4_lo            M4 (abs_tol 0) likewise: certification evidence
#     drel                    |C - M4| / M4 in MPFR, diagnostic
#
# Certification (a project rule: the numerical reference layer): every corner is
# certified by M4 (conditional integral, adaptive Gauss-Legendre;
# tools/phi2_tail_refcheck_m4.R) and Method C (Plackett from rho = -1,
# tanh-sinh; tools/phi2_plackett_c.R): different representation, algorithm
# and code. Asserted in MPFR: |C - M4| <= 1e-20 |M4|, relative, on every
# corner (no absolute floor: these corners are 1e-68 to 1e-15 and the
# relative claim is the one needed). M4 is called with abs_tol = 0, since
# its default absolute stop (1e-45) carries no relative accuracy below it.
# M1 is not used: it is absolute-only near 1e-30 and these corners are far
# below. Shared: MPFR itself and Rmpfr's elementary functions (C calls no
# pnorm, so M4's pnorm is not shared).
#
# What this does not certify: the arm-1 DECISION (see s1 above); production
# correctness (arm 0 is a characterization of current production); any
# other cell; a rate.
#
# Own manifest (MANIFEST_phi2_flush.json). Run time about 15 min.
# Run from the repository root:  Rscript fixtures/generate/99l_phi2_flush_events.R
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))

suppressMessages(library(Rmpfr))
PREC <- 256

OUT <- file.path("fixtures", "expected")
MAN <- file.path("fixtures", "MANIFEST_phi2_flush.json")
EVD <- file.path("docs", "notes", "phi2-flush-fixtures")
SRC_M1 <- file.path("tools", "make_phi2_tail_ref.R")
SRC_M4 <- file.path("tools", "phi2_tail_refcheck_m4.R")
SRC_C <- file.path("tools", "phi2_plackett_c.R")

section <- function(path, from, to) {
  s <- readLines(path); a <- grep(from, s); b <- grep(to, s) - 1
  stopifnot(length(a) == 1, length(b) == 1, a < b)
  eval(parse(text = s[a:b]), envir = globalenv())
}
section(SRC_M1, "^# ---- MPFR primitives", "^cache <- new.env")   # one, zero, MPFR helpers
section(SRC_M1, "^hilo <- function", "^# ---- grids and cell design")
section(SRC_M4, "^# ---- M4 ---", "^# ---- scope ---")
section(SRC_C, "^# ---- Method C ---", "^# ---- end Method C")

hex <- function(x) {
  m <- regmatches(x, regexec("^([+-]?)([0-9a-f.]+)X([+-])([0-9a-f]+)$", x))
  vapply(m, function(p) {
    if (length(p) != 5) stop("bad hex: ", paste(p, collapse = ""))
    e <- strtoi(p[5], 16L) * (if (p[4] == "-") -1 else 1)
    (if (p[2] == "-") -1 else 1) * as.numeric(sprintf("0x%sp0", p[3])) * 2^e
  }, 0)
}
f17 <- function(x) sprintf("%.17g", x)
rd <- function(f) read.csv(file.path(EVD, f), colClasses = "character")

# ---- inputs and their consistency -------------------------------------------
ev <- rd("events_insideD.csv")
pr <- rd("pairs.csv")
cn <- rd("corners.csv")
key <- function(d) paste(d$run, d$k)
stopifnot(nrow(ev) == 14, nrow(pr) == 14, all(ev$cell == "180"), all(ev$arm == "1"),
          all(ev$inD == "1"), setequal(key(ev), key(pr)), !anyDuplicated(key(pr)),
          length(unique(ev$run)) == 12, all(pr$match == "1"), all(ev$same == "1"))
ev <- ev[match(key(pr), key(ev)), ]
num <- function(x) as.numeric(x)
# arm 0 as rebuilt equals the study's row, field by field, as doubles
stopifnot(hex(pr$rho) == num(ev$r0), num(pr$status) == num(ev$s0),
          num(pr$cl) == num(ev$cl0), num(pr$ch) == num(ev$ch0),
          hex(pr$xl) == num(ev$xl0), hex(pr$xh) == num(ev$xh0),
          num(pr$f11) == num(ev$n11), num(pr$K1) == 2, num(pr$K2) == 2)
# every pair is a status event between arms
stopifnot(all(num(ev$s0) != num(ev$s1)))

nc <- c(outer(1:3, 1:3, function(r, s) paste0("n", r, s)))
nc <- sort(nc); fc <- sub("^n", "f", nc)
ev_out <- data.frame(run = pr$run, k = pr$k, vi = pr$vi, vj = pr$vj,
  a1 = f17(hex(pr$a1)), a2 = f17(hex(pr$a2)), b1 = f17(hex(pr$b1)), b2 = f17(hex(pr$b2)),
  pr[, nc], pr[, fc],
  s0 = pr$status, cl0 = pr$cl, ch0 = pr$ch,
  xl0 = f17(hex(pr$xl)), xh0 = f17(hex(pr$xh)), rho0 = f17(hex(pr$rho)),
  s1 = ev$s1, cl1 = ev$cl1, ch1 = ev$ch1,
  xl1 = f17(num(ev$xl1)), xh1 = f17(num(ev$xh1)), rho1 = f17(num(ev$r1)),
  replay_same = ev$same, stringsAsFactors = FALSE)

# ---- certify the corners -----------------------------------------------------
TOL_REL <- mpfr("1e-20", PREC)
stopifnot(nrow(cn) >= 1, all(key(cn) %in% key(pr)))
rows <- lapply(seq_len(nrow(cn)), function(q) {
  h <- hex(cn$h[q]); k <- hex(cn$kk[q]); rho <- hex(cn$rho[q])
  m4 <- phi2_m4(h, k, rho, abs_tol = 0, rel_tol = 1e-30)
  cc <- phi2_c(h, k, rho)$v
  d <- abs(cc - m4) / m4
  if (!(m4 > 0 && d <= TOL_REL))
    stop(sprintf("C and M4 disagree at h=%s k=%s rho=%s: rel %.3e", f17(h), f17(k), f17(rho),
                 asNumeric(d)))
  hc <- hilo(cc); h4 <- hilo(m4)
  stopifnot(hc[1] > 0, hc[1] + hc[2] == hc[1])
  data.frame(run = cn$run[q], k = cn$k[q], at = cn$at[q], i = cn$i[q], j = cn$j[q],
             h = f17(h), kc = f17(k), rho = f17(rho),
             ref_hi = f17(hc[1]), ref_lo = f17(hc[2]), m4_hi = f17(h4[1]), m4_lo = f17(h4[2]),
             drel = sprintf("%.3e", asNumeric(d)), stringsAsFactors = FALSE)
})
co_out <- do.call(rbind, rows)
cat(sprintf("corners %d; max drel %s; ref range %s to %s\n", nrow(co_out),
            max(as.numeric(co_out$drel)), min(as.numeric(co_out$ref_hi)),
            max(as.numeric(co_out$ref_hi))))

f1 <- file.path(OUT, "phi2_flush_events.csv")
f2 <- file.path(OUT, "phi2_flush_corners.csv")
write_csv_lf(ev_out, f1)
write_csv_lf(co_out, f2)

srcs <- c(SRC_M1, SRC_M4, SRC_C)
evid <- file.path(EVD, c("events_insideD.csv", "pairs.csv", "corners.csv", "extract_pairs.do"))
write_manifest(MAN, c(f1, f2),
  seeds = "study seeds 2026200000 + 1000*180 + run (Investigation D, frozen 05f34e0); none here",
  extra = list(packages = c(
    R.version.string,
    paste("Rmpfr", as.character(packageVersion("Rmpfr"))),
    paste("MPFR", as.character(mpfrVersion())),
    sprintf("precision %d bits", PREC),
    "pair inputs and arm 0: extract_pairs.do under StataMP 18 at the generator's commit; src/ identical to 05f34e0",
    "certified scope (user, 2026-11-09): the 56 corner VALUES only (C and M4 independently, max rel gap 4.4e-43, required 1e-20)",
    "NOT certified: the arm-1 decision (s1, cl1, ch1, ends, rho1). It came from the study's double-precision Mata solver with the true corners substituted; its only support is the M4 replay, about 1e-12. C certified that decision's inputs and did not reproduce the decision",
    "arm 0: a characterization of current production, not a correctness claim",
    "corner reference = Method C (Plackett from rho=-1, tanh-sinh) as double-double; certified by M4 (abs_tol 0)",
    "C-M4 agreement asserted in MPFR on every corner: |C-M4| <= 1e-20 |M4|",
    "shared: MPFR and Rmpfr elementary functions; pnorm not shared (C uses none)",
    paste("source md5", srcs, tools::md5sum(srcs)),
    paste("evidence md5", evid, tools::md5sum(evid)))))

cat(sprintf("%d pairs to %s, %d corners to %s; %s written.\n", nrow(ev_out), f1, nrow(co_out), f2, MAN))
