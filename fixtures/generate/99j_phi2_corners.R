# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- numerical reference layer: Phi2 corners
#
# The 1,293 distinct bivariate-normal corners (see "Scope" below) whose MPFR
# reference is certified by the
# agreement of two independent methods (user, 2026-10-25;
# E-063):
#   M1  arcsine-form integral, tanh-sinh at step 1/64 (tools/make_phi2_tail_ref.R,
#       "MPFR primitives" section)
#   M4  conditional integral, adaptive 20-point Gauss-Legendre
#       (tools/phi2_tail_refcheck_m4.R, "M4" section)
# Both are loaded from those files the way the M4 script loads M1 (by section
# markers), not copied, so this generator cannot drift from them. Scope is
# exactly the M4-checked set: the 569 half-step failures and the 724-corner
# fixed-seed sample (phi2tail_m4check.csv), and arm C's 3 remaining design
# corners (phi2tail_m4check_armC.csv; its 3 rows are already in the main set). Values are recomputed here at 256 bits;
# the check's CSVs supply only the corners and a double-precision cross-check.
#
#   phi2_corners.csv  one row per corner
#     h, k, rho        the corner, %.17g (exact doubles)
#     ref_hi, ref_lo   M1 as a double-double: hi = round(M1), lo = round(M1 - hi)
#     m4_hi, m4_lo     M4 likewise: certification evidence, not the reference
#     group            half_step_fail | pass_sample (as m4_check)
#     armc             1 if the corner is one of arm C's 3 supplement rows
#     grid             1 if h, k, rho are prereg §3 grid values (after the
#                      snap described under "Scope"); 0 for the 64 off-grid
#                      corners, whose doubles are the parsed check strings
#     chk_gap          |M4 - the check's recorded M4|, diagnostic; enforced
#                      only on grid corners
#     rho_band         |rho| in [0,0.9) | [0.9,0.99) | [0.99,0.999] | (0.999,1)
#                      (descriptive; 0.999 is prereg §3's B, the cert bound)
#     cert             1 if in prereg §3's certification domain (in_cert)
#     dd               full    |hi| >= 2^-969: lo is a normal double
#                      partial 0 < |hi| < 2^-969: lo subnormal or flushed
#                      zero    hi = 0: below double range; only the absolute
#                              tolerance is meaningful
#     d14              |M1 - M4| (MPFR, then rounded), diagnostic
#     tol              max(1e-30, 1e-20 |M4|), the reference tolerance
#     rel_cert         1 if |M1 - M4| <= 1e-20 |M4| (in MPFR): the pair is
#                      certified relatively. 0: certified to 1e-30 ABSOLUTE
#                      only; ref_hi/ref_lo then carry no relative accuracy and
#                      can be negative (e.g. h=-1.75, k=0.25, rho=-0.999:
#                      M1 = -4.5e-52, M4 = 2.0e-250). 642 of 1,293 corners in
#                      m4_check's own figures. Tests must use absolute
#                      comparison on these rows.
#
# Decided (user, 2026-10-28): M1 stays the primary reference (ref_hi/ref_lo);
# M4 is carried alongside as explicit certification evidence (m4_hi/m4_lo),
# not as a substitute, including where rel_cert = 0.
#
# Agreement is asserted in MPFR, before any conversion to double. M1 and M4
# share MPFR's pnorm (and MPFR itself): an error in that routine is NOT
# excluded by this certification (m4_check.md, "What this does not establish").
#
# Own manifest (MANIFEST_phi2_corners.json); fixtures/verify.R needs the
# matching optional entry. Run time roughly 15-30 min (M1 + M4, 1,293 corners;
# the check, which also ran M2, M4c and TVPACK, took 2,105 s).
#
# Run from the repository root:  Rscript fixtures/generate/99j_phi2_corners.R
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))

suppressMessages(library(Rmpfr))
PREC <- 256

OUT <- file.path("fixtures", "expected")
MAN <- file.path("fixtures", "MANIFEST_phi2_corners.json")
CHK <- file.path("docs", "notes", "phi2-tail", "m4_check")
SRC_M1 <- file.path("tools", "make_phi2_tail_ref.R")
SRC_M4 <- file.path("tools", "phi2_tail_refcheck_m4.R")
SRC_GRID <- file.path("tools", "phi2_tail_grid.R")

# ---- M1 and M4, loaded by section ------------------------------------------
section <- function(path, from, to) {
  s <- readLines(path); a <- grep(from, s); b <- grep(to, s) - 1
  stopifnot(length(a) == 1, length(b) == 1, a < b)
  eval(parse(text = s[a:b]), envir = globalenv())
}
section(SRC_M1, "^# ---- MPFR primitives", "^cache <- new.env")
section(SRC_M1, "^hilo <- function", "^# ---- grids and cell design")
section(SRC_M4, "^# ---- M4 ---", "^# ---- scope ---")
smoke <- FALSE; source(SRC_GRID)                  # L, B, G, rho_grid, in_cert

# ---- scope ------------------------------------------------------------------
# Scope accounting (user, 2026-10-28): m4_check.md counted 1,296 = 569 +
# 724 + 3, but the 3 rows of phi2tail_m4check_armC.csv are byte-identical
# corners already in phi2tail_m4check.csv's pass sample. The distinct scope
# in the committed files is 1,293. This generator writes what the files hold,
# flags arm C's rows (armc = 1), and asserts the overlap, so a corrected armC
# file stops it rather than being silently merged. This is a correction of the
# count, not a reason to add corners: the fixture set is the 1,293.
rd <- function(f) read.csv(file.path(CHK, f), colClasses = "character")
main <- rd("phi2tail_m4check.csv")
armc <- rd("phi2tail_m4check_armC.csv")
key <- function(d) paste(d$h, d$k, d$rho, sep = "|")
stopifnot(nrow(main) == 1293, nrow(armc) == 3,
          all(key(armc) %in% key(main)))
sc <- main
sc$armc <- as.integer(key(sc) %in% key(armc))
sc$h <- as.numeric(sc$h); sc$k <- as.numeric(sc$k); sc$rho <- as.numeric(sc$rho)
# The check's CSV was printed on Windows ARM, whose "%.17g" is not correctly
# rounded: some strings parse back to a NEIGHBOURING double (measured
# 2026-10-28: 139 rho values lie 1-4 ulps off a rho_grid value, none on h or
# k; e.g. "-0.99499999999999993" parses to -0.99499999999999988, not -0.995).
# Grid corners came from the grid, so values within 4 ulps of a grid value
# are snapped back to it. The 64 corners off the grid (arm C's fitted rho and
# thresholds) cannot be recovered this way; their corner is the double the
# string parses to, and `grid` = 0 marks them.
snap <- function(x, grid) vapply(x, function(z) {
  g <- grid[abs(z - grid) <= 4 * 2^-52 * max(1, abs(z))]
  if (length(g) == 1) g else z
}, 0)
on_grid <- function(d) d$h %in% G & d$k %in% G & d$rho %in% rho_grid
rho0 <- sc$rho
sc$h <- snap(sc$h, G); sc$k <- snap(sc$k, G); sc$rho <- snap(sc$rho, rho_grid)
sc$grid <- as.integer(on_grid(sc))
stopifnot(sum(sc$rho != rho0) == 139, sum(sc$grid == 0) == 64)
stopifnot(sum(sc$group == "half_step_fail") == 569,
          sum(sc$group == "pass_sample") == 724,
          !anyDuplicated(sprintf("%.17g|%.17g|%.17g", sc$h, sc$k, sc$rho)),
          all(is.finite(c(sc$h, sc$k, sc$rho))), all(abs(sc$rho) < 1),
          all(sc$agree41 == "TRUE"))
if (nzchar(Sys.getenv("PHI2C_SMOKE_N")))          # smoke: first n of each group
  sc <- do.call(rbind, lapply(split(sc, sc$group), head, as.integer(Sys.getenv("PHI2C_SMOKE_N"))))

# ---- compute and certify ----------------------------------------------------
TOL_ABS <- mpfr("1e-30", PREC); TOL_REL <- mpfr("1e-20", PREC)
DD_FULL <- 2^-969
f17 <- function(x) sprintf("%.17g", x)
rows <- lapply(seq_len(nrow(sc)), function(q) {
  h <- sc$h[q]; k <- sc$k[q]; rho <- sc$rho[q]
  m1 <- phi2_m1(h, k, rho); m4 <- phi2_m4(h, k, rho)
  d <- abs(m1 - m4); tol <- max(TOL_ABS, TOL_REL * abs(m4))
  if (!(d <= tol))
    stop(sprintf("M1 and M4 disagree at h=%s k=%s rho=%s: |d|=%.3e tol=%.3e",
                 f17(h), f17(k), f17(rho), asNumeric(d), asNumeric(tol)))
  # On grid corners the recomputed M4 must reproduce the check's recorded
  # value to the reference tolerance, widened by that value's own rounding
  # (its printing is not exact either: up to 4 ulps, as for rho above; e.g.
  # h=-3.5 k=3.5 rho=-0.95 is recorded 1.25 ulps off, measured 2026-10-28). Off the grid the check may have
  # evaluated a neighbouring corner, so the gap is reported, not enforced;
  # the certification there is the M1-M4 agreement above.
  x4 <- mpfr(sc$m4[q], PREC); g4 <- abs(m4 - x4)
  chk_ok <- sc$grid[q] == 0 || g4 <= max(tol, 4 * 2^-52 * abs(x4))
  # Violations are collected and reported together after the loop.
  hl <- hilo(m1); h4 <- hilo(m4)
  stopifnot(hl[1] + hl[2] == hl[1] || hl[1] == 0)   # lo is below hi's last bit
  a <- abs(rho)
  data.frame(h = f17(h), k = f17(k), rho = f17(rho),
             ref_hi = f17(hl[1]), ref_lo = f17(hl[2]),
             m4_hi = f17(h4[1]), m4_lo = f17(h4[2]),
             group = sc$group[q], armc = sc$armc[q], grid = sc$grid[q],
             rho_band = if (a < 0.9) "lt0.9" else if (a < 0.99) "0.9-0.99"
                        else if (a <= B) "0.99-0.999" else "gt0.999",
             cert = as.integer(in_cert(c(h, k), rho)),
             dd = if (hl[1] == 0) "zero" else if (abs(hl[1]) >= DD_FULL) "full" else "partial",
             d14 = f17(asNumeric(d)), tol = f17(asNumeric(tol)),
             rel_cert = as.integer(d <= TOL_REL * abs(m4)),
             chk_gap = f17(asNumeric(g4)), chk_ok = chk_ok)
})
out <- do.call(rbind, rows)
if (any(!out$chk_ok)) {
  print(out[!out$chk_ok, c("h", "k", "rho", "group", "m4_hi", "chk_gap", "tol")])
  stop(sprintf("%d grid corner(s): M4 moved since m4_check (listed above)", sum(!out$chk_ok)))
}
out$chk_ok <- NULL

# m4_check.md: 186 of the failures and 486 of the sample lie in the cert domain.
if (!nzchar(Sys.getenv("PHI2C_SMOKE_N")))
  stopifnot(sum(out$cert[out$group == "half_step_fail"]) == 186,
            sum(out$cert[out$group == "pass_sample"]) == 486)

for (g in unique(out$group)) {
  z <- out[out$group == g, ]
  cat(sprintf("%s: n=%d, cert=%d, rel_cert=%d, ref_hi<0: %d, dd full/partial/zero=%d/%d/%d, max|M1-M4|=%.3e\n",
              g, nrow(z), sum(z$cert), sum(z$rel_cert), sum(as.numeric(z$ref_hi) < 0),
              sum(z$dd == "full"), sum(z$dd == "partial"),
              sum(z$dd == "zero"), max(as.numeric(z$d14))))
}

if (nzchar(Sys.getenv("PHI2C_SMOKE_N"))) {
  print(out); quit(save = "no")
}

f <- file.path(OUT, "phi2_corners.csv")
write_csv_lf(out, f)

write_manifest(MAN, f,
  seeds = "none (corners from m4_check; its pass sample used seed 20261025)",
  extra = list(packages = c(
    R.version.string,
    paste("Rmpfr", as.character(packageVersion("Rmpfr"))),
    paste("MPFR", as.character(mpfrVersion())),
    sprintf("precision %d bits", PREC),
    "reference = M1 (tanh-sinh, arcsine form) as double-double; certified by M4 (adaptive Gauss-Legendre, conditional integral)",
    "M1-M4 agreement asserted in MPFR on every corner: |M1-M4| <= max(1e-30, 1e-20 |M4|)",
    "shared dependency: M1 and M4 both use Rmpfr pnorm (MPFR); an error there is not excluded",
    paste("source md5", names(tools::md5sum(c(SRC_M1, SRC_M4, SRC_GRID))),
          tools::md5sum(c(SRC_M1, SRC_M4, SRC_GRID))),
    paste("scope md5", file.path(CHK, c("phi2tail_m4check.csv", "phi2tail_m4check_armC.csv")),
          tools::md5sum(file.path(CHK, c("phi2tail_m4check.csv", "phi2tail_m4check_armC.csv")))))))

cat(sprintf("%d corners written to %s; %s written.\n", nrow(out), f, MAN))
