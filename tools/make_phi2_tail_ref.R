# make_phi2_tail_ref.R -- MPFR reference for the Phi2 tail-validation study
# (docs/notes/phi2-tail/prereg.md). Tool data, NOT a fixture.
#
# NOT YET EXECUTED on the frozen grid. Smoke-run only (PHI2T_SMOKE=1).
#
# Usage, from the repository root:
#   Rscript tools/make_phi2_tail_ref.R
# Environment: PHI2T_OUT (output directory, default "."); PHI2T_SMOKE=1 for
# the smoke grid (off-design points, synthetic weights).
#
# Output (all numbers %.17g, LF):
#   phi2tail_corners.csv  Arm A: one row per (h, k, rho) corner, h <= k
#   phi2tail_cells.csv    Arm B: one row per cell of the cell design, with the
#                         reference P, u, v, w and tau_U for C0 and C1
#   phi2tail_design.csv   Arm C: the cells of the design datasets, weights,
#                         fitted thresholds, slopes and rho
#   phi2tail_blocks.csv   Arm C: reference U_rho and U_theta per dataset x pair
#   phi2tail_uni.csv      univariate add-on: P, phi(l)/P, phi(u)/P
#
# Reference method: that of tools/make_binormal_ref.R (M1, tanh-sinh on the
# asin(rho) form, checked against M2 and against the half step), at 256 bits.
# The cell quantities evaluate @sec-condx-hess's closed forms in MPFR from
# MPFR corners: the formulas are shared on purpose, the arithmetic is not.
# Every reference is a double-double pair ref_hi + ref_lo.
#
# tau_U per cell comes from cxh_cell_err (the derivation file), in double. It
# is affine in delta_Phi, so two columns are written per quantity: tau at
# delta_Phi = 0 (t0) and the increment per 2^-52 (t1). The do-file uses
# tau = t0 + (delta_Phi / 2^-52) * t1 with the delta_Phi that Arm A certifies.
# Arm B tests against tau_U2 (sec_cols below), formed from those columns.

suppressMessages(library(Rmpfr))
PREC <- 256
smoke <- Sys.getenv("PHI2T_SMOKE") == "1"
out_dir <- Sys.getenv("PHI2T_OUT", ".")

# ---- the oracle's chunks: closed forms, tau_U, Stage 1 ---------------------
qmd <- readLines("spec/derivations/dwlsmv-derivations.qmd")
st <- grep("^```[{]r ", qmd); lab <- sub("^```[{]r ([^,}]+).*", "\\1", qmd[st])
for (w in c("cx-setup", "cx-stage1", "cxh-forms", "cxh-tauu", "cxh-tauu-2nd-def")) {
  a <- st[lab == w]; b <- a + which(qmd[(a + 1):length(qmd)] == "```")[1]
  eval(parse(text = qmd[(a + 1):(b - 1)]), envir = globalenv())
}

# ---- MPFR primitives (as tools/make_binormal_ref.R) ------------------------
two_pi <- 2 * Const("pi", PREC); one <- mpfr(1, PREC); zero <- mpfr(0, PREC)
de_rule <- function(step, tmax = 4.5) {
  t <- mpfr(seq(-tmax, tmax, by = step), PREC)
  u <- Const("pi", PREC) / 2 * sinh(t)
  list(x = tanh(u), w = step * Const("pi", PREC) / 2 * cosh(t) / cosh(u)^2)
}
DE <- de_rule(1/64); DE_CHECK <- de_rule(1/32)
g_int <- function(h, k, lo, hi, rule = DE) {
  if (lo == hi) return(zero)
  h <- mpfr(h, PREC); k <- mpfr(k, PREC)
  t <- (lo + hi) / 2 + (hi - lo) / 2 * rule$x
  g <- exp(-(h * h - 2 * h * k * sin(t) + k * k) / (2 * cos(t)^2))
  g[!is.finite(g)] <- 0
  (hi - lo) / 2 * sum(rule$w * g)
}
Pn <- function(x) pnorm(mpfr(x, PREC))
phi2_m1 <- function(h, k, rho, rule = DE) {
  Pn(h) * Pn(k) + g_int(h, k, zero, asin(mpfr(rho, PREC)), rule) / two_pi
}
phi2_m2 <- function(h, k, rho) {
  if (rho >= 0) Pn(min(h, k)) - g_int(h, k, asin(mpfr(rho, PREC)), Const("pi", PREC) / 2) / two_pi
  else Pn(h) - phi2_m2(h, -k, -rho)
}
# Agreement required of M1, M2 and the half step: absolute 1e-30, or
# relative 1e-20 (C1's small corners need relative accuracy).
ref_tol <- function(x) max(1e-30, 1e-20 * abs(asNumeric(x)))
cache <- new.env()
n_check_fail <- 0
checks <- new.env()
PHI2 <- function(h, k, rho) {           # corner, +-Inf handled, cached, checked
  if (h == -Inf || k == -Inf) return(zero)
  if (h == Inf && k == Inf) return(one)
  if (h == Inf) return(Pn(k)); if (k == Inf) return(Pn(h))
  if (h > k) { t <- h; h <- k; k <- t }
  key <- sprintf("%.17g|%.17g|%.17g", h, k, rho)
  if (!is.null(cache[[key]])) return(cache[[key]])
  m1 <- phi2_m1(h, k, rho); tol <- ref_tol(m1)
  d_m2 <- abs(asNumeric(m1 - phi2_m2(h, k, rho)))
  d_half <- abs(asNumeric(m1 - phi2_m1(h, k, rho, DE_CHECK)))
  if (d_m2 > tol || d_half > tol) {
    n_check_fail <<- n_check_fail + 1
  }
  # every self-check is written out (phi2tail_refcheck.csv), so that a failing
  # corner can be identified rather than only counted
  checks[[key]] <- c(h, k, rho, asNumeric(m1), d_m2, d_half, tol)
  if (m1 < 0) m1 <- zero                # below the reference's own floor
  cache[[key]] <- m1
  m1
}
hilo <- function(x) { hi <- asNumeric(x); c(hi, asNumeric(x - hi)) }

# ---- grids and cell design (shared) ----------------------------------------
source("tools/phi2_tail_grid.R")

# ---- Arm A: corners --------------------------------------------------------
crows <- list()
for (i in seq_along(G)) for (j in i:length(G)) for (rho in rho_grid) {
  h <- G[i]; k <- G[j]; m1 <- PHI2(h, k, rho); hl <- hilo(m1)
  crows[[length(crows) + 1]] <- c(h, k, rho, in_cert(c(h, k), rho), hl)
}
corners <- as.data.frame(do.call(rbind, crows))
names(corners) <- c("h", "k", "rho", "cert", "ref_hi", "ref_lo")

# ---- cell quantities in MPFR -----------------------------------------------
# Closed forms of @sec-condx-hess, on bounds a = (a1, a2), b = (b1, b2) and
# gradient lists A, B (exact small integers). Returns P, u, v, w (mpfr).
dnm <- function(x) exp(-mpfr(x, PREC)^2 / 2) / sqrt(two_pi)
ref_cell <- function(a, b, rho, A, B) {
  r <- mpfr(rho, PREC); s2 <- 1 - r * r; s <- sqrt(s2)
  ni <- length(A[[1]]); nj <- length(B[[1]])
  P <- Pr <- Prr <- zero
  dPi <- dPri <- rep(zero, ni); dPj <- dPrj <- rep(zero, nj)
  for (k in 1:2) for (l in 1:2) {
    sg <- (-1)^(k + l); ak <- a[k]; bl <- b[l]
    P <- P + sg * PHI2(ak, bl, rho)
    if (is.finite(ak) && is.finite(bl)) {
      am <- mpfr(ak, PREC); bm <- mpfr(bl, PREC)
      Q <- (am^2 - 2 * r * am * bm + bm^2) / s2
      f <- exp(-Q / 2) / (two_pi * s)
      Pr <- Pr + sg * f
      Prr <- Prr + sg * f * (r * (1 - Q) + am * bm) / s2
      dPri <- dPri - sg * f * (am - r * bm) / s2 * A[[k]]
      dPrj <- dPrj - sg * f * (bm - r * am) / s2 * B[[l]]
    }
    cond <- function(x, y) {           # Phi((y - rho x)/s), x finite
      if (y == Inf) return(one); if (y == -Inf) return(zero)
      pnorm((mpfr(y, PREC) - r * mpfr(x, PREC)) / s)
    }
    if (is.finite(ak)) dPi <- dPi + sg * dnm(ak) * cond(ak, bl) * A[[k]]
    if (is.finite(bl)) dPj <- dPj + sg * dnm(bl) * cond(bl, ak) * B[[l]]
  }
  u <- Pr / P
  list(P = P, u = u, v = Prr / P - u^2,
       w = c(dPri, dPrj) / P - u * c(dPi, dPj) / P)
}
# C1's reflection: an axis whose cell lies mostly above 0 (a1 + a2 > 0) is
# mapped to x -> -x, with rho -> -rho, bounds reversed and negated, and the
# gradients negated. u and w then carry the product of the signs; v does not.
reflect <- function(a, A) {
  if (a[1] + a[2] > 0) list(a = -rev(a), A = lapply(rev(A), `-`), s = -1)
  else list(a = a, A = A, s = 1)
}
tau_cols <- function(a, b, rho, A, B) {   # tau_U at delta_Phi = 0 and per 2^-52
  # in the cell's units of 2^k (cxh_cell_err's exponent factoring; k is
  # written by sec_cols), so that no column underflows
  t0 <- cxh_cell_err(a, b, rho, A, B, dPhi = 0)$hat$tau
  t1 <- cxh_cell_err(a, b, rho, A, B, dPhi = 2^-52)$hat$tau
  pad <- function(x) c(x, rep(0, 8 - length(x)))
  out <- c(t0$u, t1$u - t0$u, t0$v, t1$v - t0$v, pad(t0$w), pad(t1$w - t0$w))
  # -1: not defined, the double P underflows (far below delta_P, outside Arm B)
  out[!is.finite(out)] <- -1
  out
}
# tau_U2 (@sec-condx-hess-2nd): tau_U plus the exact remainder. It is not
# affine in delta_Phi, so its pieces are written instead, each affine (t0, t1)
# or free of delta_Phi (checked here), all divided by the double |P|:
# eta-bar; the primitive parts pu, pv, pw_k; for v: x = X_rr/P,
# y = (dP_rr + X_rr)/P, z = (dP_r/P)^2, a = |v|; for w_k the same with z_k
# = (dP_r/P)(d dP_k/P), affine. All but eta-bar are in the cell's units of
# 2^k (k is the last column; vz carries 2^k once, as cxh_second does). The
# do-file forms tau_U2 / 2^k, and scales the measured error the same way,
# with eb = eta-bar,
#   u: tau_u + pu eb/(1 - eb)
#   v: tau_v + (pv + x + eb y + z + eb^2 a)/(1 - eb)^2 - pv, and w_k alike.
sec_cols <- function(a, b, rho, A, B) {
  z0 <- cxh_cell_err(a, b, rho, A, B, dPhi = 0)
  z1 <- cxh_cell_err(a, b, rho, A, B, dPhi = 2^-52)
  h0 <- z0$hat; h1 <- z1$hat; g0 <- h0$agg; g1 <- h1$agg; P <- abs(z0$P)
  same <- function(p, q) isTRUE(all.equal(p, q, tolerance = 0))
  stopifnot(same(h0$k, h1$k), same(g0$xrr, g1$xrr), same(g0$ePrr, g1$ePrr), same(g0$ePr, g1$ePr),
            same(g0$xdr, g1$xdr), same(g0$edPr, g1$edPr), same(z0$P, z1$P))
  pad <- function(x) c(x, rep(0, 8 - length(x)))
  out <- c(g0$eP / P, (g1$eP - g0$eP) / P,
           h0$prim$u, h1$prim$u - h0$prim$u, h0$prim$v, h1$prim$v - h0$prim$v,
           g0$xrr / P, (g0$ePrr + g0$xrr) / P, h0$sk * (g0$ePr / P)^2, abs(h0$v),
           pad(h0$prim$w), pad(h1$prim$w - h0$prim$w), pad(g0$xdr / P),
           pad((g0$edPr + g0$xdr) / P), pad(g0$ePr * g0$edP / P^2),
           pad(g0$ePr * (g1$edP - g0$edP) / P^2), pad(abs(h0$w)))
  out[!is.finite(out)] <- -1   # as tau_cols: the double P underflows
  c(out, h0$k)                 # k last: every piece but eb is in units of 2^k
}
sn <- c("eb_t0", "eb_t1", "pu_t0", "pu_t1", "pv_t0", "pv_t1", "vx", "vy", "vz", "va",
        paste0("pw", 1:8, "_t0"), paste0("pw", 1:8, "_t1"), paste0("wx", 1:8),
        paste0("wy", 1:8), paste0("wz", 1:8, "_t0"), paste0("wz", 1:8, "_t1"),
        paste0("wa", 1:8), "k")
cell_bounds <- function(th, cat, x) {
  K <- length(th); tt <- c(-Inf, th[-K], Inf)
  list(a = tt[cat + 1:2] - th[K] * x,
       A = list(cxh_grad(cat, K, x), cxh_grad(cat + 1, K, x)))
}

# ---- Arm D (prereg §3, split domain; user 2026-10-28) ---------------------
# For each certification cell below Arm B's domain (P < 1e3 delta_P at the
# budget), C1's corners are computed by M4 with a relative-only tolerance and
# lower limit -60 (tools/phi2_tail_refcheck_m4.R's M4, same code), and by M1.
#   D_P     C1's reference P: the M4 corners differenced in C1's form
#           (kappa <= ~1, no cancellation); +-Inf corners exact.
#   D_dom   3 floor (envelope): D_P < 2^-1074 / 1e-8, or kappa > 1e-8 / 2^-53;
#           1 certified: every finite C1 corner has |M1 - M4| <= 1e-20 |M4|;
#           2 diagnostic: otherwise (M1 has no relative accuracy there);
#           0 not in Arm D.
# The split is fixed by this rule and the reference alone, before any Stata
# output; its counts are printed (REF|D lines).
m4src <- readLines("tools/phi2_tail_refcheck_m4.R")
eval(parse(text = m4src[grep("^# ---- M4", m4src):(grep("^# ---- scope", m4src) - 1)]))
D_dP <- 4 * 4.23 * 2^-52 + 3 * 2^-52
m4cache <- new.env()
# Each C1 corner gets a status, decided in MPFR (never after conversion to
# double, where a corner below ~1e-308 and its discrepancy both become 0):
#   exact     a bound is -Inf (0), or both are +Inf (1)
#   margin    one bound +Inf: Phi of the other, by MPFR pnorm alone; not
#             cross-checked by M1-M4 (a shared-dependency limitation, reported)
#   agree     M4 finite and > 0, M1 finite, |M1 - M4| <= 1e-20 M4
#   disagree  M4 finite and > 0, but the test fails (or M1 is NaN or clipped)
#   invalid   M4 zero, negative, NaN or infinite: no relative test is defined
c1_corner <- function(h, k, r1) {     # list(m4, status)
  if (h == -Inf || k == -Inf) return(list(m4 = zero, status = "exact"))
  if (h == Inf && k == Inf) return(list(m4 = one, status = "exact"))
  if (h == Inf) return(list(m4 = Pn(k), status = "margin"))
  if (k == Inf) return(list(m4 = Pn(h), status = "margin"))
  key <- sprintf("%.17g|%.17g|%.17g", min(h, k), max(h, k), r1)
  if (is.null(m4cache[[key]])) {
    m4 <- phi2_m4(min(h, k), max(h, k), r1, abs_tol = 0, rel_tol = 1e-28, lo = -60)
    m1 <- PHI2(h, k, r1)
    st <- if (is.na(m4) || !is.finite(m4) || !(m4 > 0)) "invalid"
          else if (!is.na(m1) && is.finite(m1) && isTRUE(abs(m1 - m4) <= D_rtol * m4)) "agree"
          else "disagree"
    m4cache[[key]] <- list(m4 = m4, status = st)
  }
  m4cache[[key]]
}
D_rtol <- mpfr(1e-20, PREC)
D_Pmin <- mpfr(2, PREC)^-1074 / mpfr(1e-8, PREC)   # underflow floor, in MPFR
D_kmax <- mpfr(1e-8, PREC) / mpfr(2, PREC)^-53     # conditioning floor
D_cov <- list()                                     # per Arm D cell, for the report
armD_cols <- function(a, b, rho, Pref) {
  if (!(in_cert(c(a, b), rho) && asNumeric(Pref) < 1e3 * D_dP)) return(c(0, 0, 0))
  ra <- reflect(a, list(0, 0)); rb <- reflect(b, list(0, 0)); r1 <- ra$s * rb$s * rho
  P <- zero; S <- zero; st <- character(0); cmin <- Inf
  for (k in 1:2) for (l in 1:2) {
    cc <- c1_corner(ra$a[k], rb$a[l], r1)
    P <- P + (-1)^(k + l) * cc$m4; S <- S + abs(cc$m4); st <- c(st, cc$status)
    if (cc$status %in% c("agree", "disagree")) cmin <- min(cmin, asNumeric(log10(cc$m4)))
  }
  dom <- if (!(P > 0) || P < D_Pmin || S / P > D_kmax) 3
         else if (all(st %in% c("exact", "margin", "agree"))) 1 else 2
  D_cov[[length(D_cov) + 1]] <<- list(dom = dom, rho = rho, log10P = if (P > 0) asNumeric(log10(P)) else -Inf,
                                      cmin = cmin, margin = any(st == "margin"))
  c(hilo(P), dom)
}

brows <- list()
for (ti in seq_along(tmpl)) for (tj in ti:length(tmpl)) for (x in 0:1) {
  bs <- slope_pairs(x)
  for (bb in bs) for (rho in rho_grid) {
    thi <- c(tmpl[[ti]], bb[1]); thj <- c(tmpl[[tj]], bb[2])
    for (ci in 0:length(tmpl[[ti]])) for (cj in 0:length(tmpl[[tj]])) {
      bi <- cell_bounds(thi, ci, x); bj <- cell_bounds(thj, cj, x)
      R <- ref_cell(bi$a, bj$a, rho, bi$A, bj$A)
      ri <- reflect(bi$a, bi$A); rj <- reflect(bj$a, bj$A)
      tau0 <- tau_cols(bi$a, bj$a, rho, bi$A, bj$A)
      tau1 <- tau_cols(ri$a, rj$a, ri$s * rj$s * rho, ri$A, rj$A)
      sec0 <- sec_cols(bi$a, bj$a, rho, bi$A, bj$A)
      sec1 <- sec_cols(ri$a, rj$a, ri$s * rj$s * rho, ri$A, rj$A)
      dcol <- armD_cols(bi$a, bj$a, rho, R$P)
      w <- c(R$w, rep(zero, 8 - length(R$w)))
      vals <- c(R$P, R$u, R$v, w)
      brows[[length(brows) + 1]] <- c(ti, tj, x, bb, rho, ci, cj,
        in_cert(c(bi$a, bj$a), rho),
        unlist(lapply(seq_along(vals), function(z) hilo(vals[z]))), tau0, tau1, sec0, sec1, dcol)
    }
  }
}
cells <- as.data.frame(do.call(rbind, brows))
qn <- c("P", "u", "v", paste0("w", 1:8))
tn <- c("u_t0", "u_t1", "v_t0", "v_t1", paste0("w", 1:8, "_t0"), paste0("w", 1:8, "_t1"))
names(cells) <- c("ti", "tj", "x", "bi", "bj", "rho", "ci", "cj", "cert",
                  as.vector(t(outer(qn, c("_hi", "_lo"), paste0))),
                  paste0("c0_", tn), paste0("c1_", tn),
                  paste0("c0_", sn), paste0("c1_", sn), "D_P_hi", "D_P_lo", "D_dom")
# Arm D's classes and the exact coverage of its certified class (prereg §3).
dcv <- do.call(rbind, lapply(D_cov, function(z) data.frame(z)))
cst <- table(factor(vapply(ls(m4cache), function(z) m4cache[[z]]$status, ""),
                    levels = c("agree", "disagree", "invalid")))
cat(sprintf("REF|D|tail_cells=%d|certified=%d|diagnostic=%d|floor=%d\n", nrow(dcv),
            sum(dcv$dom == 1), sum(dcv$dom == 2), sum(dcv$dom == 3)))
cat(sprintf("REF|D|C1_corners_M4=%d|agree=%d|disagree=%d|invalid=%d\n",
            length(ls(m4cache)), cst[["agree"]], cst[["disagree"]], cst[["invalid"]]))
ce <- dcv[dcv$dom == 1, ]
if (nrow(ce)) {
  cat(sprintf("REF|D|certified|rho=%s\n", paste(sort(unique(ce$rho)), collapse = " ")))
  cat(sprintf("REF|D|certified|log10P_range=[%.2f, %.2f]|log10_min_checked_corner=%.2f|cells_with_margin_corner=%d (pnorm only)\n",
              min(ce$log10P), max(ce$log10P), min(ce$cmin), sum(ce$margin)))
}
dg <- dcv[dcv$dom == 2, ]
if (nrow(dg)) cat(sprintf("REF|D|diagnostic|log10P_range=[%.2f, %.2f]\n", min(dg$log10P), max(dg$log10P)))
cat("REF|D|certified cells are listed in phi2tail_cells.csv as D_dom == 1\n")
tmpl_tab <- data.frame(t = seq_along(tmpl),
  tau = vapply(tmpl, function(v) paste(sprintf("%.17g", v), collapse = " "), ""))

# ---- Arm C: the design datasets --------------------------------------------
if (smoke) {
  set.seed(1)                           # synthetic weights, not a design
  pop <- cx_population(cx_true); pop$w <- prop.table(runif(length(pop$w)))
  designs <- list(synthetic = pop)
} else {
  datA <- read.csv("fixtures/data/condx_A.csv", colClasses = "character")
  datA[] <- lapply(datA, as.integer)
  keyA <- do.call(paste, datA); tabA <- table(keyA); first <- match(names(tabA), keyA)
  dA <- list(Y = as.matrix(datA[first, 1:3]), x = datA$x[first],
             w = as.numeric(tabA) / nrow(datA)); dimnames(dA$Y) <- NULL
  pop <- cx_population(cx_true)
  rpop <- function(NN) {
    d <- pop; d$w <- round(NN * pop$w) / NN; keep <- d$w > 0
    d$Y <- d$Y[keep, ]; d$x <- d$x[keep]; d$w <- d$w[keep] / sum(d$w[keep]); d
  }
  designs <- list(condx_A = dA, pop_N500 = rpop(500), pop_N10000 = rpop(1e4),
                  pop_N100000 = rpop(1e5))
}
drows <- list(); krows <- list()
for (dn in names(designs)) {
  d <- designs[[dn]]; fit <- cx_stage1(d)
  for (m in seq_len(nrow(cx_pairs))) {
    i <- cx_pairs[m, 1]; j <- cx_pairs[m, 2]
    thi <- fit$th[[i]]; thj <- fit$th[[j]]; rho <- fit$rho[m]
    cells_m <- cx_pair_cells(d$Y[, i], d$Y[, j], d$x, d$w)
    Ur <- zero; Ut <- rep(zero, length(thi) + length(thj))
    for (r in seq_along(cells_m$ci)) {
      bi <- cell_bounds(thi, cells_m$ci[r], cells_m$x[r])
      bj <- cell_bounds(thj, cells_m$cj[r], cells_m$x[r])
      R <- ref_cell(bi$a, bj$a, rho, bi$A, bj$A)
      wm <- mpfr(cells_m$w[r], PREC)
      Ur <- Ur + wm * R$v; Ut <- Ut + wm * R$w
      pad4 <- function(v) c(v, rep(NA, 4 - length(v)))
      drows[[length(drows) + 1]] <- c(dataset = match(dn, names(designs)), pair = m,
        ci = cells_m$ci[r], cj = cells_m$cj[r], x = cells_m$x[r], w = cells_m$w[r],
        Ki = length(thi), Kj = length(thj), thi = pad4(thi), thj = pad4(thj), rho = rho)
    }
    Ut <- c(Ut, rep(zero, 8 - length(Ut)))
    krows[[length(krows) + 1]] <- c(match(dn, names(designs)), m, hilo(Ur),
      unlist(lapply(seq_along(Ut), function(z) hilo(Ut[z]))))
  }
}
design <- as.data.frame(do.call(rbind, drows))
names(design) <- c("dataset", "pair", "ci", "cj", "x", "w", "Ki", "Kj",
                   paste0("thi", 1:4), paste0("thj", 1:4), "rho")
blocks <- as.data.frame(do.call(rbind, krows))
names(blocks) <- c("dataset", "pair", "Ur_hi", "Ur_lo",
                   as.vector(t(outer(paste0("Ut", 1:8), c("_hi", "_lo"), paste0))))
dnames <- data.frame(dataset = seq_along(designs), name = names(designs))

# ---- univariate add-on -----------------------------------------------------
Ginf <- c(-Inf, G, Inf); urows <- list()
for (p in seq_along(Ginf)) for (q in seq_along(Ginf)) if (p < q && !(p == 1 && q == length(Ginf))) {
  l <- Ginf[p]; u <- Ginf[q]
  P <- (if (u == Inf) one else Pn(u)) - (if (l == -Inf) zero else Pn(l))
  rl <- if (is.finite(l)) dnm(l) / P else zero
  ru <- if (is.finite(u)) dnm(u) / P else zero
  urows[[length(urows) + 1]] <- c(l, u, all(!is.finite(c(l, u)) | abs(c(l, u)) <= L),
                                  hilo(P), hilo(rl), hilo(ru))
}
uni <- as.data.frame(do.call(rbind, urows))
names(uni) <- c("l", "u", "cert", "P_hi", "P_lo", "rl_hi", "rl_lo", "ru_hi", "ru_lo")

# ---- write ---------------------------------------------------------------
fmt <- function(x) ifelse(is.na(x), "", ifelse(is.infinite(x), ifelse(x > 0, "inf", "-inf"),
                                               sprintf("%.17g", x)))
write_lf <- function(df, file) {
  for (j in names(df)) if (is.numeric(df[[j]])) df[[j]] <- fmt(df[[j]])
  con <- file(file.path(out_dir, file), "wb")
  writeLines(paste(names(df), collapse = ","), con, sep = "\n")
  writeLines(do.call(paste, c(df, sep = ",")), con, sep = "\n")
  close(con)
}
write_lf(corners, "phi2tail_corners.csv"); write_lf(cells, "phi2tail_cells.csv")
write_lf(tmpl_tab, "phi2tail_templates.csv"); write_lf(design, "phi2tail_design.csv")
write_lf(blocks, "phi2tail_blocks.csv"); write_lf(dnames, "phi2tail_datasets.csv")
write_lf(uni, "phi2tail_uni.csv")
rc <- as.data.frame(do.call(rbind, as.list(checks)))
names(rc) <- c("h", "k", "rho", "m1", "d_m2", "d_half", "tol")
write_lf(rc, "phi2tail_refcheck.csv")

cat(sprintf("REF|smoke=%s corners=%d cells=%d design_cells=%d blocks=%d uni=%d\n",
            smoke, nrow(corners), nrow(cells), nrow(design), nrow(blocks), nrow(uni)))
cat(sprintf("REF|reference self-check failures (M1/M2/half step) = %d\n", n_check_fail))
cat(sprintf("REF|R %s, Rmpfr %s, MPFR %s, PREC %d\n", getRversion(),
            packageVersion("Rmpfr"), mpfrVersion(), PREC))
