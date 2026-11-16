# ---------------------------------------------------------------------------
# Method C for the bivariate-normal corner Phi2(h, k; rho): Plackett's identity
# anchored at rho = -1 (docs/notes/phi2-flush-fixtures/design.md, P3).
#
#   Phi2(h, k; rho) = Phi2(h, k; -1) + integral_{-1}^{rho} phi2(h, k; t) dt,
#   Phi2(h, k; -1)  = max(0, Phi(h) + Phi(k) - 1).
#
# For h + k < 0 the anchor is 0, so the value is an integral of a positive
# integrand with no cancellation, and it keeps RELATIVE accuracy where the
# value is tiny (to 1e-68 and below). Restricted to h + k < 0; it stops
# otherwise.
#
# Independence from M4 (tools/phi2_tail_refcheck_m4.R): the representation
# differs (an integral of the density over rho, against M4's integral over x
# of phi(x) Phi((k - rho x)/s)); the algorithm differs (tanh-sinh on a fixed
# step halved to convergence, against adaptive 20-point Gauss-Legendre); the
# code is separate. It shares MPFR and Rmpfr's exp/sqrt/pi, but NOT pnorm:
# Method C calls no normal CDF at all. It SHARES its representation family
# with M1 (the arcsine form is Plackett's identity under t = sin(theta),
# anchored at 0), so C and M1 are not two methods.
#
# Defines phi2_c(h, k, rho) -> list(v = mpfr value, step, conv). Load by
# section, from "# ---- Method C" to "# ---- end Method C".
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek
# ---------------------------------------------------------------------------

# ---- Method C ---------------------------------------------------------------
phi2_c <- function(h, k, rho, prec = PREC, tol = 1e-40) {
  stopifnot(is.finite(h), is.finite(k), h + k < 0, rho > -1, rho < 1)
  H <- mpfr(h, prec); K <- mpfr(k, prec); R <- mpfr(rho, prec)
  pi_ <- Const("pi", prec); hw <- (R + 1) / 2
  f <- function(t) {
    q <- 1 - t * t
    exp(-(H * H - 2 * t * H * K + K * K) / (2 * q)) / (2 * pi_ * sqrt(q))
  }
  level <- function(step) {
    x <- mpfr(seq(-ceiling(4.5 / step), ceiling(4.5 / step)) * step, prec)
    u <- tanh(pi_ / 2 * sinh(x)); w <- pi_ / 2 * cosh(x) / cosh(pi_ / 2 * sinh(x))^2
    keep <- as.logical(abs(u) < 1)              # in MPFR: a double test drops nodes near rho
    t <- -1 + hw * (u[keep] + 1); w <- w[keep]
    ok <- as.logical((1 - t * t) > 0)
    hw * step * sum(w[ok] * f(t[ok]))
  }
  step <- 1 / 8; prev <- level(step)
  repeat {
    step <- step / 2; cur <- level(step)
    conv <- asNumeric(abs((cur - prev) / cur))
    if (conv < tol) break
    if (step < 1 / 1024) stop(sprintf("phi2_c: no convergence at h=%.17g k=%.17g rho=%.17g (%.2e)",
                                      h, k, rho, conv))
    prev <- cur
  }
  list(v = cur, step = step, conv = conv)
}
# ---- end Method C -----------------------------------------------------------
