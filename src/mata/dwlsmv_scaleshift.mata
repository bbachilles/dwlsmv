*! dwlsmv -- Module 9: the scale a and shift b of T3 = a T + b
*!
*! Executed under StataMP 18 on 2026-09-24; passed on the first run. tests/unit/test_inference.do checks 1, 2b and 6, and nine rc 3498 error paths.
*!
*! Spec Section 4.9; derivation @sec-gate5 (Asparouhov & Muthen 2010, eq. 4):
*!
*!     a = sqrt(D / Tr(M^2)),    b = D - a Tr(M)
*!
*! chosen so E(T3) = D and Var(T3) = 2D, with df = D unchanged. b is written
*! D - a Tr(M), which is the spec's D - sqrt(D Tr(M)^2 / Tr(M^2)) because
*! Tr(M) > 0. It is a CANCELLATION when the correction is small -- case D has
*! D = 2 and b ~ 0.001 -- so it is asserted with an absolute tolerance.
*!
*! D = 0, A JUST-IDENTIFIED MODEL (decided 2026-09-24, spec 4.9). a would be
*! sqrt(0/0): M has no nonzero eigenvalues and both traces are rounding. No
*! scale or shift exists, and (., .) is returned rather than an error, so the
*! caller can still report standard errors. This is lavaan 0.7.2's behaviour,
*! measured once in R on a three-indicator model of case D: scaling.factor NA,
*! no shift parameter, the scaled statistic equal to the unscaled one (2.0e-13),
*! SEs reported. dwls_inference() then sets T3 = T_raw.
*!
*! D > 0 needs Tr(M) > 0 and Tr(M^2) > 0; anything else can only come from a
*! degenerate Gamma, which is Module 10's open question (E-054,
*! Rank-deficient Gamma), so it aborts rather than returning a number.
*!
*! Takes scalars, not matrices, so the D = 0 convention and the guards live in
*! one place and are testable without a fit.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real rowvector dwlsmv_scaleshift(real scalar trM, real scalar trM2, real scalar df)
{
    real scalar a

    if (missing(df) | df < 0 | df != floor(df)) {
        errprintf("dwlsmv_scaleshift: df = %f; the degrees of freedom must be an integer >= 0\n", df)
        _error("dwlsmv_scaleshift: df must be an integer >= 0")
    }
    if (df == 0) {
        return((., .))
    }
    if (missing(trM) | missing(trM2) | trM <= 0 | trM2 <= 0) {
        errprintf("dwlsmv_scaleshift: Tr(M) = %g and Tr(M^2) = %g with df = %f; both must be positive.\n",
                  trM, trM2, df)
        errprintf("M = U Gamma has df nonzero, non-negative eigenvalues, so this indicates a degenerate Gamma.\n")
        _error("dwlsmv_scaleshift: Tr(M) or Tr(M^2) is not positive")
    }

    a = sqrt(df / trM2)
    return((a, df - a * trM))
}

end
