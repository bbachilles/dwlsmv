*! dwlsmv -- Module 8: the robust sandwich ACOV(theta-hat)
*!
*! Executed under StataMP 18 on 2026-09-23. tests/unit/test_acov.do: fed lavaan's own Delta and Gamma it reproduces lavaan's vcov to <= 4.5e-15 on all five fits (check 2), and ACOV_tautau = Gamma_tautau at a 1e-14 per-entry budget, measured ~5e-16 (check 1). Adversarially reviewed 2026-09-24, with 13 mutations run against the test; the claims below were corrected from that review.
*!
*! Spec Section 4.10, derived in the derivation's @sec-acov:
*!
*!     ACOV = (Delta' W Delta)^-1 Delta' W Gamma W Delta (Delta' W Delta)^-1
*!
*! with W = W_D = diag(Gamma)^-1. ACOV is the estimated asymptotic covariance
*! of sqrt(N)(theta-hat - theta_0), built from Gamma EXACTLY AS STORED --
*! Section 4.4's N-scaled (1/N) sum psi psi', with no N/(N-1) factor. It is not
*! Var(theta-hat) and it is not e(V): dwlsmv_ev() is the one step that turns
*! it into e(V), and nothing here divides by any sample size. A stray factor
*! of N between Gamma, the objective and e(V) is spec Section 6's highest-rated
*! Module 4 risk, and keeping the two steps apart is what lets each be
*! unit-tested in isolation.
*!
*! AN INTERNAL PRIMITIVE, NOT THE LOCKED ENTRY POINT. dwls_inference() calls
*! it and fills acov together with Module 9's T3/shift_a/shift_b, so a
*! partially filled inference state never exists (E-038, 2026-09-23).
*! It shares nothing with dwlsmv_traces(), deliberately (E-032,
*! 2026-09-24): this routine's arithmetic stays exactly as tested.
*! This routine takes matrices rather than a fit so the same code can be fed
*! lavaan's own Delta and Gamma -- the tight test in test_acov.do.
*!
*! THE BREAD IS THE EXPECTED INFORMATION Delta' W Delta, not the Hessian of F.
*! lavaan's default is information = "expected"; "observed" would move case
*! D's SEs by up to 3.0e-3 relative (derivation @sec-acov, Level 3). Check 2
*! pins the information type at 1e-12 on lavaan's own inputs. Until
*! 2026-09-24 the end-to-end check 3 could NOT see it on case D (0.64 of the
*! old ~0.5%-relative SE tier); the current 1e-5 tier can, by two orders.
*! It is also the Gauss-Newton matrix
*! dwlsmv_fdwls() steps with, so the matrix inverted here is the one the
*! optimizer used.
*!
*! W IS A VECTOR, re-formed here as 1 :/ diagonal(Gamma) and never stored on
*! the fit (a field would be a spec change). It is W_D itself,
*! not its inverse -- the Muthen et al. (1997) convention trap in the
*! derivation's @sec-gate5.
*!
*! THE GENERAL FORMULA, deliberately, although the derivation shows Delta's
*! identity threshold block makes ACOV_tautau = Gamma_tautau exactly. That
*! identity is test_acov.do's oracle, and an implementation built on it could
*! not be tested by it.
*!
*! GAMMA IS NOT INVERTED and its rank is not checked: rank-deficient Gamma is
*! an open decision deferred to Module 10 (E-054, spec 4.5.6). Only its
*! diagonal must be strictly positive, which W needs. The bread IS inverted,
*! by cholinv(), and what that protects against is narrow:
*!
*!   cholinv() returns missing ONLY WHEN A FLOATING-POINT PIVOT IS <= 0.
*!   Singularity in exact arithmetic does not guarantee that. Measured
*!   2026-09-24: a bread [s s; s s], from two identical Delta columns, returns
*!   missing at s = 14 and s = 4 but a finite inverse with entries near 2.3e15
*!   at s = 2, whose second pivot rounds to +4.4e-16. A bread with an
*!   identically zero column -- a parameter sigma(theta) does not depend on
*!   at all -- returns missing under any rounding. (This header first said
*!   "exactly singular aborts" and then "numerically singular does not";
*!   the first half was also too strong. The 2026-09-24 review found it.)
*!
*!   So the abort below is reliable only for an identically zero Delta
*!   column. Every other rank deficiency -- identical or collinear columns,
*!   or near-deficiency from finite-difference noise (probe matrix,
*!   eigenvalue ratio 6.8e-17: entries near 1.3e14) -- can pass through as
*!   a finite inverse and absurd SEs. Against invsym() the gain is only that
*!   invsym() would silently zero the offending rows and report SEs of
*!   exactly 0. Catching the rest needs a conditioning criterion on
*!   Delta' W_D Delta, which is Module 10's open question (E-054,
*!   Rank-deficient Gamma) and is deliberately not improvised here. Module 3
*!   checks Delta's rank at the starting values only.
*!
*! ARITHMETIC. Every product is a quadcross(), for the project's numerical rule: cross()
*! and '*' move with the processor count. With At = W Delta Bi (D x npar),
*!
*!     ACOV = At' Gamma At = quadcross(At, quadcross(Gamma, At))
*!
*! which uses Gamma' = Gamma. The result is symmetrized as (V + V')/2: the two
*! triangles are separate accumulations and can differ in the last ulp, and a
*! covariance matrix should be exactly symmetric by construction rather than
*! leave that to whichever consumer checks first. That symmetrization would
*! also HIDE an asymmetric Gamma, which is why, since 2026-09-25, Gamma's
*! symmetry is an enforced input contract. The routine aborts with rc 3498 when
*! dwlsmv_asym() puts it outside its budget, with no tolerance of its own
*! (tests/unit/test_asym.do check 4). Requires dwlsmv_asym.
*!
*! src/mata/dwlsmv_types.mata need not be loaded; no struct is touched.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_acov(real matrix Delta, real matrix Gamma)
{
    real colvector g, w
    real rowvector as
    real matrix Bi, At, V

    if (rows(Gamma) != cols(Gamma) | rows(Delta) != rows(Gamma)) {
        errprintf("dwlsmv_acov: Delta is %f x %f and Gamma %f x %f; Gamma must be D x D and Delta D x npar\n",
                  rows(Delta), cols(Delta), rows(Gamma), cols(Gamma))
        _error("dwlsmv_acov: Delta and Gamma are not conformable")
    }
    if (cols(Delta) < 1 | cols(Delta) > rows(Delta)) {
        errprintf("dwlsmv_acov: Delta has %f column(s) and %f row(s); need 1 <= npar <= D\n",
                  cols(Delta), rows(Delta))
        _error("dwlsmv_acov: Delta has no columns or more columns than rows")
    }
    if (missing(Delta) | missing(Gamma)) {
        _error("dwlsmv_acov: Delta or Gamma contains missing values")
    }

    g = diagonal(Gamma)
    if (min(g) <= 0) {
        errprintf("dwlsmv_acov: diag(Gamma) has a zero or negative entry, so W_D = diag(Gamma)^-1 does not exist.\n")
        _error("dwlsmv_acov: diag(Gamma) is not strictly positive")
    }
    as = dwlsmv_asym(Gamma)
    if (as[1] > as[2]) {
        errprintf("dwlsmv_acov: Gamma is not symmetric: max |G_ij - G_ji| / sqrt(G_ii G_jj) = %9.2e,\n", as[1])
        errprintf("  against a budget of %9.2e (4 D eps). The sandwich uses Gamma' = Gamma and its\n", as[2])
        errprintf("  result is symmetrized, which would hide the asymmetry. See derivation @sec-gamma-sym.\n")
        _error("dwlsmv_acov: Gamma is not symmetric")
    }
    w = 1 :/ g

    Bi = cholinv(quadcross(Delta, w, Delta))
    if (missing(Bi)) {
        errprintf("dwlsmv_acov: the bread Delta' W_D Delta is not positive definite, so Delta is\n")
        errprintf("rank deficient at theta-hat: some direction in the parameters leaves sigma(theta)\n")
        errprintf("unchanged to first order. Module 3 checks this at the starting values only.\n")
        _error("dwlsmv_acov: Delta' W_D Delta is not positive definite")
    }

    At = quadcross((w :* Delta)', Bi)
    V  = quadcross(At, quadcross(Gamma, At))
    return((V + V') :/ 2)
}

end
