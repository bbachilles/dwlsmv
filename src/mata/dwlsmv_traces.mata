*! dwlsmv -- Module 9: Tr(M) and Tr(M^2), M = U Gamma, without forming U or M
*!
*! Executed under StataMP 18 on 2026-09-24; passed on the first run. tests/unit/test_inference.do: fed lavaan's own Delta and Gamma it reproduces lavaan's stored traces within 4.0e-2 of a 1e-12 relative budget on all five fits (check 1), and U formed explicitly within 1.4e-1 of 1e-13 (check 2).
*!
*! Spec Section 4.9; derivation @sec-gate5-traces, verified there against the
*! U-forming definition (<= 7.2e-16) and against lavaan 0.7.2's stored traces
*! on all five tier-2 fits (<= 5.7e-15). With w = diag(W_D) = 1 :/ diag(Gamma),
*! E = Delta' W Delta (the sandwich bread), G = Gamma W Delta,
*! B = (W Delta)' G and C = G' W G:
*!
*!     Tr(M)   = w' diag(Gamma) - Tr(E^-1 B)
*!     Tr(M^2) = w' (Gamma o Gamma) w - 2 Tr(E^-1 C) + Tr((E^-1 B)^2)
*!
*! U AND M ARE NEVER FORMED (decided 2026-09-24). Forming them is O(p*^3); the
*! expansion's largest step is G, O(p*^2 q). U stays the definition and the
*! unit-test oracle: tests/unit/test_inference.do forms it on small cases.
*!
*! w' diag(Gamma) is p* in exact arithmetic. It is computed as written rather
*! than replaced by rows(Gamma), so the routine is the boxed formula and
*! nothing else; the two differ by rounding only.
*!
*! NOT SHARED WITH dwlsmv_acov(), deliberately (decided 2026-09-24). That
*! routine forms E^-1 and Gamma W Delta E^-1, never G or B, so sharing would
*! mean changing its arithmetic or re-deriving these traces in terms of its
*! intermediates. What sharing would save was measured once, in a scratch
*! do-file not committed: ~80 ms at D = 940, q = 206 (bread 14 ms, G 65 ms),
*! ~1.3% of Stage 1 on the same problem.
*!
*! A MATRIX PRIMITIVE, like dwlsmv_acov(): it takes (Delta, Gamma) rather than
*! a fit, so the same code can be fed lavaan's own Delta and Gamma -- the tight
*! test. The guards repeat dwlsmv_acov()'s for that reason: a standalone caller
*! must not reach cholinv() with inputs dwls_inference() would have refused.
*! The bread is inverted by cholinv() for the reason dwlsmv_acov()'s header
*! records, and protects against the same, narrow, set of failures.
*!
*! ARITHMETIC. Every product and every trace is a quadcross(), for the project's
*! processor-count reason. Tr(X'Y) is taken as sum(diagonal(quadcross(X, Y))),
*! which costs a q x q x q product in quad precision and keeps the reduction
*! order fixed. E^-1 is symmetric, so quadcross(Ei, B) is E^-1 B.
*!
*! GAMMA MUST BE SYMMETRIC, and since 2026-09-25 that is ENFORCED: G =
*! quadcross(Gamma, WD) is Gamma' W Delta, so an asymmetric Gamma would give
*! silently different traces. The routine aborts with rc 3498 when
*! dwlsmv_asym() puts Gamma outside its budget. It carries no tolerance of its
*! own, so it cannot disagree with dwlsmv_acov() or with Module 10's diagnostic
*! (tests/unit/test_asym.do check 4). Stage 1's Gamma reads exactly 0;
*! lavaan's fixtures read 5.8e-16 and 6.8e-16 against 1.4e-14 and 2.4e-14.
*! Requires dwlsmv_asym.
*!
*! The missing-value guard is LOAD-BEARING, not defensive: quadcross() drops
*! rows containing missing values and sum() counts missing as 0, so without it
*! a missing Delta entry returns finite traces (review mutation 20, caught by
*! test_inference.do check 8).
*!
*! Returns (Tr(M), Tr(M^2)) as a 1 x 2 rowvector.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real rowvector dwlsmv_traces(real matrix Delta, real matrix Gamma)
{
    real colvector g, w
    real rowvector as
    real matrix WD, Ei, G, B, C, EB
    real scalar trM, trM2

    if (rows(Gamma) != cols(Gamma) | rows(Delta) != rows(Gamma)) {
        errprintf("dwlsmv_traces: Delta is %f x %f and Gamma %f x %f; Gamma must be D x D and Delta D x npar\n",
                  rows(Delta), cols(Delta), rows(Gamma), cols(Gamma))
        _error("dwlsmv_traces: Delta and Gamma are not conformable")
    }
    if (cols(Delta) < 1 | cols(Delta) > rows(Delta)) {
        errprintf("dwlsmv_traces: Delta has %f column(s) and %f row(s); need 1 <= npar <= D\n",
                  cols(Delta), rows(Delta))
        _error("dwlsmv_traces: Delta has no columns or more columns than rows")
    }
    if (missing(Delta) | missing(Gamma)) {
        _error("dwlsmv_traces: Delta or Gamma contains missing values")
    }

    g = diagonal(Gamma)
    if (min(g) <= 0) {
        errprintf("dwlsmv_traces: diag(Gamma) has a zero or negative entry, so W_D = diag(Gamma)^-1 does not exist.\n")
        _error("dwlsmv_traces: diag(Gamma) is not strictly positive")
    }
    as = dwlsmv_asym(Gamma)
    if (as[1] > as[2]) {
        errprintf("dwlsmv_traces: Gamma is not symmetric: max |G_ij - G_ji| / sqrt(G_ii G_jj) = %9.2e,\n", as[1])
        errprintf("  against a budget of %9.2e (4 D eps). This routine uses Gamma' = Gamma, so it\n", as[2])
        errprintf("  would return traces of a different matrix. See derivation @sec-gamma-sym.\n")
        _error("dwlsmv_traces: Gamma is not symmetric")
    }
    w  = 1 :/ g
    WD = w :* Delta

    Ei = cholinv(quadcross(Delta, w, Delta))
    if (missing(Ei)) {
        errprintf("dwlsmv_traces: the bread Delta' W_D Delta is not positive definite, so Delta is\n")
        errprintf("rank deficient: some direction in the parameters leaves sigma(theta)\n")
        errprintf("unchanged to first order. Module 3 checks this at the starting values only.\n")
        _error("dwlsmv_traces: Delta' W_D Delta is not positive definite")
    }

    G  = quadcross(Gamma, WD)
    B  = quadcross(WD, G)
    C  = quadcross(G, w, G)
    EB = quadcross(Ei, B)

    trM  = quadcross(w, g) - sum(diagonal(EB))
    trM2 = quadcross(w, quadcross(Gamma :* Gamma, w)) -
           2 * sum(diagonal(quadcross(Ei, C))) +
           sum(diagonal(quadcross(EB', EB)))

    return((trM, trM2))
}

end
