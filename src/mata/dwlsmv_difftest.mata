*! dwlsmv -- Module 12: the DIFFTEST statistic for one nested pair, any G
*!
*! Executed under StataMP 18 on 2026-10-08: tests/unit/test_difftest.do passes.
*!
*! Spec Module 12, *The statistic, frozen* and *Architecture* item 1; derivation
*! @sec-difftest. Inputs (E-003, 2026-10-08): Delta1 (D x q1) at the H1
*! estimate, Delta0 (D x q0) at the H0 estimate, Gamma (D x D), and
*! Td = T0 - T1 with each T = (N-1)F. Returns the row vector
*!
*!     (T_d3, d, p, a_d, b_d, Tr(M1), Tr(M1^2))
*!
*! where, with W = W_D = diag(Gamma)^-1 as a vector,
*!
*!     P1   = Delta1' W Delta1
*!     V1   = dwlsmv_acov(Delta1, Gamma)      the sqrt(N)-scale ACOV, as is
*!     h    = pinv(Delta1) Delta0             h-hat, NOT the selection matrix
*!     M1   = (P1 - P1 h (h' P1 h)^-1 h' P1) V1
*!     d    = q1 - q0
*!     a_d, b_d = dwlsmv_scaleshift(Tr(M1), Tr(M1^2), d)
*!     T_d3 = a_d Td + b_d,   p = chi2tail(d, T_d3)
*!
*! NO SAMPLE SIZE. Gamma's scale cancels: Gamma -> c Gamma sends W_D -> W_D/c,
*! P1 -> P1/c and V1 -> c V1, so M1, a_d and b_d are unchanged (test_difftest.do
*! check 5, measured exactly). The one scale that matters is Td's: it must be
*! T0 - T1 with each T = (N-1)F (Section 4.9), as lavaan's.
*!
*! WHAT IT DOES NOT DO. It does not check nesting, the sign of Td, or the
*! epsilon band (spec Module 12 checks 7-10): that is validation, which holds
*! T0 and T1. It does not refuse d = 0 by reporting "equivalent"; an empty step
*! is the caller's to report and skip, and here it is an error.
*!
*! G > 1 (2026-10-14; MG7(3), mg-difftest-omega): the same routine, unchanged.
*! The caller passes the stacked Delta1 and Delta0 and, in the Gamma slot,
*! Omega = blockdiag(Gamma_g / pi_g) from dwlsmv_omega(); diag(Omega)^-1 is then
*! weighting (a) and V1 the pooled ACOV (MG6). Td must be T0 - T1 with each
*! T = (N - G)F. h-hat stays the unweighted pinv on stacked rows (closed
*! 2026-10-09, E-019). The caller that applies the
*! status rule is dwlsmv_dtpair.
*!
*! h-hat uses pinv(), the Moore-Penrose inverse, as lavaan's MASS::ginv(). For
*! a full-column-rank Delta1 it is (Delta1' Delta1)^-1 Delta1'. h' P1 h is
*! inverted by cholinv(); a missing result (H0's columns do not span a q0-
*! dimensional subspace of H1's tangent space to first order) aborts. As in
*! dwlsmv_acov(), cholinv() catches only a floating-point pivot <= 0.
*!
*! ARITHMETIC. Every product is a quadcross() (the project's numerical rule). M1 is not
*! symmetric; Tr(M1^2) = sum(M1 :* M1').
*!
*! Requires dwlsmv_acov, dwlsmv_asym, dwlsmv_scaleshift. No struct is touched.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real rowvector dwlsmv_difftest(real matrix Delta1, real matrix Delta0,
                               real matrix Gamma, real scalar Td)
{
    real colvector w
    real matrix V1, P1, h, Ph, K, U, M
    real scalar d, trM, trM2, T
    real rowvector ab

    if (rows(Delta1) != rows(Gamma) | rows(Delta0) != rows(Gamma)) {
        errprintf("dwlsmv_difftest: Delta1 is %f x %f, Delta0 %f x %f, Gamma %f x %f; all need D rows\n",
                  rows(Delta1), cols(Delta1), rows(Delta0), cols(Delta0), rows(Gamma), cols(Gamma))
        _error("dwlsmv_difftest: Delta1, Delta0 and Gamma are not conformable")
    }
    d = cols(Delta1) - cols(Delta0)
    if (cols(Delta0) < 1 | d < 1) {
        errprintf("dwlsmv_difftest: H1 has %f and H0 %f free parameter(s); need 1 <= q0 < q1.\n",
                  cols(Delta1), cols(Delta0))
        errprintf("An empty step (d = 0) is the caller's to report and skip.\n")
        _error("dwlsmv_difftest: need 1 <= q0 < q1")
    }
    if (missing(Delta0) | missing(Td)) {
        _error("dwlsmv_difftest: Delta0 or Td is missing")
    }

    V1 = dwlsmv_acov(Delta1, Gamma)

    w  = 1 :/ diagonal(Gamma)
    P1 = quadcross(Delta1, w, Delta1)
    h  = quadcross(pinv(Delta1)', Delta0)
    Ph = quadcross(P1, h)
    K  = cholinv(quadcross(h, Ph))
    if (missing(K)) {
        errprintf("dwlsmv_difftest: h' P1 h is not positive definite, with h = pinv(Delta1) Delta0.\n")
        errprintf("Delta0's columns do not map to %f independent directions in H1's parameters.\n",
                  cols(Delta0))
        _error("dwlsmv_difftest: h' P1 h is not positive definite")
    }
    U = P1 - quadcross(Ph', quadcross(K, Ph'))
    M = quadcross(U, V1)

    trM  = sum(diagonal(M))
    trM2 = sum(M :* M')
    ab   = dwlsmv_scaleshift(trM, trM2, d)
    T    = ab[1] * Td + ab[2]
    return((T, d, chi2tail(d, T), ab[1], ab[2], trM, trM2))
}

end
