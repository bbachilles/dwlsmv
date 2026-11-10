*! dwlsmv -- the model-implied latent correlation matrix R_y* (Module 5)
*!
*! Executed under StataMP 18 on 2026-09-22. tests/unit/test_model.do check 4
*! asserts the unit diagonal at exactly 0.0 on four configurations, and the
*! standardization gap at 1e-15 -- a tolerance, not zero, and deliberately so:
*! Sigma[j,j] is formed as fl(a + fl(1 - a)), which is not guaranteed to be
*! bitwise 1.0, so D^-1/2 Sigma D^-1/2 can differ from Sigma in the last ulp
*! even though writing 1 down the diagonal does not. The measured gap is
*! ~1e-16. An earlier version of this line claimed both were asserted at 0.0.
*!
*! Spec Section 4.6, and Gate 2's Finding 1 in the derivation file:
*!
*!     Phi_eta = (I - B)^-1 Psi (I - B)^-T
*!     Sigma   = Lambda Phi_eta Lambda' + Theta
*!     R_y*    = D^-1/2 Sigma D^-1/2,   D = diag(Sigma)
*!
*! UNDER THE DELTA PARAMETERIZATION THOSE ARE NOT THREE STEPS. Delta fixes
*! Var(y_j*) = 1, which makes Theta DERIVED,
*! theta_jj = 1 - (Lambda Phi Lambda')_jj, so Sigma already has unit diagonal
*! and R_y* = Sigma identically. The standardization is a no-op -- verified to
*! machine zero in the derivation, and asserted again from Stata in
*! tests/unit/test_model.do.
*!
*! So this routine forms Lambda Phi Lambda' and writes 1 down the diagonal, and
*! that IS the standardization, performed exactly rather than to within a
*! division. Computing D^-1/2 Sigma D^-1/2 numerically instead would give the
*! same matrix to a few ulps and cost a square root per indicator, and the
*! ulps would then be indistinguishable from a real disagreement with lavaan.
*!
*! Theta is not returned, and there is nowhere in the API to put it. It is
*! derived, never free (spec Section 4.6): a free residual variance is
*! unidentified from correlations alone, since scaling lambda_j by c and
*! theta_jj by c^2 leaves R_y* untouched. test_model.do asserts exactly that
*! invariance, so the claim is measured here and not only in the derivation.
*!
*! (I - B) SINGULARITY IS CHECKED, not assumed. For a recursive model B is
*! strictly triangular under some ordering and the inverse always exists, but
*! the canonical form of Section 4.8 can express a non-recursive B -- 'beta F2
*! F1' and 'beta F1 F2' together -- and luinv() does not abort on a singular
*! argument, it returns a matrix of missing values. Unchecked, that propagates
*! into sigma(theta) as missing and surfaces as an optimizer failure several
*! modules downstream.
*!
*! WHAT THAT CHECK DOES NOT CATCH, so nobody reads it as more than it is: a
*! cyclic B whose (I - B) is NONSINGULAR -- 'beta F1 F2' and 'beta F2 F1' at
*! values like 0.3 and 0.4. Phi is then perfectly computable. (Whether B itself
*! is invertible is irrelevant; only I - B enters.) This used to be a gap, and
*! is now closed upstream: since 2026-09-23 dwls_identify() rejects any cycle
*! in B, because v1 is recursive-only (spec Section 4.8). So in v1 this
*! routine only ever sees an acyclic B, whose (I - B) is unit triangular under
*! some ordering and cannot be singular. The singularity check stays as the
*! backstop for anyone calling it on a hand-built table.
*!
*! NOR IS ADMISSIBILITY CHECKED, and the asymmetry with the R oracle is
*! deliberate. This routine writes 1 down the diagonal unconditionally, so a
*! specification implying (Lambda Phi Lambda')_jj > 1 yields a derived
*! Theta_jj < 0 and |R_jk| > 1 with no complaint, where the tier-2 generator's
*! own sigma_spec() stops with "stated matrices are inadmissible under delta".
*! The generator can be stricter because it is handed final values; this
*! routine is called at TRIAL values by Module 7's optimizer, which must be
*! free to step through an inadmissible region on its way somewhere better.
*! Spec Section 5 assigns the admissibility check to Module 7, before a fitted
*! theta-hat is accepted, and that is the right place for it.
*!
*! src/mata/dwlsmv_types.mata must be loaded first.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_rystar(real matrix Lam, real matrix Psi, real matrix B)
{
    real matrix A, Phi, R
    real scalar nfac, j

    nfac = rows(Psi)
    A    = luinv(I(nfac) - B)

    if (missing(A)) {
        errprintf("dwlsmv_rystar: (I - B) is singular, so Phi = (I-B)^-1 Psi (I-B)^-T does not exist.\n")
        errprintf("Check the beta directives for a feedback loop among the factors.\n")
        _error("dwlsmv_rystar: (I - B) is singular")
    }

    Phi = A * Psi * A'
    R   = Lam * Phi * Lam'

    // Theta = diag(1 - diag(Lambda Phi Lambda') is what makes this the
    // correlation matrix; adding it touches the diagonal and nothing else.
    for (j = 1; j <= rows(R); j++) {
        R[j, j] = 1
    }

    return(R)
}

end
