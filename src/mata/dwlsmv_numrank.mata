*! dwlsmv -- the one numerical-rank criterion (Modules 3 and 10)
*!
*! Executed under StataMP 18 on 2026-09-26. tests/unit/test_breaddiag.do check 1:
*! the criterion edge from both sides (2e-8 kept, 5e-9 dropped), scale
*! invariance, V as rows; and every test that runs dwls_identify().
*!
*!     real scalar dwlsmv_numrank(real matrix M, real colvector sv,
*!                                real matrix V)
*!
*! Returns the numerical rank of M: the number of singular values greater than
*! 1e-8 times the largest. Fills sv with the singular values, descending, and V
*! with the right singular vectors AS ROWS -- svd()'s thin convention,
*! M = U diag(sv) V, pinned by measurement (tests/mata_allowlist.txt, svd) --
*! so a caller can name the null direction without a second decomposition.
*! An identically zero M has rank 0.
*!
*! THE ONE DEFINITION. dwls_identify() applies it to Delta at the starting
*! values and aborts on a deficiency (Module 3); dwlsmv_breaddiag() applies it
*! to W_D^1/2 Delta at theta-hat and reports (Module 10). Until 2026-09-26 the
*! constant lived inline in dwls_identify(); a second copy in the diagnostic is
*! how a guard and a report come to disagree.
*!
*! WHY 1e-8 (relative). Delta's structural block is a central difference at
*! h = 1e-5 max(1, |v|) (dwlsmv_delta_cols); every correlation is at most
*! quadratic in one loading or variance, so its error is rounding alone,
*! about eps / h per entry. Module 3 sized the criterion on the unweighted Delta:
*! noise ~1e-10, the smallest genuine ratio 0.226 across the tier-2 oracle's
*! identified configurations -- two decades over the noise, seven under the
*! signal. Derivation @sec-bread-rank (2026-09-26) shows the row weights of
*! W_D^1/2 cancel from the noise bound, nu <= sqrt(P) max_k ||E_k|| / ||Delta_k||,
*! and measures, on the Monte Carlo reference's stress designs, noise 9e-13 to
*! 4e-12, genuine ratios 0.005 to 0.12 and exact deficiencies <= 2e-14.
*! test_breaddiag.do holds production's own Delta to the same margins.
*!
*! A NUMERICAL criterion, not an algebraic one: an exact deficiency reads at
*! the noise floor, not at zero. Structural deficiency (npar > D, a zero
*! column) is detected from its cause by the caller, never by this cutoff.
*!
*! CONTRACT (rc 3498): M non-empty, no missing value.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_numrank(real matrix M, real colvector sv, real matrix V)
{
    real matrix U

    if (rows(M) < 1 | cols(M) < 1) {
        errprintf("dwlsmv_numrank: M is %f x %f; it must be non-empty.\n", rows(M), cols(M))
        _error("dwlsmv_numrank: M is empty")
    }
    if (missing(M)) {
        _error("dwlsmv_numrank: M contains missing values")
    }

    U  = J(0, 0, 0)
    V  = J(0, 0, 0)
    sv = J(0, 1, 0)
    svd(M, U, sv, V)
    if (sv[1] <= 0) {
        return(0)
    }
    return(sum(sv :> 1e-8 * sv[1]))
}

end
