*! dwlsmv -- the small-count family: category counts (Module 10, spec Section 5)
*!
*! Executed under StataMP 18 on 2026-09-26. tests/unit/test_countdiag.do: a
*! hand-built Y, ties, all-binary, an empty category, and case D against
*! Stata's own -count- (boundary 706 at y1's category 1, interior 1500 at
*! y2's category 3).
*!
*!     real matrix dwlsmv_countdiag(real matrix Y, real colvector K)
*!
*! Spec Section 5, Module 10: small counts are a SEPARATE family from Gamma's
*! conditioning. A boundary category with one expected observation leaves that
*! threshold's mean SE at about 2x its sampling SD (derivation
*! @sec-gamma-smallcount) while Gamma-hat can be perfectly well conditioned.
*! This routine is never reported as a Gamma-conditioning failure, and it
*! routes to the sparse-cell policy, not the Gamma policy.
*!
*! Returns 2 x 4, one row per quantity, always in this order, in the
*! [flag, value, i, j] layout of dwlsmv_admissible():
*!
*!     row  quantity                            value     i, j
*!      1   smallest BOUNDARY-category count    count     variable, category
*!      2   smallest INTERIOR-category count    count     variable, category
*!
*! Boundary means category 0 or K_j - 1; interior means 1 .. K_j - 2. They are
*! apart because they fail differently: a small boundary category drives one
*! threshold to the tail, which has an exact binomial reference
*! (@sec-gamma-smallcount); a small interior category squeezes two thresholds
*! together and has none (E-060, 2026-09-25, b1 against b2).
*!
*! The category is its POSITION, 1 .. K_j (code + 1), so that 0 can mean "no
*! location" as it does throughout Module 10. A tie goes to the first
*! (variable, category) in s's order.
*!
*! VALUES ONLY. Every flag is .n: evaluated, no flag defined in v1. That is
*! not `.`, which means not evaluated (decided 2026-09-26). A threshold waits
*! on @sec-gamma-smallcount and the S partition of the Monte Carlo reference.
*! Row 2 on an all-binary Y has no interior category to look at: it is
*! (.n, ., 0, 0), evaluated with nothing to report.
*!
*! WHAT IS DELIBERATELY NOT HERE: pair-cell counts. An empty cell is not by
*! itself damage -- 11,978 of the 14,000 P1 Monte Carlo replicates have one,
*! a structurally empty corner at high rho, and their SEs track -- and (Y, K)
*! alone cannot tell a structural empty cell from a damaging one. Whether
*! cell rows exist at all waits on the S results (decided 2026-09-26).
*! dwlsmv_cellsparse() already classifies a single pair's table for Module 2.
*!
*! AN EMPTY DECLARED CATEGORY READS 0. Module 1 aborts on one; this routine
*! diagnoses and never acts, so it reports the 0 rather than refusing.
*!
*! Same (Y, K) and same per-column validation as dwls_sample_stats(),
*! dwlsmv_npatterns() and dwlsmv_thralias(): dwlsmv_catcounts(), rc 3498 on a
*! missing, non-integer or out-of-range code.
*!
*! Requires dwlsmv_catcounts.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_countdiag(real matrix Y, real colvector K)
{
    real matrix A
    real colvector n
    real scalar j, c

    if (cols(Y) != rows(K)) {
        errprintf("dwlsmv_countdiag: Y has %f columns but K names %f variables.\n",
                  cols(Y), rows(K))
        _error("dwlsmv_countdiag: Y and K disagree on the number of variables")
    }

    A = (.n, ., 0, 0 \ .n, ., 0, 0)
    for (j = 1; j <= cols(Y); j++) {
        n = dwlsmv_catcounts(Y[., j], K[j])
        for (c = 1; c <= K[j]; c++) {
            if (c == 1 | c == K[j]) {
                if (n[c] < A[1, 2]) {
                    A[1, .] = (.n, n[c], j, c)
                }
            }
            else if (n[c] < A[2, 2]) {
                A[2, .] = (.n, n[c], j, c)
            }
        }
    }
    return(A)
}

end
