*! dwlsmv -- the number of distinct response patterns (Module 10, spec 4.5.6)
*!
*! Executed under StataMP 18 on 2026-09-25. tests/unit/test_structural.do:
*! case D gives 139 and the aliased design 82, both the derivation's counts.
*!
*!     real scalar dwlsmv_npatterns(real matrix Y, real colvector K)
*!
*! THE DEFINITION IS SPEC 4.5.6's, and it is narrow on purpose. Y is the matrix
*! exactly as dwls_sample_stats() receives it: after listwise deletion,
*! restricted to the p indicators whose thresholds and polychorics make up s,
*! in K's column order, coded 0 .. K_j - 1. A pattern is one row of that
*! matrix. Variables outside the p indicators and rows removed by listwise
*! deletion do not count, because psi_i is a function of exactly this row and
*! nothing else (derivation @sec-gamma-rank, Lemma 1). The count feeds the
*! structural bound rank(Gamma-hat) <= min(D, n_pat - 1).
*!
*! It takes the same (Y, K) as dwls_sample_stats() and applies the same
*! per-column validation, dwlsmv_catcounts(): a missing value, a non-integer
*! code or a code outside 0 .. K_j - 1 aborts with rc 3498. Counting `.` as a
*! code would inflate n_pat and make the structural check pass wrongly.
*!
*! v1 ONLY. With a continuous indicator the pattern set is not finite; the v2
*! mixed-type extension must re-derive the bound, not reuse this.
*!
*! Requires dwlsmv_catcounts.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_npatterns(real matrix Y, real colvector K)
{
    real scalar j

    if (cols(Y) != rows(K)) {
        errprintf("dwlsmv_npatterns: Y has %f columns but K names %f variables.\n",
                  cols(Y), rows(K))
        _error("dwlsmv_npatterns: Y and K disagree on the number of variables")
    }
    for (j = 1; j <= cols(Y); j++) {
        (void) dwlsmv_catcounts(Y[., j], K[j])
    }
    return(rows(uniqrows(Y)))
}

end
