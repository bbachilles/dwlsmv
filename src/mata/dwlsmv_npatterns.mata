*! dwlsmv -- the number of distinct response patterns (Module 10, spec 4.5.6)
*!
*! Executed under StataMP 18 on 2026-09-25. tests/unit/test_structural.do:
*! case D gives 139 and the aliased design 82, both the derivation's counts.
*!
*!     real scalar dwlsmv_npatterns(real matrix Y, real colvector K | real colvector x)
*!
*! WITH A COVARIATE (decided 2026-11-15, the owner's: Module 10 in conditional
*! mode). The optional third argument x is the N x 1 covariate of a conditional
*! fit (one binary 0/1 x). Then psi_i is a function of (y_i, x_i) alone
*! (derivation @sec-condx (b), @sec-condx-m10 M1), so a pattern is one row of
*! (Y, x) and the count is taken over those N x (p + 1) rows. Counting on Y
*! alone is wrong there: a design in which every y-pattern occurs in both
*! groups has n_pat = 2 n_pat(y), and the y-only count can flag a deficiency
*! that is not there. x must have one row per row of Y, no missing value and
*! only the codes 0 and 1 (rc 3498). Without x the function is unchanged.
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

real scalar dwlsmv_npatterns(real matrix Y, real colvector K, | real colvector x)
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
    if (args() >= 3) {
        if (rows(x) != rows(Y)) {
            errprintf("dwlsmv_npatterns: x has %f rows but Y has %f.\n", rows(x), rows(Y))
            _error("dwlsmv_npatterns: x and Y disagree on the number of rows")
        }
        if (missing(x) > 0) {
            errprintf("dwlsmv_npatterns: x has %f missing values; apply listwise deletion first.\n",
                      missing(x))
            _error("dwlsmv_npatterns: x has missing values")
        }
        if (sum((x :!= 0) :& (x :!= 1)) > 0) {
            errprintf("dwlsmv_npatterns: x must be coded 0 and 1; %f values are neither.\n",
                      sum((x :!= 0) :& (x :!= 1)))
            _error("dwlsmv_npatterns: x must be coded 0 and 1")
        }
        return(rows(uniqrows((Y, x))))
    }
    return(rows(uniqrows(Y)))
}

end
