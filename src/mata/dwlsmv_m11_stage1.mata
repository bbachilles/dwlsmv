*! dwlsmv -- Module 11: D9's category recoding, then Stage 1
*!
*! Executed under StataMP 18 on 2026-10-04: tests/unit/test_dwlsmv.do.
*! Each indicator's distinct observed values, sorted, become 0..K-1 (spec D9). Called under capture noisily (D8a).
*! With a fourth argument naming the covariate (a conditional fit, one binary 0/1 x, rows already
*! marked out jointly with the indicators), Stage 1 is dwls_sample_stats_cx(Y, K, x), which has no
*! keep_bound; keepbound is then ignored. The covariate is NOT recoded (spec 4.8, Coding): its
*! coding check is dwlsmv_cx_check's, rc 3498, and the command's own check (rc 198) runs first.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_m11_stage1(string scalar touse, string scalar indlist, real scalar keepbound,
                       | string scalar covar)
{
    string rowvector ind
    real matrix X, Y, Xc
    real colvector K, u, x, r
    real scalar j, k
    pointer() scalar pst, pY, px

    ind = tokens(indlist)
    X   = st_data(., ind, touse)
    Y   = J(rows(X), cols(X), 0)
    K   = J(cols(X), 1, 0)
    for (j = 1; j <= cols(X); j++) {
        x = X[., j]
        u = uniqrows(x)
        r = J(rows(x), 1, 0)
        for (k = 2; k <= rows(u); k++) {
            r = r + (x :>= u[k])
        }
        Y[., j] = r
        K[j]    = rows(u)
    }
    pst = crexternal("dwlsmv_m11_st")
    if (pst == NULL) {
        pst = findexternal("dwlsmv_m11_st")
    }
    pY = crexternal("dwlsmv_m11_Y")
    if (pY == NULL) {
        pY = findexternal("dwlsmv_m11_Y")
    }
    *pY  = Y
    if (args() >= 4) {
        if (covar != "") {
            Xc   = st_data(., covar, touse)
            // kept for Module 10 (2026-11-15): the pattern count and the per-group counts need x
            px = crexternal("dwlsmv_m11_x")
            if (px == NULL) {
                px = findexternal("dwlsmv_m11_x")
            }
            *px  = Xc
            *pst = dwls_sample_stats_cx(Y, K, Xc)
            return
        }
    }
    *pst = dwls_sample_stats(Y, K, keepbound)
}

end
