*! dwlsmv -- Module 11: D9's category recoding, then Stage 1
*!
*! Executed under StataMP 18 on 2026-10-04: tests/unit/test_dwlsmv.do.
*! Each indicator's distinct observed values, sorted, become 0..K-1 (spec D9). Called under capture noisily (D8a).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_m11_stage1(string scalar touse, string scalar indlist, real scalar keepbound)
{
    string rowvector ind
    real matrix X, Y
    real colvector K, u, x, r
    real scalar j, k
    pointer() scalar pst, pY

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
    *pst = dwls_sample_stats(Y, K, keepbound)
}

end
