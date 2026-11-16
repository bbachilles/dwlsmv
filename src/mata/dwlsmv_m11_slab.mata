*! dwlsmv -- Module 11: Module 10's statistic labels with variable names
*!
*! Executed under StataMP 18 on 2026-10-04: tests/unit/test_dwlsmv.do.
*! Used by dwlsmv_m11_fit() for dwlsmv_diagadvice().
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// a pair's, threshold's or slope's map label with variable names in place of yj. A conditional map's
// slope rows (stat_type 3, 2026-11-15) read "y1~x", the covariate's name being the optional third
// argument (default "x"); before this the slope rows fell into the pair branch and named a
// second indicator that is not there (var_j is the covariate's index).
string colvector dwlsmv_m11_slab(struct dwlsmv_statmap scalar map, string rowvector ind,
                                 | string scalar cvn)
{
    string colvector L
    string scalar cv
    real scalar k

    cv = "x"
    if (args() >= 3) {
        cv = cvn
    }
    L = J(map.D, 1, "")
    for (k = 1; k <= map.D; k++) {
        if (map.stat_type[k] == 1) {
            L[k] = ind[map.var_i[k]] + "|t" + strofreal(map.thr_k[k])
        }
        else if (map.stat_type[k] == 3) {
            L[k] = ind[map.var_i[k]] + "~" + cv
        }
        else {
            L[k] = ind[map.var_i[k]] + "~~" + ind[map.var_j[k]]
        }
    }
    return(L)
}

end
