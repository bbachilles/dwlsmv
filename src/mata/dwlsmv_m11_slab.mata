*! dwlsmv -- Module 11: Module 10's statistic labels with variable names
*!
*! Executed under StataMP 18 on 2026-10-04: tests/unit/test_dwlsmv.do.
*! Used by dwlsmv_m11_fit() for dwlsmv_diagadvice().
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// a pair's or threshold's map label with variable names in place of yj
string colvector dwlsmv_m11_slab(struct dwlsmv_statmap scalar map, string rowvector ind)
{
    string colvector L
    real scalar k

    L = J(map.D, 1, "")
    for (k = 1; k <= map.D; k++) {
        if (map.stat_type[k] == 1) {
            L[k] = ind[map.var_i[k]] + "|t" + strofreal(map.thr_k[k])
        }
        else {
            L[k] = ind[map.var_i[k]] + "~~" + ind[map.var_j[k]]
        }
    }
    return(L)
}

end
