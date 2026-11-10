*! dwlsmv -- Module 11: sem's structural names for D6's check
*!
*! Executed under StataMP 18 on 2026-10-04: tests/unit/test_translate.do.
*! Called by _dwlsmv_translate.ado.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// sem's structural names, less _cons and indicator error variances
string colvector dwlsmv_semnames()
{
    string matrix st
    string colvector ind, out
    real scalar j
    string scalar co

    st  = st_matrixcolstripe("e(b)")
    ind = tokens(st_global("e(oyvars)"))'
    out = J(0, 1, "")
    for (j = 1; j <= rows(st); j++) {
        co = st[j, 2]
        if (co == "_cons") {
            continue
        }
        if (substr(co, 1, 6) == "var(e.") {
            if (dwlsmv_tr_pos(ind, substr(co, 7, strlen(co) - 7))) {
                continue
            }
        }
        // as colfullnames shows it: y1:F1, F2:F1, /var(F1)
        out = out \ (st[j, 1] == "/" ? "/" + co : st[j, 1] + ":" + co)
    }
    return(out)
}

end
