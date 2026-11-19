*! dwlsmv -- Module 11: the threshold rows of an efa() fit, for display (release A)
*!
*! Run 2026-11-18 in this form (suite; M13): the cut that kept every column but the
*! thresholds, and the cut of a matrix by rows, went to dwlsmv_m11_ctab.mata's planner and are
*! removed; only the efa() threshold table (_dwlsmv_efathr) still calls this.
*!
*!     real scalar dwlsmv_m11_dispsub(string scalar bn, string scalar vn, string scalar obn,
*!                                    string scalar ovn)
*!
*! For display only; e() is never touched. A threshold column is one whose equation starts with "/"
*! and is not "/" alone: "/y1:cut1" at G = 1, "/y1:1bn.g#c.cut1" at G > 1 (spec 4.8, D6). Variances,
*! covariances, omega and alpha have the equation "/" and are not thresholds.
*! bn 1 x n b; vn n x n V ("" for none); obn, ovn receive b and V cut to the threshold columns,
*! with the stripes of the kept columns. Returns the number of threshold columns (0: nothing filled).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_m11_dispsub(string scalar bn, string scalar vn, string scalar obn,
                               string scalar ovn)
{
    real rowvector b
    real matrix V
    string matrix stp
    real colvector thr, sel
    real scalar k

    b   = st_matrix(bn)
    stp = st_matrixcolstripe(bn)
    thr = (substr(stp[., 1], 1, 1) :== "/") :& (stp[., 1] :!= "/")
    sel = selectindex(thr)
    k   = length(sel)
    if (k == 0) {
        return(0)
    }
    st_matrix(obn, b[1, sel'])
    st_matrixcolstripe(obn, stp[sel, .])
    if (vn != "") {
        V = st_matrix(vn)
        st_matrix(ovn, V[sel, sel])
        st_matrixcolstripe(ovn, stp[sel, .])
        st_matrixrowstripe(ovn, stp[sel, .])
    }
    return(k)
}

end
