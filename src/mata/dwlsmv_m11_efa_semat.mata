*! dwlsmv -- Module 11: standard-error matrices of the rotated estimates, for the efa() display (release A)
*!
*!     void dwlsmv_m11_efa_semat(string scalar bn, string scalar vn, string scalar ln,
*!                               string scalar sln, string scalar spn)
*!
*! bn, vn  copies of the posted rotated e(b) and e(V) (e(V) is posted over every e(b) column);
*! ln      a copy of e(Lambda), whose stripes name the indicators (rows) and factors (columns);
*! sln     filled with the p x m matrix of sqrt(V[i,i]) at each loading, e(Lambda)'s stripes;
*! spn     filled with the m x m matrix of sqrt(V[i,i]) at each factor correlation, both triangles,
*!         the diagonal missing (fixed at 1), the factor names on both stripes.
*! Entries are found by NAME from e(b)'s stripe (loading: equation = indicator, name = factor;
*! correlation: name cov(a,b)), never by position. A negative variance (rounding) gives missing.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_m11_efa_semat(string scalar bn, string scalar vn, string scalar ln,
                          string scalar sln, string scalar spn)
{
    real matrix V, SL, SP
    string matrix stp, lrow, lcol
    string rowvector ind, lab
    string scalar nm, a, b2
    real scalar i, c, ja, jb, j, k, m

    V    = st_matrix(vn)
    stp  = st_matrixcolstripe(bn)
    lrow = st_matrixrowstripe(ln)
    lcol = st_matrixcolstripe(ln)
    ind  = lrow[., 2]'
    lab  = lcol[., 2]'
    m    = cols(lab)
    SL   = J(cols(ind), m, .)
    SP   = J(m, m, .)
    for (i = 1; i <= rows(stp); i++) {
        nm = stp[i, 2]
        c  = (V[i, i] >= 0 ? sqrt(V[i, i]) : .)
        j  = selectindex(ind :== stp[i, 1])
        k  = selectindex(lab :== nm)
        if (length(j) == 1 & length(k) == 1) {
            SL[j, k] = c
        }
        else if (substr(nm, 1, 4) == "cov(") {
            a  = substr(nm, 5, strpos(nm, ",") - 5)
            b2 = substr(nm, strpos(nm, ",") + 1, strlen(nm) - strpos(nm, ",") - 1)
            ja = selectindex(lab :== a)
            jb = selectindex(lab :== b2)
            if (length(ja) == 1 & length(jb) == 1) {
                SP[ja, jb] = c
                SP[jb, ja] = c
            }
        }
    }
    st_matrix(sln, SL)
    st_matrixrowstripe(sln, lrow)
    st_matrixcolstripe(sln, lcol)
    st_matrix(spn, SP)
    st_matrixrowstripe(spn, lcol)
    st_matrixcolstripe(spn, lcol)
}

end
