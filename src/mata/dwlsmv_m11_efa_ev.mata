*! dwlsmv -- Module 11: e(Ev), the eigenvalues of the sample correlation matrix (release A, 2026-11-18)
*!
*!     void dwlsmv_m11_efa_ev(string scalar evn)
*!
*! Called by dwlsmv.ado for an efa() fit, after dwlsmv_m11_fit() and before dwlsmv_m11_clear(),
*! while the Stage 1 external dwlsmv_m11_st is alive. Builds the p x p sample correlation matrix
*! R from Stage 1's s (the polychoric / tetrachoric entries, stat_type 2, at (var_i, var_j); unit
*! diagonal; the one-empty-cell repair and keepbound pairs enter as Stage 1 has them) and fills
*! the Stata matrix evn with its eigenvalues, descending, as a 1 x p row vector named Ev1 .. Evp:
*! the layout Stata's screeplot reads after pca and factor (e(Ev), with "eigen" in e(properties)).
*! R can be indefinite (a pairwise polychoric matrix need not be positive semidefinite), so
*! negative eigenvalues are possible; they are not clipped. The eigenvalues sum to p.
*!
*! symeigenvalues(), sort(), I(), st_matrix() and st_matrixcolstripe() all predate Stata 16.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_m11_efa_ev(string scalar evn)
{
    struct dwlsmv_stage1 scalar st
    real matrix R
    real rowvector L
    real scalar d, p

    st = *findexternal("dwlsmv_m11_st")
    p  = st.map.p
    R  = I(p)
    for (d = 1; d <= st.map.D; d++) {
        if (st.map.stat_type[d] == 2) {
            R[st.map.var_i[d], st.map.var_j[d]] = st.s[d]
            R[st.map.var_j[d], st.map.var_i[d]] = st.s[d]
        }
    }
    L = sort(symeigenvalues(R)', -1)'
    st_matrix(evn, L)
    st_matrixcolstripe(evn, (J(p, 1, ""), "Ev" :+ strofreal(1::p)))
}

end
