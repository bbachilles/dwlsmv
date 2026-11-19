*! dwlsmv -- Module 11: standard errors of the rotated EFA solution (release A, unit A5, posting layer)
*!
*! Drafted 2026-11-17 by a read-only agent; compiled and executed by the main session under
*! StataMP 18 (test_rotse_post.do, test_efa_surface.do, the suite), 2026-11-17.
*!
*!     real scalar dwlsmv_m11_efa_se(string scalar bn, string scalar vn, string scalar canon,
*!                                   string scalar rbn, string scalar lamn, string scalar phin,
*!                                   string scalar vrn, string scalar efalist)
*!
*! Called by dwlsmv.ado after dwlsmv_m11_efa_rotate(), while the Stage 1 external
*! (dwlsmv_m11_st) still exists. Inputs, all by name:
*!     bn, vn   the UNROTATED echelon fit's e(b) and e(V) matrices (dwlsmv_m11_fit(): every
*!              table row, fixed rows with exact zeros in V)
*!     canon    the echelon canonical text the fit used
*!     rbn      the rotated e(b) (dwlsmv_m11_efa_rotate(): loadings, factor correlations, thresholds)
*!     lamn, phin   the rotated Lambda and Phi, the GPA stall point (A3.15)
*!     efalist  efa()'s namelist (labels the rotated columns)
*! and the local rot_eps. Fills vrn with the rotated covariance, in rbn's column order and with
*! rbn's stripes on rows and columns, and the local rot_cond (the condition number of the
*! bordered matrix [Delta*'W Delta*, R'; R, 0], A3.2: stored, no cutoff, no note).
*!
*! Returns a status:
*!     0  vrn filled (A3.5, A3.13: J ACOV_u J' at the stall point)
*!     1  A3.2: the section 4.8 rank cut fails on W^1/2 Delta* Z, or R lacks full row rank
*!     2  A3.4: R G_T is singular (degenerate critical point), or Delta*/R hold a missing value
*!     3  A3.4: the unrotated covariance of (tau, a) holds a missing value
*!
*! Nothing is posted here and nothing is withheld here: dwlsmv.ado turns the status into e()
*! and the message. The status 1 to 3 paths post the rotated e(b) without e(V) (A3.2, A3.4).
*!
*! ORDERS. theta_u = (tau in the fit's table order, then the free echelon loadings a = A[i,j],
*! i >= j, j outer) is built by NAME from the fit's stripe, so the order of the fit's table is
*! not assumed. theta* = (tau, vec Lambda column-major, phi strict lower triangle column-major)
*! is dwlsmv_rotse.mata's. The permutation from theta* to the posted rotated e(b) is built by
*! name from rbn's stripe. Delta* is route 1 of RS7: dwlsmv_delta_cols() on the echelon table's
*! own rows, loadings and factor covariances set to the rotated values, no boundary change (A2b).
*!
*! W = diag(Gamma)^-1 (spec W_D), the fit's weighting; e(V) scale, so ACOV_u is the fit's e(V)
*! (the product J ACOV_u J' is linear in it, dwlsmv_rotse_cov()).
*!
*! Load order: dwlsmv_types, dwlsmv_rotse, dwlsmv_geomin (via rotse), dwlsmv_delta_cols,
*! dwls_identify, dwlsmv_numrank. All load from src/mata in the library build.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_m11_efa_se(string scalar bn, string scalar vn, string scalar canon,
                              string scalar rbn, string scalar lamn, string scalar phin,
                              string scalar vrn, string scalar efalist)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    string rowvector lines, fac, ind, lab
    string matrix stp, rst
    string scalar nm, inner
    real rowvector b
    real matrix A, L, Phi, Vu, ACOVu, Vs, Vr, Ds, R, J0
    real colvector tcol, acol, pi, sel, vals, w, ix, tauSel
    real scalar i, j, k, c, r, p, m, nth, nu, nstar, eps, ok, ratio, cnd, q, t, a1, a2, cm

    st    = *findexternal("dwlsmv_m11_st")
    lines = select(tokens(canon, ";"), tokens(canon, ";") :!= ";")
    fac   = tokens(substr(lines[1], 8, .))
    ind   = tokens(substr(lines[2], 5, .))
    lab   = tokens(efalist)
    p     = cols(ind)
    m     = cols(fac)
    b     = st_matrix(bn)
    stp   = st_matrixcolstripe(bn)
    rst   = st_matrixcolstripe(rbn)
    eps   = strtoreal(st_local("rot_eps"))
    L     = st_matrix(lamn)
    Phi   = st_matrix(phin)
    if (rows(L) != p | cols(L) != m | rows(Phi) != m | cols(Phi) != m | missing(eps)) {
        _error("dwlsmv: internal error, the rotated solution does not conform (efa_se)")
    }

    // the unrotated A (Phi = I) from bn's loadings, as dwlsmv_m11_efa_rotate() reads it
    A = J(p, m, .)
    for (i = 1; i <= cols(b); i++) {
        j = selectindex(ind :== stp[i, 1])
        k = selectindex(fac :== stp[i, 2])
        if (length(j) == 1 & length(k) == 1) {
            A[j, k] = b[i]
        }
    }
    if (missing(A)) {
        _error("dwlsmv: internal error, a loading of the EFA block is not in e(b) (efa_se)")
    }

    // ---- theta_u = (tau, a) by name, and ACOV_u from the fit's e(V) ----------------
    tcol = J(0, 1, .)
    for (i = 1; i <= cols(b); i++) {
        if (substr(stp[i, 1], 1, 1) == "/" & substr(stp[i, 2], 1, 3) == "cut") {
            tcol = tcol \ i
        }
    }
    nth  = rows(tcol)
    acol = J(0, 1, .)
    for (j = 1; j <= m; j++) {
        for (i = j; i <= p; i++) {
            ix = selectindex((stp[., 1] :== ind[i]) :& (stp[., 2] :== fac[j]))
            if (length(ix) != 1) {
                _error("dwlsmv: internal error, a free loading is not unique in e(b) (efa_se)")
            }
            acol = acol \ ix
        }
    }
    nu    = nth + rows(acol)
    nstar = nth + p * m + m * (m - 1) / 2
    Vu    = st_matrix(vn)
    if (rows(Vu) != cols(b) | cols(Vu) != cols(b)) {
        _error("dwlsmv: internal error, the unrotated e(V) does not match e(b) (efa_se)")
    }
    ACOVu = Vu[(tcol \ acol), (tcol \ acol)]
    // A3.4: a missing value in the unrotated covariance (a singular bread) means no SEs
    if (missing(ACOVu)) {
        return(3)
    }

    // ---- the permutation from theta* to the posted rotated e(b) ---------------------
    if (cols(rst) != 2 | rows(rst) != nstar) {
        _error("dwlsmv: internal error, the rotated e(b) has the wrong size (efa_se)")
    }
    pi = J(nstar, 1, .)
    for (c = 1; c <= nstar; c++) {
        if (substr(rst[c, 1], 1, 1) == "/" & substr(rst[c, 2], 1, 3) == "cut") {
            ix = selectindex((stp[tcol, 1] :== rst[c, 1]) :& (stp[tcol, 2] :== rst[c, 2]))
            if (length(ix) != 1) {
                _error("dwlsmv: internal error, a threshold of the rotated e(b) is not in the fit (efa_se)")
            }
            pi[c] = ix
        }
        else if (substr(rst[c, 2], 1, 4) == "cov(") {
            nm    = rst[c, 2]
            inner = substr(nm, 5, strlen(nm) - 5)
            cm    = strpos(inner, ",")
            a1    = selectindex(lab :== substr(inner, 1, cm - 1))
            a2    = selectindex(lab :== substr(inner, cm + 1, .))
            if (length(a1) != 1 | length(a2) != 1) {
                _error("dwlsmv: internal error, a factor correlation is not in efa() (efa_se)")
            }
            pi[c] = nth + p * m + dwlsmv_rotse_phidx(a1, a2, m)
        }
        else {
            j = selectindex(ind :== rst[c, 1])
            k = selectindex(lab :== rst[c, 2])
            if (length(j) != 1 | length(k) != 1) {
                _error("dwlsmv: internal error, a rotated loading is not in the block (efa_se)")
            }
            pi[c] = nth + (k - 1) * p + j
        }
    }
    if (any(sort(pi, 1) :!= (1::nstar))) {
        _error("dwlsmv: internal error, the rotated e(b) is not a permutation of theta* (efa_se)")
    }

    // ---- Delta* at the rotated point, RS7 route 1 -----------------------------------
    par    = dwls_identify(subinstr(canon, ";", char(10)), st.map)
    tauSel = selectindex(par.ptype :== 4)
    if (rows(tauSel) != nth) {
        _error("dwlsmv: internal error, the table's thresholds differ from e(b)'s (efa_se)")
    }
    sel  = tauSel
    vals = par.value
    for (k = 1; k <= m; k++) {
        for (j = 1; j <= p; j++) {
            ix = selectindex((par.ptype :== 1) :& (par.row :== j) :& (par.col :== k))
            if (length(ix) != 1) {
                _error("dwlsmv: internal error, a loading row is not unique in the table (efa_se)")
            }
            sel     = sel \ ix
            vals[ix] = L[j, k]
        }
    }
    for (c = 1; c <= m - 1; c++) {
        for (r = c + 1; r <= m; r++) {
            ix = selectindex((par.ptype :== 2) :& (((par.row :== r) :& (par.col :== c)) :|
                                                   ((par.row :== c) :& (par.col :== r))))
            if (length(ix) != 1) {
                _error("dwlsmv: internal error, a factor covariance row is not unique in the table (efa_se)")
            }
            sel      = sel \ ix
            vals[ix] = Phi[r, c]
        }
    }
    Ds = dwlsmv_delta_cols(vals, par, st.map, sel)
    R  = (J(m * (m - 1), nth, 0), dwlsmv_rotse_R(L, Phi, eps))
    if (missing(Ds) | missing(R)) {
        return(2)
    }

    // ---- A3.2: the rank cut on W^1/2 Delta* Z, and the condition number ------------
    w  = 1 :/ diagonal(st.Gamma)
    ok = dwlsmv_rotse_rank(Ds, w, R, ratio, cnd)
    st_local("rot_cond", strofreal(cnd, "%21.0g"))
    if (!ok) {
        return(1)
    }

    // ---- A3.13: J ACOV_u J' at the stall point, then onto the posted order ---------
    Vs = dwlsmv_rotse(A, L, Phi, eps, nth, ACOVu, J0)
    if (missing(Vs)) {
        return(2)
    }
    Vr = Vs[pi, pi]
    st_matrix(vrn, Vr)
    st_matrixcolstripe(vrn, rst)
    st_matrixrowstripe(vrn, rst)
    return(0)
}

end
