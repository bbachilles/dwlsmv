*! dwlsmv -- Module 11: the standardized solution (dwlsmv_m11_std)
*!
*! Executed under StataMP 18 on 2026-10-05: tests/unit/test_standardized.do.
*!
*!   dwlsmv_m11_std(bn, vn, cn)
*!       D7 step (3). From e(b), e(V) and e(indicators) alone, so replay needs
*!       nothing stored (D8d): the standardized values g(theta) and their
*!       delta-method covariance J_g V J_g' (derivation @sec-standardize), in
*!       e(b)'s order, followed by one Theta row /var(e.y) per indicator.
*!       Writes the Stata matrices bn (1 x m), vn (m x m) and, for the rows
*!       whose standardized value is a constant (SE exactly 0), cn; the local
*!       nfixs holds that row count.
*!
*! The parameter table is rebuilt from e(b)'s stripe (D6): y:F a loading,
*! T:F a path, /var(F) /var(e.F) /cov(.,.) a psi entry, /y:cutk a threshold.
*! Every factor has a psi_ff entry, which fixes the factor order. Fixed
*! parameters keep their e(b) value and have zero rows in e(V) (D5), so the
*! full-table J_g gives the free-block result. v1 is recursive (spec 4.8):
*! B is nilpotent, and (I-B)^-1 = I + B + ... + B^(n-1) is computed exactly,
*! so a structural zero of J_g stays zero and an exogenous fixed psi_ff = 1
*! standardizes to 1 with SE exactly 0.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_m11_std(string scalar bn, string scalar vn, string scalar cn)
{
    real rowvector b, g, kind, r1, r2
    real matrix V, L, B, P, A, Bk, Phi, dL, dB, dP, dPhi, Jg, Vs, Cns, LPhi
    real colvector s, ds
    string matrix stp, out
    string rowvector ind, fac
    string scalar eq, nm, x, y
    real scalar n, q, nI, nf, i, j, k, m, f, t, ss, dss, nfix, c

    b   = st_matrix("e(b)")
    V   = st_matrix("e(V)")
    stp = st_matrixcolstripe("e(b)")
    ind = tokens(st_global("e(indicators)"))
    q   = cols(b)
    nI  = cols(ind)

    // factors, in the order of their psi_ff entries
    fac = J(1, 0, "")
    for (i = 1; i <= q; i++) {
        if (stp[i, 1] == "/" & substr(stp[i, 2], 1, 4) == "var(") {
            x = substr(stp[i, 2], 5, strlen(stp[i, 2]) - 5)
            if (substr(x, 1, 2) == "e.") {
                x = substr(x, 3, .)
            }
            fac = fac, x
        }
    }
    nf = cols(fac)

    // kind: 1 loading, 2 path, 3 psi, 4 threshold; r1, r2 its row and column
    kind = J(1, q, 0)
    r1   = J(1, q, 0)
    r2   = J(1, q, 0)
    for (i = 1; i <= q; i++) {
        eq = stp[i, 1]
        nm = stp[i, 2]
        if (eq == "/") {
            kind[i] = 3
            if (substr(nm, 1, 4) == "var(") {
                x = subinstr(substr(nm, 5, strlen(nm) - 5), "e.", "")
                y = x
            }
            else {
                x = substr(nm, 5, strlen(nm) - 5)
                y = substr(x, strpos(x, ",") + 1, .)
                x = substr(x, 1, strpos(x, ",") - 1)
                x = subinstr(x, "e.", "")
                y = subinstr(y, "e.", "")
            }
            r1[i] = selectindex(fac :== x)
            r2[i] = selectindex(fac :== y)
        }
        else if (substr(eq, 1, 1) == "/") {
            kind[i] = 4
        }
        else if (anyof(ind, eq)) {
            kind[i] = 1
            r1[i] = selectindex(ind :== eq)
            r2[i] = selectindex(fac :== nm)
        }
        else {
            kind[i] = 2
            r1[i] = selectindex(fac :== eq)
            r2[i] = selectindex(fac :== nm)
        }
    }
    if (anyof(kind, 0)) {
        _error("dwlsmv: internal error, e(b) stripe not recognized")
    }

    L = J(nI, nf, 0)
    B = J(nf, nf, 0)
    P = J(nf, nf, 0)
    for (i = 1; i <= q; i++) {
        if (kind[i] == 1) {
            L[r1[i], r2[i]] = b[i]
        }
        else if (kind[i] == 2) {
            B[r1[i], r2[i]] = b[i]
        }
        else if (kind[i] == 3) {
            P[r1[i], r2[i]] = b[i]
            P[r2[i], r1[i]] = b[i]
        }
    }
    // (I-B)^-1 by the finite series; B is nilpotent (v1 recursive)
    A  = I(nf)
    Bk = I(nf)
    for (k = 1; k < nf; k++) {
        Bk = Bk * B
        A  = A + Bk
    }
    if (nf > 1) {
        if (any(Bk * B)) {
            _error("dwlsmv: internal error, B is not nilpotent")
        }
    }
    Phi  = A * P * A'
    s    = sqrt(diagonal(Phi))
    LPhi = L * Phi

    m = q + nI
    g = J(1, m, .)
    for (i = 1; i <= q; i++) {
        if (kind[i] == 1) {
            g[i] = b[i] * s[r2[i]]
        }
        else if (kind[i] == 2) {
            g[i] = b[i] * s[r2[i]] / s[r1[i]]
        }
        else if (kind[i] == 3) {
            g[i] = b[i] / (s[r1[i]] * s[r2[i]])
        }
        else {
            g[i] = b[i]
        }
    }
    for (j = 1; j <= nI; j++) {
        g[q + j] = 1 - LPhi[j, .] * L[j, .]'
    }

    // J_g, one column per e(b) entry (derivation @sec-standardize)
    Jg = J(m, q, 0)
    for (c = 1; c <= q; c++) {
        if (kind[c] == 4) {
            Jg[c, c] = 1
            continue
        }
        dL = J(nI, nf, 0)
        dB = J(nf, nf, 0)
        dP = J(nf, nf, 0)
        if (kind[c] == 1) {
            dL[r1[c], r2[c]] = 1
        }
        else if (kind[c] == 2) {
            dB[r1[c], r2[c]] = 1
        }
        else {
            dP[r1[c], r2[c]] = 1
            dP[r2[c], r1[c]] = 1
        }
        dPhi = A * dB * Phi + Phi * dB' * A' + A * dP * A'
        ds   = diagonal(dPhi) :/ (2 :* s)
        for (i = 1; i <= q; i++) {
            if (kind[i] == 1) {
                j = r1[i]
                f = r2[i]
                Jg[i, c] = dL[j, f] * s[f] + b[i] * ds[f]
            }
            else if (kind[i] == 2) {
                t = r1[i]
                f = r2[i]
                Jg[i, c] = dB[t, f] * s[f] / s[t] +
                           b[i] * (ds[f] / s[t] - s[f] * ds[t] / s[t]^2)
            }
            else if (kind[i] == 3) {
                f   = r1[i]
                t   = r2[i]
                ss  = s[f] * s[t]
                dss = s[t] * ds[f] + s[f] * ds[t]
                Jg[i, c] = dP[f, t] / ss - b[i] * dss / ss^2
            }
        }
        for (j = 1; j <= nI; j++) {
            Jg[q + j, c] = -2 * LPhi[j, .] * dL[j, .]' - L[j, .] * dPhi * L[j, .]'
        }
    }
    // an exogenous factor's psi_ff / Phi_ff is identically 1 (Phi_ff = psi_ff):
    // its row is zero by construction, not by rounding of 1/psi - psi/psi^2
    for (i = 1; i <= q; i++) {
        if (kind[i] == 3) {
            if (r1[i] == r2[i]) {
                if (!any(B[r1[i], .])) {
                    Jg[i, .] = J(1, q, 0)
                    g[i]     = 1
                }
            }
        }
    }
    Vs = Jg * V * Jg'
    Vs = (Vs + Vs') / 2

    out = stp \ J(nI, 2, "/")
    for (j = 1; j <= nI; j++) {
        out[q + j, 2] = "var(e." + ind[j] + ")"
    }
    // rows with a constant standardized value (SE exactly 0) are posted as
    // constraints, as the raw table's fixed rows are
    nfix = sum(diagonal(Vs) :== 0)
    if (nfix) {
        Cns = J(nfix, m + 1, 0)
        k   = 0
        for (i = 1; i <= m; i++) {
            if (Vs[i, i] == 0) {
                k++
                Cns[k, i]     = 1
                Cns[k, m + 1] = g[i]
            }
        }
        st_matrix(cn, Cns)
        st_matrixcolstripe(cn, out \ ("", "_r"))
    }
    st_matrix(bn, g)
    st_matrixcolstripe(bn, out)
    st_matrix(vn, Vs)
    st_matrixcolstripe(vn, out)
    st_matrixrowstripe(vn, out)
    st_local("nfixs", strofreal(nfix))
}

end
