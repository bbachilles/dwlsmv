*! dwlsmv -- Module 11: the standardized solution at G > 1 (dwlsmv_m11_std_g)
*!
*! Executed under StataMP 18 on 2026-11-01: tests/unit/test_standardized_g.do.
*!
*!   dwlsmv_m11_std_g(bn, vn, cn)
*!       D7 at G > 1 (spec 4.8 O17, derivation @sec-standardize-mg). From e()
*!       alone, so replay needs nothing stored (D8d). Writes the Stata matrices
*!       bn, vn and, for constant rows, cn; the local nfixs holds their count.
*!
*!   dwlsmv_m11_std_gc(b, V, cs, rs, ind, K, gv, gvar, canon, rank, g, Vs, out)
*!       The same from its arguments, filling g, Vs and out. Separate so a test
*!       can feed it another estimate or malformed metadata.
*!
*! The map is the frozen one (SG2): loading omega*lambda*s, psi psi/(s s),
*! alpha alpha/s, tau omega*tau, nu omega*nu, omega the constant 1, and one
*! Theta row 1 - omega^2 V_jj per group and indicator, after e(b)'s rows.
*! Its Jacobian (SG3) is block-structured: an output row of group g reads only
*! group g's columns. Tied members are not merged: e(V) repeats their rows
*! (spec 4.8 G2), so J_rows V_rows J_rows' is the result (SG4). An exogenous
*! factor's psi_ff is the constant 1 even when free (SG5): an invariant, its J row
*! exactly zero rather than cancelling to zero (mutant C3, decisions 2026-11-01), by the same
*! rule as dwlsmv_m11_std.
*!
*! Fails closed (user, 2026-11-01): nothing is displayed unless every row has
*! one group and one kind, the rows make the complete table of e(K), e(canon)'s
*! ties agree with e(b), e(V) and e(rank), and e(V) fits e(b). beta is an
*! internal error: invariance() refuses structural paths (spec 4.8 G5). A
*! factor variance at or below 0 has no scale (SG7 assumes Phi_ff > 0), so that
*! is rc 3498 with its own message.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// the index of s in v, or 0 unless s appears there exactly once
real scalar dwlsmv_sg_at(string scalar s, string rowvector v)
{
    real scalar i, at

    at = 0
    for (i = 1; i <= cols(v); i++) {
        if (v[i] == s) {
            if (at) {
                return(0)
            }
            at = i
        }
    }
    return(at)
}

// a name into the part beside its group tag, and that group; 0 if not exactly one tag
real scalar dwlsmv_sg_tag(string scalar nm, string rowvector tags, string scalar head)
{
    real scalar p, a, b

    p = strpos(nm, "#")
    if (p == 0) {
        head = ""
        return(dwlsmv_sg_at(nm, tags))
    }
    // tag#c.x for tau and lambda, x#tag for the "/" rows: exactly one side is a tag
    a = dwlsmv_sg_at(substr(nm, 1, p - 1), tags)
    b = dwlsmv_sg_at(substr(nm, p + 1, .), tags)
    if ((a > 0) + (b > 0) != 1) {
        head = ""
        return(0)
    }
    head = (a ? substr(nm, p + 1, .) : substr(nm, 1, p - 1))
    return(a + b)
}

// the table position of a row, or 0
real scalar dwlsmv_sg_find(real matrix tab, real scalar kd, real scalar gr,
                           real scalar r, real scalar c)
{
    real scalar i

    for (i = 1; i <= rows(tab); i++) {
        if (tab[i, 1] == kd) {
            if (tab[i, 2] == gr & tab[i, 3] == r & tab[i, 4] == c) {
                return(i)
            }
        }
    }
    return(0)
}

void dwlsmv_m11_std_g(string scalar bn, string scalar vn, string scalar cn)
{
    real rowvector g
    real matrix Vs, Cns
    string matrix out, K_s
    real scalar m, nfix, i, k

    K_s = st_matrixrowstripe("e(K)")
    dwlsmv_m11_std_gc(st_matrix("e(b)"), st_matrix("e(V)"), st_matrixcolstripe("e(b)"),
                      st_matrixrowstripe("e(V)"), tokens(st_global("e(indicators)")),
                      st_matrix("e(K)"), K_s[., 2]', st_global("e(groupvar)"),
                      st_global("e(canon)"), st_numscalar("e(rank)"), g, Vs, out)
    m    = cols(g)
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

void dwlsmv_m11_std_gc(real rowvector b, real matrix V, string matrix cs, string matrix rs,
                       string rowvector ind, real matrix K, string rowvector gv,
                       string scalar gvar, string scalar canon, real scalar rank,
                       real rowvector g, real matrix Vs, string matrix out)
{
    // tab: kind (1 tau, 2 lambda, 3 psi, 4 nu, 5 omega, 6 alpha), group, r, c
    real matrix tab, L, P, Jg, dL, dP, Lk, Pk
    real colvector om, al, s, ds, dom
    real rowvector seen
    string rowvector fac, tags, w, lines
    string colvector ln
    string scalar eq, nm, head, x, y
    real scalar q, nI, nf, G, i, j, f, h, c, gr, m, kd, r, cc, a, nfree, nextra, first, pos, gc, gi, Vjj
    pointer(real matrix) rowvector Lp, Pp
    pointer(real colvector) rowvector omp, alp

    // ---- dimensions --------------------------------------------------------
    q  = cols(b)
    nI = cols(ind)
    G  = cols(gv)
    if (gvar == "" | nI == 0 | G < 2) {
        _error("dwlsmv: internal error, group metadata missing from e()")
    }
    if (rows(b) != 1 | rows(V) != q | cols(V) != q | rows(cs) != q | rows(rs) != q) {
        _error("dwlsmv: internal error, e(V) does not fit e(b)")
    }
    if (cs != rs) {
        _error("dwlsmv: internal error, e(V) stripe differs from e(b)")
    }
    if (rows(K) != G | cols(K) != nI) {
        _error("dwlsmv: internal error, e(K) does not fit the groups and indicators")
    }
    if (missing(b) | missing(V) | missing(K)) {
        _error("dwlsmv: internal error, missing value in e(b), e(V) or e(K)")
    }
    for (i = 1; i <= nI; i++) {
        if (dwlsmv_sg_at(ind[i], ind) == 0) {
            _error("dwlsmv: internal error, indicator named twice")
        }
    }
    tags = J(1, G, "")
    for (gr = 1; gr <= G; gr++) {
        if (gv[gr] == "" | dwlsmv_sg_at(gv[gr], gv) == 0) {
            _error("dwlsmv: internal error, group values not unique")
        }
        tags[gr] = gv[gr] + (gr == 1 ? "bn." : ".") + gvar
    }
    // a tolerance for a guard, not a comparison: e(V) is symmetrized on posting
    if (mreldif(V, V') > 1e-8) {
        _error("dwlsmv: internal error, e(V) is not symmetric")
    }

    // ---- factors, from the psi_ff rows --------------------------------------
    fac = J(1, 0, "")
    for (i = 1; i <= q; i++) {
        if (cs[i, 1] == "/" & substr(cs[i, 2], 1, 4) == "var(") {
            gr = dwlsmv_sg_tag(cs[i, 2], tags, head)
            if (gr == 1) {
                x = substr(head, 5, strlen(head) - 5)
                fac = fac, x
            }
        }
    }
    nf = cols(fac)
    if (nf == 0) {
        _error("dwlsmv: internal error, no factor variance row in e(b)")
    }
    for (f = 1; f <= nf; f++) {
        if (dwlsmv_sg_at(fac[f], fac) == 0 | anyof(ind, fac[f])) {
            _error("dwlsmv: internal error, ambiguous factor name in e(b)")
        }
    }

    // ---- one group and one kind per row ------------------------------------
    tab = J(q, 4, 0)
    for (i = 1; i <= q; i++) {
        eq = cs[i, 1]
        nm = cs[i, 2]
        gr = dwlsmv_sg_tag(nm, tags, head)
        if (gr == 0) {
            _error("dwlsmv: internal error, e(b) row has no group")
        }
        if (eq == "/") {
            if (substr(head, 1, 4) == "var(" | substr(head, 1, 4) == "cov(") {
                x = substr(head, 5, strlen(head) - 5)
                y = x
                if (substr(head, 1, 4) == "cov(") {
                    a = strpos(x, ",")
                    if (a == 0) {
                        _error("dwlsmv: internal error, e(b) cov row not recognized")
                    }
                    y = substr(x, a + 1, .)
                    x = substr(x, 1, a - 1)
                }
                r  = dwlsmv_sg_at(x, fac)
                cc = dwlsmv_sg_at(y, fac)
                if (r == 0 | cc == 0) {
                    _error("dwlsmv: internal error, e(b) psi row not recognized")
                }
                tab[i, .] = (3, gr, max((r, cc)), min((r, cc)))
            }
            else if (substr(head, 1, 6) == "scale(") {
                r = dwlsmv_sg_at(substr(head, 7, strlen(head) - 7), ind)
                if (r == 0) {
                    _error("dwlsmv: internal error, e(b) omega row not recognized")
                }
                tab[i, .] = (5, gr, r, 1)
            }
            else if (substr(head, 1, 5) == "mean(") {
                r = dwlsmv_sg_at(substr(head, 6, strlen(head) - 6), fac)
                if (r == 0) {
                    _error("dwlsmv: internal error, e(b) alpha row not recognized")
                }
                tab[i, .] = (6, gr, r, 1)
            }
            else {
                _error("dwlsmv: internal error, e(b) row not recognized")
            }
        }
        else if (substr(eq, 1, 1) == "/") {
            r = dwlsmv_sg_at(substr(eq, 2, .), ind)
            if (r == 0 | substr(head, 1, 5) != "c.cut") {
                _error("dwlsmv: internal error, e(b) threshold row not recognized")
            }
            cc = strtoreal(substr(head, 6, .))
            if (cc >= . | cc != floor(cc) | cc < 1) {
                _error("dwlsmv: internal error, e(b) threshold row not recognized")
            }
            tab[i, .] = (1, gr, r, cc)
        }
        else if (anyof(fac, eq)) {
            _error("dwlsmv: internal error, a structural path reached standardization")
        }
        else {
            r = dwlsmv_sg_at(eq, ind)
            if (r == 0) {
                _error("dwlsmv: internal error, e(b) row not recognized")
            }
            if (head == "") {
                tab[i, .] = (4, gr, r, 1)
            }
            else {
                cc = (substr(head, 1, 2) == "c." ? dwlsmv_sg_at(substr(head, 3, .), fac) : 0)
                if (cc == 0) {
                    _error("dwlsmv: internal error, e(b) loading row not recognized")
                }
                tab[i, .] = (2, gr, r, cc)
            }
        }
        if (i > 1) {
            if (dwlsmv_sg_find(tab[1..i - 1, .], tab[i, 1], tab[i, 2], tab[i, 3], tab[i, 4])) {
                _error("dwlsmv: internal error, e(b) has a row twice")
            }
        }
    }

    // ---- the rows make the complete table -----------------------------------
    m = 0
    for (gr = 1; gr <= G; gr++) {
        for (j = 1; j <= nI; j++) {
            if (K[gr, j] != floor(K[gr, j]) | K[gr, j] < 2) {
                _error("dwlsmv: internal error, e(K) entry is not a category count")
            }
            m = m + K[gr, j] - 1 + 2
            for (c = 1; c <= K[gr, j] - 1; c++) {
                if (!dwlsmv_sg_find(tab, 1, gr, j, c)) {
                    _error("dwlsmv: internal error, e(b) lacks a threshold of e(K)")
                }
            }
            if (!dwlsmv_sg_find(tab, 4, gr, j, 1) | !dwlsmv_sg_find(tab, 5, gr, j, 1)) {
                _error("dwlsmv: internal error, e(b) lacks a nu or omega row")
            }
        }
        for (f = 1; f <= nf; f++) {
            if (!dwlsmv_sg_find(tab, 3, gr, f, f) | !dwlsmv_sg_find(tab, 6, gr, f, 1)) {
                _error("dwlsmv: internal error, e(b) lacks a psi or alpha row")
            }
        }
    }
    if (sum(tab[., 1] :== 1) + sum(tab[., 1] :== 4) + sum(tab[., 1] :== 5) != m) {
        _error("dwlsmv: internal error, e(b) has rows e(K) does not explain")
    }

    // ---- e(V)'s fixed rows, the ties and e(rank) ----------------------------
    nfree = 0
    for (i = 1; i <= q; i++) {
        if (V[i, i] < 0) {
            _error("dwlsmv: internal error, negative diagonal in e(V)")
        }
        if (V[i, i] == 0) {
            if (any(V[i, .]) | any(V[., i])) {
                _error("dwlsmv: internal error, a fixed row of e(V) is not zero")
            }
        }
        else {
            nfree++
        }
    }
    nextra = 0
    seen   = J(1, q, 0)
    ln = tokens(canon, char(10))'
    lines = select(ln, ln :!= char(10))'
    for (i = 1; i <= cols(lines); i++) {
        w = tokens(lines[i])
        if (cols(w) == 0) {
            continue
        }
        if (w[1] != "tie") {
            continue
        }
        if (mod(cols(w) - 1, 4) | cols(w) < 9) {
            _error("dwlsmv: internal error, e(canon) tie line not recognized")
        }
        first = 0
        for (a = 2; a <= cols(w); a = a + 4) {
            gr = strtoreal(w[a + 1])
            if (gr >= . | gr < 1 | gr > G) {
                _error("dwlsmv: internal error, e(canon) tie names no group")
            }
            pos = 0
            if (w[a] == "tau") {
                r = dwlsmv_sg_at(w[a + 2], ind)
                pos = (r ? dwlsmv_sg_find(tab, 1, gr, r, strtoreal(w[a + 3])) : 0)
            }
            else if (w[a] == "lambda") {
                f = dwlsmv_sg_at(w[a + 2], fac)
                r = dwlsmv_sg_at(w[a + 3], ind)
                if (f & r) {
                    pos = dwlsmv_sg_find(tab, 2, gr, r, f)
                }
            }
            if (pos == 0) {
                _error("dwlsmv: internal error, e(canon) tie names no e(b) row")
            }
            if (seen[pos]) {
                _error("dwlsmv: internal error, an e(b) row is in two ties")
            }
            seen[pos] = 1
            if (first == 0) {
                first = pos
            }
            else {
                if (b[pos] != b[first] | V[pos, .] != V[first, .]) {
                    _error("dwlsmv: internal error, tied rows differ in e(b) or e(V)")
                }
                nextra++
            }
        }
    }
    if (nfree - nextra != rank) {
        _error("dwlsmv: internal error, e(rank) does not fit e(V) and the ties")
    }

    // ---- each group's table --------------------------------------------------
    Lp = J(1, G, NULL)
    Pp = J(1, G, NULL)
    omp = J(1, G, NULL)
    alp = J(1, G, NULL)
    for (gr = 1; gr <= G; gr++) {
        Lk = J(nI, nf, 0)
        Pk = J(nf, nf, 0)
        om = J(nI, 1, 1)
        al = J(nf, 1, 0)
        for (i = 1; i <= q; i++) {
            if (tab[i, 2] != gr) {
                continue
            }
            r  = tab[i, 3]
            cc = tab[i, 4]
            kd = tab[i, 1]
            if (kd == 2) {
                Lk[r, cc] = b[i]
            }
            else if (kd == 3) {
                Pk[r, cc] = b[i]
                Pk[cc, r] = b[i]
            }
            else if (kd == 5) {
                om[r] = b[i]
            }
            else if (kd == 6) {
                al[r] = b[i]
            }
        }
        Lp[gr]  = &(Lk * 1)
        Pp[gr]  = &(Pk * 1)
        omp[gr] = &(om * 1)
        alp[gr] = &(al * 1)
        if (any(diagonal(Pk) :<= 0)) {
            errprintf("dwlsmv: a factor variance is at or below 0 in group %f (%s = %s), so the\n",
                      gr, gvar, gv[gr])
            errprintf("standardized solution is not defined for this fit.\n")
            _error("dwlsmv: standardized needs positive factor variances")
        }
    }

    // ---- the map (SG2) -------------------------------------------------------
    m  = q + G * nI
    g  = J(1, m, .)
    Jg = J(m, q, 0)
    for (i = 1; i <= q; i++) {
        kd = tab[i, 1]
        gr = tab[i, 2]
        r  = tab[i, 3]
        cc = tab[i, 4]
        s  = sqrt(diagonal(*Pp[gr]))
        if (kd == 1 | kd == 4) {
            g[i] = (*omp[gr])[r] * b[i]
        }
        else if (kd == 2) {
            g[i] = (*omp[gr])[r] * b[i] * s[cc]
        }
        else if (kd == 3) {
            g[i] = (r == cc ? 1 : b[i] / (s[r] * s[cc]))
        }
        else if (kd == 5) {
            g[i] = 1
        }
        else {
            g[i] = b[i] / s[r]
        }
    }
    for (gr = 1; gr <= G; gr++) {
        L = *Lp[gr]
        P = *Pp[gr]
        for (j = 1; j <= nI; j++) {
            g[q + (gr - 1) * nI + j] = 1 - (*omp[gr])[j]^2 * (L[j, .] * P * L[j, .]')
        }
    }

    // ---- J (SG3), one column per e(b) row, group by group ---------------------
    for (c = 1; c <= q; c++) {
        gc = tab[c, 2]
        L  = *Lp[gc]
        P  = *Pp[gc]
        om = *omp[gc]
        al = *alp[gc]
        s  = sqrt(diagonal(P))
        dL = J(nI, nf, 0)
        dP = J(nf, nf, 0)
        dom = J(nI, 1, 0)
        kd = tab[c, 1]
        r  = tab[c, 3]
        cc = tab[c, 4]
        if (kd == 2) {
            dL[r, cc] = 1
        }
        else if (kd == 3) {
            dP[r, cc] = 1
            dP[cc, r] = 1
        }
        else if (kd == 5) {
            dom[r] = 1
        }
        ds = diagonal(dP) :/ (2 :* s)
        for (i = 1; i <= q; i++) {
            gi = tab[i, 2]
            if (gi != gc) {
                continue
            }
            j = tab[i, 3]
            h = tab[i, 4]
            if (tab[i, 1] == 1 | tab[i, 1] == 4) {
                Jg[i, c] = b[i] * dom[j] + om[j] * (i == c)
            }
            else if (tab[i, 1] == 2) {
                Jg[i, c] = s[h] * (om[j] * dL[j, h] + L[j, h] * dom[j]) + om[j] * L[j, h] * ds[h]
            }
            else if (tab[i, 1] == 3) {
                if (j != h) {
                    Jg[i, c] = dP[j, h] / (s[j] * s[h]) -
                               P[j, h] * (s[h] * ds[j] + s[j] * ds[h]) / (s[j] * s[h])^2
                }
            }
            else if (tab[i, 1] == 6) {
                Jg[i, c] = (i == c) / s[j] - al[j] * ds[j] / s[j]^2
            }
        }
        for (j = 1; j <= nI; j++) {
            Vjj = L[j, .] * P * L[j, .]'
            Jg[q + (gc - 1) * nI + j, c] = -2 * om[j] * Vjj * dom[j] -
                om[j]^2 * (2 * (L[j, .] * P) * dL[j, .]' + L[j, .] * dP * L[j, .]')
        }
    }
    Vs = Jg * V * Jg'
    Vs = (Vs + Vs') / 2

    out = cs \ J(G * nI, 2, "/")
    for (gr = 1; gr <= G; gr++) {
        for (j = 1; j <= nI; j++) {
            out[q + (gr - 1) * nI + j, 2] = "var(e." + ind[j] + ")#" + tags[gr]
        }
    }
}

end
