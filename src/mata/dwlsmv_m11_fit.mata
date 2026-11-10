*! dwlsmv -- Module 11: the estimation driver behind dwlsmv.ado (dwlsmv_m11_fit)
*!
*! Executed under StataMP 18 on 2026-10-04: tests/unit/test_dwlsmv.do.
*!
*! One file per function (test_loading.do's contract): dwlsmv_m11_stage1,
*! _legend, _fit, _slab, _notes, _clear. Two steps, called by dwlsmv.ado in
*! two separate mata: lines, because the two can fail differently (D8a): a Stage 1 abort is followed by the legend
*! and the pointer line (D8b, D8f) before rc 3498 is rethrown; an abort after
*! it is rethrown as it is.
*!
*!   dwlsmv_m11_stage1(touse, indlist, keepbound)
*!       D9's category recoding, then dwls_sample_stats(Y, K, keepbound).
*!       Each indicator's distinct observed values, sorted, become 0..K-1:
*!       recoding the observed support, never collapsing a declared category
*!       (nothing declares one in v1a; value labels are not read).
*!   dwlsmv_m11_legend(indlist)
*!       D8b: one line per indicator, yj = name, in e(oyvars) order.
*!   dwlsmv_m11_fit(canon, indlist, semnames, bn, vn, cn, boundn, metric)
*!       dwls_identify, dwls_fit, dwls_inference; under metric "marker",
*!       D11's recovery (dwlsmv_m11_recover()), whose result is the local
*!       recovery and whose note opens the notes block; e(b), e(V), e(Cns) and
*!       e(bound) into the Stata matrices named bn, vn, cn, boundn (D5, D6,
*!       D8d), scalars into locals, and the notes block (D8e) kept for
*!       dwlsmv_m11_notes(). Names are built from the parameter map and
*!       checked against sem's stripe in both directions (D6); a mismatch is
*!       an internal error.
*!   dwlsmv_m11_notes()
*!       prints the notes block: the keepbound footnote, the one-empty-cell repair, the near-bound note,
*!       admissibility, Module 10's notes and values.
*!
*! State between the steps lives in three Mata externals, removed by
*! dwlsmv_m11_clear(): the Stage 1 struct, the recoded Y, and the notes text.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_m11_fit(string scalar canon, string scalar indlist, string scalar semnames,
                    string scalar bn, string scalar vn, string scalar cn, string scalar boundn,
                    string scalar metric)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit    scalar fit
    string scalar rec, rnote
    string rowvector ind, fac, sem, built
    string colvector nm, slab, pname
    string matrix stp, bnm
    real colvector endo, bvec, rho, rep
    real matrix V, Ev, Cns, A, Gd, Cd, Bd, bnd
    real scalar i, k, f, t, nfix, P, T, N
    string scalar txt, nl, a, z, line
    pointer() scalar p

    nl  = char(10)
    st  = *findexternal("dwlsmv_m11_st")
    ind = tokens(indlist)
    fac = tokens(substr(tokens(canon, ";")[1], 8, .))

    // the fit and D11's recovery: dwlsmv_m11_estimate (shared with difftest's refit)
    dwlsmv_m11_estimate(canon, metric, st, par, fit, A, Bd, rec, rnote)

    // ---- D6: names from the parameter map ------------------------------
    // endogenous = has a beta row, or sem names its error term: a path fixed
    // at 0 is no table row (and no e(b) column, as after sem), yet sem still
    // reads the factor as endogenous and posts /var(e.F)
    endo = J(cols(fac), 1, 0)
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 3) {
            endo[par.row[i]] = 1
        }
    }
    for (i = 1; i <= cols(fac); i++) {
        if (anyof(tokens(semnames), "/var(e." + fac[i] + ")")) {
            endo[i] = 1
        }
    }
    // stripe = (equation, name); nm is the full name as colfullnames shows
    // it: y1:F1, F2:F1, /var(F1), /cov(F1,F2), /y1:cut1
    stp = J(par.nall, 2, "")
    t   = 0
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 4) {
            t++
            stp[i, .] = ("/" + ind[st.map.var_i[t]], "cut" + strofreal(st.map.thr_k[t]))
        }
        else if (par.ptype[i] == 1) {
            stp[i, .] = (ind[par.row[i]], fac[par.col[i]])
        }
        else if (par.ptype[i] == 3) {
            stp[i, .] = (fac[par.row[i]], fac[par.col[i]])
        }
        else {
            a = (endo[par.col[i]] ? "e." : "") + fac[par.col[i]]
            z = (endo[par.row[i]] ? "e." : "") + fac[par.row[i]]
            stp[i, .] = ("/", (par.row[i] == par.col[i] ? "var(" + z + ")" : "cov(" + a + "," + z + ")"))
        }
    }
    nm = J(par.nall, 1, "")
    for (i = 1; i <= par.nall; i++) {
        nm[i] = (stp[i, 1] == "/" ? "/" + stp[i, 2] : stp[i, 1] + ":" + stp[i, 2])
    }
    built = select(nm, par.ptype :!= 4)'
    sem   = tokens(semnames)
    for (i = 1; i <= cols(built); i++) {
        if (!anyof(sem, built[i])) {
            errprintf("dwlsmv: internal error: built name %s is not in sem's stripe.\n", built[i])
            errprintf("    built: %s\n    sem:   %s\n", invtokens(built), semnames)
            _error("dwlsmv: internal error, parameter names differ from sem's")
        }
    }
    for (i = 1; i <= cols(sem); i++) {
        if (!anyof(built, sem[i])) {
            errprintf("dwlsmv: internal error: sem's name %s was not built.\n", sem[i])
            errprintf("    built: %s\n    sem:   %s\n", invtokens(built), semnames)
            _error("dwlsmv: internal error, parameter names differ from sem's")
        }
    }

    // ---- D5: e(b), e(V), e(Cns) ----------------------------------------
    N    = st.N
    bvec = J(par.nall, 1, .)
    for (i = 1; i <= par.nall; i++) {
        bvec[i] = (par.free[i] ? fit.theta[par.fidx[i]] : par.value[i])
    }
    Ev = dwlsmv_ev(fit.acov, N)
    V  = J(par.nall, par.nall, 0)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            for (k = 1; k <= par.nall; k++) {
                if (par.free[k]) {
                    V[i, k] = Ev[par.fidx[i], par.fidx[k]]
                }
            }
        }
    }
    nfix = par.nall - par.npar
    Cns  = J(nfix, par.nall + 1, 0)
    k    = 0
    for (i = 1; i <= par.nall; i++) {
        if (!par.free[i]) {
            k++
            Cns[k, i]            = 1
            Cns[k, par.nall + 1] = par.value[i]
        }
    }
    st_matrix(bn, bvec')
    st_matrixcolstripe(bn, stp)
    st_matrix(vn, V)
    st_matrixcolstripe(vn, stp)
    st_matrixrowstripe(vn, stp)
    if (nfix) {
        st_matrix(cn, Cns)
        st_matrixcolstripe(cn, stp \ ("", "_r"))
    }

    // ---- D8d: e(bound), stripe varA:varB ---------------------------------
    T   = st.map.T
    P   = st.map.P
    bnd = st.bound'
    bnm = J(P, 2, "")
    for (k = 1; k <= P; k++) {
        bnm[k, 1] = ind[st.map.var_i[T + k]]
        bnm[k, 2] = ind[st.map.var_j[T + k]]
    }
    st_matrix(boundn, bnd)
    st_matrixcolstripe(boundn, bnm)

    st_local("N",          strofreal(N, "%21.0g"))
    st_local("T3",         strofreal(fit.T3, "%21.0g"))
    st_local("df",         strofreal(fit.T3_df, "%21.0g"))
    dwlsmv_m11_fitidx(st, fit.T3, fit.T3_df)
    st_local("converged",  strofreal(fit.converged))
    st_local("iterations", strofreal(fit.iterations))
    st_local("rank",       strofreal(par.npar))
    st_local("nfix",       strofreal(nfix))

    st_local("recovery",   rec)

    // ---- D8e: the notes block; D11's recovery line first -----------------
    txt = rnote
    if (any(st.bound)) {
        txt = txt + "Pairs included at the bound under keepbound:" + nl
        for (k = 1; k <= P; k++) {
            if (st.bound[k]) {
                txt = txt + "    " + bnm[k, 1] + "~~" + bnm[k, 2] + nl
            }
        }
    }
    // one-empty-cell pairs (O6, decided 2026-10-06): from Y, as Stage 1 and
    // Module 10 derive them (dwlsmv_emptycell_pairs())
    rep  = dwlsmv_emptycell_pairs(*findexternal("dwlsmv_m11_Y"), st.map)
    if (any(rep)) {
        txt = txt + "One-empty-cell repair applied to:" + nl
        for (k = 1; k <= P; k++) {
            if (rep[k]) {
                txt = txt + "    " + bnm[k, 1] + "~~" + bnm[k, 2] + "   rho = " +
                      strtrim(strofreal(st.s[T + k], "%9.4f")) + nl
            }
        }
        txt = txt + "0.5 was added to the empty cell and its diagonal partner and taken from" + nl +
              "the other two; margins and thresholds are unchanged. lavaan and Mplus apply" + nl +
              "the same repair and the same Gamma. Its effect on standard errors and T3 at" + nl +
              "small N is not established." + nl
    }
    rho  = st.s[(T + 1)::(T + P)]
    line = ""
    for (k = 1; k <= P; k++) {
        if (!st.bound[k]) {
            if (abs(rho[k]) >= 0.99) {
                line = line + "    " + bnm[k, 1] + "~~" + bnm[k, 2] + "   rho = " +
                       strtrim(strofreal(rho[k], "%9.4f")) + nl
            }
        }
    }
    if (line != "") {
        txt = txt + "Polychoric correlations at or beyond 0.99 in absolute value:" + nl + line +
              "A correlation this close to the bound has several possible readings: near-" + nl +
              "duplicate (redundant) items, strong skew collapsing the table, or too few" + nl +
              "observations. The data cannot tell them apart; this is not a diagnosis." + nl
    }
    txt = txt + dwlsmv_admissadvice(A, ind', fac')
    slab  = dwlsmv_m11_slab(st.map, ind)
    pname = select(nm, par.free)
    Gd  = dwlsmv_gammadiag(st, *findexternal("dwlsmv_m11_Y"))
    Cd  = dwlsmv_countdiag(*findexternal("dwlsmv_m11_Y"), st.map.K)
    txt = txt + dwlsmv_diagadvice(Gd, Cd, Bd, ind', slab, pname)
    // Module 10's singular-value ratios: stored, not printed (2026-10-07, the
    // user's). Missing when the row was not evaluated (flag `.`) or there is
    // no supported rank.
    st_local("gratio",  strofreal(Gd[7, 1] == . ? . : Gd[7, 2], "%21.0g"))
    st_local("gratior", strofreal(Gd[7, 1] == . ? . : Gd[6, 2], "%21.0g"))
    st_local("grank",   strofreal(Gd[7, 1] == . ? . : Gd[7, 3], "%21.0g"))
    st_local("bratio",  strofreal(Bd[4, 1] == . ? . : Bd[4, 2], "%21.0g"))

    p = crexternal("dwlsmv_m11_txt")
    if (p == NULL) {
        p = findexternal("dwlsmv_m11_txt")
    }
    *p = txt
}

end
