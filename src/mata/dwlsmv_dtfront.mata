*! dwlsmv -- difftest's Mata side: refit, checks 6-11, the statistic, output
*! (Module 12, front end)
*!
*! Written 2026-10-23 (user: build difftest for single-group comparisons;
*! design E-020, U1-U8). difftest.ado
*! runs checks 1-5 on the stored results, then:
*!
*!   dwlsmv_dtfront_init()          empties the model list
*!   dwlsmv_dtfront_add(name)       appends the ACTIVE e() results: e(b), e(V),
*!                                  e(bound), e(T3), e(df_T3), e(converged),
*!                                  e(canon), e(metric), e(recovery)
*!   dwlsmv_dtfront(touse, ind, kb) the rest, below
*!
*! dwlsmv_dtfront: one shared Stage 1 (dwlsmv_m11_stage1 on the common
*! estimation sample), then for each model the refit through
*! dwlsmv_m11_estimate (the code dwlsmv itself runs) with full inference;
*! check 7 on the stored fit; check 6 through dwlsmv_dt_check6 (its only
*! entry), whose refusal gains the Stata version/edition line (decided
*! 2026-10-23); check 7 on the refit (converged, same e(recovery)); checks 8-9
*! (dwlsmv_dtchain); per adjacent pair checks 10-11 and the statistic
*! (dwlsmv_dtpair, c = N - 1, nrtol 1e-12 as dwlsmv_minimize sets it); the
*! fit indices of each refit (dwlsmv_fitidx) and their differences
*! (dwlsmv_dtfront_delta: restricted minus less restricted, F6). Prints the models table and the steps table
*! (U5 (b)), status 2 as "within eps; not tested" (U6 (a)), the routines'
*! notes verbatim, and posts r(models), r(steps) (U7 (a)) and r(N).
*! Every refusal is rc 3498.
*!
*! G > 1 (admitted 2026-10-31, user; narrow claim, E-021):
*! gvar and gvals name the group variable and its values in the sample; Stage
*! 1 is one per group (dwlsmv_m11_gsupport, dwlsmv_m11_stage1_g); each model
*! is refitted from e(canon) by dwls_identify, dwls_fit and dwls_inference
*! (dwlsmv_m11_fit_g's path), e(b) and e(V) by dwlsmv_m11_gbv (the posting
*! code); e(bound) is not posted there, so check 6 compares 0 x 0; check 8
*! gets e(canon_configural), e(canon), e(invariance) and e(K) for the
*! generated-chain exemption; c = N - G, Omega and the fit indices over G.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

struct dwlsmv_dtm {
    string scalar  name, canon, metric, rec, inv, cc
    real rowvector b
    real matrix    V, bound, K
    real scalar    T3, df, conv
}

void dwlsmv_dtfront_init()
{
    pointer() scalar p

    p = crexternal("dwlsmv_dtfront_M")
    if (p == NULL) {
        p = findexternal("dwlsmv_dtfront_M")
    }
    *p = dwlsmv_dtm(0, 1)
}

void dwlsmv_dtfront_add(string scalar name)
{
    struct dwlsmv_dtm scalar m
    struct dwlsmv_dtm colvector M, Mo
    pointer() scalar p
    real scalar k

    m.name   = name
    m.b      = st_matrix("e(b)")
    m.V      = st_matrix("e(V)")
    m.bound  = st_matrix("e(bound)")
    m.T3     = st_numscalar("e(T3)")
    m.df     = st_numscalar("e(df_T3)")
    m.conv   = st_numscalar("e(converged)")
    m.canon  = st_global("e(canon)")
    m.metric = st_global("e(metric)")
    m.rec    = st_global("e(recovery)")
    m.inv    = st_global("e(invariance)")
    m.cc     = st_global("e(canon_configural)")
    m.K      = st_matrix("e(K)")
    p  = findexternal("dwlsmv_dtfront_M")
    Mo = *p
    M  = dwlsmv_dtm(rows(Mo) + 1, 1)
    for (k = 1; k <= rows(Mo); k++) {
        M[k] = Mo[k]
    }
    M[rows(M)] = m
    *p = M
}

// e(b) and e(V) as dwlsmv_m11_fit posts them (D5)
void dwlsmv_dtfront_post(struct dwlsmv_parmap scalar par, struct dwlsmv_fit scalar fit,
                         real scalar N, real rowvector b, real matrix V)
{
    real matrix Ev
    real scalar i, k

    b = J(1, par.nall, .)
    for (i = 1; i <= par.nall; i++) {
        b[i] = (par.free[i] ? fit.theta[par.fidx[i]] : par.value[i])
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
}

// check 6 for model k, run under _stata() so its rc can be caught
void dwlsmv_dtfront_c6(real scalar k)
{
    struct dwlsmv_dtm colvector M, Mo
    pointer() matrix R
    struct dwlsmv_dtm scalar m

    M = *findexternal("dwlsmv_dtfront_M")
    R = *findexternal("dwlsmv_dtfront_R")
    m = M[k]
    dwlsmv_dt_check6(m.b, m.V, m.T3, m.bound, m.df,
                     *R[1, k], *R[2, k], *R[3, k], *R[4, k], *R[5, k], m.name)
}

void dwlsmv_dtfront(string scalar touse, string scalar indlist, real scalar kb,
                    string scalar gvar, string scalar gvals)
{
    struct dwlsmv_dtm colvector M, Mo
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_statmap vector maps
    real colvector bv
    real matrix P
    real scalar G, g, nfac
    string colvector cl, ccs, cans, invs
    struct dwlsmv_parmap vector pars
    struct dwlsmv_fit vector fits
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar fit
    pointer() matrix R
    pointer() scalar p
    pointer(real matrix) colvector Kp
    string colvector names, notes, blank
    string scalar rec, rnote, note, lab, cell
    real colvector df, ord
    real rowvector b1, base, fi, r
    real matrix V1, A, Bd, info, FI, steps, models, Om
    real scalar match, n, k, a, c, N, rc, w

    M = *findexternal("dwlsmv_dtfront_M")
    n = rows(M)
    if (n < 2) {
        _error("difftest: need at least two models")
    }

    // ---- one shared Stage 1 ------------------------------------------------
    // (G > 1: one per group, as dwlsmv's invariance() flow builds it)
    if (gvar == "") {
        dwlsmv_m11_stage1(touse, indlist, kb)
        st  = *findexternal("dwlsmv_m11_st")
        sts = dwlsmv_stage1(1, 1)
        sts[1] = st
    }
    else {
        dwlsmv_m11_gsupport(touse, indlist, gvar, gvals)
        dwlsmv_m11_stage1_g(kb)
        sts = *findexternal("dwlsmv_m11_st")
    }
    G    = rows(sts)
    maps = dwlsmv_statmap(G, 1)
    N    = 0
    for (g = 1; g <= G; g++) {
        maps[g] = sts[g].map
        N       = N + sts[g].N
    }

    // ---- e(K) against this Stage 1 (spec 4.8 G8) ----------------------------
    if (G > 1) {
        for (k = 1; k <= n; k++) {
            if (M[k].inv != "") {
                match = (rows(M[k].K) == G)
                for (g = 1; match & g <= G; g++) {
                    match = (cols(M[k].K) == rows(sts[g].map.K))
                    if (match) {
                        match = (M[k].K[g, .]' == sts[g].map.K)
                    }
                }
                if (!match) {
                    errprintf("difftest (%s): e(K) differs from the category counts of the refit's own Stage 1.\n", M[k].name)
                    exit(3498)
                }
            }
        }
    }

    // ---- check 7 on the stored fits ------------------------------------------
    for (k = 1; k <= n; k++) {
        if (M[k].conv != 1) {
            errprintf("difftest check 7 (%s): the stored fit did not converge; no comparison is reported.\n", M[k].name)
            exit(3498)
        }
    }

    // ---- the refits, with full inference ---------------------------------------
    pars = dwlsmv_parmap(n, 1)
    fits = dwlsmv_fit(n, 1)
    R    = J(5, n, NULL)
    for (k = 1; k <= n; k++) {
        if (G == 1) {
            dwlsmv_m11_estimate(M[k].canon, M[k].metric, st, par, fit, A, Bd, rec, rnote)
            dwlsmv_dtfront_post(par, fit, N, b1, V1)
        }
        else {
            // the invariance() fit's own path (dwlsmv_m11_fit_g): identify, fit,
            // inference; e(b) and e(V) by the posting code (dwlsmv_m11_gbv)
            par  = dwls_identify(subinstr(M[k].canon, ";", char(10)), maps, sts)
            fit  = dwls_fit(sts, par)
            dwls_inference(fit)
            rec  = "none"
            cl   = dwlsmv_canon_lines(subinstr(M[k].canon, ";", char(10)))
            nfac = 0
            for (g = 1; g <= rows(cl); g++) {
                if (cols(tokens(cl[g]))) {
                    if (tokens(cl[g])[1] == "factor") {
                        nfac = cols(tokens(cl[g])) - 1
                        break
                    }
                }
            }
            dwlsmv_m11_gbv(par, fit, N, G, cols(tokens(indlist)), nfac, P, bv, V1)
            b1 = bv'
        }
        pars[k] = par
        fits[k] = fit
        R[1, k] = &(b1 * 1)
        R[2, k] = &(V1 * 1)
        R[3, k] = &(fit.T3 * 1)
        // e(bound) is not posted at G > 1
        R[4, k] = (G == 1 ? &(st.bound') : &J(0, 0, .))
        R[5, k] = &(fit.T3_df * 1)
        if (fit.converged != 1) {
            errprintf("difftest check 7 (%s): the refit did not converge; no comparison is reported.\n", M[k].name)
            exit(3498)
        }
        if (rec != M[k].rec) {
            errprintf("difftest check 7 (%s): the stored fit's recovery was %s, the refit's %s.\n",
                      M[k].name, M[k].rec, rec)
            errprintf("The refit did not reproduce the recovery; no comparison is made.\n")
            exit(3498)
        }
    }
    p = crexternal("dwlsmv_dtfront_R")
    if (p == NULL) {
        p = findexternal("dwlsmv_dtfront_R")
    }
    *p = R

    // ---- check 6, through its one wrapper ----------------------------------------
    for (k = 1; k <= n; k++) {
        rc = _stata("mata: dwlsmv_dtfront_c6(" + strofreal(k) + ")")
        if (rc) {
            errprintf("If %s was fitted in a different version or edition of Stata, refit it in\n", M[k].name)
            errprintf("this one and store it again.\n")
            exit(3498)
        }
    }

    // ---- checks 8-9 ------------------------------------------------------------------
    names = J(n, 1, "")
    df    = J(n, 1, .)
    Kp    = J(n, 1, NULL)
    for (k = 1; k <= n; k++) {
        names[k] = M[k].name
        df[k]    = fits[k].D_stat
        Kp[k]    = &J(0, 0, .)
    }
    blank = J(n, 1, "")
    // at G > 1 each model's e(canon_configural), e(canon), e(invariance) and
    // e(K), for check 8's generated-chain exemption
    ccs  = blank
    cans = blank
    invs = blank
    if (G > 1) {
        for (k = 1; k <= n; k++) {
            ccs[k]  = M[k].cc
            cans[k] = M[k].canon
            invs[k] = M[k].inv
            Kp[k]   = &(M[k].K * 1)
        }
    }
    ord   = dwlsmv_dtchain(pars, df, names, ccs, cans, invs, Kp, info)

    // ---- checks 10-11 and the statistic, pair by pair ------------------------
    Om    = dwlsmv_omega(sts)
    c     = N - G
    steps = J(n - 1, 6, .)
    notes = J(n - 1, 1, "")
    for (k = 1; k < n; k++) {
        a = ord[k]
        w = ord[k + 1]
        note = ""
        r = dwlsmv_dtpair(pars[a], fits[a].theta, info[k, 1], info[k, 2],
                          fits[a].T_raw, fits[w].T_raw, c, 1e-12, fits[a].D_stat,
                          ccs[a], *Kp[a], invs[a], invs[w], names[a], names[w],
                          fits[w].Delta, fits[a].Delta, Om, note)
        steps[k, 1] = r[2]
        steps[k, 2] = r[5]
        steps[k, 3] = r[6]
        steps[k, 6] = r[1]
        notes[k]    = note
    }

    // ---- fit indices and their differences ----------------------------------------
    base   = dwlsmv_fitidx_base(sts)
    models = J(n, 6, .)
    for (k = 1; k <= n; k++) {
        a  = ord[k]
        fi = dwlsmv_fitidx(fits[a].T3, fits[a].T3_df, N - G, G, base[1], base[2])
        models[k, .] = (pars[a].npar, fits[a].T3, fits[a].T3_df, fi)
    }
    steps[., 4..5] = dwlsmv_dtfront_delta((models[., 5], models[., 4]), steps[., 6])

    // ---- output -----------------------------------------------------------------------
    printf("\n{txt}Difference test (DIFFTEST) between nested dwlsmv models\n")
    printf("{txt}Number of obs = {res}%g{txt}; models refitted over one shared Stage 1\n", N)
    printf("\n{txt}Models, most restricted first\n")
    printf("{hline 78}\n")
    printf("{txt}%-20s %6s %12s %6s %10s %10s %10s\n", "Model", "npar", "T3", "df",
           "RMSEA", "CFI", "TLI")
    printf("{hline 78}\n")
    for (k = 1; k <= n; k++) {
        printf("{txt}%-20s {res}%6.0f %12.3f %6.0f %10.4f %10.4f %10.4f\n",
               dwlsmv_dtfront_lab(names[ord[k]], 20), models[k, 1], models[k, 2], models[k, 3],
               models[k, 4], models[k, 5], models[k, 6])
    }
    printf("{hline 78}\n")
    printf("\n{txt}Steps (H0 against H1)\n")
    printf("{hline 78}\n")
    printf("{txt}%-30s %4s %12s %10s %8s %9s\n", "H0 -> H1", "d", "T_d3", "Prob>chi2",
           "dCFI", "dRMSEA")
    printf("{hline 78}\n")
    for (k = 1; k < n; k++) {
        lab = dwlsmv_dtfront_lab(names[ord[k]], 14) + " -> " + dwlsmv_dtfront_lab(names[ord[k + 1]], 12)
        if (steps[k, 6] == 1) {
            cell = sprintf("{res}%4.0f %12.3f %10.4f", steps[k, 1], steps[k, 2], steps[k, 3])
        }
        else if (steps[k, 6] == 2) {
            cell = sprintf("{res}%4.0f {txt}%23s", steps[k, 1], "within eps; not tested")
        }
        else {
            cell = sprintf("{res}%4.0f {txt}%23s", steps[k, 1], "skipped (d = 0)")
        }
        printf("{txt}%-30s %s {res}%8.4f %9.4f\n", lab, cell, steps[k, 4], steps[k, 5])
    }
    printf("{hline 78}\n")
    for (k = 1; k < n; k++) {
        if (notes[k] != "") {
            dwlsmv_dtfront_print(notes[k])
        }
    }
    printf("{txt}No published dCFI or dRMSEA cutoff has been validated for DWLS estimation with\n")
    printf("{txt}the T3 statistic; cutoffs developed for maximum likelihood should not be\n")
    printf("{txt}assumed to carry over.\n")

    // ---- r() --------------------------------------------------------------------------------
    st_matrix("r(models)", models)
    st_matrixcolstripe("r(models)", (J(6, 1, ""), ("npar" \ "T3" \ "df" \ "rmsea" \ "cfi" \ "tli")))
    // "." (the active results) is not a valid stripe name
    names = subinstr(names, ".", "_active")
    st_matrixrowstripe("r(models)", (J(n, 1, ""), names[ord]))
    st_matrix("r(steps)", steps)
    st_matrixcolstripe("r(steps)", (J(6, 1, ""), ("d" \ "T_d3" \ "p" \ "dCFI" \ "dRMSEA" \ "status")))
    st_matrixrowstripe("r(steps)", (names[ord[1..n - 1]], names[ord[2..n]]))
    st_numscalar("r(N)", N)
}

// one printf per line: one %s holding char(10)s lets SMCL wrap mid-word
void dwlsmv_dtfront_print(string scalar s)
{
    real scalar k

    while ((k = strpos(s, char(10))) > 0) {
        printf("{txt}%s\n", substr(s, 1, k - 1))
        s = substr(s, k + 1, .)
    }
    if (s != "") {
        printf("{txt}%s\n", s)
    }
}

end

mata:

// dCFI, dRMSEA per step, the RESTRICTED model's value minus the less
// restricted one's (H0 minus H1; FI3, user 2026-10-31), whatever its sign.
// X: (CFI, RMSEA), one row per model, most restricted first, as r(models);
// status: one per step (0 skipped, d = 0). dwlsmv_fitidx_delta takes the
// chain least restricted first (F6's order): a skipped step's restricted
// model is the skipped level, the next difference against the last fitted.
real matrix dwlsmv_dtfront_delta(real matrix X, real colvector status)
{
    real scalar n
    real colvector rv
    real matrix dX

    n = rows(X)
    if (n < 2 | cols(X) != 2 | rows(status) != n - 1) {
        _error("dwlsmv_dtfront_delta: X needs two columns and one status per step")
    }
    rv = (n..1)'
    dX = dwlsmv_fitidx_delta(X[rv, .], 0 \ (status[(n - 1..1)'] :== 0))
    return(dX[rv[1..n - 1], .])
}

// a model name for display; abbrev() returns "" for "."
string scalar dwlsmv_dtfront_lab(string scalar s, real scalar w)
{
    return(s == "." ? "." : abbrev(s, w))
}

end
