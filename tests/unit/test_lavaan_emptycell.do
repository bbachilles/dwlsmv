*! O6 against lavaan, end to end: data with a one-empty-cell pair
*! (fixtures/generate/90_lavaan_emptycell.R; lavaan 0.7.2 default
*! zero.add = c(0.5, 0), zero.keep.margins = TRUE, which is the repair and the
*! Gamma construction (a) that dwlsmv adopted, spec Module 2,
*! *One-empty-cell repair*). Every tolerance is spec 7.2's tier for that
*! quantity; none is set here.
*!
*!   1. s: thresholds and every correlation, repaired y1~~y2 included, at the
*!      7.2 tier 1e-8 abs + rel. lavaan's repaired rho differs from the exact
*!      keep-margins root by 4.5e-9 (the generator's self-check), inside it.
*!   2. Gamma at LAVAAN'S OWN rho-hat (the matched-rho tier, 1e-12 abs + rel):
*!      production's assembly with every pair's information at lavaan's rho
*!      and W = n/N, the raw table, for the repaired pair too. This is the
*!      construction question, free of convergence.
*!   3. Gamma end to end (the rho tier, 1e-8 abs + rel), production's own:
*!      PRINTED, NOT ASSERTED (ratio 7.9, from lavaan's 4.5e-9 rho gap; see the
*!      code), as test_gamma.do 5d for case E.
*!   4. theta-hat (1e-5 abs, 1e-4 rel), SEs (1e-8 abs, 1e-5 rel), T3 (1e-6
*!      abs, 1e-5 rel) and df (exact), marker identification, as lavaan.
*!   5. Negative control: Gamma at lavaan's rho but with the repaired pair's
*!      information at W = n+/N (the rejected hybrid) must FAIL check 2's
*!      tier, so check 2 can tell the constructions apart.
*!
*! Executed under StataMP 18 on 2026-10-06 and it passes.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
local load tests/dwlsmv_assert.mata src/mata/dwlsmv_types.mata
foreach f of local srcs {
    local load `load' src/mata/`f'
}
foreach f of local load {
    capture noisily do "`f'"
    if _rc {
        display as error "test_lavaan_emptycell: `f' failed to load"
        exit 9
    }
}

* reference fixtures: stringcols(_all) + strtoreal(), delimiters(",") always
import delimited using "fixtures/expected/lavaan_emptycell_gamma.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
qui ds label, not
local gcols `r(varlist)'
mata: G_LAV = strtoreal(st_sdata(., "`gcols'"))
mata: G_LAB = st_sdata(., "label")

import delimited using "fixtures/expected/lavaan_emptycell_sampstat.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: S_LAB = st_sdata(., "label")
mata: S_LAV = strtoreal(st_sdata(., "value"))

import delimited using "fixtures/expected/lavaan_emptycell_parameters.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: P_LAB = st_sdata(., "label")
mata: P_EST = strtoreal(st_sdata(., "est"))
mata: P_SE  = strtoreal(st_sdata(., "se"))

import delimited using "fixtures/expected/lavaan_emptycell_test.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: T_LAB = st_sdata(., "label")
mata: T_VAL = strtoreal(st_sdata(., "value"))

import delimited using "fixtures/data/emptycell.csv", clear varnames(1) delimiters(",")
mata: Y = st_data(., .)

mata:

real scalar tle_find(string colvector L, string scalar s)
{
    real scalar i

    for (i = 1; i <= rows(L); i++) {
        if (L[i] == s) return(i)
    }
    return(0)
}

// Gamma assembled by production's dwls_gamma() from pair fits whose
// information is taken at the given rhos; rep_plus = 1 puts the repaired
// pairs at W = n+/N (the hybrid, check 5's negative control).
real matrix tle_gamma(real matrix Y, struct dwlsmv_statmap scalar m,
                      real colvector s, real colvector rho, real scalar rep_plus)
{
    struct dwlsmv_pairfit colvector pairs
    struct dwlsmv_pairfit scalar f
    real colvector ti, tj, rep
    real matrix n, W
    real scalar k, pos, vi, vj, N

    N     = rows(Y)
    rep   = dwlsmv_emptycell_pairs(Y, m)
    pairs = dwlsmv_pairfit(m.P, 1)
    for (k = 1; k <= m.P; k++) {
        pos = m.T + k
        vi  = m.var_i[pos]
        vj  = m.var_j[pos]
        ti  = s[|m.Coff[vi] + 1 \ m.Coff[vi] + m.C[vi]|]
        tj  = s[|m.Coff[vj] + 1 \ m.Coff[vj] + m.C[vj]|]
        n   = dwlsmv_crosstab(Y[., vi], Y[., vj], m.K[vi], m.K[vj])
        W   = n :/ N
        if (rep_plus & rep[k]) {
            W = dwlsmv_emptycell_repair(n) :/ N
        }
        // Start from a full dwls_polychoric() fit, as test_gamma.do's
        // rebuild does: building on an empty element of dwlsmv_pairfit(P, 1)
        // made dwls_gamma() crash Stata 18 with a segmentation fault (found
        // 2026-10-06 writing this test; not chased further). Only rho and the
        // information below are used from it.
        f     = dwls_polychoric(Y[., vi], Y[., vj], ti, tj)
        f.rho = rho[k]
        dwlsmv_pairinfo(f, ti, tj, rho[k], W)
        f.i = vi
        f.j = vj
        pairs[k] = f
    }
    return(dwls_gamma(Y, m, pairs, dwlsmv_psi_tau_stack(Y, m)))
}

// The worst gap ratio between our Gamma and lavaan's, matched by label.
real scalar tle_gratio(real matrix G, struct dwlsmv_statmap scalar m,
                       real matrix GL, string colvector GLAB, real scalar tol)
{
    real scalar a, b, i, j, w

    w = 0
    for (a = 1; a <= m.D; a++) {
        i = tle_find(GLAB, m.label[a])
        if (i == 0) return(.)
        for (b = 1; b <= m.D; b++) {
            j = tle_find(GLAB, m.label[b])
            if (j == 0) return(.)
            w = max((w, dwlsmv_gapratio(G[a, b], GL[i, j], tol, tol)))
        }
    }
    return(w)
}

real scalar tle_main(real matrix Y,
                     real matrix G_LAV, string colvector G_LAB,
                     string colvector S_LAB, real colvector S_LAV,
                     string colvector P_LAB, real colvector P_EST, real colvector P_SE,
                     string colvector T_LAB, real colvector T_VAL)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit    scalar fit
    string colvector FAIL
    real colvector rhoL, sL, se
    real matrix G2, V
    real scalar k, i, r1, r2, r3, r5, rt, rs, re, tT3, j, N
    string scalar canon, lab

    FAIL = J(0, 1, "")
    N    = rows(Y)
    st   = dwls_sample_stats(Y, (2 \ 2 \ 2 \ 2))
    if (dwlsmv_emptycell_pairs(Y, st.map) != (1 \ 0 \ 0 \ 0 \ 0 \ 0)) {
        FAIL = FAIL \ "0: the fixture's one-empty-cell pair is not exactly y1~~y2"
    }

    // ---- 1. s ------------------------------------------------------------
    sL = J(st.map.D, 1, .)
    for (k = 1; k <= st.map.D; k++) {
        i = tle_find(S_LAB, st.map.label[k])
        if (i) sL[k] = S_LAV[i]
    }
    if (missing(sL)) {
        FAIL = FAIL \ "1: a label of s is missing from lavaan's wls.obs"
        sL = editmissing(sL, 0)
    }
    r1 = dwlsmv_gapratio(st.s, sL, 1e-8, 1e-8)
    printf("{txt}  1. s vs lavaan           : ratio %g (1e-8 abs and rel); repaired rho ours %17.15f, lavaan %17.15f\n",
           r1, st.s[st.map.T + 1], sL[st.map.T + 1])
    if (!(r1 <= 1)) FAIL = FAIL \ ("1: s ratio " + strofreal(r1))

    // ---- 2. Gamma at lavaan's rho ------------------------------------------
    rhoL = sL[(st.map.T + 1)::st.map.D]
    G2 = tle_gamma(Y, st.map, sL, rhoL, 0)
    r2 = tle_gratio(G2, st.map, G_LAV, G_LAB, 1e-12)
    printf("{txt}  2. Gamma at lavaan's rho : ratio %g (1e-12 abs and rel)\n", r2)
    if (!(r2 <= 1)) FAIL = FAIL \ ("2: Gamma at lavaan's rho, ratio " + strofreal(r2))

    // ---- 3. Gamma end to end -----------------------------------------------
    // PRINTED, NOT ASSERTED, as test_gamma.do's 5d (case E) for the same
    // cause: lavaan's repaired rho-hat sits 4.5e-9 from the exact keep-margins
    // root (the generator's self-check), ours 2.6e-15 from it
    // (test_emptycell.do 3), and Gamma carries that gap past the 1e-8 tier:
    // measured ratio 7.9 (2026-10-06). Check 2 shows it is the rho gap and
    // nothing else: at lavaan's own rho the two agree to 1e-12. The tier is
    // not loosened.
    r3 = tle_gratio(st.Gamma, st.map, G_LAV, G_LAB, 1e-8)
    printf("{txt}  3. Gamma end to end      : ratio %g (1e-8 abs and rel; printed, not asserted)\n", r3)

    // ---- 5. negative control: the hybrid must fail check 2's tier ------------
    r5 = tle_gratio(tle_gamma(Y, st.map, sL, rhoL, 1), st.map, G_LAV, G_LAB, 1e-12)
    printf("{txt}  5. hybrid at lavaan's rho: ratio %g (must exceed 1)\n", r5)
    if (!(r5 > 1)) FAIL = FAIL \ "5: the hybrid passes check 2's tier; check 2 has no teeth"

    // ---- 4. theta-hat, SEs, T3 ---------------------------------------------
    canon = "factor f" + char(10) + "ind y1 y2 y3 y4" + char(10) +
            "lambda f y2 free" + char(10) + "lambda f y3 free" + char(10) +
            "lambda f y4 free" + char(10) + "psi f f free" + char(10) +
            "metric marker f y1" + char(10)
    par = dwls_identify(canon, st.map)
    fit = dwls_fit(st, par)
    dwls_inference(fit)
    V  = dwlsmv_ev(fit.acov, N)
    se = sqrt(diagonal(V))
    rt = 0
    rs = 0
    for (i = 1; i <= par.nall; i++) {
        if (!par.free[i]) continue
        if (par.ptype[i] == 4)      lab = par.name[i]
        else if (par.ptype[i] == 1) lab = "f=~y" + strofreal(par.row[i])
        else                        lab = "f~~f"
        j = tle_find(P_LAB, lab)
        if (j == 0) {
            FAIL = FAIL \ ("4: no lavaan row " + lab)
            continue
        }
        rt = max((rt, dwlsmv_gapratio(fit.theta[par.fidx[i]], P_EST[j], 1e-5, 1e-4)))
        rs = max((rs, dwlsmv_gapratio(se[par.fidx[i]], P_SE[j], 1e-8, 1e-5)))
    }
    tT3 = T_VAL[tle_find(T_LAB, "chisq.scaled")]
    re  = dwlsmv_gapratio(fit.T3, tT3, 1e-6, 1e-5)
    printf("{txt}  4. theta %g (1e-5, 1e-4); SE %g (1e-8, 1e-5); T3 %g (1e-6, 1e-5), %17.12f vs %17.12f; df %g vs %g\n",
           rt, rs, re, fit.T3, tT3, fit.T3_df, T_VAL[tle_find(T_LAB, "df.scaled")])
    if (!(rt <= 1)) FAIL = FAIL \ ("4: theta ratio " + strofreal(rt))
    if (!(rs <= 1)) FAIL = FAIL \ ("4: SE ratio " + strofreal(rs))
    if (!(re <= 1)) FAIL = FAIL \ ("4: T3 ratio " + strofreal(re))
    if (fit.T3_df != T_VAL[tle_find(T_LAB, "df.scaled")]) FAIL = FAIL \ "4: df differs"

    for (i = 1; i <= rows(FAIL); i++) {
        printf("{err}%s\n", FAIL[i])
    }
    return(rows(FAIL))
}

end

mata: st_local("nfail", strofreal(tle_main(Y, G_LAV, G_LAB, S_LAB, S_LAV, P_LAB, P_EST, P_SE, T_LAB, T_VAL)))

if `nfail' {
    display as error _n "test_lavaan_emptycell: `nfail' failure(s)"
    exit 9
}
display as result "test_lavaan_emptycell: s, Gamma at lavaan's rho, theta, SEs and T3 agree with lavaan on a one-empty-cell pair; the hybrid does not"
exit 0
