*! Executed under StataMP 18 on 2026-10-14: passes (first run).
*!
*! dwlsmv_dtpair (gate step 3c; E-027, 2026-10-14, option (a)): the
*! statistic is called only at dwlsmv_dtstep status 1. Inputs as in
*! test_mg_difftest.do: one invariance(configural) fit per design posts
*! e(canon_configural) and e(K); each level fitted inside Mata from
*! dwlsmv_levelgen's text over one shared Stage 1; d and exempt from
*! dwlsmv_dtnest with the generated inputs, as test_dtstep.do.
*!
*!   p1  E thresholds vs configural, the fits' own T: status 1, and
*!       (T_d3, p, a_d, b_d) bitwise those of dwlsmv_difftest called directly
*!       on (Delta1, Delta0, Omega, T0 - T1)
*!   p2  E metric vs thresholds, T0 = T1 + 1e-12 (inside the band): status 2,
*!       the four missing
*!   p3  B thresholds vs configural (d = 0), the fits' own T: status 0, rc 0,
*!       the four missing, note dwlsmv_dtstep's (non-empty, passed through)
*!   p4  p1 with metric's Jacobian as Delta0: rc 3498 (not this pair's)
*!   p5  p1 with T0 = T1 - 1 (check 10's abort passes through): rc 3498
*!
*! Mutants (2026-10-14, scratch worktree, diff printed before each run,
*! %TEMP%\exgate\mutate3c.sh), predicted | observed, all red:
*!   M1  the statistic called at d = 0 (owed from step 2)   p3 | p3 (rc 3498)
*!   M2  the statistic called at status 2                   p2 | p2
*!   M6  the pair swapped in the call                       p1 | p1 (rc 3498)
*!   M8  the d guard removed                                p4 | p4
*! Not run here: Td in lavaan's convention and Omega without pi are inputs
*! the caller builds (the front end, not yet written); test_mg_difftest.do
*! checks 4 and 6 are their negative controls. Delta0 at theta1 has no
*! meaning (H0's columns are not H1's), so no mutant was written for it.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

local nfail 0

mata:

struct tdp_lev {
    struct dwlsmv_parmap scalar par
    real colvector theta
    real matrix    Delta
    real scalar    T, df
}

void tdp_keep(string scalar X)
{
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_statmap vector maps
    struct dwlsmv_fit scalar fit
    struct tdp_lev vector L
    string rowvector lev
    real matrix K, Y
    real colvector gv
    real scalar G, g, k, D
    pointer() scalar p

    K    = st_matrix("e(K)")
    G    = rows(K)
    gv   = uniqrows(st_data(., "g"))
    sts  = dwlsmv_stage1(G, 1)
    maps = dwlsmv_statmap(G, 1)
    D    = 0
    for (g = 1; g <= G; g++) {
        Y       = select(st_data(., ("y1", "y2", "y3", "y4", "y5")), st_data(., "g") :== gv[g])
        sts[g]  = dwls_sample_stats(Y, K[g, .]')
        maps[g] = sts[g].map
        D = D + sts[g].map.D
    }
    lev = ("configural", "thresholds", "metric", "scalar")
    L   = tdp_lev(4, 1)
    for (k = 1; k <= 4; k++) {
        L[k].par   = dwls_identify(dwlsmv_levelgen(st_global("e(canon_configural)"), K, lev[k]), maps, sts)
        fit        = dwls_fit(sts, L[k].par)
        L[k].theta = fit.theta
        L[k].Delta = fit.Delta
        L[k].T     = fit.T_raw
        L[k].df    = D - L[k].par.npar
    }
    p = crexternal("TDP_L_" + X)
    if (p == NULL) p = findexternal("TDP_L_" + X)
    *p = L
    p = crexternal("TDP_S_" + X)
    if (p == NULL) p = findexternal("TDP_S_" + X)
    *p = sts
    p = crexternal("TDP_cc_" + X)
    if (p == NULL) p = findexternal("TDP_cc_" + X)
    *p = st_global("e(canon_configural)")
    p = crexternal("TDP_K_" + X)
    if (p == NULL) p = findexternal("TDP_K_" + X)
    *p = K
}

real scalar tdp_k(string scalar lev)
{
    if (lev == "configural") return(1)
    if (lev == "thresholds") return(2)
    if (lev == "metric")     return(3)
    return(4)
}

// One pair H0 = lev0, H1 = lev1 of design X through dwlsmv_dtnest, then
// dwlsmv_dtpair. mod: "band" T0 = T1 + 1e-12; "neg" T0 = T1 - 1; "wrong"
// Delta0 from the metric level. Sets locals st, t3, miss, same, note.
void tdp_pair(string scalar X, string scalar lev0, string scalar lev1, string scalar mod)
{
    struct tdp_lev vector L
    struct dwlsmv_stage1 vector sts
    string scalar cc, can0, can1, nt
    real matrix K, Om, D0
    real rowvector b, r
    real scalar d, ex, k0, k1, T0, T1, N, G

    L   = *findexternal("TDP_L_" + X)
    sts = *findexternal("TDP_S_" + X)
    cc  = *findexternal("TDP_cc_" + X)
    K   = *findexternal("TDP_K_" + X)
    k0  = tdp_k(lev0)
    k1  = tdp_k(lev1)
    can0 = dwlsmv_levelgen(cc, K, lev0)
    can1 = dwlsmv_levelgen(cc, K, lev1)
    ex  = 0
    d   = dwlsmv_dtnest(L[k0].par, L[k1].par, L[k0].df, L[k1].df, X + ":" + lev0, X + ":" + lev1,
                        cc, can0, lev0, K, cc, can1, lev1, K, ex)
    G   = rows(sts)
    N   = sts[1].N + sts[2].N
    Om  = dwlsmv_omega(sts)
    T1  = L[k1].T
    T0  = L[k0].T
    if (mod == "band") T0 = T1 + 1e-12
    if (mod == "neg")  T0 = T1 - 1
    D0  = L[k0].Delta
    if (mod == "wrong") D0 = L[3].Delta
    nt  = ""
    b   = dwlsmv_dtpair(L[k0].par, L[k0].theta, d, ex, T0, T1, N - G, 1e-12, L[k0].df,
                        cc, K, lev0, lev1, X + ":" + lev0, X + ":" + lev1,
                        L[k1].Delta, D0, Om, nt)
    st_local("st", strofreal(b[1]))
    st_local("miss", strofreal(missing(b[5..8]) == 4))
    st_local("note", strofreal(strlen(nt) > 0))
    st_local("same", "0")
    if (b[1] == 1) {
        r = dwlsmv_difftest(L[k1].Delta, L[k0].Delta, Om, T0 - T1)
        st_local("same", strofreal(b[5..8] == (r[1], r[3], r[4], r[5])))
        printf("{txt}    T_d3 %14.10f p %12.10f a %12.10f b %12.10f\n", b[5], b[6], b[7], b[8])
    }
}
end

foreach X in E B {
    import delimited using "fixtures/data/mg_`X'.csv", clear varnames(1) delimiters(",")
    quietly dwlsmv (F -> y1 y2 y3 y4 y5), group(g) invariance(configural) metric(unitvariance)
    mata: tdp_keep("`X'")
    display as text "    design `X': four levels fitted"
}

* p1
capture noisily mata: tdp_pair("E", "thresholds", "configural", "")
if _rc | "`st'" != "1" | "`same'" != "1" {
    display as error "    p1 FAILED: rc `=_rc' status `st' same `same'"
    local ++nfail
}
else display as text "    p1 ok"

* p2
capture noisily mata: tdp_pair("E", "metric", "thresholds", "band")
if _rc | "`st'" != "2" | "`miss'" != "1" {
    display as error "    p2 FAILED: rc `=_rc' status `st' missing `miss'"
    local ++nfail
}
else display as text "    p2 ok"

* p3
capture noisily mata: tdp_pair("B", "thresholds", "configural", "")
if _rc | "`st'" != "0" | "`miss'" != "1" | "`note'" != "1" {
    display as error "    p3 FAILED: rc `=_rc' status `st' missing `miss' note `note'"
    local ++nfail
}
else display as text "    p3 ok"

* p4
capture noisily mata: tdp_pair("E", "thresholds", "configural", "wrong")
if _rc != 3498 {
    display as error "    p4 FAILED: rc `=_rc', expected 3498"
    local ++nfail
}
else display as text "    p4 ok"

* p5
capture noisily mata: tdp_pair("E", "thresholds", "configural", "neg")
if _rc != 3498 {
    display as error "    p5 FAILED: rc `=_rc', expected 3498"
    local ++nfail
}
else display as text "    p5 ok"

if `nfail' {
    display as error "test_dtpair: `nfail' failure(s)"
    exit 9
}
display as result "test_dtpair: the statistic runs only at dtstep status 1 and equals dwlsmv_difftest's there; status 0 and 2 leave it missing; the note and dtstep's aborts pass through; another pair's Jacobian is refused"
