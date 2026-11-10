*! Module 12: the invariance level generator, dwlsmv_levelgen() (spec 4.8 G8;
*! input contract L1 (b), level names L2 (a), decided 2026-10-12).
*!
*! Written and executed 2026-10-12 by the main session (StataMP 18): green on
*! the second run (the first failed e2 only: the log wrapped the message at 80
*! columns, so the search substring was moved within one line; no code change).
*!
*! Reference texts: design E (spec 4.8 G8's worked example; one factor F,
*! y1-y3 with 4 categories, y4-y5 binary, two groups), written out here by
*! tlg_ref(), a copy of test_fit_groups.do's tfg_text() with the level names
*! of L2 (metric = loadings, scalar = intercepts). That text is the one the
*! end-to-end test fits against lavaan, and it follows G8's emission order.
*!
*! Cases:
*!   r1-r4  generated text == reference, byte for byte, at each level
*!   r5     round trip: the configural output regenerates itself exactly, and
*!          equals the configural input when that input is already canonical
*!   d1     dwls_identify() on each generated text: free counts 32, 29, 27, 23
*!          (df 10, 13, 15, 19 with D = 42; spec G8, Result I4)
*!   q1     Q-5: a fixed nonzero cross-loading is copied into every group
*!          unchanged at every level, and no tie names it
*!   q2     Q-5: a fixed nonzero factor covariance likewise
*!   g3     G = 3: ties carry three operands; group 3 gets G8's later-group rows
*!   c1     configural accepts unequal K (no thresholds are tied there)
*!   Refusals, each rc 3498 with its message:
*!   e1     group 2 differs from group 1 (a loading start)
*!   e2     unequal K at the thresholds level
*!   e3     an unknown level name (the derivation's "loadings")
*!   e4     K of the wrong shape
*!   e5     a beta line in the configural input
*!
*! Mutants (each must turn this test red). Predicted | observed 2026-10-12, all red:
*!   M1  lambda ties emitted before tau ties            r3 r4    | r3 r4
*!   M2  omega for binary items at thresholds           r2       | r2 d1
*!   M3  metric-level psi rule inverted                 r3       | r3 d1
*!   M4  mirror validation removed                      e1       | e1
*!   M5  tau ties for k = 1..K, not 1..K-1              r2 r3 r4 | r2 r3 r4 d1 g3
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

clear all
set more off

local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
local load src/mata/dwlsmv_types.mata
foreach f of local srcs {
    local load `load' src/mata/`f'
}
foreach f of local load {
    capture noisily do "`f'"
    if _rc {
        display as error "test_levelgen: `f' failed to load"
        exit 9
    }
}

import delimited using "fixtures/data/mg_E.csv", clear varnames(1) delimiters(",")
mata: YALL = st_data(., .)

mata:

string scalar tlg_join(string rowvector lines)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= cols(lines); i++) {
        s = s + lines[i] + char(10)
    }
    return(s)
}

// Design E at a level: test_fit_groups.do's tfg_text(), level names per L2.
string scalar tlg_ref(string scalar level)
{
    string rowvector L
    real rowvector C
    real scalar j, k

    L = ("groups 2", "group 1", "factor F", "ind y1 y2 y3 y4 y5")
    for (j = 1; j <= 5; j++) {
        L = L, ("lambda F y" + strofreal(j) + " free")
    }
    L = L, ("metric unitvariance F", "group 2", "factor F", "ind y1 y2 y3 y4 y5")
    for (j = 1; j <= 5; j++) {
        L = L, ("lambda F y" + strofreal(j) + " free")
    }
    if (level == "configural" | level == "thresholds") {
        L = L, "psi F F fixed 1"
    }
    else {
        L = L, "psi F F free"
    }
    if (level == "thresholds" | level == "metric") {
        for (j = 1; j <= 5; j++) {
            L = L, ("nu y" + strofreal(j) + " free")
        }
    }
    if (level == "thresholds") {
        for (j = 1; j <= 3; j++) {
            L = L, ("omega y" + strofreal(j) + " free")
        }
    }
    if (level == "metric" | level == "scalar") {
        for (j = 1; j <= 5; j++) {
            L = L, ("omega y" + strofreal(j) + " free")
        }
    }
    if (level == "scalar") {
        L = L, "alpha F free"
    }
    if (level != "configural") {
        C = (3, 3, 3, 1, 1)
        for (j = 1; j <= 5; j++) {
            for (k = 1; k <= C[j]; k++) {
                L = L, ("tie tau 1 y" + strofreal(j) + " " + strofreal(k) +
                        " tau 2 y" + strofreal(j) + " " + strofreal(k))
            }
        }
    }
    if (level == "metric" | level == "scalar") {
        for (j = 1; j <= 5; j++) {
            L = L, ("tie lambda 1 F y" + strofreal(j) + " lambda 2 F y" + strofreal(j))
        }
    }
    return(tlg_join(L))
}

real matrix tlg_KE()
{
    return((4, 4, 4, 2, 2) \ (4, 4, 4, 2, 2))
}

// Two factors, six items, a fixed nonzero cross-loading and covariance, G groups
string scalar tlg_two(real scalar G)
{
    string rowvector L
    real scalar g

    L = "groups " + strofreal(G)
    for (g = 1; g <= G; g++) {
        L = L, ("group " + strofreal(g)), "factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
            "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
            "lambda F2 y3 fixed 0.3", "lambda F2 y4 free", "lambda F2 y5 free",
            "lambda F2 y6 free"
        if (g == 1) {
            L = L, "metric unitvariance F1", "metric unitvariance F2"
        }
        else {
            L = L, "psi F1 F1 fixed 1", "psi F2 F2 fixed 1"
        }
        L = L, "psi F2 F1 fixed 0.4"
    }
    return(tlg_join(L))
}

// Count of lines of t equal to s; count of tie lines containing u.
real scalar tlg_nline(string scalar t, string scalar s)
{
    string colvector v
    return(sum(dwlsmv_canon_lines(t) :== s))
}

real scalar tlg_ntie(string scalar t, string scalar u)
{
    string colvector v
    real scalar i, n

    v = dwlsmv_canon_lines(t)
    n = 0
    for (i = 1; i <= rows(v); i++) {
        if (substr(v[i], 1, 4) == "tie ") {
            if (strpos(v[i], u) > 0) n++
        }
    }
    return(n)
}

// 1 if the case passes, 0 if not (prints why)
real scalar tlg_case(string scalar id)
{
    external real matrix YALL
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_statmap vector maps
    struct dwlsmv_parmap scalar par
    string rowvector lev
    string scalar c, t, t2
    real colvector KE, sel
    real rowvector want
    real scalar g, i, ok
    real matrix K

    lev = ("configural", "thresholds", "metric", "scalar")
    c   = tlg_ref("configural")
    ok  = 1
    if (id == "r") {
        for (i = 1; i <= 4; i++) {
            t = dwlsmv_levelgen(c, tlg_KE(), lev[i])
            if (t != tlg_ref(lev[i])) {
                printf("    r%g (%s): generated text differs from the reference\n", i, lev[i])
                printf("%s\n--- reference ---\n%s\n", t, tlg_ref(lev[i]))
                ok = 0
            }
        }
    }
    else if (id == "r5") {
        t  = dwlsmv_levelgen(c, tlg_KE(), "configural")
        t2 = dwlsmv_levelgen(t, tlg_KE(), "configural")
        if (t != c | t2 != t) {
            printf("    r5: configural round trip is not exact\n")
            ok = 0
        }
    }
    else if (id == "d1") {
        KE   = (4 \ 4 \ 4 \ 2 \ 2)
        sts  = dwlsmv_stage1(2, 1)
        maps = dwlsmv_statmap(2, 1)
        for (g = 1; g <= 2; g++) {
            sel     = selectindex(YALL[., 1] :== g)
            sts[g]  = dwls_sample_stats(YALL[sel, 2..6], KE)
            maps[g] = sts[g].map
        }
        want = (32, 29, 27, 23)
        for (i = 1; i <= 4; i++) {
            par = dwls_identify(dwlsmv_levelgen(c, tlg_KE(), lev[i]), maps, sts)
            if (par.npar != want[i] | 42 - par.npar != (10, 13, 15, 19)[i]) {
                printf("    d1 (%s): npar %g, wanted %g\n", lev[i], par.npar, want[i])
                ok = 0
            }
        }
    }
    else if (id == "q1" | id == "q2") {
        K = J(2, 6, 4)
        for (i = 1; i <= 4; i++) {
            t = dwlsmv_levelgen(tlg_two(2), K, lev[i])
            if (id == "q1") {
                if (tlg_nline(t, "lambda F2 y3 fixed 0.3") != 2 | tlg_ntie(t, "F2 y3") != 0) {
                    printf("    q1 (%s): the fixed cross-loading is not copied unchanged, or is tied\n", lev[i])
                    ok = 0
                }
            }
            else {
                if (tlg_nline(t, "psi F2 F1 fixed 0.4") != 2 | tlg_ntie(t, "psi") != 0) {
                    printf("    q2 (%s): the fixed covariance is not copied unchanged, or is tied\n", lev[i])
                    ok = 0
                }
            }
        }
    }
    else if (id == "g3") {
        K = J(3, 6, 4)
        K[., 6] = J(3, 1, 2)
        t = dwlsmv_levelgen(tlg_two(3), K, "scalar")
        if (tlg_nline(t, "tie tau 1 y1 1 tau 2 y1 1 tau 3 y1 1") != 1 |
            tlg_nline(t, "tie lambda 1 F1 y1 lambda 2 F1 y1 lambda 3 F1 y1") != 1 |
            tlg_nline(t, "alpha F2 free") != 2 | tlg_nline(t, "omega y6 free") != 2 |
            tlg_nline(t, "psi F1 F1 free") != 2 | tlg_ntie(t, "y6 2") != 0) {
            printf("    g3: the three-group scalar text lacks a required line\n%s\n", t)
            ok = 0
        }
    }
    else if (id == "c1") {
        K = tlg_KE()
        K[2, 4] = 3
        t = dwlsmv_levelgen(c, K, "configural")
        if (t != c) {
            printf("    c1: configural with unequal K changed the text\n")
            ok = 0
        }
    }
    return(ok)
}

void tlg_err(string scalar id)
{
    string scalar c
    real matrix K

    c = tlg_ref("configural")
    if (id == "e1") {
        c = subinstr(c, "group 2" + char(10) + "factor F" + char(10) + "ind y1 y2 y3 y4 y5" + char(10) +
                     "lambda F y1 free", "group 2" + char(10) + "factor F" + char(10) +
                     "ind y1 y2 y3 y4 y5" + char(10) + "lambda F y1 free 0.5")
        (void) dwlsmv_levelgen(c, tlg_KE(), "metric")
    }
    else if (id == "e2") {
        K = tlg_KE()
        K[2, 4] = 3
        (void) dwlsmv_levelgen(c, K, "thresholds")
    }
    else if (id == "e3") {
        (void) dwlsmv_levelgen(c, tlg_KE(), "loadings")
    }
    else if (id == "e4") {
        (void) dwlsmv_levelgen(c, (4, 4, 4, 2, 2), "metric")
    }
    else {
        c = c + "beta F F fixed 0" + char(10)
        (void) dwlsmv_levelgen(c, tlg_KE(), "metric")
    }
}
end

local nfail 0
foreach id in r r5 d1 q1 q2 g3 c1 {
    capture noisily mata: st_local("ok", strofreal(tlg_case("`id'")))
    if _rc {
        display as error "    `id' FAILED: aborted, rc " _rc
        local ++nfail
    }
    else if !`ok' {
        display as error "    `id' FAILED"
        local ++nfail
    }
}

tempfile lg
foreach c in "e1|does not mirror group 1" "e2|has 2 categories in group 1 and 3" ///
             "e3|unknown level loadings" "e4|K must be 2 x 5" "e5|a configural input has no beta" {
    gettoken id msg : c, parse("|")
    local msg = substr("`msg'", 2, .)
    capture log close tlg
    log using "`lg'", text replace name(tlg)
    capture noisily mata: tlg_err("`id'")
    local rc = _rc
    log close tlg
    mata: st_local("hit", strofreal(strpos(invtokens(cat(st_local("lg"))'), st_local("msg")) > 0))
    if `rc' != 3498 | !`hit' {
        display as error "    `id' FAILED: rc `rc' (wanted 3498), message found `hit'"
        local ++nfail
    }
}

if `nfail' {
    display as error "test_levelgen: `nfail' failure(s)"
    exit 9
}
display as result "test_levelgen: the generator reproduces design E's four levels byte for byte, regenerates configural exactly, copies fixed rows untied (Q-5), handles G = 3, and refuses invalid inputs with rc 3498"
exit 0
