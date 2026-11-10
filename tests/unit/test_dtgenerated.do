*! Module 12, check 8's generated-chain exemption (dwlsmv_dtgenerated) and its
*! use by checks 8 and 9 (dwlsmv_dtnest, dwlsmv_dtchain); exemption gate step
*! 2b. Spec Module 12 checks 8 and 9 and 4.8 G8, decided 2026-10-08 and
*! 2026-10-13 (non-adjacent pairs (a), condition (4) option (ii), a failed
*! generated claim (a), check 9 level-rank tie-break (a)).
*!
*! Written and executed by the main session 2026-10-13 under StataMP 18; it
*! passes (green on the third run: the first two failed in the harness only,
*! tg_want reading a caller's local and Mata's results set in tg_log's scope,
*! then a compound-quoted permutation list split into seven; no routine changed).
*!
*! Inputs are posted and generated, not hand-typed: for each design (E, B, D of
*! fixtures/generate/99e_lavaan_mg_difftest.R; F -> y1-y5, group(g),
*! metric(unitvariance)) one invariance(configural) fit posts
*! e(canon_configural) and e(K); each level's e(canon) is then
*! dwlsmv_levelgen's output for it, and its table dwls_identify's, with Stage 1
*! computed here per group (the starts only). df = D - npar.
*!
*! Cases (want):
*!   p0 precondition: design E npar 32 29 27 23 (Result I4, test_levelgen d1)
*!   g1 E thresholds in configural, exempt                     rc 0, d 3
*!   g2 E scalar in configural (non-adjacent), exempt          rc 0, d 9
*!   g3 E, H1's e(canon_configural) differs                    3498, note (1) + tie
*!   g4 E, one element of H1's e(K) differs                    3498, note (2) + tie
*!   g5 E, H0's e(canon) has one extra byte                    3498, note (3) + tie
*!   g6 E, H0 posts no level (hand-written claim)              3498, no note, tie
*!   g7 E, a beta row added to H0: the beta guard runs first   3498, beta, no note
*!   g8 E, a direct exempt call with H0 = configural, H1 = scalar  3498, less restricted
*!   g9 E, an exempt pair whose df disagree with d             3498, df message
*!   c1 E chain, four levels named in a shuffled order: ord = scalar, metric,
*!      thresholds, configural; d = 4 2 3, every pair exempt
*!   z1 B chain (thresholds step d = 0) named in all six orders of
*!      configural, thresholds, metric: one ordering (metric, thresholds,
*!      configural); info d = (2, 0) [B's counts read in p0b], both exempt
*!   z2 B, configural (generated) tied with thresholds that posts no level:
*!      input order kept, refused by check 8                3498, tie message
*!   z3 D, thresholds and metric (the I2 step, d = 0): exempt, d = 0
*!
*! Mutants (each must turn this test red). Predicted | observed 2026-10-13
*! (scratch worktree, diff printed before each run, %TEMP%\exgate\mutate2b.sh),
*! all red:
*!   X1  condition (1) removed                     g3          | g3 (the edited
*!       base then fails in the generator, rc 3498 without the note)
*!   X2  condition (2) removed                     g4          | g4
*!   X3  condition (3) removed                     g5          | g5
*!   X4  exemption after the tie refusal (emulated: a tied exempt pair falls
*!       through)                                  g1 g2 c1 z1 z3 | those and g8 g9
*!   X5  beta guard after the exemption (emulated: skipped when generated
*!       inputs are passed)                        g7          | g7
*!   X6  d >= 1 kept for exempt pairs              z1 z3       | z1 z3
*!   X7  level-rank tie-break dropped              z1          | z1 (orders 1 2 5)
*!   X8  the failed-claim note removed             g3 g4 g5    | g3 g4 g5
*!   X9  the note printed for a single-level pair  g6          | g6
*!   X10 conditions (1) and (2) notes swapped      g3 g4       | g3 g4
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

* run a command under a log; rc into tg_rc, the log text into Mata TGLOG
capture program drop tg_log
program tg_log
    args cmd
    tempfile lg
    capture log close tgl
    quietly log using "`lg'", text replace name(tgl)
    capture noisily `cmd'
    local rc = _rc
    quietly log close tgl
    mata: TGLOG = invtokens(cat(st_local("lg"))')
    c_local tg_rc `rc'
    * results the Mata helpers set in this scope, passed to the caller
    foreach v in d ex ord dd ee {
        c_local `v' `"``v''"'
    }
end

* 1 if the text is in the last log
capture program drop tg_has
program tg_has
    args s
    mata: st_local("h", strofreal(strpos(TGLOG, st_local("s")) > 0))
    c_local hit `h'
end

* rc check: tg_want id got want
capture program drop tg_want
program tg_want
    args id got rc
    if `got' != `rc' {
        display as error "    `id' FAILED: rc `got', wanted `rc'"
        c_local bad 1
    }
end

mata:
// Keep design X's posted configural text and K, and its Stage 1 per group
// (computed here from the data in memory, groups in sorted g order).
void tg_keep(string scalar X)
{
    struct dwlsmv_stage1 vector sts
    real matrix K, Y
    real colvector gv
    real scalar G, g
    pointer() scalar p

    K  = st_matrix("e(K)")
    G  = rows(K)
    gv = uniqrows(st_data(., "g"))
    sts = dwlsmv_stage1(G, 1)
    for (g = 1; g <= G; g++) {
        Y      = select(st_data(., ("y1", "y2", "y3", "y4", "y5")), st_data(., "g") :== gv[g])
        sts[g] = dwls_sample_stats(Y, K[g, .]')
    }
    p = crexternal("TG_sts_" + X)
    if (p == NULL) p = findexternal("TG_sts_" + X)
    *p = sts
    p = crexternal("TG_cc_" + X)
    if (p == NULL) p = findexternal("TG_cc_" + X)
    *p = st_global("e(canon_configural)")
    p = crexternal("TG_K_" + X)
    if (p == NULL) p = findexternal("TG_K_" + X)
    *p = K
}

string scalar tg_cc(string scalar X)
{
    return(*findexternal("TG_cc_" + X))
}

real matrix tg_K(string scalar X)
{
    return(*findexternal("TG_K_" + X))
}

struct dwlsmv_parmap scalar tg_par(string scalar X, string scalar lev)
{
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_statmap vector maps
    real scalar g

    sts  = *findexternal("TG_sts_" + X)
    maps = dwlsmv_statmap(rows(sts), 1)
    for (g = 1; g <= rows(sts); g++) {
        maps[g] = sts[g].map
    }
    return(dwls_identify(dwlsmv_levelgen(tg_cc(X), tg_K(X), lev), maps, sts))
}

real scalar tg_D(string scalar X)
{
    struct dwlsmv_stage1 vector sts
    real scalar g, D

    sts = *findexternal("TG_sts_" + X)
    D   = 0
    for (g = 1; g <= rows(sts); g++) {
        D = D + sts[g].map.D
    }
    return(D)
}

// One pair through check 8, H0 = lev0, H1 = lev1, design X. mod alters the
// posted inputs: "cc" H1's base, "K" one element of H1's K, "can" one byte
// of H0's canon, "single" H0 posts no level, "beta" a beta row in H0, "df"
// H0's df off by one. Sets locals d and ex.
void tg_pair(string scalar X, string scalar lev0, string scalar lev1,
             string scalar mod)
{
    struct dwlsmv_parmap scalar p0, p1
    string scalar cc0, cc1, can0, can1, inv0
    real matrix K0, K1
    real scalar df0, df1, d, ex

    p0   = tg_par(X, lev0)
    p1   = tg_par(X, lev1)
    cc0  = tg_cc(X)
    cc1  = cc0
    K0   = tg_K(X)
    K1   = K0
    can0 = dwlsmv_levelgen(cc0, K0, lev0)
    can1 = dwlsmv_levelgen(cc0, K0, lev1)
    inv0 = lev0
    df0  = tg_D(X) - p0.npar
    df1  = tg_D(X) - p1.npar
    if (mod == "cc") {
        cc1 = cc1 + "groups 2" + char(10)
    }
    if (mod == "K") {
        K1[2, 1] = K1[2, 1] + 1
    }
    if (mod == "can") {
        can0 = can0 + " "
    }
    if (mod == "single") {
        inv0 = ""
    }
    if (mod == "beta") {
        p0.nall  = p0.nall + 1
        p0.ptype = p0.ptype \ 3
        p0.grp   = p0.grp \ 2
        p0.row   = p0.row \ 1
        p0.col   = p0.col \ 1
        p0.free  = p0.free \ 0
        p0.value = p0.value \ .5
        p0.fidx  = p0.fidx \ 0
        p0.name  = p0.name \ "g2:beta[F,F]"
    }
    if (mod == "df") {
        df0 = df0 + 1
    }
    ex = 0
    d  = dwlsmv_dtnest(p0, p1, df0, df1, X + ":" + lev0, X + ":" + lev1,
                       cc0, can0, inv0, K0, cc1, can1, lev1, K1, ex)
    st_local("d", strofreal(d))
    st_local("ex", strofreal(ex))
}

// A chain of design X with levels in the given (user's) order. Sets locals
// ord (level names in the returned order), dd (info d), ee (info exempt).
// inv0 = "" for every position listed in nolev.
void tg_chain(string scalar X, string rowvector levs, real rowvector nolev)
{
    struct dwlsmv_parmap vector pars
    real colvector df, ord
    string colvector names, cc, can, inv
    pointer(real matrix) colvector Kp
    real matrix info, K
    real scalar n, k
    string scalar s

    n     = cols(levs)
    pars  = dwlsmv_parmap(n, 1)
    df    = J(n, 1, .)
    names = levs'
    cc    = J(n, 1, tg_cc(X))
    can   = J(n, 1, "")
    inv   = levs'
    Kp    = J(n, 1, NULL)
    K     = tg_K(X)
    for (k = 1; k <= n; k++) {
        pars[k] = tg_par(X, levs[k])
        df[k]   = tg_D(X) - pars[k].npar
        can[k]  = dwlsmv_levelgen(cc[k], K, levs[k])
        Kp[k]   = &K
    }
    for (k = 1; k <= cols(nolev); k++) {
        inv[nolev[k]] = ""
    }
    ord = dwlsmv_dtchain(pars, df, names, cc, can, inv, Kp, info)
    s = ""
    for (k = 1; k <= n; k++) {
        s = s + " " + levs[ord[k]]
    }
    st_local("ord", strtrim(s))
    st_local("dd", invtokens(strofreal(info[., 1]')))
    st_local("ee", invtokens(strofreal(info[., 2]')))
}

// Permutation k (1-6) of configural, thresholds, metric, as the user's order.
void tg_chain_perm(string scalar X, real scalar k)
{
    string matrix P

    P = ("configural", "thresholds", "metric") \
        ("configural", "metric", "thresholds") \
        ("thresholds", "configural", "metric") \
        ("thresholds", "metric", "configural") \
        ("metric", "configural", "thresholds") \
        ("metric", "thresholds", "configural")
    tg_chain(X, P[k, .], J(1, 0, .))
}

void tg_npar(string scalar X)
{
    struct dwlsmv_parmap scalar p
    string rowvector L
    real scalar k
    string scalar s

    L = ("configural", "thresholds", "metric", "scalar")
    s = ""
    for (k = 1; k <= 4; k++) {
        p = tg_par(X, L[k])
        s = s + " " + strofreal(p.npar)
    }
    st_local("np", strtrim(s))
}
end

* ---------------------------------------------------------------- inputs
foreach X in E B D {
    import delimited using "fixtures/data/mg_`X'.csv", clear varnames(1) delimiters(",")
    quietly dwlsmv (F -> y1 y2 y3 y4 y5), group(g) invariance(configural) metric(unitvariance)
    if "`e(invariance)'" != "configural" {
        display as error "    inputs FAILED: design `X' posted no configural fit"
        exit 9
    }
    mata: tg_keep("`X'")
    mata: tg_npar("`X'")
    local np`X' `np'
    display as text "    design `X': npar (configural thresholds metric scalar) = `np'"
}
if "`npE'" != "32 29 27 23" {
    display as error "    p0 FAILED: design E npar `npE', wanted 32 29 27 23"
    exit 9
}
* B: the thresholds step d = 0; D: the metric step d = 0 (99e's designs)
local nB : word 1 of `npB'
local tB : word 2 of `npB'
local mB : word 3 of `npB'
local tD : word 2 of `npD'
local mD : word 3 of `npD'
if `nB' != `tB' | `tD' != `mD' {
    display as error "    p0 FAILED: B (`npB') or D (`npD') lacks its d = 0 step"
    exit 9
}
local dB2 = `tB' - `mB'
display as text "    p0 ok"

local tie "a model has equality constraints (tied parameters)"
local nte "the pair is not treated as a generated invariance"

* ---------------------------------------------------------------- pairs
local bad 0
tg_log `"mata: tg_pair("E", "thresholds", "configural", "")"'
tg_want g1 `tg_rc' 0
if `tg_rc' == 0 & ("`d'" != "3" | "`ex'" != "1") {
    display as error "    g1 FAILED: d `d' exempt `ex', wanted 3 1"
    local bad 1
}
tg_has "`nte'"
if `hit' {
    display as error "    g1 FAILED: a note on an exempt pair"
    local bad 1
}
if `bad' local ++nfail
else display as text "    g1 ok"

local bad 0
tg_log `"mata: tg_pair("E", "scalar", "configural", "")"'
tg_want g2 `tg_rc' 0
if `tg_rc' == 0 & ("`d'" != "9" | "`ex'" != "1") {
    display as error "    g2 FAILED: d `d' exempt `ex', wanted 9 1"
    local bad 1
}
if `bad' local ++nfail
else display as text "    g2 ok"

local i 0
foreach m in cc K can {
    local ++i
    local id = "g" + string(`i' + 2)
    local why1 "different e(canon_configural)"
    local why2 "different e(K)"
    local why3 "e(canon) of E:metric is not the generator's text"
    local bad 0
    tg_log `"mata: tg_pair("E", "metric", "configural", "`m'")"'
    tg_want `id' `tg_rc' 3498
    tg_has "`nte'"
    local h1 `hit'
    tg_has "`why`i''"
    local h2 `hit'
    tg_has "`tie'"
    if !`h1' | !`h2' | !`hit' {
        display as error "    `id' FAILED: note `h1', condition `h2', tie message `hit'"
        local bad 1
    }
    if `bad' local ++nfail
    else display as text "    `id' ok"
}

local bad 0
tg_log `"mata: tg_pair("E", "metric", "configural", "single")"'
tg_want g6 `tg_rc' 3498
tg_has "`nte'"
local h1 `hit'
tg_has "`tie'"
if `h1' | !`hit' {
    display as error "    g6 FAILED: note `h1' (want 0), tie message `hit'"
    local bad 1
}
if `bad' local ++nfail
else display as text "    g6 ok"

local bad 0
tg_log `"mata: tg_pair("E", "metric", "configural", "beta")"'
tg_want g7 `tg_rc' 3498
tg_has "groups and a structural path (beta)"
local h1 `hit'
tg_has "`nte'"
if !`h1' | `hit' {
    display as error "    g7 FAILED: beta message `h1', note `hit' (want 1, 0)"
    local bad 1
}
if `bad' local ++nfail
else display as text "    g7 ok"

local bad 0
tg_log `"mata: tg_pair("E", "configural", "scalar", "")"'
tg_want g8 `tg_rc' 3498
tg_has "is the less restricted model"
if !`hit' {
    display as error "    g8 FAILED: not the d < 0 message"
    local bad 1
}
if `bad' local ++nfail
else display as text "    g8 ok"

local bad 0
tg_log `"mata: tg_pair("E", "metric", "configural", "df")"'
tg_want g9 `tg_rc' 3498
tg_has "d disagrees with the difference in e(df_T3)"
if !`hit' {
    display as error "    g9 FAILED: not the df message"
    local bad 1
}
if `bad' local ++nfail
else display as text "    g9 ok"

* ---------------------------------------------------------------- chains
local bad 0
tg_log `"mata: tg_chain("E", ("metric", "configural", "scalar", "thresholds"), J(1, 0, .))"'
tg_want c1 `tg_rc' 0
if `tg_rc' == 0 & ("`ord'" != "scalar metric thresholds configural" | "`dd'" != "4 2 3" | "`ee'" != "1 1 1") {
    display as error "    c1 FAILED: ord `ord', d `dd', exempt `ee'"
    local bad 1
}
if `bad' local ++nfail
else display as text "    c1 ok"

local bad 0
forvalues k = 1/6 {
    tg_log `"mata: tg_chain_perm("B", `k')"'
    if `tg_rc' != 0 | "`ord'" != "metric thresholds configural" | "`dd'" != "`dB2' 0" | "`ee'" != "1 1" {
        display as error "    z1 order `k' FAILED: rc `tg_rc', ord `ord', d `dd', exempt `ee'"
        local bad 1
    }
}
if `bad' local ++nfail
else display as text "    z1 ok (six orders)"

local bad 0
tg_log `"mata: tg_chain("B", ("configural", "thresholds"), 2)"'
tg_want z2 `tg_rc' 3498
tg_has "`tie'"
if !`hit' {
    display as error "    z2 FAILED: not the tie message"
    local bad 1
}
if `bad' local ++nfail
else display as text "    z2 ok"

local bad 0
tg_log `"mata: tg_pair("D", "metric", "thresholds", "")"'
tg_want z3 `tg_rc' 0
if `tg_rc' == 0 & ("`d'" != "0" | "`ex'" != "1") {
    display as error "    z3 FAILED: d `d' exempt `ex', wanted 0 1"
    local bad 1
}
if `bad' local ++nfail
else display as text "    z3 ok"

if `nfail' {
    display as error "test_dtgenerated: `nfail' failure(s)"
    exit 9
}
display as result "test_dtgenerated: all passed"
