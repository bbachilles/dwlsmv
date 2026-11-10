*! Module 12, one DIFFTEST pair after checks 8 and 9: dwlsmv_dtstep (check 11's
*! comparison and backstop callers, check 8's zero-df pairs) and check 10's
*! negative-phi branch in dwlsmv_dtband. Exemption gate step 2c; spec Module 12
*! checks 8, 10 and 11 (decided 2026-10-13: zero-df option (c), non-adjacent
*! wording (a), negative H0 phi option (c), messages and option (ii)).
*!
*! Written and executed by the main session 2026-10-13 under StataMP 18; it
*! passes (green on the fourth run: the first three failed to build or run in
*! the harness, a missing helper, a struct scalar declared optional in
*! dwlsmv_dtband (Mata refuses it; declared without scalar), and the note
*! external; no case or wanted value changed).
*!
*! Inputs as in test_dtgenerated.do: one invariance(configural) fit per design
*! posts e(canon_configural) and e(K); levels are dwlsmv_levelgen's output, tables
*! dwls_identify's, d and exempt from dwlsmv_dtnest with the generated inputs.
*! theta0 is H0's start vector, edited where a case needs an estimate. T0, T1
*! are constructed; c = 1000, nrtol = 1e-12, so eps = 1e-9. Designs E, B, D of
*! 99e, and F: mg_D.csv with y1 collapsed to 3 categories (every item K <= 3,
*! one polytomous item), so configural, thresholds and metric have equal npar.
*!
*! Cases (want):
*!   p0  preconditions: B thresholds = configural npar; D thresholds = metric;
*!       F configural = thresholds = metric
*!  zero-df (check 8):
*!   z1  D metric vs thresholds (I2), |Td| = 1e-12          status 0, Td 0, prefix + verified
*!   z2  the same pair, Td = 1                              3498, I2 message
*!   z3  B thresholds vs configural, |Td| = 1e-12           status 0, equal-df, no info line
*!   z4  B, Td = 1                                          rc 0, status 0, equal-df + info line
*!   z5  D chain configural, thresholds, metric             (metric, thresholds) status 0;
*!                                                          (thresholds, configural) status 1
*!   z6  F metric vs configural (non-adjacent, I2 holds), Td = 1  rc 0, equal-df + info, not verified
*!  negative H0 phi (check 10), E metric vs thresholds, d = 2:
*!   n1  exempt, phi_g2 = -.2, T0 - T1 = -1                 3498, negative-variance message
*!   n2  the same values user-written (exempt 0)            3498, old message
*!   n3  exempt, phi_g2 = 0                                 3498, old message
*!   n4  exempt, phi_g2 = -.2, T0 - T1 = +1                 status 1
*!  check 11's callers, E metric with psi_g2[F,F] fixed at 0:
*!   b1  exempt                                             3498, backstop sentence
*!   b2  user-written                                       3498, comparison sentence
*!
*! Mutants (each must turn this test red). Predicted | observed 2026-10-13
*! (scratch worktree, diff printed before each run, %TEMP%\exgate\mutate2c.sh),
*! all red:
*!   Z1  the statistic path taken at d = 0           z1 z3 z4 z5 z6 | those and z2
*!   Z2  the I2 abort removed                        z2             | z2
*!   Z3  the abort applied to unverified pairs       z4 z6          | z4 z6
*!   Z4  the band replaced by exact zero             z1 z3          | z1 z3 z5
*!   Z5  verified wording for any pair with metric   z6             | z6 (aborts)
*!   N1  the negative-phi branch removed             n1             | n1
*!   N2  the branch not limited to exempt pairs (a user-written pair handed
*!       to dtband as exempt)                        n2             | n2
*!   N3  <= for <                                    n3             | n3
*!   C1  backstop and comparison sentences swapped   b1 b2          | b1 b2
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

capture program drop ts_log
program ts_log
    args cmd
    tempfile lg
    capture log close tsl
    quietly log using "`lg'", text replace name(tsl)
    capture noisily `cmd'
    local rc = _rc
    quietly log close tsl
    mata: TSLOG = invtokens(cat(st_local("lg"))')
    c_local ts_rc `rc'
    foreach v in st td d ex st2 {
        c_local `v' `"``v''"'
    }
end

* 1 if the text is in the last log (where = LOG) or in the returned note (NOTE)
capture program drop ts_has
program ts_has
    args where s
    mata: st_local("h", strofreal(strpos(TS`where', st_local("s")) > 0))
    c_local hit `h'
end

mata:
TSNOTE = ""

void ts_keep(string scalar X)
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
    p = crexternal("TS_sts_" + X)
    if (p == NULL) p = findexternal("TS_sts_" + X)
    *p = sts
    p = crexternal("TS_cc_" + X)
    if (p == NULL) p = findexternal("TS_cc_" + X)
    *p = st_global("e(canon_configural)")
    p = crexternal("TS_K_" + X)
    if (p == NULL) p = findexternal("TS_K_" + X)
    *p = K
}

struct dwlsmv_parmap scalar ts_par(string scalar X, string scalar lev)
{
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_statmap vector maps
    real scalar g

    sts  = *findexternal("TS_sts_" + X)
    maps = dwlsmv_statmap(rows(sts), 1)
    for (g = 1; g <= rows(sts); g++) {
        maps[g] = sts[g].map
    }
    return(dwls_identify(dwlsmv_levelgen(*findexternal("TS_cc_" + X),
                                         *findexternal("TS_K_" + X), lev), maps, sts))
}

real scalar ts_D(string scalar X)
{
    struct dwlsmv_stage1 vector sts
    real scalar g, D

    sts = *findexternal("TS_sts_" + X)
    D   = 0
    for (g = 1; g <= rows(sts); g++) {
        D = D + sts[g].map.D
    }
    return(D)
}

// Position of group 2's psi[F,F] row in par (0 if none).
real scalar ts_phi2(struct dwlsmv_parmap scalar par)
{
    real scalar i

    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 2 & par.grp[i] == 2 & par.row[i] == par.col[i]) {
            return(i)
        }
    }
    return(0)
}

// One pair H0 = lev0, H1 = lev1 of design X through check 8 (generated
// inputs), then dwlsmv_dtstep with T0, T1. mod: "neg" phi_g2 = -.2 in theta0,
// "zero" phi_g2 = 0, "user" exempt passed as 0, "bound" psi_g2[F,F] fixed at
// 0 in par0. Sets locals st, td, d, ex; TSNOTE the returned note.
void ts_step(string scalar X, string scalar lev0, string scalar lev1,
             real scalar T0, real scalar T1, string scalar mod)
{
    struct dwlsmv_parmap scalar p0, p1
    string scalar cc, can0, can1, nt
    real matrix K
    real colvector th
    real rowvector b
    real scalar d, ex, i, df0, df1

    p0  = ts_par(X, lev0)
    p1  = ts_par(X, lev1)
    cc  = *findexternal("TS_cc_" + X)
    K   = *findexternal("TS_K_" + X)
    can0 = dwlsmv_levelgen(cc, K, lev0)
    can1 = dwlsmv_levelgen(cc, K, lev1)
    df0 = ts_D(X) - p0.npar
    df1 = ts_D(X) - p1.npar
    ex  = 0
    d   = dwlsmv_dtnest(p0, p1, df0, df1, X + ":" + lev0, X + ":" + lev1,
                        cc, can0, lev0, K, cc, can1, lev1, K, ex)
    th = J(p0.npar, 1, .)
    for (i = 1; i <= p0.nall; i++) {
        if (p0.free[i]) {
            th[p0.fidx[i]] = p0.value[i]
        }
    }
    i = ts_phi2(p0)
    if (strpos(mod, "neg") > 0) {
        th[p0.fidx[i]] = -.2
    }
    if (strpos(mod, "zero") > 0) {
        th[p0.fidx[i]] = 0
    }
    if (strpos(mod, "user") > 0) {
        ex = 0
    }
    if (strpos(mod, "bound") > 0) {
        p0.free[i]  = 0
        p0.value[i] = 0
        p0.fidx[i]  = 0
    }
    nt = ""
    b  = dwlsmv_dtstep(p0, th, d, ex, T0, T1, 1000, 1e-12, df0, cc, K, lev0, lev1,
                       X + ":" + lev0, X + ":" + lev1, nt)
    *findexternal("TSNOTE") = nt
    st_local("st", strofreal(b[1]))
    st_local("td", strofreal(b[3]))
    st_local("d", strofreal(d))
    st_local("ex", strofreal(ex))
}

// D's chain named configural, thresholds, metric; dwlsmv_dtstep on each pair
// in dwlsmv_dtchain's order with T(metric) = 10, T(thresholds) = 10 - 1e-12,
// T(configural) = 5. Sets locals st (first pair) and st2 (second).
void ts_chain()
{
    struct dwlsmv_parmap vector pars
    string rowvector levs
    string colvector cc, can, inv, nm
    pointer(real matrix) colvector Kp
    real colvector df, ord, T, th
    real matrix info, K
    real rowvector b
    real scalar k, a, i
    string scalar nt

    levs = ("configural", "thresholds", "metric")
    T    = (5 \ 10 - 1e-12 \ 10)
    pars = dwlsmv_parmap(3, 1)
    df   = J(3, 1, .)
    cc   = J(3, 1, *findexternal("TS_cc_D"))
    can  = J(3, 1, "")
    Kp   = J(3, 1, NULL)
    K    = *findexternal("TS_K_D")
    for (k = 1; k <= 3; k++) {
        pars[k] = ts_par("D", levs[k])
        df[k]   = ts_D("D") - pars[k].npar
        can[k]  = dwlsmv_levelgen(cc[k], K, levs[k])
        Kp[k]   = &K
    }
    inv = levs'
    nm  = levs'
    ord = dwlsmv_dtchain(pars, df, nm, cc, can, inv, Kp, info)
    for (k = 1; k <= 2; k++) {
        a  = ord[k]
        th = J(pars[a].npar, 1, .)
        for (i = 1; i <= pars[a].nall; i++) {
            if (pars[a].free[i]) {
                th[pars[a].fidx[i]] = pars[a].value[i]
            }
        }
        b = dwlsmv_dtstep(pars[a], th, info[k, 1], info[k, 2], T[a], T[ord[k + 1]],
                          1000, 1e-12, df[a], cc[a], K, levs[a], levs[ord[k + 1]],
                          levs[a], levs[ord[k + 1]], nt)
        st_local(k == 1 ? "st" : "st2", strofreal(b[1]))
    }
}

void ts_npar(string scalar X)
{
    struct dwlsmv_parmap scalar p
    string rowvector L
    real scalar k
    string scalar s

    L = ("configural", "thresholds", "metric", "scalar")
    s = ""
    for (k = 1; k <= 4; k++) {
        p = ts_par(X, L[k])
        s = s + " " + strofreal(p.npar)
    }
    st_local("np", strtrim(s))
}
end

* ---------------------------------------------------------------- inputs
foreach X in E B D F {
    local src `X'
    if "`X'" == "F" local src D
    import delimited using "fixtures/data/mg_`src'.csv", clear varnames(1) delimiters(",")
    if "`X'" == "F" {
        quietly replace y1 = min(y1, 2)
    }
    quietly dwlsmv (F -> y1 y2 y3 y4 y5), group(g) invariance(configural) metric(unitvariance)
    mata: ts_keep("`X'")
    mata: ts_npar("`X'")
    local np`X' `np'
    display as text "    design `X': npar (configural thresholds metric scalar) = `np'"
}
tokenize `npB'
local okB = (`1' == `2')
tokenize `npD'
local okD = (`2' == `3')
tokenize `npF'
local okF = (`1' == `2' & `2' == `3')
if !`okB' | !`okD' | !`okF' {
    display as error "    p0 FAILED: B `npB', D `npD', F `npF'"
    exit 9
}
display as text "    p0 ok"

local pre "No test for"
local ver "is equivalent to invariance(thresholds) for this model"
local edf "has the same degrees of freedom as invariance("
local inf "so equivalence of the two models could not be established"
local i2m "is mathematically the"
local neg "estimates include a factor variance below 0"
local old "did not reach its"
local bks "The level generator introduced this fixed value"
local cmp "where the chi-squared reference of the difference test"

* ---------------------------------------------------------------- zero-df
local bad 0
ts_log `"mata: ts_step("D", "metric", "thresholds", 10, 10 - 1e-12, "")"'
ts_has NOTE "`pre'"
local h1 `hit'
ts_has NOTE "`ver'"
if `ts_rc' != 0 | "`st'" != "0" | "`td'" != "0" | "`ex'" != "1" | !`h1' | !`hit' {
    display as error "    z1 FAILED: rc `ts_rc' status `st' Td `td' exempt `ex' prefix `h1' verified `hit'"
    local bad 1
}
if `bad' local ++nfail
else display as text "    z1 ok"

local bad 0
ts_log `"mata: ts_step("D", "metric", "thresholds", 11, 10, "")"'
ts_has LOG "`i2m'"
if `ts_rc' != 3498 | !`hit' {
    display as error "    z2 FAILED: rc `ts_rc', I2 message `hit'"
    local bad 1
}
if `bad' local ++nfail
else display as text "    z2 ok"

local bad 0
ts_log `"mata: ts_step("B", "thresholds", "configural", 10, 10 + 1e-12, "")"'
ts_has NOTE "`edf'"
local h1 `hit'
ts_has NOTE "`inf'"
local h2 `hit'
ts_has NOTE "`ver'"
if `ts_rc' != 0 | "`st'" != "0" | "`td'" != "0" | !`h1' | `h2' | `hit' {
    display as error "    z3 FAILED: rc `ts_rc' status `st' Td `td' equal-df `h1' info `h2' verified `hit'"
    local bad 1
}
if `bad' local ++nfail
else display as text "    z3 ok"

local bad 0
ts_log `"mata: ts_step("B", "thresholds", "configural", 11, 10, "")"'
ts_has NOTE "`edf'"
local h1 `hit'
ts_has NOTE "`inf'"
if `ts_rc' != 0 | "`st'" != "0" | "`td'" != "1" | !`h1' | !`hit' {
    display as error "    z4 FAILED: rc `ts_rc' status `st' Td `td' equal-df `h1' info `hit'"
    local bad 1
}
if `bad' local ++nfail
else display as text "    z4 ok"

local bad 0
ts_log `"mata: ts_chain()"'
if `ts_rc' != 0 | "`st'" != "0" | "`st2'" != "1" {
    display as error "    z5 FAILED: rc `ts_rc', first pair status `st' (want 0), second `st2' (want 1)"
    local bad 1
}
if `bad' local ++nfail
else display as text "    z5 ok"

local bad 0
ts_log `"mata: ts_step("F", "metric", "configural", 11, 10, "")"'
ts_has NOTE "`edf'"
local h1 `hit'
ts_has NOTE "`inf'"
local h2 `hit'
ts_has NOTE "`ver'"
if `ts_rc' != 0 | "`st'" != "0" | "`ex'" != "1" | !`h1' | !`h2' | `hit' {
    display as error "    z6 FAILED: rc `ts_rc' status `st' exempt `ex' equal-df `h1' info `h2' verified `hit'"
    local bad 1
}
if `bad' local ++nfail
else display as text "    z6 ok"

* ---------------------------------------------------------------- check 10
* n1 "neg": new message; n2 "neg user" and n3 "zero": the old message
local i 0
foreach m in neg "neg user" zero {
    local ++i
    local id n`i'
    local want  = cond(`i' == 1, "`neg'", "`old'")
    local nwant = cond(`i' == 1, "`old'", "`neg'")
    local bad 0
    ts_log `"mata: ts_step("E", "metric", "thresholds", 1, 2, "`m'")"'
    ts_has LOG "`want'"
    local h1 `hit'
    ts_has LOG "`nwant'"
    if `ts_rc' != 3498 | !`h1' | `hit' {
        display as error "    `id' FAILED: rc `ts_rc', wanted message `h1', other message `hit'"
        local bad 1
    }
    if `bad' local ++nfail
    else display as text "    `id' ok"
}

local bad 0
ts_log `"mata: ts_step("E", "metric", "thresholds", 2, 1, "neg")"'
if `ts_rc' != 0 | "`st'" != "1" {
    display as error "    n4 FAILED: rc `ts_rc' status `st'"
    local bad 1
}
if `bad' local ++nfail
else display as text "    n4 ok"

* ---------------------------------------------------------------- check 11
local bad 0
ts_log `"mata: ts_step("E", "metric", "thresholds", 2, 1, "bound")"'
ts_has LOG "`bks'"
local h1 `hit'
ts_has LOG "`cmp'"
if `ts_rc' != 3498 | !`h1' | `hit' {
    display as error "    b1 FAILED: rc `ts_rc', backstop `h1', comparison `hit'"
    local bad 1
}
if `bad' local ++nfail
else display as text "    b1 ok"

local bad 0
ts_log `"mata: ts_step("E", "metric", "thresholds", 2, 1, "bound user")"'
ts_has LOG "`cmp'"
local h1 `hit'
ts_has LOG "`bks'"
if `ts_rc' != 3498 | !`h1' | `hit' {
    display as error "    b2 FAILED: rc `ts_rc', comparison `h1', backstop `hit'"
    local bad 1
}
if `bad' local ++nfail
else display as text "    b2 ok"

if `nfail' {
    display as error "test_dtstep: `nfail' failure(s)"
    exit 9
}
display as result "test_dtstep: all passed"
