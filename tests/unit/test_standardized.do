*! Module 11 D7 step (3): dwlsmv, standardized (dwlsmv_m11_std).
*!
*! Executed under StataMP 18 on 2026-10-05 and it passes; six mutations of
*! dwlsmv_m11_std and dwlsmv.ado (J_g terms, Phi vs psi, the estimates hold)
*! each fail it.
*!
*! Builds build/ldwlsmv.mlib first and runs against the library and src/ado,
*! as a user would; both adopath entries are removed at the end.
*!
*! Every tolerance below was set before the first run.
*!
*! Checks:
*!   1  AGAINST LAVAAN (fixtures/expected/lavaan_case_E_standardized.csv,
*!      tier 2): Case E's structural (F2 <- F1) and correlated (F1 ~~ F2)
*!      models, marker and unit variance. Every one of the 27 rows per fit
*!      is matched by name (loadings, path, factor ~~ rows, Theta rows
*!      /var(e.y), thresholds) and nothing is left over. Values at spec
*!      7.2's theta-hat tier (1e-5 abs + 1e-4 rel), SEs at its SE tier
*!      (1e-8 abs + 1e-5 rel): lavaan's forward-difference SE error, ~1.7e-8
*!      relative (derivation @sec-standardize), is far inside it. A row
*!      lavaan gives SE 0 must have SE exactly 0 here, and only those.
*!   2  J_g AGAINST CENTRAL DIFFERENCES of dwlsmv_m11_std itself, h = 1e-6,
*!      over e(b)'s free entries (e(b) reposted): J_fd V J_fd' against the
*!      posted standardized V, max abs gap over max |V_std| < 1e-8 (the
*!      derivation's claim-2 budget: J error ~1e-10 on O(1) entries).
*!   3  THE DISPLAY ROWS: structural, var(e.F2) = 1 - (F2:F1)^2 (1 - R^2,
*!      one predictor), 1e-12; correlated, cov(F1,F2) is r_12 and var(F1),
*!      var(F2) are 1 with SE 0.
*!   4  NOTHING STORED (D8d): after the fit with standardized and after
*!      replay with standardized, e(b), e(V), e(Cns) are bit-identical to the
*!      plain fit's, no e() name is added, and every scalar, macro (but
*!      e(cmdline)) and e(bound) entry is unchanged; replay with standardized
*!      after another command is rc 301. Names are compared sorted: the
*!      display's _estimates hold/restore returns the scalars in reverse
*!      listing order (seen on the first run), which is not a stored result.
*!   5  THE nlcom STOPGAP (owed by D7): the help file's two expressions,
*!      marker fit, against lavaan's standardized y2 and y5 at the tiers of
*!      check 1.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 120

quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

local nfail 0

* ---- the fixture, into Mata --------------------------------------------------
import delimited using "fixtures/expected/lavaan_case_E_standardized.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: tst_fx = st_sdata(., ("fit", "lhs", "op", "rhs", "est_std", "se"))

mata:
// lavaan's row -> dwlsmv's full name (D6), for Case E
string scalar tst_name(string rowvector r)
{
    string scalar l, o, h
    l = r[2]
    o = r[3]
    h = r[4]
    if (o == "=~") {
        return(h + ":" + strupper(l))
    }
    if (o == "~") {
        return(strupper(l) + ":" + strupper(h))
    }
    if (o == "|") {
        return("/" + l + ":cut" + substr(h, 2, .))
    }
    if (l != h) {
        return("/cov(" + strupper(l) + "," + strupper(h) + ")")
    }
    if (substr(l, 1, 1) == "y") {
        return("/var(e." + l + ")")
    }
    if (substr(r[1], 1, 3) == "str" & l == "f2") {
        return("/var(e.F2)")
    }
    return("/var(" + strupper(l) + ")")
}

// check 1 for one fit: dwlsmv_m11_std on the current e() against the fixture
void tst_vs_lavaan(string scalar fitn)
{
    external string matrix tst_fx
    string matrix stp, fx
    string colvector nm
    real rowvector g
    real matrix V
    real scalar i, k, n, a, b, se, lse, re, rs, bad, used

    dwlsmv_m11_std("tst_bs", "tst_Vs", "tst_Cs")
    g   = st_matrix("tst_bs")
    V   = st_matrix("tst_Vs")
    stp = st_matrixcolstripe("tst_bs")
    nm  = J(rows(stp), 1, "")
    for (i = 1; i <= rows(stp); i++) {
        nm[i] = (stp[i, 1] == "/" ? "/" + stp[i, 2] : stp[i, 1] + ":" + stp[i, 2])
    }
    fx   = select(tst_fx, tst_fx[., 1] :== fitn)
    n    = rows(fx)
    used = J(rows(nm), 1, 0)
    bad  = 0
    re   = 0
    rs   = 0
    for (i = 1; i <= n; i++) {
        k = selectindex(nm :== tst_name(fx[i, .]))
        if (length(k) != 1) {
            printf("{err}    %s: lavaan row %s %s %s has no dwlsmv row\n", fitn, fx[i, 2], fx[i, 3], fx[i, 4])
            bad++
            continue
        }
        used[k] = 1
        a   = g[k]
        b   = strtoreal(fx[i, 5])
        se  = sqrt(V[k, k])
        lse = strtoreal(fx[i, 6])
        re  = max((re, abs(a - b) / (1e-5 + 1e-4 * abs(b))))
        if (lse == 0) {
            if (V[k, k] != 0) {
                printf("{err}    %s: %s has SE 0 in lavaan, %g here\n", fitn, nm[k], se)
                bad++
            }
        }
        else {
            if (V[k, k] == 0) {
                printf("{err}    %s: %s has SE 0 here, %g in lavaan\n", fitn, nm[k], lse)
                bad++
            }
            rs = max((rs, abs(se - lse) / (1e-8 + 1e-5 * abs(lse))))
        }
    }
    printf("{txt} 1. %-11s %f rows matched of %f; value gap %8.2e of tier, SE gap %8.2e of tier\n",
           fitn, sum(used), rows(nm), re, rs)
    if (bad | sum(used) != rows(nm) | n != rows(nm) | re > 1 | rs > 1) {
        st_local("nfail", strofreal(strtoreal(st_local("nfail")) + 1))
    }
}

// check 2: collect one central-difference column
void tst_fdcol(real scalar k, real scalar h)
{
    external real matrix tst_Jfd
    tst_Jfd[., k] = (st_matrix("tst_gp") - st_matrix("tst_gm"))' / (2 * h)
}
end

* e()'s names, sorted, and every scalar's and macro's value (except
* e(cmdline), which records the call); e(bound) by value. Sorted because
* _estimates hold/restore returns e()'s scalars in another order: the
* listing order is not a stored result.
capture program drop tst_snap
program tst_snap, rclass
    local s : e(scalars)
    local s : list sort s
    local m : e(macros)
    local m : list sort m
    local x : e(matrices)
    local x : list sort x
    local out "S `s' M `m' X `x' |"
    foreach k of local s {
        local out "`out' `k'=`=strofreal(e(`k'), "%21.0g")'"
    }
    foreach k of local m {
        if "`k'" != "cmdline" local out `"`out' `k'=`e(`k')'"'
    }
    tempname bd
    matrix `bd' = e(bound)
    forvalues k = 1/`=colsof(`bd')' {
        local out "`out' `=`bd'[1, `k']'"
    }
    return local snap `"`out'"'
end

capture program drop tst_repost
program tst_repost, eclass
    ereturn repost b = `1'
end

* ---- the data ------------------------------------------------------------------
import delimited using "fixtures/data/case_E.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
quietly destring y*, replace
local mS "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)"
local mC "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6)"

* ---- 1. against lavaan -----------------------------------------------------------
quietly dwlsmv `mS'
mata: tst_vs_lavaan("str_marker")
quietly dwlsmv `mS', metric(unitvariance)
mata: tst_vs_lavaan("str_unitvar")
quietly dwlsmv `mC'
mata: tst_vs_lavaan("cor_marker")
quietly dwlsmv `mC', metric(unitvariance)
mata: tst_vs_lavaan("cor_unitvar")

* ---- 2. J_g against central differences ------------------------------------------
quietly dwlsmv `mS'
tempname b0 bp V0 Vs
matrix `b0' = e(b)
matrix `V0' = e(V)
mata: dwlsmv_m11_std("tst_bs", "tst_Vs", "tst_Cs")
matrix `Vs' = tst_Vs
local q = colsof(`b0')
local h = 1e-6
mata: tst_Jfd = J(cols(st_matrix("tst_bs")), `q', 0)
forvalues k = 1/`q' {
    if el(`V0', `k', `k') == 0 continue
    matrix `bp' = `b0'
    matrix `bp'[1, `k'] = `b0'[1, `k'] + `h'
    tst_repost `bp'
    mata: dwlsmv_m11_std("tst_gp", "tst_x", "tst_y")
    matrix `bp' = `b0'
    matrix `bp'[1, `k'] = `b0'[1, `k'] - `h'
    tst_repost `bp'
    mata: dwlsmv_m11_std("tst_gm", "tst_x", "tst_y")
    mata: tst_fdcol(`k', `h')
}
matrix `bp' = `b0'
tst_repost `bp'
mata: st_local("r2", strofreal(max(abs(tst_Jfd * st_matrix("`V0'") * tst_Jfd' - st_matrix("`Vs'"))) / max(abs(st_matrix("`Vs'")))))
display as text " 2. J_g vs central differences: max |J_fd V J_fd' - V_std| / max |V_std| = " %9.2e `r2' " (1e-8)"
if `r2' >= 1e-8 | missing(`r2') {
    local ++nfail
}

* ---- 3. the display rows ---------------------------------------------------------
quietly dwlsmv `mS', metric(unitvariance)
mata: dwlsmv_m11_std("tst_bs", "tst_Vs", "tst_Cs")
local d3 = abs(tst_bs[1, colnumb(tst_bs, "/var(e.F2)")] - (1 - tst_bs[1, colnumb(tst_bs, "F2:F1")]^2))
quietly dwlsmv `mC'
mata: dwlsmv_m11_std("tst_bs", "tst_Vs", "tst_Cs")
local i1 = colnumb(tst_bs, "/var(F1)")
local i2 = colnumb(tst_bs, "/var(F2)")
local ok3 = tst_bs[1, `i1'] == 1 & tst_bs[1, `i2'] == 1 & tst_Vs[`i1', `i1'] == 0 & tst_Vs[`i2', `i2'] == 0
local ic = colnumb(tst_bs, "/cov(F1,F2)")
display as text " 3. |var(e.F2) - (1 - beta_std^2)| = " %9.2e `d3' " (1e-12); correlated var(F1), var(F2) = 1 with SE 0: `ok3'; r_12 = " tst_bs[1, `ic']
if `d3' >= 1e-12 | !`ok3' | missing(`ic') {
    local ++nfail
}

* ---- 4. nothing stored ------------------------------------------------------------
tempname B V C
quietly dwlsmv `mS'
matrix `B' = e(b)
matrix `V' = e(V)
matrix `C' = e(Cns)
tst_snap
local snap0 `"`r(snap)'"'
quietly dwlsmv `mS', standardized
tst_snap
local s1 = mreldif(e(b), `B') == 0 & mreldif(e(V), `V') == 0 & mreldif(e(Cns), `C') == 0
local s1 = `s1' & `"`r(snap)'"' == `"`snap0'"'
capture quietly dwlsmv, standardized
local r4 = _rc
tst_snap
local s2 = mreldif(e(b), `B') == 0 & mreldif(e(V), `V') == 0 & mreldif(e(Cns), `C') == 0
local s2 = `s2' & `"`r(snap)'"' == `"`snap0'"'
quietly regress y1 y2
capture dwlsmv, standardized
local r301 = _rc
display as text " 4. e() unchanged after fit, standardized: `s1'; after replay (rc `r4'): `s2'; replay after regress rc `r301'"
if !`s1' | !`s2' | `r4' | `r301' != 301 {
    local ++nfail
}

* ---- 5. the nlcom stopgap ----------------------------------------------------------
quietly dwlsmv `mS'
quietly nlcom (a: _b[y2:F1]*sqrt(_b[/var(F1)])) ///
              (b: _b[y5:F2]*sqrt(_b[F2:F1]^2*_b[/var(F1)] + _b[/var(e.F2)]))
tempname nb nV
matrix `nb' = r(b)
matrix `nV' = r(V)
mata:
void tst_nlcom()
{
    external string matrix tst_fx
    string matrix fx
    real rowvector lv, ls, nb
    real colvector ns
    real scalar re, rs
    fx = select(tst_fx, tst_fx[., 1] :== "str_marker")
    lv = strtoreal(select(fx[., 5], (fx[., 2] :== "f1" :& fx[., 4] :== "y2") :| (fx[., 2] :== "f2" :& fx[., 4] :== "y5")))'
    ls = strtoreal(select(fx[., 6], (fx[., 2] :== "f1" :& fx[., 4] :== "y2") :| (fx[., 2] :== "f2" :& fx[., 4] :== "y5")))'
    nb = st_matrix(st_local("nb"))
    ns = sqrt(diagonal(st_matrix(st_local("nV"))))
    re = max(abs(nb - lv) :/ (1e-5 :+ 1e-4 :* abs(lv)))
    rs = max(abs(ns' - ls) :/ (1e-8 :+ 1e-5 :* abs(ls)))
    printf("{txt} 5. nlcom vs lavaan, y2 and y5: value gap %8.2e of tier, SE gap %8.2e of tier\n", re, rs)
    if (cols(lv) != 2 | re > 1 | rs > 1) {
        st_local("nfail", strofreal(strtoreal(st_local("nfail")) + 1))
    }
}
end
mata: tst_nlcom()

adopath - "`here'/src/ado"
adopath - "`here'/build"

if `nfail' {
    display as error _n "test_standardized: `nfail' failure(s)"
    exit 9
}
display as result "test_standardized: the standardized solution matches lavaan on Case E in four fits, J_g matches central differences, the 1 - R^2 and r_12 rows hold, nothing is stored, and the nlcom stopgap matches lavaan"
exit 0
