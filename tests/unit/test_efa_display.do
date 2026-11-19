*! Executed under StataMP 18 on 2026-11-18; green, and the suite 102/0/0 at 6 and at 4 processors.
*! REVISED 2026-11-18 (compact table; drafted by a read-only agent): check 7 now expects "(marker)", "Estimate" and "Std. err." instead of
*! "(constrained)" and "Coefficient"; checks 1 and 2 add "P>|z|" absent. Run 2026-11-18 in this form.
*! What broke on the first run, and the fixes:
*!   - nothresholds was ignored everywhere: with THResholds and NOTHResholds in one syntax
*!     statement, syntax takes nothresholds alone as THResholds' negation and records neither.
*!     Every parse now reads NOTHResholds in a first pass. Check 8 caught it (and now also
*!     checks nothresholds at estimation, where the default shows the rows).
*!   - check 1's phrase "Wald interval" matched the default's pointer line; now "% Wald interval)".
*!   - the long ci form padded its brackets; now compact, as the wide form.
*!   - _dwlsmv_efaci is private to dwlsmv.ado; the test runs the file to call it, and asserts rc 0.
*!
*! Release A display (owner 2026-11-17; spec Module 12 A3.6 revised; threshold rows; e(Ev); D8e):
*!   the efa() default output, options se and ci, the threshold rows (hidden by default, shown by
*!   default only when invariance() constrains them; thresholds / nothresholds), and e(Ev) with
*!   "eigen" in e(properties) (screeplot).
*!
*! End to end through build/ldwlsmv.mlib and src/ado, as test_efa_surface.do. Display checks read the
*! command's own log text (tdp_run keeps it); a phrase test says what is in the output, never what a
*! number is, except where the expected string is built here from e() with the ado's own format.
*!
*! Tolerances, fixed before the first run (2026-11-18):
*!   e(Ev) against the eigenvalues of R assembled in this file from dwls_sample_stats() (labels
*!     parsed here, symeigensystem()): max abs difference 1e-12 (the DECIDING criterion for e(Ev));
*!     sum of e(Ev) = p, 1e-12.
*!   e(Ev) against lavaan's polychoric matrix (lavaan_case_E_sampstat.csv): printed, a diagnostic,
*!     not asserted.
*!   standard-error matrices against sqrt(diag e(V)) read by colnumb(): 1e-14 relative.
*!
*! Checks:
*!   1  default efa() output: matrices, one line on e(V), thresholds note; no threshold rows, no
*!      coefficient table, no caveat
*!   2  se: SE matrices (values), caveat, no coefficient table, no threshold rows
*!   3  ci: bracketed intervals (strings built here from e(Lambda), e(V), invnormal), wide form at
*!      linesize 120; long grouped form at linesize 60; Wald phi note when an interval crosses 1
*!   4  thresholds / nothresholds on efa(); both together rc 198; se or ci without efa() rc 198
*!   5  replay with each option; level(90); e(b), e(V) untouched by display
*!   6  e(Ev): 1 x p, descending, sums to p, equals the independent eigenvalues; e(properties)
*!      is b V eigen; screeplot runs after an efa() fit and is refused after a CFA fit (rc 321);
*!      a CFA fit has no e(Ev) (efa() fits only)
*!   7  ordinary fit: thresholds hidden with the note, the marker shows as "(marker)", shown with
*!      thresholds; standardized likewise; e(b) unchanged
*!   8  invariance(): thresholds shown by default at thresholds / metric / scalar, hidden at
*!      configural and under nothresholds
*!
*! Not covered (needs a fit this file cannot construct cheaply): an efa() fit without e(V)
*! (e(rot_nose)): its "b eigen" properties and the "need e(V)" line; a covariate fit that did not
*! converge (the last-iterate table's hidden threshold rows).
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

global TDP_NF 0
global TDP_LG "`c(tmpdir)'/tdp_a7_log.txt"

* runs a command with its output in the log file $TDP_LG; returns the command's rc
capture program drop tdp_run
program tdp_run, rclass
    capture quietly log close tdplog
    quietly log using "$TDP_LG", text replace name(tdplog)
    capture noisily `0'
    local rc = _rc
    quietly log close tdplog
    return scalar rc = `rc'
end

* want (0/1): is the phrase in the log of the last tdp_run
capture program drop tdp_chk
program tdp_chk
    args want phrase label
    mata: st_local("got", strofreal(tdp_has(st_local("phrase"))))
    if `got' != `want' {
        display as error "test_efa_display: `label' (phrase [`phrase'] found `got', expected `want')"
        global TDP_NF = $TDP_NF + 1
    }
end

capture program drop tdp_ok
program tdp_ok
    args cond label
    if !(`cond') {
        display as error "test_efa_display: `label'"
        global TDP_NF = $TDP_NF + 1
    }
end

mata:
real scalar tdp_has(string scalar phrase)
{
    string colvector L

    L = cat(st_global("TDP_LG"))
    return(strpos(invtokens(L'), phrase) > 0)
}

// the p x p polychoric matrix assembled from dwls_sample_stats() with the pair labels parsed here
// (Stage 1's labels are positional, "y1~~y2"), its eigenvalues by symeigensystem(), descending
void tdp_ev(string scalar vars, string scalar outn, string scalar rn)
{
    struct dwlsmv_stage1 scalar st
    string rowvector ind
    real matrix X, Y, R, XX
    real colvector K, u, x, r
    real rowvector L
    real scalar j, k, d, p, pos, a, b
    string scalar lab

    ind = tokens(vars)
    p   = cols(ind)
    X   = st_data(., ind)
    Y   = J(rows(X), p, 0)
    K   = J(p, 1, 0)
    for (j = 1; j <= p; j++) {
        x = X[., j]
        u = uniqrows(x)
        r = J(rows(x), 1, 0)
        for (k = 2; k <= rows(u); k++) {
            r = r + (x :>= u[k])
        }
        Y[., j] = r
        K[j]    = rows(u)
    }
    st = dwls_sample_stats(Y, K)
    R  = I(p)
    for (d = 1; d <= st.map.D; d++) {
        if (st.map.stat_type[d] == 2) {
            lab = st.map.label[d]
            pos = strpos(lab, "~~")
            a   = strtoreal(substr(lab, 2, pos - 2))
            b   = strtoreal(substr(lab, pos + 3, .))
            R[a, b] = st.s[d]
            R[b, a] = st.s[d]
        }
    }
    symeigensystem(R, XX, L)
    L = sort(L', -1)'
    st_matrix(outn, L)
    st_matrix(rn, R)
}

// eigenvalues of the matrix read from a lavaan fixture (rows = the table's rows)
void tdp_lav(string scalar outn)
{
    string matrix S
    real matrix R, XX
    real rowvector L
    real scalar i, j

    S = st_sdata(., ("y1", "y2", "y3", "y4", "y5", "y6"))
    R = J(6, 6, .)
    for (i = 1; i <= 6; i++) {
        for (j = 1; j <= 6; j++) {
            R[i, j] = strtoreal(S[i, j])
        }
    }
    symeigensystem(R, XX, L)
    st_matrix(outn, sort(L', -1)')
}
end

import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")

* ---- 1  the default efa() output ----
tdp_run dwlsmv (F1 F2 -> y1-y6), efa(F1 F2)
tdp_ok "r(rc) == 0" "1: efa() fit rc"
tdp_chk 1 "Rotated loadings"        "1: loadings matrix title"
tdp_chk 1 "Factor correlations"     "1: factor correlation matrix title"
tdp_chk 1 "is posted in e(V)"       "1: the e(V) line"
tdp_chk 0 "Thresholds are not shown" "1: no thresholds-hidden note on efa() output (owner veto, 2026-11-18)"
tdp_chk 1 "They are approximate, least accurate near zero" "1: the caveat pointer on the default display (owner 2026-11-18)"
tdp_chk 0 "cut1"                    "1: no threshold row"
tdp_chk 0 "Coefficient"             "1: no coefficient table (no SE, z, p, interval columns)"
tdp_chk 0 "P>|z|"                  "1: no z or p column either (the compact table is not printed for efa())"
tdp_chk 0 "Standard errors of"      "1: no SE matrices"
tdp_chk 0 "Standard errors of the rotated estimates are a local first-order" "1: the caveat belongs to se and ci"
tdp_chk 0 "% Wald interval)"        "1: no interval table (the pointer line names the option only)"
tempname b0 V0 L0 P0
matrix `b0' = e(b)
matrix `V0' = e(V)
matrix `L0' = e(Lambda)
matrix `P0' = e(Phi)

* ---- 2  se ----
tdp_run dwlsmv, se
tdp_ok "r(rc) == 0" "2: replay with se rc"
tdp_chk 1 "Standard errors of the rotated loadings"        "2: SE matrix of the loadings"
tdp_chk 1 "Standard errors of the factor correlations"     "2: SE matrix of the correlations"
tdp_chk 1 "Standard errors of the rotated estimates are a local first-order" "2: the A3.1 caveat"
tdp_chk 0 "Coefficient"                                    "2: no coefficient table"
tdp_chk 0 "P>|z|"                                         "2: no z or p column"
tdp_chk 0 "cut1"                                           "2: no threshold row"
tdp_chk 0 "is posted in e(V)"                              "2: the e(V) line is the default's"
tdp_chk 0 "% Wald interval)"                               "2: no interval without ci"

tempname SL SP
mata: dwlsmv_m11_efa_semat("`b0'", "`V0'", "`L0'", "`SL'", "`SP'")
local worst 0
forvalues j = 1/6 {
    forvalues k = 1/2 {
        local c = colnumb(`b0', "y`j':F`k'")
        local x = sqrt(`V0'[`c', `c'])
        local worst = max(`worst', abs(`SL'[`j', `k'] - `x') / `x')
    }
}
local c = colnumb(`b0', "/:cov(F1,F2)")
local x = sqrt(`V0'[`c', `c'])
local worst = max(`worst', abs(`SP'[1, 2] - `x') / `x', abs(`SP'[2, 1] - `x') / `x')
display as text "  se: worst relative difference to sqrt(diag e(V)) " %9.2e `worst'
tdp_ok "`worst' <= 1e-14" "2: SE matrices equal sqrt(diag e(V)) read by name"
local sa : colfullnames `SL'
local sb : colfullnames `L0'
local same = ("`sa'" == "`sb'")
tdp_ok "`same'" "2: SE loadings stripes are e(Lambda)'s"
tdp_ok "missing(`SP'[1, 1]) & missing(`SP'[2, 2])" "2: the diagonal of the correlation SEs is missing (fixed at 1)"

* ---- 3  ci ----
local z = invnormal(.975)
tdp_run dwlsmv, ci
tdp_ok "r(rc) == 0" "3: replay with ci rc"
tdp_chk 1 "95% Wald interval)"   "3: the interval title"
tdp_chk 1 "Standard errors of the rotated estimates are a local first-order" "3: the A3.1 caveat"
tdp_chk 0 "Standard errors of the rotated loadings" "3: ci alone prints no SE matrices"
tdp_chk 0 "Loadings on F1"       "3: wide form at linesize 120"
local ok 1
foreach jk in "1 1" "3 2" "6 1" {
    local j : word 1 of `jk'
    local k : word 2 of `jk'
    local est = `L0'[`j', `k']
    local sev = `SL'[`j', `k']
    local cell = strtrim(strofreal(`est', "%7.3f")) + " [" + strofreal(`est' - `z' * `sev', "%7.3f") + ", " + strofreal(`est' + `z' * `sev', "%7.3f") + "]"
    tdp_chk 1 "`cell'" "3: wide cell, loading (`j', `k')"
}
local est = `P0'[2, 1]
local sev = `SP'[2, 1]
local cell = strtrim(strofreal(`est', "%7.3f")) + " [" + strofreal(`est' - `z' * `sev', "%7.3f") + ", " + strofreal(`est' + `z' * `sev', "%7.3f") + "]"
tdp_chk 1 "`cell'" "3: wide cell, factor correlation (2, 1)"

* level(90) uses the 90% quantile
local z9 = invnormal(.95)
tdp_run dwlsmv, ci level(90)
tdp_chk 1 "90% Wald interval)" "3: level(90) in the title"
local est = `L0'[1, 1]
local sev = `SL'[1, 1]
local cell = strtrim(strofreal(`est', "%7.3f")) + " [" + strofreal(`est' - `z9' * `sev', "%7.3f") + ", " + strofreal(`est' + `z9' * `sev', "%7.3f") + "]"
tdp_chk 1 "`cell'" "3: level(90) cell"

* the long grouped form when the columns do not fit the line (m = 2 at linesize 60)
set linesize 60
tdp_run dwlsmv, ci
set linesize 120
tdp_chk 1 "Loadings on F1"       "3: long form at linesize 60, group F1"
tdp_chk 1 "Loadings on F2"       "3: long form, group F2"
tdp_chk 1 "Correlations with F1" "3: long form, correlations with F1"
tdp_chk 0 "Correlations with F2" "3: long form, no group for the last factor's correlations"
local est = `L0'[2, 1]
local sev = `SL'[2, 1]
local cell = strtrim(strofreal(`est', "%9.3f")) + "   [" + strofreal(`est' - `z' * `sev', "%8.3f") + ", " + strofreal(`est' + `z' * `sev', "%8.3f") + "]"
tdp_chk 1 "`cell'" "3: long form cell, loading (2, 1)"

* the Wald interval of a correlation is printed as computed; a note when it crosses 1
tempname E S
matrix `E' = (1, . \ .95, 1)
matrix `S' = (., . \ .05, .)
matrix rownames `E' = F1 F2
matrix colnames `E' = F1 F2
matrix rownames `S' = F1 F2
matrix colnames `S' = F1 F2
* _dwlsmv_efaci is private to dwlsmv.ado; running the file defines its subprograms globally
* (dropping the autoloaded dwlsmv first; the rest of the test then runs the same file's programs)
program drop dwlsmv
quietly run src/ado/dwlsmv.ado
tdp_run _dwlsmv_efaci `E', se(`S') level(95) title("Factor correlations") phi
tdp_ok "r(rc) == 0" "3: _dwlsmv_efaci rc (crossing)"
tdp_chk 1 "crosses -1 or 1"     "3: the note when a Wald interval of a correlation crosses 1"
tdp_chk 1 "1.048"               "3: the crossing interval is printed unclipped"
matrix `S' = (., . \ .01, .)
tdp_run _dwlsmv_efaci `E', se(`S') level(95) title("Factor correlations") phi
tdp_ok "r(rc) == 0" "3: _dwlsmv_efaci rc (inside)"
tdp_chk 1 "0.950"               "3: the inside case printed its table (the negative control is not vacuous)"
tdp_chk 0 "crosses -1 or 1"     "3: no note when the interval stays inside"

* ---- 4  thresholds / nothresholds on efa(); the refusals ----
tdp_run dwlsmv, thresholds
tdp_ok "r(rc) == 0" "4: replay with thresholds rc"
tdp_chk 1 "cut1"                       "4: thresholds shows the rows"
* Stata's own table (ereturn display): "Coefficient" from Stata 17, "Coef." on 16.1 (the "Std. Err." case of
* test_condx_marker.do; changed before 0.5.0's milestone, 2026-11-18). "Coef" matches both.
tdp_chk 1 "Coef"                       "4: and as a coefficient table with e(V)"
tdp_chk 0 "Thresholds are not shown"   "4: no hidden-note when shown"
tdp_run dwlsmv, nothresholds
tdp_chk 0 "cut1"                       "4: nothresholds hides them"
tdp_chk 0 "Thresholds are not shown"   "4: no note on efa() output, even under nothresholds (owner veto)"
tdp_chk 1 "Rotated loadings"           "4: nothresholds printed the matrices"
capture dwlsmv, thresholds nothresholds
tdp_ok "_rc == 198" "4: replay thresholds with nothresholds is rc 198"
capture dwlsmv (F1 F2 -> y1-y6), efa(F1 F2) thresholds nothresholds
tdp_ok "_rc == 198" "4: fit with thresholds and nothresholds is rc 198"
capture dwlsmv (F1 -> y1-y3) (F2 -> y4-y6), se
tdp_ok "_rc == 198" "4: se without efa() is rc 198"
capture dwlsmv (F1 -> y1-y3) (F2 -> y4-y6), ci
tdp_ok "_rc == 198" "4: ci without efa() is rc 198"
* the refusals ran before any fit: the efa() result is still the one posted
local same = ("`e(efa)'" == "F1 F2")
tdp_ok "`same'" "4: nothing was fitted by a refused call (the efa() result is still posted)"

* ---- 5  replay: every option, and e() untouched ----
tdp_run dwlsmv
tdp_chk 1 "is posted in e(V)"   "5: plain replay"
tdp_chk 0 "cut1"                "5: plain replay hides thresholds"
tdp_run dwlsmv, se ci thresholds level(90)
tdp_ok "r(rc) == 0" "5: replay with every option"
tdp_chk 1 "90% Wald interval)"  "5: se ci thresholds level(90)"
tdp_chk 1 "Standard errors of the rotated loadings" "5: se with ci"
tdp_chk 1 "cut1"                "5: thresholds with ci"
tempname b1 V1
matrix `b1' = e(b)
matrix `V1' = e(V)
tdp_ok "mreldif(`b1', `b0') == 0 & mreldif(`V1', `V0') == 0" "5: e(b) and e(V) after replay equal the fit's"
local same = ("`e(cmd)'" == "dwlsmv")
tdp_ok "`same'" "5: e(cmd)"

* ---- 6  e(Ev) and screeplot ----
capture confirm matrix e(Ev)
tdp_ok "_rc == 0" "6: e(Ev) is posted"
tempname Ev Evi Ri
matrix `Ev' = e(Ev)
tdp_ok "rowsof(`Ev') == 1 & colsof(`Ev') == 6" "6: e(Ev) is 1 x p"
local cnames : colnames `Ev'
local same = ("`cnames'" == "Ev1 Ev2 Ev3 Ev4 Ev5 Ev6")
tdp_ok "`same'" "6: e(Ev) column names"
local desc 1
forvalues j = 2/6 {
    if `Ev'[1, `j'] > `Ev'[1, `j' - 1] local desc 0
}
tdp_ok "`desc'" "6: e(Ev) is descending"
local sm = 0
forvalues j = 1/6 {
    local sm = `sm' + `Ev'[1, `j']
}
display as text "  e(Ev) sum - p = " %9.2e (`sm' - 6)
tdp_ok "abs(`sm' - 6) <= 1e-12" "6: e(Ev) sums to p"
mata: tdp_ev("y1 y2 y3 y4 y5 y6", "`Evi'", "`Ri'")
local d = 0
forvalues j = 1/6 {
    local d = max(`d', abs(`Ev'[1, `j'] - `Evi'[1, `j']))
}
display as text "  e(Ev) vs independently assembled eigenvalues: max abs " %9.2e `d'
tdp_ok "`d' <= 1e-12" "6: e(Ev) equals the eigenvalues of R assembled from dwls_sample_stats() (deciding criterion)"
local same = ("`e(properties)'" == "b V eigen")
tdp_ok "`same'" "6: e(properties) is b V eigen"
capture noisily screeplot, nodraw
tdp_ok "_rc == 0" "6: screeplot runs after an efa() fit"
capture graph drop _all
* diagnostic against lavaan's polychoric matrix on Case E (printed, not asserted)
tempname Evl
preserve
import delimited using "fixtures/expected/lavaan_case_E_sampstat.csv", clear varnames(1) delimiters(",") stringcols(_all)
mata: tdp_lav("`Evl'")
restore
local dl = 0
forvalues j = 1/6 {
    local dl = max(`dl', abs(`Ev'[1, `j'] - `Evl'[1, `j']))
}
display as text "  diagnostic: e(Ev) vs lavaan polychoric eigenvalues, max abs " %9.2e `dl'

* a CFA fit: no e(Ev), screeplot refused (e(Ev) belongs to efa() fits)
quietly dwlsmv (F1 -> y1-y3) (F2 -> y4-y6)
capture confirm matrix e(Ev)
tdp_ok "_rc == 111" "6: an ordinary fit posts no e(Ev)"
capture screeplot, nodraw
tdp_ok "_rc == 321" "6: screeplot after an ordinary fit is refused by screeplot itself"
local same = ("`e(properties)'" == "b V")
tdp_ok "`same'" "6: an ordinary fit keeps b V"

* ---- 7  an ordinary fit: thresholds hidden by default ----
tempname bc Vc
matrix `bc' = e(b)
matrix `Vc' = e(V)
tdp_run dwlsmv
tdp_chk 0 "cut1"                       "7: ordinary replay hides the thresholds"
tdp_chk 1 "Thresholds are not shown"   "7: and says so"
tdp_chk 1 "(marker)"                   "7: the fixed marker loadings are labelled (marker)"
tdp_chk 0 "(constrained)"              "7: no user constraint here, so no (constrained) label"
tdp_chk 1 "Estimate"                   "7: the coefficient table is printed"
tdp_chk 1 "Std. err."                  "7: with its standard errors"
tdp_run dwlsmv, thresholds
tdp_chk 1 "cut1"                       "7: thresholds shows them"
tdp_chk 0 "Thresholds are not shown"   "7: no note when shown"
tdp_run dwlsmv, standardized
tdp_ok "r(rc) == 0" "7: standardized rc"
tdp_chk 0 "cut1"                       "7: standardized hides the thresholds"
tdp_chk 1 "Thresholds are not shown"   "7: standardized says so"
tdp_chk 1 "Standardized (std.all)"     "7: standardized heading is kept"
tdp_run dwlsmv, standardized thresholds
tdp_chk 1 "cut1"                       "7: standardized thresholds"
tdp_ok "mreldif(e(b), `bc') == 0 & mreldif(e(V), `Vc') == 0" "7: e(b) and e(V) are not changed by the display"
tdp_run dwlsmv (F1 -> y1-y3) (F2 -> y4-y6), nothresholds
tdp_chk 0 "cut1"                       "7: estimation with nothresholds"
tdp_run dwlsmv (F1 -> y1-y3) (F2 -> y4-y6), thresholds
tdp_chk 1 "cut1"                       "7: estimation with thresholds"

* ---- 8  invariance(): thresholds shown by default only when constrained ----
generate g = mod(_n, 2)
local k 0
foreach lev in configural thresholds metric scalar {
    tdp_run dwlsmv (F1 -> y1-y3) (F2 -> y4-y6), group(g) invariance(`lev') metric(unitvariance)
    tdp_ok "r(rc) == 0" "8: invariance(`lev') rc"
    local want = ("`lev'" != "configural")
    tdp_chk `want'  "cut1" "8: invariance(`lev') threshold rows by default"
    tdp_chk `=1 - `want'' "Thresholds are not shown" "8: invariance(`lev') hidden note by default"
    tdp_run dwlsmv, nothresholds
    tdp_chk 0 "cut1" "8: invariance(`lev') then nothresholds (replay)"
    tdp_run dwlsmv, thresholds
    tdp_chk 1 "cut1" "8: invariance(`lev') then thresholds (replay)"
    * nothresholds at estimation, where the default shows them (not vacuous, unlike check 7's)
    if "`lev'" != "configural" {
        tdp_run dwlsmv (F1 -> y1-y3) (F2 -> y4-y6), group(g) invariance(`lev') metric(unitvariance) nothresholds
        tdp_ok "r(rc) == 0" "8: invariance(`lev') nothresholds at estimation rc"
        tdp_chk 0 "cut1" "8: invariance(`lev') nothresholds at estimation"
        tdp_chk 1 "Thresholds are not shown" "8: invariance(`lev') nothresholds at estimation, the note"
    }
}

* ---- 9  the column-order tie note (A7 wording, owner 2026-11-18), on the exact-tie fixture ----
import delimited using "fixtures/data/efatie_T.csv", clear varnames(1) delimiters(",")
tdp_run dwlsmv (F1 F2 -> y1-y8), efa(F1 F2)
tdp_ok "r(rc) == 0" "9: Case T fit rc"
tdp_chk 1 "Factors F1 and F2 tie on the statistic that orders the factors." "9: the tie note"
tdp_chk 1 "They are ordered by the sum of squared loadings, larger first." "9: the deciding key (key 2 on Case T)"
tdp_chk 1 "Remarks on the order of the factors" "9: the pointer to the help"
import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
tdp_run dwlsmv (F1 F2 -> y1-y6), efa(F1 F2)
tdp_chk 1 "Rotated loadings" "9: Case E printed"
tdp_chk 0 "tie on the statistic" "9: no tie note without a tie (Case E)"

if $TDP_NF {
    display as error _n "test_efa_display: $TDP_NF failure(s)"
    exit 9
}
display as result "test_efa_display: efa() default output, se, ci (wide, long, unclipped phi), thresholds / nothresholds, replay, e(Ev) with screeplot, and the invariance threshold default hold"
exit 0
