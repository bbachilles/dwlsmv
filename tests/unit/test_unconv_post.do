*! Executed under StataMP 18 on 2026-11-18: green on the first run (the real fit: e(converged) 0 after
*! 301 iterations, e(recovery) ineligible). Check 4's rc was then fixed at the observed 321.
*!
*! Non-convergence is an inference boundary (owner 2026-11-18; spec Module 11, the paragraph of that
*! name): a fit that does not converge is posted for inspection, with no inference, whatever the fit.
*!
*! A REAL non-converged fit (deciding): the marker-trap fixture with a second loading fixed at 1, the
*! model of test_marker_recovery.do check 6. It is ineligible for D11's recovery, and its direct fit does
*! not converge, so the posted fit is the last iterate. A premise that no longer holds (the fit converges)
*! is a FAILURE of this test (the data no longer provoke the case), never a skip.
*!   1  posted: e(b), e(Cns_b), e(N), e(rank), e(df_T3), e(converged) = 0, e(iterations), the Module 10
*!      ratios, e(bound), e(canon), e(cmdline), e(recovery) ineligible. Not posted: e(V), e(Cns), e(T3),
*!      e(p_T3), e(rmsea_s), e(cfi_s), e(tli_s), e(T3_b), e(df_T3_b).
*!   2  the output: the error-colour lines under the title, the T3 line "not reported (not a solution)",
*!      the D11 remedy, the last-iterate table; no Std. err., no Prob > chi2, no fit index
*!   3  replay and estimates store / restore say the same
*!   4  lincom and test cannot run (no e(V)); rc as observed, asserted exactly
*!   5  standardized is refused, rc 198, at estimation (the fit stays posted) and on replay
*!   6  thresholds shows the threshold rows of the last-iterate table
*!   7  CONTROL: the same data under metric(unitvariance) converge and post e(V), T3, p and the indices
*! Hand-set (display only; no real input is known, spec: no iteration-limit option):
*!   8  a G > 1 fit and an efa() fit with e(converged) set to 0: the header withholds T3 and the indices,
*!      the G > 1 table is the last-iterate table, the efa() display shows its matrices
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

global TUP_NF 0
global TUP_LG "`c(tmpdir)'/tup_log.txt"

capture program drop tup_run
program tup_run, rclass
    capture quietly log close tuplog
    quietly log using "$TUP_LG", text replace name(tuplog)
    capture noisily `0'
    local rc = _rc
    quietly log close tuplog
    return scalar rc = `rc'
end

* want (0/1): is the phrase in the log of the last tup_run
capture program drop tup_chk
program tup_chk
    args want phrase label
    mata: st_local("got", strofreal(strpos(invtokens(cat(st_global("TUP_LG"))', " "), st_local("phrase")) > 0))
    if `got' != `want' {
        display as error "test_unconv_post: `label' (phrase [`phrase'] found `got', expected `want')"
        global TUP_NF = $TUP_NF + 1
    }
end

capture program drop tup_ok
program tup_ok
    args cond label
    if !(`cond') {
        display as error "test_unconv_post: `label'"
        global TUP_NF = $TUP_NF + 1
    }
end

* is e(name) posted, as a scalar or macro (1) or not (0)
capture program drop tup_has
program tup_has, rclass
    args nm
    return scalar has = ("`e(`nm')'" != "")
end

* sets e(converged) to 0 on the posted results, as test_dwlsmv_gaps.do G10 does
capture program drop tup_unconv
program tup_unconv, eclass
    ereturn scalar converged = 0
end

* the checks of 2 and 3 on the log of the last tup_run
capture program drop tup_out
program tup_out
    args tag
    tup_chk 1 "convergence not achieved: the optimizer stopped after" "`tag': the first error line"
    tup_chk 1 "The estimates below are the last iterate, not a solution." "`tag': last iterate"
    tup_chk 1 "No standard errors, T3 or p-value are reported for them." "`tag': nothing reported"
    tup_chk 1 "Scaled-and-shifted test (T3): not reported (not a solution); df =" "`tag': the T3 line"
    tup_chk 1 "Last iterate, not a solution; no standard errors" "`tag': the last-iterate table"
    tup_chk 0 "Std. err."     "`tag': no standard-error column"
    tup_chk 0 "Prob > chi2"   "`tag': no p-value"
    tup_chk 0 "RMSEA (scaled)" "`tag': no fit index"
end

import delimited using "fixtures/data/marker_trap.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
quietly destring q* att*, replace
local mI "(MathAb -> q1@1 q2@1 q3-q8) (MathAtt -> att1-att5) (MathAtt -> MathAb)"

* ---- 1  the real non-converged fit: what is posted ----
tup_run dwlsmv `mI'
tup_ok "r(rc) == 0" "1: the fit runs, rc 0"
display as text "  1: e(converged) " e(converged) ", e(iterations) " e(iterations) ", e(recovery) `e(recovery)'"
tup_ok "e(converged) == 0" "1: PREMISE: the fit does not converge (the data no longer provoke the case)"
local same = ("`e(recovery)'" == "ineligible")
tup_ok "`same'" "1: e(recovery) ineligible"
foreach m in b Cns_b bound {
    capture confirm matrix e(`m')
    tup_ok "_rc == 0" "1: e(`m') is posted"
}
foreach m in V Cns {
    capture confirm matrix e(`m')
    tup_ok "_rc == 111" "1: e(`m') is not posted"
}
foreach s in N rank df_T3 converged iterations gamma_ratio gamma_ratio_raw gamma_rank bread_ratio {
    tup_has `s'
    tup_ok "r(has) == 1" "1: e(`s') is posted"
}
foreach s in T3 p_T3 rmsea_s cfi_s tli_s T3_b df_T3_b {
    tup_has `s'
    tup_ok "r(has) == 0" "1: e(`s') is not posted"
}
foreach s in canon cmdline indicators metric cmd {
    tup_has `s'
    tup_ok "r(has) == 1" "1: e(`s') is posted"
}
tempname b0
matrix `b0' = e(b)

* ---- 2  the output ----
tup_out "2"
tup_chk 1 "Refit with metric(unitvariance). No recovery refit was attempted:" "2: the D11 remedy"
tup_chk 0 "dwlsmv has no option for other starting values" "2: the unit-variance remedy is not the marker fit's"

* ---- 3  replay, and estimates store / restore ----
tup_run dwlsmv
tup_ok "r(rc) == 0" "3: replay rc"
tup_out "3 replay"
estimates store tupI
quietly dwlsmv (MathAb -> q1-q8) (MathAtt -> att1-att5) (MathAtt -> MathAb), metric(unitvariance)
estimates restore tupI
tup_run dwlsmv
tup_out "3 after estimates restore"
capture confirm matrix e(V)
tup_ok "_rc == 111" "3: no e(V) after estimates restore"
tup_ok "mreldif(e(b), `b0') == 0" "3: e(b) after estimates restore"

* ---- 4  lincom and test cannot run ----
local c1 : word 1 of `: colfullnames e(b)'
capture lincom _b[`c1']
local rcl = _rc
capture test _b[`c1'] = 0
local rct = _rc
display as text "  4: lincom rc `rcl', test rc `rct' (no e(V))"
* rc 321 (Stata's own, without e(V)): observed on the first run, 2026-11-18, then asserted
tup_ok "`rcl' == 321" "4: lincom refused, rc `rcl'"
tup_ok "`rct' == 321" "4: test refused, rc `rct'"

* ---- 5  standardized is refused ----
tup_run dwlsmv `mI', standardized
tup_ok "r(rc) == 198" "5: standardized at estimation is rc 198"
tup_chk 1 "standardized is not available for a fit that did not converge" "5: the message at estimation"
local same = ("`e(cmd)'" == "dwlsmv")
tup_ok "`same' & e(converged) == 0" "5: the fit stays posted after the refusal"
tup_run dwlsmv, standardized
tup_ok "r(rc) == 198" "5: standardized on replay is rc 198"
tup_chk 1 "standardized is not available for a fit that did not converge" "5: the message on replay"
tup_chk 0 "DWLS estimation" "5: the replay refusal comes before the header"

* ---- 6  thresholds on the last-iterate table ----
tup_run dwlsmv
tup_chk 1 "Last iterate, not a solution" "6: the last-iterate table is printed"
tup_chk 0 "cut1" "6: no threshold rows by default"
tup_run dwlsmv, thresholds
tup_chk 1 "cut1" "6: thresholds shows them"

* ---- 7  CONTROL: a converged fit of the same data ----
tup_run dwlsmv (MathAb -> q1-q8) (MathAtt -> att1-att5) (MathAtt -> MathAb), metric(unitvariance)
tup_ok "r(rc) == 0 & e(converged) == 1" "7: CONTROL converges"
capture confirm matrix e(V)
tup_ok "_rc == 0" "7: CONTROL posts e(V)"
foreach s in T3 p_T3 rmsea_s cfi_s tli_s T3_b df_T3_b {
    tup_has `s'
    tup_ok "r(has) == 1" "7: CONTROL posts e(`s')"
}
tup_chk 1 "Prob > chi2" "7: CONTROL prints p"
tup_chk 0 "not reported (not a solution)" "7: CONTROL has no non-convergence line"

* ---- 8  hand-set: G > 1 and efa() ----
import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
generate g = mod(_n, 2)
quietly dwlsmv (F1 -> y1-y3) (F2 -> y4-y6), group(g) invariance(metric) metric(unitvariance)
tup_unconv
tup_run dwlsmv
tup_ok "r(rc) == 0" "8: G > 1 hand-set replay rc"
tup_out "8 G > 1"
tup_chk 1 "cut1" "8: G > 1 last-iterate table shows the thresholds by default under invariance(metric)"
quietly dwlsmv (F1 F2 -> y1-y6), efa(F1 F2)
tup_unconv
tup_run dwlsmv
tup_ok "r(rc) == 0" "8: efa() hand-set replay rc"
tup_chk 1 "convergence not achieved: the optimizer stopped after" "8 efa(): the first error line"
tup_chk 1 "not reported (not a solution)" "8 efa(): the T3 line"
tup_chk 0 "Prob > chi2" "8 efa(): no p-value"
tup_chk 1 "Rotated loadings" "8 efa(): the matrices"
tup_chk 0 "is posted in e(V)" "8 efa(): the display does not offer e(V) for a fit that did not converge"

if $TUP_NF {
    display as error _n "test_unconv_post: $TUP_NF failure(s)"
    exit 9
}
display as result "test_unconv_post: a fit that did not converge is posted for inspection with no e(V), T3, p or fit indices; display, replay, restore, lincom/test and standardized follow"
