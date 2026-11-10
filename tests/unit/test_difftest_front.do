*! difftest, the front end (Module 12; single-group comparisons, 2026-10-23).
*! Executed 2026-10-23 under StataMP 16.1, 4 processors: passes (step 1
*! T_d3 506.93462, step 2 0.04993383). Mutants D1-D14 (applied in place, each
*! checked applied, restored; difftest.ado and dwlsmv_dtfront.mata): check 6
*! skipped, the version/edition line removed, checks 1, 2, 3, 4, 5, the stored
*! check 7, the G > 1 refusal, the U6 branch, CFI/RMSEA swapped in the delta,
*! r() not returned, duplicate names allowed: all red. D12 (active results not
*! restored) survived the first version of check 3, whose last restored model
*! was the active one; the second difftest call in check 3 was added and it is
*! red. Refusals assert their own message, not only rc 3498, so a check removed
*! cannot pass by a later check refusing.
*!
*! Data fixtures/data/case_E.csv; reference fixtures/expected/
*! lavaan_case_E_difftest.csv (97_lavaan_difftest.R, lavaan 0.7.2, rows
*! variant "reference"). Models, metric(unitvariance) throughout (the
*! cross-loading model does not converge under the marker metric):
*!   h0  (F1 -> y1 y2 y3) (F2 -> y4 y5 y6), cov(F1*F2@0)   lavaan beta0 h0
*!   h1  (F1 -> y1 y2 y3) (F2 -> y4 y5 y6)                 beta0 h1 = xload2 h0
*!   h2  (F1 -> y1 y2 y3 y4) (F2 -> y4 y5 y6 y3)           xload2 h1
*! (beta0 is written with a covariance where lavaan has a path: the same
*! model, and the statistic does not depend on the parameterization.)
*!
*! Checks:
*!   1  pair: difftest h0 h1: r(steps) d = 1, T_d3 and p against lavaan's
*!      beta0 statistic at spec 7.2's test-statistic tier (1e-8 abs + 1e-5
*!      rel, chosen before the first run), status 1; r(models) T3 and df
*!      equal the stored e(T3), e(df_T3); the printed table rows appear.
*!   2  chain: difftest h2 h0 h1 (named out of order): order h0, h1, h2;
*!      step 2 against lavaan's xload2 statistic at the same tier; the dCFI
*!      column equals the CFI differences of r(models), H0 minus H1.
*!   3  one name: difftest h0 with h2 active: d = 3; the active results are
*!      h2's afterwards (restored).
*!   4  U6: H0 fixes cov(F1,F2) at H1's estimate (21 digits): status 2,
*!      T_d3 and p missing, "within eps; not tested" printed.
*!   7  G > 1 (2026-10-31): (a) e(groupvar), e(datasignature), e(keepbound)
*!      posted; (b) the four generated levels of 99e's design E as a chain,
*!      named out of order: d, T_d3 and p per step against the fixture's
*!      (N - G) statistic at check 1's tier, status 1, r(models) order;
*!      (c) same-session gap (user): every refit bitwise its stored e(b),
*!      e(V), e(T3), e(df_T3); (d) refusals, rc 3498: G > 1 against G = 1,
*!      another group variable, check 1, check 6 (e(b) altered; hint line),
*!      check 2 (one group value swapped).
*!      Mutants 2026-10-31 (StataMP 18, %TEMP% copies): exemption inputs
*!      dropped (7b red); check 4's group clause removed (7d red); c = N - 1
*!      SURVIVES, equivalent: c enters only check 10's eps (DT-M1).
*!   8  sign orientation (user, 2026-10-31): dwlsmv_dtfront_delta on chosen
*!      indices returns the restricted minus the less restricted value,
*!      whatever its sign: CFI .940 vs .960, RMSEA .060 vs .045 give -.020,
*!      +.015; CFI .955 vs .950, RMSEA .040 vs .050 give +.005, -.010; a
*!      three-model chain with its middle step skipped (F6); malformed
*!      input rc 3498.
*!   5  refusals, each _rc == 3498: a name twice; an unknown name; a regress
*!      result; check 1, e(version) altered;
*!      check 1, e(canon) removed; check 2, a data value changed; check 3, a
*!      different sample; check 4, other indicators, other metric, keepbound
*!      on one model; check 5, the same indicators in another order; check 6,
*!      e(b) altered (and the Stata version/edition line printed); check 7,
*!      e(converged) set to 0; check 8, a non-nested pair; check 9, three
*!      models that do not form a chain.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

quietly do tools/build_mlib.do
quietly do tests/dwlsmv_assert.mata
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

import delimited using "fixtures/expected/lavaan_case_E_difftest.csv", ///
    clear varnames(1) stringcols(_all) delimiters(",")
mata: DTX = st_sdata(., .)
mata: DTV = st_varname(1..st_nvar())

local nfail 0

mata:
real scalar tdf_ref(string scalar pair, string scalar col)
{
    external string matrix DTX
    external string rowvector DTV
    real colvector k
    real scalar j

    k = selectindex((DTX[., 1] :== pair) :& (DTX[., 2] :== "reference"))
    j = selectindex(DTV :== col)
    if (rows(k) != 1 | cols(j) != 1) {
        _error("test_difftest_front: fixture cell not found")
    }
    return(strtoreal(DTX[k, j]))
}
end

capture program drop tdf_set
program tdf_set, eclass
    args what val
    if "`what'" == "version"   ereturn local version "`val'"
    if "`what'" == "canon"     ereturn local canon ""
    if "`what'" == "converged" ereturn scalar converged = `val'
    if "`what'" == "sig"       ereturn local datasignature ""
    if "`what'" == "recovery"  ereturn local recovery "`val'"
    if "`what'" == "b" {
        tempname b
        matrix `b' = e(b)
        matrix `b'[1, 1] = `b'[1, 1] * (1 + `val')
        ereturn repost b = `b'
    }
end

capture program drop tdf_log
program tdf_log
    args cmd
    tempfile lg
    capture log close tdf
    quietly log using "`lg'", text replace name(tdf)
    capture noisily `cmd'
    local rc = _rc
    capture matrix drop TDS TDM
    capture matrix TDS = r(steps)
    capture matrix TDM = r(models)
    quietly log close tdf
    mata: TDFLOG = invtokens(cat(st_local("lg"))')
    c_local tdf_rc `rc'
end

capture program drop tdf_has
program tdf_has
    args s
    mata: st_local("h", strofreal(strpos(TDFLOG, st_local("s")) > 0))
    c_local hit `h'
end

capture program drop tdf_refuse
program tdf_refuse
    args lbl cmd msg
    tdf_log "`cmd'"
    tdf_has "`msg'"
    if `tdf_rc' != 3498 | !`hit' {
        display as error "check 5 FAILED (`lbl'): rc `tdf_rc'; message '`msg'' printed: `hit'"
        c_local bad 1
    }
    else c_local bad 0
end

import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
local mA (F1 -> y1 y2 y3) (F2 -> y4 y5 y6)
quietly dwlsmv `mA', cov(F1*F2@0) metric(unitvariance)
estimates store h0
quietly dwlsmv `mA', metric(unitvariance)
estimates store h1
local v2 = string(_b[/cov(F1,F2)], "%21.0g")
quietly dwlsmv (F1 -> y1 y2 y3 y4) (F2 -> y4 y5 y6 y3), metric(unitvariance)
estimates store h2

* ---- check 1 ---------------------------------------------------------------
tdf_log "difftest h0 h1"
matrix S = TDS
matrix Mo = TDM
mata: st_local("r1", strofreal(dwlsmv_gapratio(st_matrix("S")[1, 2], tdf_ref("beta0", "stat"), 1e-8, 1e-5)))
mata: st_local("r2", strofreal(dwlsmv_gapratio(st_matrix("S")[1, 3], tdf_ref("beta0", "p"), 1e-12, 1e-5)))
if `tdf_rc' | S[1, 1] != 1 | S[1, 6] != 1 | !(`r1' <= 1) | !(`r2' <= 1) {
    display as error "check 1 FAILED: rc `tdf_rc', d " S[1, 1] ", T_d3 " S[1, 2] " (ratio `r1'), p ratio `r2'"
    local ++nfail
}
quietly estimates restore h0
local e0T = e(T3)
local e0d = e(df_T3)
if reldif(Mo[1, 2], `e0T') > 1e-12 | Mo[1, 3] != `e0d' {
    display as error "check 1 FAILED: r(models) row 1 is not h0's stored T3/df"
    local ++nfail
}
foreach s in "Models, most restricted first" "h0 -> h1" "No published dCFI or dRMSEA cutoff" {
    tdf_has "`s'"
    if !`hit' {
        display as error "check 1 FAILED: '`s'' not printed"
        local ++nfail
    }
}
display "check 1 done: T_d3 " S[1, 2]

* ---- check 1b: lavaan's own beta0 form, the path fixed at 0 (possible since
* the 2026-10-23 fix; the active results are h2's again afterwards) ----------
quietly dwlsmv `mA' (F2 <- F1@0), metric(unitvariance)
estimates store p0
quietly dwlsmv `mA' (F2 <- F1), metric(unitvariance)
estimates store p1
quietly difftest p0 p1
matrix S = r(steps)
mata: st_local("r1", strofreal(dwlsmv_gapratio(st_matrix("S")[1, 2], tdf_ref("beta0", "stat"), 1e-8, 1e-5)))
if S[1, 1] != 1 | !(`r1' <= 1) {
    display as error "check 1b FAILED: d " S[1, 1] ", T_d3 " S[1, 2] " (ratio `r1')"
    local ++nfail
}
quietly estimates restore h2
display "check 1b done: T_d3 " S[1, 2]

* ---- check 2 ---------------------------------------------------------------
quietly difftest h2 h0 h1
matrix S = r(steps)
matrix Mo = r(models)
local rn : rownames Mo
* row names are documented (help, Stored results, 2026-10-31): models by
* stored name, steps as H0:H1
local rs : rowfullnames S
mata: st_local("r1", strofreal(dwlsmv_gapratio(st_matrix("S")[2, 2], tdf_ref("xload2", "stat"), 1e-8, 1e-5)))
if "`rn'" != "h0 h1 h2" | "`rs'" != "h0:h1 h1:h2" | S[2, 1] != 2 | !(`r1' <= 1) {
    display as error "check 2 FAILED: order `rn', steps `rs', d " S[2, 1] ", T_d3 " S[2, 2] " (ratio `r1')"
    local ++nfail
}
if reldif(S[1, 4], Mo[1, 5] - Mo[2, 5]) > 1e-12 | reldif(S[2, 5], Mo[2, 4] - Mo[3, 4]) > 1e-12 {
    display as error "check 2 FAILED: dCFI/dRMSEA are not the differences of r(models)"
    local ++nfail
}
display "check 2 done: T_d3 " S[2, 2]

* ---- check 3 ---------------------------------------------------------------
quietly estimates restore h2
scalar tdf_bh2 = _b[y4:F1]
quietly difftest h0
matrix S = r(steps)
if S[1, 1] != 3 | _b[y4:F1] != tdf_bh2 {
    display as error "check 3 FAILED: d " S[1, 1] " or the active results were not restored"
    local ++nfail
}
* the last model restored inside (h1) is not the active one (h2): mutant
* D12 (no restore) survived the line above, where it was
quietly difftest h0 h1
scalar tdf_y4 = .
capture scalar tdf_y4 = _b[y4:F1]
if tdf_y4 != tdf_bh2 {
    display as error "check 3 FAILED: the active results were not restored after difftest h0 h1"
    local ++nfail
}
display "check 3 done"

* ---- check 4: U6 ------------------------------------------------------------
quietly dwlsmv `mA', cov(F1*F2@`v2') metric(unitvariance)
estimates store hfix
tdf_log "difftest hfix h1"
matrix S = TDS
tdf_has "within eps; not tested"
if `tdf_rc' | S[1, 6] != 2 | !missing(S[1, 2]) | !missing(S[1, 3]) | !`hit' {
    display as error "check 4 FAILED: rc `tdf_rc', status " S[1, 6] ", T_d3 " S[1, 2]
    local ++nfail
}
* spec F6: only a skipped d = 0 row has its deltas blank; this row has d = 1
if missing(S[1, 4]) | missing(S[1, 5]) {
    display as error "check 4 FAILED: dCFI/dRMSEA blank on a status-2 row (DT-M11)"
    local ++nfail
}
display "check 4 done"

* ---- check 6: survivor-ledger coverage (E-076) ---------
* r(N) and the models table's indices against the stored fit (DT-M14, DT-M2)
quietly difftest h0 h1
matrix Mo = r(models)
local rN = r(N)
quietly estimates restore h1
if `rN' != e(N) {
    display as error "check 6 FAILED: r(N) `rN' is not e(N) " e(N) " (DT-M14)"
    local ++nfail
}
* h1 sits at the clamps (RMSEA 0, CFI 1), so h0's row is the one that can move
local i 0
foreach m in h0 h1 {
    local ++i
    quietly estimates restore `m'
    local j 4
    foreach e in rmsea_s cfi_s tli_s {
        if reldif(Mo[`i', `j'], e(`e')) > 1e-8 {
            display as error "check 6 FAILED: r(models) `m' col `j' " Mo[`i', `j'] " vs e(`e') " e(`e') " (DT-M2)"
            local ++nfail
        }
        local ++j
    }
}
* the active results' name is unchanged, stored and unstored (DT-M22)
quietly estimates restore h2
quietly difftest h0
if "`e(_estimates_name)'" != "h2" {
    display as error "check 6 FAILED: e(_estimates_name) '`e(_estimates_name)'' after difftest h0, not h2 (DT-M22)"
    local ++nfail
}
quietly dwlsmv `mA', metric(unitvariance)
quietly difftest h0
if "`e(_estimates_name)'" != "" {
    display as error "check 6 FAILED: unstored active results now named '`e(_estimates_name)'' (DT-M22)"
    local ++nfail
}
* keepbound with a real boundary pair (DT-M15): y7 = y1 puts y1~~y7 at the
* bound, which Stage 1 refuses without keepbound; check 6 must reproduce
* both stored fits over the shared keepbound Stage 1 (no lavaan reference:
* lavaan has no keepbound counterpart, so nothing beyond check 6 is claimed)
preserve
quietly generate byte y7 = y1
local mK (F1 -> y1 y2 y3 y7) (F2 -> y4 y5 y6)
quietly dwlsmv `mK', cov(F1*F2@0) metric(unitvariance) keepbound
estimates store hk0
quietly dwlsmv `mK', metric(unitvariance) keepbound
estimates store hk1
tdf_log "difftest hk0 hk1"
matrix S = TDS
if `tdf_rc' | S[1, 1] != 1 | S[1, 6] != 1 {
    display as error "check 6 FAILED: keepbound boundary pair, rc `tdf_rc' (DT-M15)"
    local ++nfail
}
restore
quietly estimates restore h2
display "check 6 done"

* ---- check 5: refusals -------------------------------------------------------
tdf_refuse "name twice" "difftest h0 h0" "is named twice"
local nfail = `nfail' + `bad'
tdf_refuse "unknown name" "difftest h0 nosuch" "no stored estimation results named"
local nfail = `nfail' + `bad'

preserve
quietly regress y1 y2 y3
estimates store rg
restore
tdf_refuse "regress" "difftest h0 rg" "was not fitted by dwlsmv"
local nfail = `nfail' + `bad'

* check 1: version, canon
foreach w in version canon {
    quietly estimates restore h1
    tdf_set `w' 0.0.0
    estimates store h1`w'
    tdf_refuse "check 1 `w'" "difftest h0 h1`w'" "difftest check 1"
    local nfail = `nfail' + `bad'
}
* check 7: converged
quietly estimates restore h1
tdf_set converged 0
estimates store h1nc
tdf_refuse "check 7" "difftest h0 h1nc" "difftest check 7"
local nfail = `nfail' + `bad'
* check 1: no data signature (DT-M16); check 7: recovery mismatch (DT-M6)
quietly estimates restore h1
tdf_set sig
estimates store h1sig
tdf_refuse "check 1 sig" "difftest h0 h1sig" "difftest check 1"
local nfail = `nfail' + `bad'
quietly estimates restore h1
tdf_set recovery recovered
estimates store h1rec
tdf_refuse "check 7 recovery" "difftest h0 h1rec" "did not reproduce the recovery"
local nfail = `nfail' + `bad'
* check 6: e(b), and the version/edition line
quietly estimates restore h1
tdf_set b 1e-6
estimates store h1b
tdf_log "difftest h0 h1b"
tdf_has "different version or edition of Stata"
if `tdf_rc' != 3498 | !`hit' {
    display as error "check 5 FAILED (check 6): rc `tdf_rc', hint line `hit'"
    local ++nfail
}
* check 3: another sample
quietly dwlsmv `mA' if _n <= 4000, metric(unitvariance)
estimates store hs
tdf_refuse "check 3" "difftest h0 hs" "difftest check 3"
local nfail = `nfail' + `bad'
* check 4: indicators, metric, keepbound
quietly dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5), metric(unitvariance)
estimates store hi
tdf_refuse "check 4 indicators" "difftest h0 hi" "different indicators"
local nfail = `nfail' + `bad'
quietly dwlsmv `mA'
estimates store hm
tdf_refuse "check 4 metric" "difftest h0 hm" "different metrics"
local nfail = `nfail' + `bad'
quietly dwlsmv `mA', metric(unitvariance) keepbound
estimates store hk
tdf_refuse "check 4 keepbound" "difftest h0 hk" "keepbound was specified"
local nfail = `nfail' + `bad'
* check 5: order
quietly dwlsmv (F2 -> y4 y5 y6) (F1 -> y1 y2 y3), metric(unitvariance)
estimates store ho
tdf_refuse "check 5 order" "difftest h0 ho" "difftest check 5"
local nfail = `nfail' + `bad'
* check 8: not nested
quietly dwlsmv (F1 -> y1 y2 y3 y4) (F2 -> y5 y6), metric(unitvariance)
estimates store hn
tdf_refuse "check 8" "difftest h1 hn" "dwlsmv check 8"
local nfail = `nfail' + `bad'
* check 9: not a chain. dwlsmv_dtchain orders by npar and refuses through
* check 8 on the adjacent step, so the message names check 8 and the step
* (h1, hn), which a two-model call cannot produce; its own _error()s guard
* inputs the front end always builds complete
tdf_refuse "check 9" "difftest h0 h1 hn" "dwlsmv check 8 (h1 nested in hn)"
local nfail = `nfail' + `bad'
* check 2: a data value changed (last: it alters the data)
replace y1 = 1 - y1 in 1
tdf_refuse "check 2" "difftest h0 h1" "difftest check 2"
local nfail = `nfail' + `bad'

display "check 5 done"

* ---- check 7: G > 1 (admitted 2026-10-31, narrow claim; E-095) --------
* Design E of fixture 99e (lavaan 0.7.2, two groups, N 1200/2800), its four
* generated levels. Tolerances are check 1's, spec 7.2's test-statistic tier,
* fixed before the first G > 1 run: T_d3 1e-8 abs + 1e-5 rel against the
* fixture's (N - G) statistic T_d3_NG, p 1e-12 abs + 1e-5 rel against p_NG.
import delimited using "fixtures/expected/lavaan_mg_difftest.csv", ///
    clear varnames(1) case(lower) delimiters(",") stringcols(_all)
mata: MGK = st_sdata(., ("design", "h1", "h0"))
mata: MGV = strtoreal(st_sdata(., ("d", "t_d3_ng", "p_ng")))
import delimited using "fixtures/data/mg_E.csv", clear varnames(1) delimiters(",")
foreach L in configural thresholds metric scalar {
    quietly dwlsmv (F -> y1 y2 y3 y4 y5), group(g) invariance(`L') metric(unitvariance)
    estimates store g_`L'
}
* 7a  what a G > 1 fit posts for difftest's refit
quietly estimates restore g_metric
if "`e(groupvar)'" != "g" | "`e(datasignature)'" == "" | "`e(keepbound)'" != "" {
    display as error "check 7a FAILED: e(groupvar) '`e(groupvar)'', e(datasignature) '`e(datasignature)''"
    local ++nfail
}
* 7b  the chain, named out of order: steps scalar-metric, metric-thresholds,
*     thresholds-configural, each against 99e's E row
tdf_log "difftest g_metric g_configural g_scalar g_thresholds"
if `tdf_rc' {
    display as error "check 7b FAILED: the G > 1 chain was refused, rc `tdf_rc'"
    local ++nfail
    matrix TDS = J(3, 6, .)
    matrix TDM = J(4, 6, .)
}
matrix S = TDS
matrix GM = TDM
local h0s scalar metric thresholds
local h1s metric thresholds configural
forvalues i = 1/3 {
    local a : word `i' of `h0s'
    local w : word `i' of `h1s'
    mata: k = selectindex(MGK[., 1] :== "E" :& MGK[., 2] :== "`w'" :& MGK[., 3] :== "`a'")
    mata: st_local("dref", strofreal(MGV[k, 1]))
    mata: st_local("r1", strofreal(dwlsmv_gapratio(st_matrix("S")[`i', 2], MGV[k, 2], 1e-8, 1e-5)))
    mata: st_local("r2", strofreal(dwlsmv_gapratio(st_matrix("S")[`i', 3], MGV[k, 3], 1e-12, 1e-5)))
    if `tdf_rc' | S[`i', 1] != `dref' | S[`i', 6] != 1 | !(`r1' <= 1) | !(`r2' <= 1) {
        display as error "check 7b FAILED (`a' -> `w'): rc `tdf_rc', d " S[`i', 1] " (fixture `dref'), T_d3 ratio `r1', p ratio `r2'"
        local ++nfail
    }
}
if `"`: rowfullnames TDM'"' != "g_scalar g_metric g_thresholds g_configural" {
    display as error "check 7b FAILED: r(models) rows `: rowfullnames TDM'"
    local ++nfail
}
* 7c  same-session gap (user, 2026-10-31): each refit in this session
*     reproduces its stored fit bit for bit, e(b), e(V), e(T3), e(df_T3).
*     c6-repeat certified separate processes; this is the implementation path
*     difftest itself takes
mata:
void tdf_gap(real scalar nexp)
{
    struct dwlsmv_dtm colvector M
    pointer() matrix R
    real scalar k, bad

    M   = *findexternal("dwlsmv_dtfront_M")
    R   = *findexternal("dwlsmv_dtfront_R")
    bad = 0
    for (k = 1; k <= rows(M); k++) {
        if (M[k].b != *R[1, k] | M[k].V != *R[2, k] | M[k].T3 != *R[3, k] | M[k].df != *R[5, k]) {
            printf("{err}check 7c FAILED: %s's same-session refit is not bitwise its stored fit (e(b) %g)\n",
                   M[k].name, mreldif(M[k].b, *R[1, k]))
            bad = bad + 1
        }
    }
    if (rows(M) != nexp) {
        bad = bad + 1
    }
    st_local("gapbad", strofreal(bad))
}
end
mata: tdf_gap(4)
local nfail = `nfail' + `gapbad'
* 7d  refusals at G > 1, each _rc == 3498
quietly dwlsmv (F -> y1 y2 y3 y4 y5), metric(unitvariance)
estimates store g_single
tdf_refuse "G > 1 against G = 1" "difftest g_single g_configural" "one model is a multiple-group fit"
local nfail = `nfail' + `bad'
generate g2 = g
quietly dwlsmv (F -> y1 y2 y3 y4 y5), group(g2) invariance(configural) metric(unitvariance)
estimates store g_other
tdf_refuse "G > 1 group variable" "difftest g_metric g_other" "different group variables"
local nfail = `nfail' + `bad'
quietly estimates restore g_metric
tdf_set version 0.0.0
estimates store g_mv
tdf_refuse "G > 1 check 1" "difftest g_mv g_configural" "difftest check 1"
local nfail = `nfail' + `bad'
quietly estimates restore g_metric
tdf_set b 1e-6
estimates store g_mb
tdf_log "difftest g_mb g_configural"
tdf_has "different version or edition of Stata"
if `tdf_rc' != 3498 | !`hit' {
    display as error "check 7d FAILED (G > 1 check 6): rc `tdf_rc', hint line `hit'"
    local ++nfail
}
* check 2: one group value swapped keeps the sample and changes the data
* signature, which covers the group variable at G > 1 (last: alters the data)
replace g = 3 - g in 1
tdf_refuse "G > 1 check 2" "difftest g_metric g_configural" "difftest check 2"
local nfail = `nfail' + `bad'
* 7e  (review gap 1) r(models)' fit indices over all groups equal each stored
*     fit's e(rmsea_s), e(cfi_s), e(tli_s); columns 4-6, rows most restricted first
local i 0
foreach L in scalar metric thresholds configural {
    local ++i
    quietly estimates restore g_`L'
    foreach j in 4 5 6 {
        local e : word `= `j' - 3' of rmsea_s cfi_s tli_s
        if !((missing(GM[`i', `j']) & missing(e(`e'))) | reldif(GM[`i', `j'], e(`e')) <= 1e-12) {
            display as error "check 7e FAILED (`L'): r(models) col `j' " GM[`i', `j'] " vs e(`e') " e(`e')
            local ++nfail
        }
    }
}
* 7f  (review gap 2) a group value outside the estimation sample: difftest
*     reads the group values in e(sample) only
import delimited using "fixtures/data/mg_E.csv", clear varnames(1) delimiters(",")
generate gx = g
replace gx = 9 in 1/40
foreach L in configural metric {
    quietly dwlsmv (F -> y1 y2 y3 y4 y5) if gx != 9, group(gx) invariance(`L') metric(unitvariance)
    estimates store gx_`L'
}
tdf_log "difftest gx_metric gx_configural"
if `tdf_rc' | TDS[1, 6] != 1 {
    display as error "check 7f FAILED: a group value outside e(sample), rc `tdf_rc'"
    local ++nfail
}
* 7g  (review gap 3) two factors: the factor count read from e(canon) sets
*     the alpha rows of the refit's e(b); every refit bitwise its stored fit.
*     mg_beta (6 items, two factors) as a CFA; on mg_E two factors did not
*     converge (check 7 refused, correctly; first run 2026-10-31)
preserve
import delimited using "fixtures/data/mg_beta.csv", clear varnames(1) delimiters(",")
foreach L in configural metric {
    quietly dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6), group(g) invariance(`L') metric(unitvariance)
    estimates store g2f_`L'
}
tdf_log "difftest g2f_metric g2f_configural"
if `tdf_rc' {
    display as error "check 7g FAILED: two-factor G > 1 pair, rc `tdf_rc'"
    local ++nfail
}
else {
    mata: tdf_gap(2)
    local nfail = `nfail' + `gapbad'
}
restore
* 7h  (review gap 4) keepbound at G > 1: y5 = y1 puts the pair y1~~y5 at the
*     bound in both groups; both keepbound fits are reproduced (as DT-M15)
replace y5 = y1
foreach L in configural thresholds {
    quietly dwlsmv (F -> y1 y2 y3 y4 y5), group(g) invariance(`L') metric(unitvariance) keepbound
    estimates store gk_`L'
}
tdf_log "difftest gk_thresholds gk_configural"
if `tdf_rc' | TDS[1, 6] != 1 {
    display as error "check 7h FAILED: keepbound G > 1 pair, rc `tdf_rc'"
    local ++nfail
}
display "check 7 done"

* ---- check 8: sign orientation of dCFI, dRMSEA (user, 2026-10-31) ----------
* rows most restricted first, columns (CFI, RMSEA); expected H0 minus H1
mata:
void tdf_sign()
{
    real matrix d
    real scalar bad

    bad = 0
    // restriction worsens fit: CFI .960 -> .940, RMSEA .045 -> .060
    d = dwlsmv_dtfront_delta((.940, .060 \ .960, .045), 1)
    if (rows(d) != 1 | abs(d[1, 1] - (-.020)) > 1e-12 | abs(d[1, 2] - .015) > 1e-12) {
        printf("{err}check 8 FAILED: worsening step gives (%g, %g), expected (-.020, +.015)\n", d[1, 1], d[1, 2])
        bad++
    }
    // restriction improves fit: the signed difference, not a forced sign
    d = dwlsmv_dtfront_delta((.955, .040 \ .950, .050), 1)
    if (abs(d[1, 1] - .005) > 1e-12 | abs(d[1, 2] - (-.010)) > 1e-12) {
        printf("{err}check 8 FAILED: improving step gives (%g, %g), expected (+.005, -.010)\n", d[1, 1], d[1, 2])
        bad++
    }
    // three models, the first step (restricted pair) skipped (d = 0):
    // that row blank, the next against the less restricted model
    d = dwlsmv_dtfront_delta((.90, .08 \ .90, .08 \ .95, .05), (0 \ 1))
    if (!missing(d[1, 1]) | !missing(d[1, 2]) | abs(d[2, 1] - (-.05)) > 1e-12 | abs(d[2, 2] - .03) > 1e-12) {
        printf("{err}check 8 FAILED: skipped chain gives (%g, %g / %g, %g)\n", d[1, 1], d[1, 2], d[2, 1], d[2, 2])
        bad++
    }
    st_local("sgnbad", strofreal(bad))
}
end
mata: tdf_sign()
local nfail = `nfail' + `sgnbad'
capture mata: (void) dwlsmv_dtfront_delta((.9, .05 \ .91, .04), (1 \ 1))
if _rc != 3498 {
    display as error "check 8 FAILED: wrong status length gave rc " _rc
    local ++nfail
}
display "check 8 done"

if `nfail' {
    display as error _n "test_difftest_front: `nfail' failure(s)"
    exit 9
}
display _n "test_difftest_front: all checks passed"
