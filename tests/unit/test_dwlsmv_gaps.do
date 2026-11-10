*! Module 11: dwlsmv.ado gaps (spec D6, D8a, D9, D5 e(sample), Section 3 listwise).
*!
*! Executed under StataMP 18 on 2026-10-05 and it passes. Drafted by a
*! read-only agent from the Module 11 test-gap audit (2026-10-04); run by the
*! main session, with 20 mutations of dwlsmv.ado and dwlsmv_m11_{stage1,
*! legend,fit}.mata, each of which fails it. Seventeen are reported by their
*! own section. Two stop the file at an uncaught rc in an earlier fit before
*! their section runs: reversing the cov name order (at G1's fit) and
*! sorting inside the wrapper (at G3's). One, sorting Stage 1's columns, is
*! reported by G1, which then stops at r(111) on a cut name that no longer
*! exists. The first run caught one wrapper defect (G8: the message lost its
*! space). Builds build/ldwlsmv.mlib first
*! and runs against the library and src/ado; both adopath entries are removed
*! at the end.
*!
*! Checks:
*!   G1  D9 (owed): with path order y3 u1 y2 (u1 binary), the cut names,
*!       e(indicators), the canonical ind line, Stage 1's thresholds (each
*!       invnormal of its own variable's cumulative share, 1e-6 set before
*!       measuring) and the abort legend all follow the paths.
*!   G2  listwise deletion of missing indicators (. and .a): e(N) and
*!       e(sample) are the complete cases, e(b) equals the complete-case fit
*!       exactly; a missing non-model variable drops nothing.
*!   G3  if and in: e(N), e(sample) and e(b) equal the subset fit's.
*!   G4  D6's internal error both ways (a sem name not built, a built name not
*!       in sem's stripe): rc 3498 and its message; the true stripe passes.
*!   G5  D6 names on more models: correlated factors, unit variance, two
*!       endogenous factors with a disturbance covariance.
*!   G6  state through dwlsmv on a Stage 1 abort, a Module 3 abort, an rc 198
*!       refusal and a success: data, sort order, rngstate, e(cmd), and the
*!       Mata externals removed.
*!   G7  replay: header and unchanged e(b); rc 301 after another estimator.
*!   G8  a non-numeric indicator is Module 11's refusal: rc 198 (decided by
*!       the user 2026-10-05) and a message naming the variable. Its mutation
*!       targets the wrapper's conversion itself (exit 198 -> exit with
*!       confirm's rc), so the test shows the public boundary is enforced
*!       rather than recording whatever rc Stata returns.
*!   G9  an empty estimation sample: rc 2000.
*!   G10 the header under e(converged) = 0, set by hand on a posted fit.
*!       NOT TESTED: the unforced path, a fit that actually fails to converge
*!       and is posted with e(converged) = 0. No option limits iterations and
*!       none is added for testing (decided by the user 2026-10-05).
*!   G11 sem's estimation options (vce(), method(), iterate(), technique(),
*!       from(), difficult, ...) are refused at rc 198 with a message; vce()'s
*!       says clustering and weights are planned for version 2; nolog passes.
*!   G12 level() on the fit, on replay and with standardized; e(b) unchanged.
*!   G13 build/dwlsmv_example.dta (tools/build_example.do) gives exactly the
*!       fit on fixtures/data/case_E.csv.
*!   G11-G13 added 2026-10-05 (decided by the user).
*!   G14 e(version) equals the version on src/ado/dwlsmv.ado's first line, and
*!       replay leaves it. Added 2026-10-07 (decided by the user).
*!   G15 Module 10's singular-value ratios are stored (e(gamma_ratio),
*!       e(gamma_ratio_raw), e(gamma_rank), e(bread_ratio)) and equal Module 10
*!       computed directly; the output prints no ratio and does print the
*!       category counts. Added 2026-10-07 (decided by the user).
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

capture program drop tdw_run
program tdw_run, rclass
    gettoken phrase 0 : 0
    tempfile lg
    quietly log using "`lg'", text replace name(tdwlog)
    capture noisily `0'
    local rc = _rc
    quietly log close tdwlog
    mata: st_local("found", strofreal(strpos(invtokens(cat(st_local("lg"))', " "), st_local("phrase")) > 0))
    return scalar rc    = `rc'
    return scalar found = `found'
end

import delimited using "fixtures/data/case_E.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
quietly destring y*, replace
local mE "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)"

* ---- G1. D9 owed: Stage 1 columns, ind line, /y:cutk, e(indicators), legend agree ----
* u1 is binary (K=2), y3 has K=3, so a column/name misalignment changes the cut count.
* Path order y3 u1 y2 differs from dataset order (y1..y6 then u1).
generate byte u1 = y1 > 0
local mo "(F1 -> y3 u1 y2) (F2 -> y4 y5 y6)"
quietly dwlsmv `mo'
local want "/y3:cut1 /y3:cut2 /u1:cut1 /y2:cut1 /y2:cut2 /y4:cut1 /y4:cut2 /y5:cut1 /y5:cut2 /y6:cut1 /y6:cut2"
local names : colfullnames e(b)
local cuts
foreach n of local names {
    if strpos("`n'", ":cut") local cuts `cuts' `n'
}
if "`cuts'" != "`want'" | "`e(indicators)'" != "y3 u1 y2 y4 y5 y6" {
    display as error " G1. cut names `cuts'; e(indicators) `e(indicators)'"
    local ++nfail
}
quietly _dwlsmv_translate , paths(`mo')
if strpos(`"`r(canon)'"', ";ind y3 u1 y2 y4 y5 y6;") == 0 {
    display as error " G1. ind line: `r(canon)'"
    local ++nfail
}
* Stage 1's columns: saturated thresholds equal invnormal(cumulative share) of
* the variable the name says. Tolerance 1e-6 chosen BEFORE measuring (optimizer
* stopping, not exact equality); the misaligned values differ by >> 1e-6.
quietly count if e(sample)
local n = r(N)
foreach v in y3 y2 {
    quietly count if `v' == 0 & e(sample)
    local t1 = invnormal(r(N)/`n')
    quietly count if `v' <= 1 & e(sample)
    local t2 = invnormal(r(N)/`n')
    if abs(_b[/`v':cut1] - `t1') > 1e-6 | abs(_b[/`v':cut2] - `t2') > 1e-6 {
        display as error " G1. thresholds of `v' are not `v''s"
        local ++nfail
    }
}
quietly count if u1 == 0 & e(sample)
if abs(_b[/u1:cut1] - invnormal(r(N)/`n')) > 1e-6 {
    display as error " G1. threshold of u1 is not u1's"
    local ++nfail
}
* legend under the same reordering: q7 duplicates y3 and sits 4th
generate q7 = y3
tdw_run "y1 = y3" dwlsmv (F1 -> y3 u1 y2) (F2 -> q7 y5 y6)
local l1 = r(found)
local rcg = r(rc)
tdw_run "y4 = q7" dwlsmv (F1 -> y3 u1 y2) (F2 -> q7 y5 y6)
if `rcg' != 3498 | !`l1' | !r(found) {
    display as error " G1. legend under path order: rc `rcg', y1=y3 `l1', y4=q7 " r(found)
    local ++nfail
}
drop q7 u1

* ---- G2. listwise deletion of missing indicators; non-model missing ignored ----
preserve
quietly replace y2 = .  in 1/25
quietly replace y5 = .a in 30/40
set seed 7
quietly generate double x = cond(_n <= 60, ., rnormal())
quietly count if !missing(y1, y2, y3, y4, y5, y6)
local nc = r(N)
quietly dwlsmv `mE'
tempname bm bd
matrix `bm' = e(b)
local Nm = e(N)
quietly count if e(sample)
local Sm = r(N)
quietly count if e(sample) & missing(y2, y5)
local Sbad = r(N)
quietly drop if missing(y1, y2, y3, y4, y5, y6)
quietly dwlsmv `mE'
matrix `bd' = e(b)
if `Nm' != `nc' | `Sm' != `nc' | `Sbad' != 0 | colsof(`bm') != colsof(`bd') {
    display as error " G2. listwise: e(N) `Nm', e(sample) `Sm', want `nc', missing in sample `Sbad'"
    local ++nfail
}
else if mreldif(`bm', `bd') != 0 {
    display as error " G2. listwise fit differs from the fit on complete cases"
    local ++nfail
}
restore

* ---- G3. if / in and e(sample) ----------------------------------------------------
quietly dwlsmv `mE' if _n <= 600
tempname bi bk bn
matrix `bi' = e(b)
quietly count if e(sample)
local si = r(N)
quietly count if e(sample) & _n > 600
local sout = r(N)
local Ni = e(N)
quietly dwlsmv `mE' in 1/600
matrix `bn' = e(b)
preserve
quietly keep in 1/600
quietly dwlsmv `mE'
matrix `bk' = e(b)
restore
if `Ni' != 600 | `si' != 600 | `sout' != 0 | mreldif(`bi', `bk') != 0 | mreldif(`bn', `bk') != 0 {
    display as error " G3. if/in: e(N) `Ni', e(sample) `si', outside `sout'"
    local ++nfail
}

* ---- G4. D6 internal error, both directions: rc 3498 and the message -------------
quietly _dwlsmv_translate , paths(`mE')
local canon    `"`r(canon)'"'
local ind      `r(ind)'
local semnames `r(semnames)'
local sem1 : subinstr local semnames "/var(F1)" "/var(F9)"
local sem2 `"`semnames' /var(F9)"'
generate byte tu = 1
tempname b V C bnd
mata: dwlsmv_m11_stage1("tu", "`ind'", 0)
tdw_run "is not in sem's stripe" mata: dwlsmv_m11_fit(`"`canon'"', "`ind'", "`sem1'", "`b'", "`V'", "`C'", "`bnd'", "marker")
local r1 = r(rc)
local f1 = r(found)
mata: dwlsmv_m11_stage1("tu", "`ind'", 0)
tdw_run "was not built" mata: dwlsmv_m11_fit(`"`canon'"', "`ind'", "`sem2'", "`b'", "`V'", "`C'", "`bnd'", "marker")
local r2 = r(rc)
local f2 = r(found)
* negative control: the true stripe does not abort
mata: dwlsmv_m11_stage1("tu", "`ind'", 0)
capture noisily quietly mata: dwlsmv_m11_fit(`"`canon'"', "`ind'", "`semnames'", "`b'", "`V'", "`C'", "`bnd'", "marker")
local r0 = _rc
mata: dwlsmv_m11_clear()
drop tu
if `r1' != 3498 | !`f1' | `r2' != 3498 | !`f2' | `r0' != 0 {
    display as error " G4. internal error: renamed rc `r1' msg `f1'; extra rc `r2' msg `f2'; true stripe rc `r0'"
    local ++nfail
}

* ---- G5. D6 names on more models (sets compared; order is Module 3's) ------------
* expected sets: confirm once against tools/probe_d6_names.do's records
capture noisily quietly dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6)
local rc5 = _rc
local got : colfullnames e(b)
local w5a "y1:F1 y2:F1 y3:F1 y4:F2 y5:F2 y6:F2 /var(F1) /var(F2) /cov(F1,F2)"
local miss : list w5a - got
if `rc5' != 0 | "`miss'" != "" {
    display as error " G5a. correlated: rc `rc5', missing `miss'"
    local ++nfail
}
capture noisily quietly dwlsmv `mE', metric(unitvariance)
local rc5 = _rc
local got : colfullnames e(b)
local w5b "y1:F1 y4:F2 F2:F1 /var(F1) /var(e.F2)"
local miss : list w5b - got
if `rc5' != 0 {
    display as error " G5b. unitvariance: rc `rc5'"
    local ++nfail
}
else if "`miss'" != "" | _b[/var(e.F2)] != 1 | _b[/var(F1)] != 1 | strpos("`got'", "/var(F2)") {
    display as error " G5b. unitvariance names/values: `got'"
    local ++nfail
}
capture noisily quietly dwlsmv (F1 -> y1 y2) (F2 -> y3 y4) (F3 -> y5 y6) (F2 <- F1) (F3 <- F1), cov(e.F2*e.F3)
local rc5 = _rc
local got : colfullnames e(b)
local w5c "F2:F1 F3:F1 /var(F1) /var(e.F2) /var(e.F3) /cov(e.F2,e.F3)"
local miss : list w5c - got
if `rc5' != 0 | "`miss'" != "" {
    display as error " G5c. two endogenous: rc `rc5', missing `miss'"
    local ++nfail
}

* ---- G6. state through dwlsmv: data, sort, rngstate, e(), externals --------------
set seed 4711
quietly generate double x = rnormal()
sort x
quietly regress y1 y2
quietly _datasignature
local dsig `r(datasignature)'
local rng `c(rngstate)'
capture program drop tdw_state
program tdw_state
    args label ecmd dsig rng
    local bad 0
    if "`e(cmd)'" != "`ecmd'" local bad 1
    quietly _datasignature
    if "`r(datasignature)'" != "`dsig'" local bad 1
    if "`c(rngstate)'" != "`rng'" local bad 1
    if "`: sortedby'" != "x" local bad 1
    mata: st_local("ext", strofreal((findexternal("dwlsmv_m11_st") != NULL) + (findexternal("dwlsmv_m11_Y") != NULL)))
    if `ext' local bad 1
    if `bad' display as error "  G6 `label': state NOT restored"
    c_local tdw_bad `bad'
end
generate q7 = y1
capture quietly dwlsmv (F1 -> y1 y2 y3 q7) (F2 -> y4 y5 y6)
local rca = _rc
drop q7
quietly _datasignature
local dsig `r(datasignature)'
tdw_state "Stage 1 abort" "regress" "`dsig'" "`rng'"
local nfail = `nfail' + `tdw_bad'
capture quietly dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (G -> y1@1 y2@1 y3@1 y4@1 y5@1 y6@1), cov(G*F1@0) cov(G*F2@0)
local rcb = _rc
tdw_state "Module 3 abort" "regress" "`dsig'" "`rng'"
local nfail = `nfail' + `tdw_bad'
capture quietly dwlsmv `mE', metric(total)
local rcc = _rc
tdw_state "rc 198 refusal" "regress" "`dsig'" "`rng'"
local nfail = `nfail' + `tdw_bad'
quietly dwlsmv `mE'
tdw_state "success" "dwlsmv" "`dsig'" "`rng'"
local nfail = `nfail' + `tdw_bad'
if `rca' != 3498 | `rcb' != 3498 | `rcc' != 198 {
    display as error " G6. rcs `rca' `rcb' `rcc'"
    local ++nfail
}

* ---- G7. replay: content, same table, rc 301 after another estimator --------------
quietly dwlsmv `mE', metric(unitvariance)
tempname b0
matrix `b0' = e(b)
tdw_run "Metric: factor variances fixed at 1 (for a factor predicted by other factors, its residual variance)" dwlsmv
local r7 = r(rc)
local f7 = r(found)
if `r7' != 0 | !`f7' | mreldif(`b0', e(b)) != 0 {
    display as error " G7. replay rc `r7', header `f7'"
    local ++nfail
}
quietly regress y1 y2
capture noisily dwlsmv
if _rc != 301 {
    display as error " G7. replay after regress rc " _rc
    local ++nfail
}

* ---- G8. non-numeric indicator: Module 11's refusal, rc 198 ------------------------
generate str2 s1 = string(y1)
tdw_run "not numeric: s1" dwlsmv (F1 -> s1 y2 y3) (F2 -> y4 y5 y6)
display as text " G8. string indicator rc " r(rc) ", variable named " r(found)
if r(rc) != 198 | !r(found) {
    local ++nfail
}
drop s1

* ---- G9. empty estimation sample: rc 2000 -----------------------------------------
capture noisily dwlsmv `mE' if y1 > 99
if _rc != 2000 {
    display as error " G9. empty sample rc " _rc
    local ++nfail
}

* ---- G10. header under e(converged) = 0 (header only; the unforced path is untested) -
capture program drop tdw_unconv
program tdw_unconv, eclass
    ereturn scalar converged = 0
end
quietly dwlsmv `mE'
tdw_unconv
tdw_run "convergence not achieved" dwlsmv
if r(rc) != 0 | !r(found) {
    display as error " G10. no 'convergence not achieved' under e(converged)=0"
    local ++nfail
}

* ---- G11. sem's estimation options are refused, rc 198 ------------------------------
* each would otherwise reach sem noestimate and change nothing (probed
* 2026-10-05: vce(cluster id) gave the plain fit's e(V) at rc 0)
generate long cid = ceil(_n / 5)
local r11 ""
local b11 0
foreach o in "vce(cluster cid)" "vce(robust)" "method(adf)" "iterate(3)" "technique(nr)" "from(b0)" "difficult" {
    tdw_run "is not available with dwlsmv" dwlsmv `mE', `o'
    local r11 "`r11' `=r(rc)'"
    if r(rc) != 198 | !r(found) local b11 1
}
tdw_run "clustering and sampling weights are planned" dwlsmv `mE', vce(cluster cid)
local f11 = r(found)
drop cid
capture quietly dwlsmv `mE', nolog
local n11 = _rc
display as text " G11. estimation options rc:`r11'; vce line `f11'; nolog rc `n11'"
if `b11' | !`f11' | `n11' {
    local ++nfail
}

* ---- G12. level() on the fit and on replay, both tables; nothing stored ------------
quietly dwlsmv `mE'
tempname b12
matrix `b12' = e(b)
tdw_run "[90%" dwlsmv `mE', level(90)
local f12a = r(found)
tdw_run "[80%" dwlsmv, level(80)
local f12b = r(found)
tdw_run "[99%" dwlsmv, standardized level(99)
local f12c = r(found)
tdw_run "[95%" dwlsmv
local f12d = r(found)
display as text " G12. level(90) fit `f12a', replay level(80) `f12b', standardized level(99) `f12c', default 95 `f12d'"
if !`f12a' | !`f12b' | !`f12c' | !`f12d' | mreldif(e(b), `b12') != 0 {
    local ++nfail
}

* ---- G13. the distributed example dataset is the fixture's data ---------------------
* both in file order: the data here are sorted by x since G6, and row order
* moves Gamma's sums by an ulp (seen on the first run: 5.6e-17)
quietly do tools/build_example.do
tempname b13
preserve
import delimited using "fixtures/data/case_E.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
quietly destring y*, replace
quietly dwlsmv `mE'
matrix `b13' = e(b)
quietly use "build/dwlsmv_example.dta", clear
quietly dwlsmv `mE'
local d13 = mreldif(e(b), `b13')
local n13 = e(N)
restore
display as text " G13. example data: e(N) `n13', mreldif(e(b), fixture fit) `d13'"
if `d13' != 0 | `n13' != 5000 {
    local ++nfail
}

* ---- G14. e(version) is the version on dwlsmv.ado's first line, and replay keeps it --
* Read from the ado file itself, so a release that bumps line 1 and not the
* posted value (or the reverse) fails here.
file open fh14 using "`here'/src/ado/dwlsmv.ado", read text
file read fh14 line14
file close fh14
local v14 : word 3 of `line14'
quietly dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6)
local e14 "`e(version)'"
quietly dwlsmv
local r14 "`e(version)'"
display as text " G14. line 1 says `v14'; e(version) `e14'; after replay `r14'"
if "`e14'" != "`v14'" | "`r14'" != "`v14'" | !regexm("`v14'", "^[0-9]+\.[0-9]+\.[0-9]+$") {
    local ++nfail
}

* ---- G15. Module 10's ratios: stored in e(), never printed (2026-10-07) ------------
* The stored values against Module 10 computed directly on the same data and
* model. They travel through a %21.0g local, as e(T3) does, so the comparison
* allows 1e-15 relative (the format's precision, not a fit tolerance).
mata:
real rowvector tdw_ratios(string scalar canon, string scalar ind)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit    scalar fit
    real matrix Y, G, B

    Y   = st_data(., tokens(ind))
    st  = dwls_sample_stats(Y, J(cols(Y), 1, 3))
    par = dwls_identify(subinstr(canon, ";", char(10)), st.map)
    fit = dwls_fit(st, par)
    G   = dwlsmv_gammadiag(st, Y)
    B   = dwlsmv_breaddiag(fit.Delta, st.Gamma)
    return((G[7, 2], G[6, 2], G[7, 3], B[4, 2]))
}
end
local m15 "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)"
tdw_run "singular value" dwlsmv `m15'
local sv15 = r(found)
tdw_run "Gamma" dwlsmv `m15'
local sv15 = `sv15' | r(found)
tdw_run "Category counts (the fewest observations in any one category):" dwlsmv `m15'
local cc15 = r(found)
quietly dwlsmv `m15'
* printed at 4 significant digits: e(bread_ratio) moves at ~7e-13 relative
* between processor counts (measured 2026-10-07, chosen after measuring); the
* comparison below uses the full values
local e15 : display %9.3e e(gamma_ratio) " " %9.3e e(gamma_ratio_raw) " " e(gamma_rank) " " %9.3e e(bread_ratio)
mata: TDW_e15 = (st_numscalar("e(gamma_ratio)"), st_numscalar("e(gamma_ratio_raw)"), ///
                 st_numscalar("e(gamma_rank)"), st_numscalar("e(bread_ratio)"))
quietly _dwlsmv_translate , paths(`m15')
mata: TDW_d15 = tdw_ratios(st_global("r(canon)"), "y1 y2 y3 y4 y5 y6")
mata: st_local("ok15", strofreal(!missing(TDW_e15) & all(abs(TDW_e15 :- TDW_d15) :<= 1e-15 :* abs(TDW_d15)) & TDW_e15[3] == TDW_d15[3]))
display as text " G15. stored `e15'; printed a ratio `sv15'; counts heading `cc15'; match direct Module 10 `ok15'"
if `sv15' | !`cc15' | !`ok15' {
    local ++nfail
}


* G16-G18 (2026-11-03): a first-stage message on a fit that continues gets the
* legend; difftest's r() holds only what its help lists; a refusal prints no
* Mata traceback.
capture program drop tdw_logged
program define tdw_logged, sclass
    * runs the command under a log; s(rc) its return code, s(txt) the file holding its output
    tempfile lf
    quietly log using `"`lf'"', text replace name(tdw16)
    capture noisily `0'
    local rc = _rc
    quietly log close tdw16
    local out "`c(tmpdir)'/tdw_g16.txt"
    quietly copy `"`lf'"' `"`out'"', replace
    sreturn local rc "`rc'"
    sreturn local txt `"`out'"'
end
capture program drop tdw_has
program define tdw_has, rclass
    args file text
    local s = fileread(`"`file'"')
    return scalar has = strpos(`"`macval(s)'"', `"`text'"') > 0
end

clear
set seed 20261103
set obs 1000
generate double f = rnormal()
forvalues j = 1/5 {
    generate double u`j' = 0.7*f + sqrt(0.51)*rnormal()
    generate byte v`j' = (u`j' > -0.4) + (u`j' > 0.4)
}
generate byte g = 1 + (_n > 500)
order v1 v2 v3 v4 v5

* G16: no repair or bound pair, no legend; a repair and a bound pair, the legend naming them
tdw_logged dwlsmv (F1 -> v1-v5), metric(unitvariance)
local f = s(txt)
tdw_has `"`f'"' "Positional labels above"
local none16 = !r(has)
generate byte b1 = u1 > 0
generate byte b2 = u1 > 0.5
generate byte b3 = u2 > 0
generate byte b4 = u3 > 0
generate byte b5 = u4 > 0
tdw_logged dwlsmv (F1 -> b1 b3 b2 b4 b5), metric(unitvariance)
local rc16a = s(rc)
local f = s(txt)
tdw_has `"`f'"' "y3 = b2"
local rep16 = r(has)
generate byte d1 = v1
quietly replace d1 = 1 if _n == 1 & v1 == 0
tdw_logged dwlsmv (F1 -> v1 d1 v2 v3 v4), metric(unitvariance) keepbound
local rc16b = s(rc)
local f = s(txt)
tdw_has `"`f'"' "y2 = d1"
local kb16 = r(has)
tdw_logged dwlsmv (F1 -> v1 d1 v2 v3 v4), group(g) invariance(configural) metric(unitvariance) keepbound
local rc16c = s(rc)
local f = s(txt)
tdw_has `"`f'"' "y2 = d1"
local gkb16 = r(has)
display as text " G16. none `none16'; repair `rep16' (rc `rc16a'); keepbound `kb16' (rc `rc16b'); groups `gkb16' (rc `rc16c')"
if !`none16' | !`rep16' | !`kb16' | !`gkb16' | `rc16a' | `rc16b' | `rc16c' {
    local ++nfail
}

* G17: r() at G > 1
quietly dwlsmv (F1 -> v1-v5), group(g) invariance(configural) metric(unitvariance)
estimates store g17c
quietly dwlsmv (F1 -> v1-v5), group(g) invariance(scalar) metric(unitvariance)
estimates store g17s
quietly difftest g17c g17s
local ok17 = (`"`r(levels)'"' == "") & missing(r(r)) & r(N) == 1000 & rowsof(r(steps)) == 1
display as text " G17. r() holds no levelsof leftovers: `ok17'"
if !`ok17' {
    local ++nfail
}

* G18: three refusals at their rc, with the message and no traceback
quietly generate byte v6 = (0.7*f + sqrt(0.51)*rnormal() > 0)
quietly generate byte a1 = v1 > 0
quietly generate byte a2 = v2 > 0
quietly generate byte a3 = v3 > 0
quietly generate byte a4 = v4 > 0
tdw_logged dwlsmv (F1 -> a1-a4), group(g) invariance(scalar) metric(unitvariance)
local rc18a = s(rc)
local f = s(txt)
tdw_has `"`f'"' "function returned error"
local tb18a = r(has)
tdw_has `"`f'"' "is not available for an all-binary factor"
local ms18a = r(has)
quietly dwlsmv (F1 -> v1-v3) (F2 -> v4 v5 v6), cov(F1*F2@.3) metric(unitvariance)
estimates store g18a
quietly dwlsmv (F1 -> v1-v3) (F2 -> v4 v5 v6), cov(F1*F2@0) metric(unitvariance)
estimates store g18b
tdw_logged difftest g18a g18b
local rc18b = s(rc)
local f = s(txt)
tdw_has `"`f'"' "function returned error"
local tb18b = r(has)
tdw_has `"`f'"' "the models are not nested"
local ms18b = r(has)
tdw_logged dwlsmv (F1 -> v1 d1 v2 v3 v4), metric(unitvariance)
local rc18c = s(rc)
local f = s(txt)
tdw_has `"`f'"' "function returned error"
local tb18c = r(has)
tdw_has `"`f'"' "y2 = d1"
local ms18c = r(has)
display as text " G18. rc `rc18a'/`rc18b'/`rc18c'; traceback `tb18a'/`tb18b'/`tb18c'; message `ms18a'/`ms18b'/`ms18c'"
if `rc18a' != 198 | `rc18b' != 3498 | `rc18c' != 3498 | `tb18a' | `tb18b' | `tb18c' | !`ms18a' | !`ms18b' | !`ms18c' {
    local ++nfail
}
capture erase "`c(tmpdir)'/tdw_g16.txt"

adopath - "`here'/src/ado"
adopath - "`here'/build"
if `nfail' {
    display as error _n "test_dwlsmv_gaps: `nfail' failure(s)"
    exit 9
}
display as result "test_dwlsmv_gaps: D9 order holds from Stage 1 to names and legend; listwise, if/in and e(sample) agree with the subset fit; D6's internal error fires both ways at rc 3498; state and externals are restored on every exit; replay is guarded"
exit 0
