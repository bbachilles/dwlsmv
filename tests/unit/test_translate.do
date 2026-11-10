*! Module 11: _dwlsmv_translate -- sem's reading of the model as the canonical
*! form (spec Module 11 D1, D2, D3, D9).
*!
*! Executed under StataMP 18 on 2026-10-04 and it passes; 1e and 4.9 added and
*! executed 2026-10-05 (suite green).
*!
*! The command is distributed, so this test builds build/ldwlsmv.mlib from
*! src/mata first (tools/build_mlib.do) and runs against the library and
*! src/ado, as a user would; a stale library cannot pass. Both adopath
*! entries are removed again at the end, so later tests see the source files.
*!
*! Checks:
*!   1  translation, marker metric: the expected canonical text, exactly, for
*!      D1's models (structural, two correlated factors, a fixed covariance,
*!      D9's path order y2 y1 y3), and 1e, D1's token scan: an option keyword
*!      that abbreviates a data variable (var -> var1) leaves the text as 1b's.
*!   2  metric(unitvariance), D2: the injected options and canonical text for
*!      D2's models; an endogenous factor gets var(e.F@1), never var(F@1); a
*!      factor the user normalised gets nothing.
*!   3  the D3 random-intercept case translates and leaves F1 with no metric
*!      line (Module 3's rc 3498 reports it, not Module 11).
*!   4  Module 11's refusals exit rc 198 (D8a): a non-unit loading and a
*!      non-unit variance as scale (D3), a second-order factor (D3), a
*!      residual covariance among indicators, a constraint tying two
*!      loadings, a fixed indicator residual variance, an observed exogenous
*!      variable, a path between observed variables, and more than one group
*!      (4.9, group(g), added 2026-10-05). 4.9 also requires Module 11's own
*!      message: with the earlier continuous dummy, sem refused first ("not
*!      nonnegative integer valued", also rc 198), so rc alone could not tell.
*!   5  sem's refusals keep sem's rc 498 (D1, D4): cov(F1*e.F2) beside F2 <- F1,
*!      and var(F2@1) on an endogenous factor.
*!   6  state restoration (D1): data signature, sort order, rngstate and the
*!      user's e() are unchanged after a success, a Module 11 refusal and a sem
*!      refusal.
*!   7  dwlsmv_tr_num(): fixed values from 1 to 1e-9 read back as the
*!      identical double, with no space (2026-10-07: %21.0g alone did not).
*!
*! Expected texts are D1's and D2's probe records (tools/probe_d1_sem.do,
*! tools/probe_d2_metric.do; StataMP 18 and 16.1 identical).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off

quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

import delimited using "fixtures/data/case_E.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
quietly destring y*, replace
set seed 4711
quietly generate double x = rnormal()
quietly generate byte g = mod(_n, 2)
quietly generate double var1 = rnormal()
set varabbrev on

local nfail 0
capture program drop tt_expect
program tt_expect
    args label model opts metric expect
    capture noisily _dwlsmv_translate , paths(`model') semopts(`opts') metric(`metric')
    if _rc {
        display as error "  `label': rc " _rc
        c_local tt_bad 1
        exit
    }
    if `"`r(canon)'"' != `"`expect'"' {
        display as error "  `label': canonical text differs"
        display as error "    got:  `r(canon)'"
        display as error "    want: `expect'"
        c_local tt_bad 1
        exit
    }
    display as text "  `label': ok (`r(injected)')"
    c_local tt_bad 0
    c_local tt_inj `"`r(injected)'"'
end

display as text "1. translation, marker metric"
tt_expect "1a structural" "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)" "" "" ///
    "factor F1 F2;ind y1 y2 y3 y4 y5 y6;beta F2 F1 free;lambda F1 y2 free;lambda F1 y3 free;lambda F2 y5 free;lambda F2 y6 free;metric marker F1 y1;metric marker F2 y4;psi F1 F1 free;psi F2 F2 free"
local nfail = `nfail' + `tt_bad'
tt_expect "1b correlated" "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6)" "" "marker" ///
    "factor F1 F2;ind y1 y2 y3 y4 y5 y6;lambda F1 y2 free;lambda F1 y3 free;lambda F2 y5 free;lambda F2 y6 free;metric marker F1 y1;metric marker F2 y4;psi F1 F1 free;psi F2 F1 free;psi F2 F2 free"
local nfail = `nfail' + `tt_bad'
tt_expect "1c fixed cov" "(F1 -> y1-y3) (F2 -> y4-y6)" "cov(F1*F2@.3)" "" ///
    "factor F1 F2;ind y1 y2 y3 y4 y5 y6;lambda F1 y2 free;lambda F1 y3 free;lambda F2 y5 free;lambda F2 y6 free;metric marker F1 y1;metric marker F2 y4;psi F1 F1 free;psi F2 F1 fixed .3;psi F2 F2 free"
local nfail = `nfail' + `tt_bad'
tt_expect "1d D9 order" "(F1 -> y2 y1 y3) (F2 -> y4 y5 y6)" "" "" ///
    "factor F1 F2;ind y2 y1 y3 y4 y5 y6;lambda F1 y1 free;lambda F1 y3 free;lambda F2 y5 free;lambda F2 y6 free;metric marker F1 y2;metric marker F2 y4;psi F1 F1 free;psi F2 F1 free;psi F2 F2 free"
local nfail = `nfail' + `tt_bad'
* D1's token scan: with varabbrev on, the option keyword var abbreviates the
* data's var1, which the scan adds to the dummy; the canonical text is 1b's.
tt_expect "1e token scan" "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6)" "var(F1)" "" ///
    "factor F1 F2;ind y1 y2 y3 y4 y5 y6;lambda F1 y2 free;lambda F1 y3 free;lambda F2 y5 free;lambda F2 y6 free;metric marker F1 y1;metric marker F2 y4;psi F1 F1 free;psi F2 F1 free;psi F2 F2 free"
local nfail = `nfail' + `tt_bad'

display as text "2. metric(unitvariance), D2"
tt_expect "2a structural" "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)" "" "unitvariance" ///
    "factor F1 F2;ind y1 y2 y3 y4 y5 y6;beta F2 F1 free;lambda F1 y1 free;lambda F1 y2 free;lambda F1 y3 free;lambda F2 y4 free;lambda F2 y5 free;lambda F2 y6 free;metric unitvariance F1;metric unitvariance F2"
local nfail = `nfail' + `tt_bad'
if `"`tt_inj'"' != "var(F1@1) var(e.F2@1)" {
    display as error "  2a: injected `tt_inj'"
    local ++nfail
}
tt_expect "2b correlated" "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6)" "" "unitvariance" ///
    "factor F1 F2;ind y1 y2 y3 y4 y5 y6;lambda F1 y1 free;lambda F1 y2 free;lambda F1 y3 free;lambda F2 y4 free;lambda F2 y5 free;lambda F2 y6 free;metric unitvariance F1;metric unitvariance F2;psi F2 F1 free"
local nfail = `nfail' + `tt_bad'
tt_expect "2c user marker kept" "(F1 -> y1@1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)" "" "unitvariance" ///
    "factor F1 F2;ind y1 y2 y3 y4 y5 y6;beta F2 F1 free;lambda F1 y2 free;lambda F1 y3 free;lambda F2 y4 free;lambda F2 y5 free;lambda F2 y6 free;metric marker F1 y1;metric unitvariance F2;psi F1 F1 free"
local nfail = `nfail' + `tt_bad'
if `"`tt_inj'"' != "var(e.F2@1)" {
    display as error "  2c: injected `tt_inj'"
    local ++nfail
}

display as text "3. the unanchored random-intercept case translates"
tt_expect "3 random intercept" "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (G -> y1@1 y2@1 y3@1 y4@1 y5@1 y6@1)" "cov(G*F1@0) cov(G*F2@0)" "" ///
    "factor F1 F2 G;ind y1 y2 y3 y4 y5 y6;lambda F1 y1 free;lambda F1 y2 free;lambda F1 y3 free;lambda F2 y4 free;lambda F2 y5 free;lambda F2 y6 free;lambda G y2 fixed 1;lambda G y3 fixed 1;lambda G y4 fixed 1;lambda G y5 fixed 1;lambda G y6 fixed 1;metric marker G y1;psi F1 F1 free;psi F2 F1 free;psi F2 F2 free;psi G G free"
local nfail = `nfail' + `tt_bad'

* ---- state captured once; checked after every call below -----------------
sort y3 y1
quietly regress y1 y2
local ecmd `e(cmd)'
local esig = e(N)
quietly _datasignature
local dsig `r(datasignature)'
local rng `c(rngstate)'
capture program drop tt_state
program tt_state
    args label ecmd esig dsig rng sorted
    local bad 0
    if "`e(cmd)'" != "`ecmd'" | e(N) != `esig' local bad 1
    quietly _datasignature
    if "`r(datasignature)'" != "`dsig'" local bad 1
    if "`c(rngstate)'" != "`rng'" local bad 1
    if "`: sortedby'" != "`sorted'" local bad 1
    if `bad' display as error "  `label': state NOT restored"
    else display as text "  `label': state restored"
    c_local tt_bad `bad'
end

display as text "4. Module 11 refusals, rc 198"
local r1 `"(F1 -> y1 y2@.5 y3) (F2 -> y4 y5 y6)"'
local o1 ""
local r2 `"(F1 -> y1 y2 y3) (F2 -> y4 y5 y6)"'
local o2 "var(F1@2)"
local r3 `"(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (G -> F1 F2)"'
local o3 ""
local r4 `"(F1 -> y1 y2 y3) (F2 -> y4 y5 y6)"'
local o4 "cov(e.y1*e.y2)"
local r5 `"(F1 -> y1 y2@a y3@a) (F2 -> y4 y5 y6)"'
local o5 ""
local r6 `"(F1 -> y1 y2 y3) (F2 -> y4 y5 y6)"'
local o6 "var(e.y1@.5)"
local r7 `"(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F1 <- x)"'
local o7 ""
local r8 `"(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (y6 <- y1)"'
local o8 ""
local w1 "non-unit loading as scale (D3)"
local w2 "non-unit variance as scale (D3)"
local w3 "second-order factor (D3)"
local w4 "residual covariance among indicators (v1b)"
local w5 "a constraint tying two loadings"
local w6 "a fixed indicator residual variance"
local w7 "an observed exogenous variable"
local w8 "a path between observed variables"
local r9 `"(F1 -> y1 y2 y3) (F2 -> y4 y5 y6)"'
local o9 "group(g)"
local w9 "more than one group"
forvalues i = 1/9 {
    capture noisily _dwlsmv_translate , paths(`r`i'') semopts(`o`i'')
    local rc = _rc
    if `rc' == 198 display as text "  4.`i' rc 198: `w`i''"
    else {
        display as error "  4.`i' FAILED with rc `rc': `w`i''"
        local ++nfail
    }
    tt_state "4.`i'" "`ecmd'" "`esig'" "`dsig'" "`rng'" "y3 y1"
    local nfail = `nfail' + `tt_bad'
}

* 4.9's message must be Module 11's, not sem's (both are rc 198)
tempfile m49
quietly log using "`m49'", text name(m49) replace
capture noisily _dwlsmv_translate , paths(`r9') semopts(`o9')
quietly log close m49
local hit 0
tempname fh
file open `fh' using "`m49'", read text
file read `fh' line
while r(eof) == 0 {
    if strpos(`"`macval(line)'"', "more than one group; dwlsmv fits a single group") local hit 1
    file read `fh' line
}
file close `fh'
if `hit' display as text "  4.9 message is Module 11's"
else {
    display as error "  4.9 FAILED: Module 11's multi-group message not shown"
    local ++nfail
}

display as text "5. sem's refusals keep rc 498"
capture noisily _dwlsmv_translate , paths((F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)) semopts(cov(F1*e.F2))
local rc = _rc
if `rc' == 498 display as text "  5.1 rc 498: cov(F1*e.F2) (D4)"
else {
    display as error "  5.1 FAILED with rc `rc': cov(F1*e.F2)"
    local ++nfail
}
tt_state "5.1" "`ecmd'" "`esig'" "`dsig'" "`rng'" "y3 y1"
local nfail = `nfail' + `tt_bad'
capture noisily _dwlsmv_translate , paths((F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)) semopts(var(F2@1))
local rc = _rc
if `rc' == 498 display as text "  5.2 rc 498: var(F2@1) on an endogenous factor"
else {
    display as error "  5.2 FAILED with rc `rc': var(F2@1) on an endogenous factor"
    local ++nfail
}

display as text "6. state after a success"
quietly _dwlsmv_translate , paths((F1 -> y1 y2 y3) (F2 -> y4 y5 y6)) metric(unitvariance)
tt_state "6" "`ecmd'" "`esig'" "`dsig'" "`rng'" "y3 y1"
local nfail = `nfail' + `tt_bad'

* ---- 7. a fixed value written as canonical text reads back identically ----
* %21.0g alone does not below about 1e-3 (measured 2026-10-07)
mata:
real scalar tt_num()
{
    real colvector x
    real scalar i, bad
    string scalar s

    x   = (1 \ 0.5 \ -0.25 \ 1/3 \ (1/3) * 1e-3 \ (1/3) * 1e-4 \ -(2/3) * 1e-9 \ 123456.789012345678)
    bad = 0
    for (i = 1; i <= rows(x); i++) {
        s = dwlsmv_tr_num(x[i])
        if (strtoreal(s) != x[i]) {
            printf("{err}  7. %21x written as %s reads back as %21x\n", x[i], s, strtoreal(s))
            bad++
        }
        if (strpos(s, " ")) {
            printf("{err}  7. %s carries a space\n", s)
            bad++
        }
    }
    printf("{txt}  7. dwlsmv_tr_num(): %f of %f values read back identically\n", rows(x) - bad, rows(x))
    return(bad)
}
end
mata: st_local("t7", strofreal(tt_num()))
local nfail = `nfail' + `t7'

adopath - "`here'/src/ado"
adopath - "`here'/build"

if `nfail' {
    display as error _n "test_translate: `nfail' failure(s)"
    exit 9
}
display as result "test_translate: sem's reading translates to the expected canonical text under both metrics, Module 11 refuses what v1a does not fit at rc 198, sem's refusals keep rc 498, and data, sort, rngstate and e() are restored"
exit 0
