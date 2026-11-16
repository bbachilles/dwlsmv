*! Executed 2026-11-14 by the main session (Stata 18 MP), mutant-tested.
*! Module 11, conditional mode: _dwlsmv_translate reads one binary covariate x as gamma and kappa
*! lines. Drafted by a read-only agent 2026-11-14; target path tests/unit/test_translate_cx.do.
*! Remove this banner only after the file has been run and what broke has been fixed.
*!
*! The command is distributed, so this builds build/ldwlsmv.mlib from src/mata first (as
*! test_translate.do does) and runs against the library and src/ado.
*!
*! The translator never reads the user's values (it runs sem on dummy data), so the data here
*! carry names only: y1-y4, x, z (a second covariate), g (a group variable).
*!
*! Checks
*!   1  translation, exact canonical text, r(cov) and r(injected):
*!      1a (F <- x)                         gamma F x free
*!      1b (F <- x) (y4 <- x)               gamma and kappa
*!      1c both, metric(unitvariance)       injected var(e.F@1), metric unitvariance F
*!      1d (F <- x@.4) (y4 <- x@-.3)        gamma F x fixed .4, kappa y4 x fixed -.3 (P3)
*!      1e (y4 <- x) with cov(x*F@0)        kappa only; no gamma row for F (P2)   [needs probe p1]
*!      1f no covariate: r(cov) is empty, the canonical text has no covariate line
*!   2  refusals, each rc 198 exactly, each with its own message:
*!      2a two covariates (P1)              2b (y4 <- x) alone: sem's cov(x,F) (P2)
*!      2c (y4 <- y1) beside a covariate    2d covariate with group()
*!      2e var(x@1)                         2f an equality constraint through the covariate (P3)
*!      2g cov(x*F@.2): a covariance fixed at a value other than 0   [needs probe p2]
*!      2h two factors, x predicts F2 only: cov(x,F1) is listed for F1 (P2)
*!   3  state restoration: the data signature and the rngstate are unchanged after a covariate
*!      success and after a covariate refusal.
*!
*! Provisional choices under test (the owner has not confirmed them): P1 one covariate, P2 the
*! covariance refusal, P3 fixed values literal and equality refused. If the owner changes one,
*! this is the file that changes with it.
*!
*! Cases that depend on a sem fact the probe settles (probe_cov.do): 1e, 2e, 2g. If sem
*! rejects the syntax the case reports rc and the message says so; adjust the case, not the code.
*!
*! Return codes: 0 pass, 9 fail.
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

set seed 4711
quietly set obs 60
quietly generate byte x = mod(_n, 2)
quietly generate byte z = mod(_n, 3) > 0
quietly generate byte g = _n > 30
forvalues j = 1/4 {
    quietly generate byte y`j' = mod(_n + `j', 3)
}

local nfail 0

* a translation that must succeed: exact text, r(cov), r(injected)
capture program drop tx_ok
program tx_ok
    args label model opts metric expect expcov expinj
    capture noisily _dwlsmv_translate , paths(`model') semopts(`opts') metric(`metric')
    if _rc {
        display as error "  `label': rc " _rc
        c_local tx_bad 1
        exit
    }
    if `"`r(canon)'"' != `"`expect'"' | `"`r(cov)'"' != `"`expcov'"' | `"`r(injected)'"' != `"`expinj'"' {
        display as error "  `label': differs"
        display as error "    got:  `r(canon)'  | cov `r(cov)' | injected `r(injected)'"
        display as error "    want: `expect'  | cov `expcov' | injected `expinj'"
        c_local tx_bad 1
        exit
    }
    display as text "  `label': ok"
    c_local tx_bad 0
end

* a refusal: rc 198 exactly, and the message contains the phrase (rc alone cannot tell sem's
* own 198 from Module 11's)
capture program drop tx_refuse
program tx_refuse
    args label phrase model opts extra
    tempfile lg
    quietly log using "`lg'", text replace name(txlog)
    capture noisily _dwlsmv_translate , paths(`model') semopts(`opts') `extra'
    local rc = _rc
    quietly log close txlog
    mata: st_local("txt", invtokens(cat(st_local("lg"))', " "))
    mata: st_local("found", strofreal(strpos(st_local("txt"), st_local("phrase")) > 0))
    if `rc' != 198 | !`found' {
        display as error "  `label': rc `rc' (want 198); phrase found: `found' (`phrase')"
        c_local tx_bad 1
        exit
    }
    display as text "  `label': rc 198, message names it"
    c_local tx_bad 0
end

display as text "1. translation"
tx_ok "1a factor on x" "(F -> y1 y2 y3 y4) (F <- x)" "" "marker" ///
    "factor F;ind y1 y2 y3 y4;covariate x;gamma F x free;lambda F y2 free;lambda F y3 free;lambda F y4 free;metric marker F y1;psi F F free" ///
    "x" ""
local nfail = `nfail' + `tx_bad'
tx_ok "1b factor and indicator on x" "(F -> y1 y2 y3 y4) (F <- x) (y4 <- x)" "" "marker" ///
    "factor F;ind y1 y2 y3 y4;covariate x;gamma F x free;kappa y4 x free;lambda F y2 free;lambda F y3 free;lambda F y4 free;metric marker F y1;psi F F free" ///
    "x" ""
local nfail = `nfail' + `tx_bad'
tx_ok "1c both, unitvariance" "(F -> y1 y2 y3 y4) (F <- x) (y4 <- x)" "" "unitvariance" ///
    "factor F;ind y1 y2 y3 y4;covariate x;gamma F x free;kappa y4 x free;lambda F y1 free;lambda F y2 free;lambda F y3 free;lambda F y4 free;metric unitvariance F" ///
    "x" "var(e.F@1)"
local nfail = `nfail' + `tx_bad'
tx_ok "1d fixed values" "(F -> y1 y2 y3 y4) (F <- x@.4) (y4 <- x@-.3)" "" "marker" ///
    "factor F;ind y1 y2 y3 y4;covariate x;gamma F x fixed .4;kappa y4 x fixed -.3;lambda F y2 free;lambda F y3 free;lambda F y4 free;metric marker F y1;psi F F free" ///
    "x" ""
local nfail = `nfail' + `tx_bad'
tx_ok "1e kappa only, cov(x*F@0)" "(F -> y1 y2 y3 y4) (y4 <- x)" "cov(x*F@0)" "marker" ///
    "factor F;ind y1 y2 y3 y4;covariate x;kappa y4 x free;lambda F y2 free;lambda F y3 free;lambda F y4 free;metric marker F y1;psi F F free" ///
    "x" ""
local nfail = `nfail' + `tx_bad'
tx_ok "1f no covariate" "(F -> y1 y2 y3 y4)" "" "marker" ///
    "factor F;ind y1 y2 y3 y4;lambda F y2 free;lambda F y3 free;lambda F y4 free;metric marker F y1;psi F F free" ///
    "" ""
local nfail = `nfail' + `tx_bad'

display as text "2. refusals, rc 198"
tx_refuse "2a two covariates" "more than one observed exogenous" "(F -> y1 y2 y3 y4) (F <- x) (y4 <- z)" ""
local nfail = `nfail' + `tx_bad'
tx_refuse "2b sem's cov(x,F)" "(F <- x)" "(F -> y1 y2 y3 y4) (y4 <- x)" ""
local nfail = `nfail' + `tx_bad'
tx_refuse "2b2 names the zero form" "cov(x*F@0)" "(F -> y1 y2 y3 y4) (y4 <- x)" ""
local nfail = `nfail' + `tx_bad'
tx_refuse "2c indicator on indicator" "a path between observed variables (y4 <- y1)" ///
    "(F -> y1 y2 y3 y4) (F <- x) (y4 <- y1)" ""
local nfail = `nfail' + `tx_bad'
tx_refuse "2d covariate with group()" "group()" "(F -> y1 y2 y3 y4) (F <- x)" "" ///
    "metric(unitvariance) group(g) gvals(0 1)"
local nfail = `nfail' + `tx_bad'
tx_refuse "2e var(x@1)" "variance or mean of the covariate" "(F -> y1 y2 y3 y4) (F <- x)" "var(x@1)"
local nfail = `nfail' + `tx_bad'
tx_refuse "2f equality constraint" "ties several parameters" "(F -> y1 y2 y3 y4) (F <- x@a) (y4 <- x@a)" ""
local nfail = `nfail' + `tx_bad'
tx_refuse "2g cov(x*F@.2)" "fixed at zero" "(F -> y1 y2 y3 y4) (y4 <- x)" "cov(x*F@.2)"
local nfail = `nfail' + `tx_bad'
tx_refuse "2h unpredicted exogenous factor" "F1" "(F1 -> y1 y2) (F2 -> y3 y4) (F2 <- x)" ""
local nfail = `nfail' + `tx_bad'

display as text "3. state restoration"
quietly _datasignature
local dsig `r(datasignature)'
local rng `c(rngstate)'
quietly _dwlsmv_translate , paths((F -> y1 y2 y3 y4) (F <- x) (y4 <- x)) metric(unitvariance)
quietly _datasignature
local bad = ("`r(datasignature)'" != "`dsig'") | ("`c(rngstate)'" != "`rng'")
capture quietly _dwlsmv_translate , paths((F -> y1 y2 y3 y4) (F <- x) (y4 <- z))
local rc = _rc
quietly _datasignature
local bad = `bad' | ("`r(datasignature)'" != "`dsig'") | ("`c(rngstate)'" != "`rng'")
if `rc' != 198 | `bad' {
    display as error "  3. covariate refusal rc `rc'; state not restored: `bad'"
    local ++nfail
}
else display as text "  3. data signature and rngstate restored after a covariate success and a refusal"

adopath - "`here'/src/ado"
adopath - "`here'/build"

if `nfail' {
    display as error _n "test_translate_cx: `nfail' failure(s)"
    exit 9
}
display as result "test_translate_cx: one binary covariate is translated to covariate, gamma and kappa lines, fixed values literally; two covariates, sem's free cov(x,F), a constrained variance, an equality constraint and a covariate with group() are refused at rc 198"
exit 0
