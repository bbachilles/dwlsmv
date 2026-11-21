*! Drafted by an agent; executed and revised 2026-11-20 by the main session (StataMP 18).
*! Module 11, release B (B3.1, B3.2, B3.5): second-order factors through _dwlsmv_translate.
*!
*! Drives _dwlsmv_translate on the data in memory (ordinal y1-y12, a binary x, a group
*! variable g), as test_translate.do does. The command is distributed, so the test builds
*! build/ldwlsmv.mlib first and runs against it and src/ado.
*!
*! Checks:
*!   1  G over three and over four first-order factors, metric(marker): the canonical text
*!      exactly ('metric unitvariance G', beta lines, no lambda for G, first-order factors
*!      anchored on their first items), r(injected) var(G@1), r(sofac) G, and the one printed
*!      line (B3.2).
*!   2  the same under metric(unitvariance): no printed line (D2 already injects), every
*!      factor 'metric unitvariance'.
*!   3  a user-written var(G@1) under marker: the text of 1, nothing injected, no line.
*!   4  refusals at rc 198 with their text and no printed Note line (B3.1, B3.5): a nonzero
*!      fixed path out of G (naming var(G@1)), with and without a user var(G@1) and under
*!      unitvariance; no indicator and no path; group(); a covariate; efa() (through dwlsmv).
*!   5  a fixed-zero path is no path (G's other paths keep it a second-order factor).
*!   6  a first-order model: canonical text unchanged, r(sofac) empty, no printed line;
*!      6b the same under unitvariance (no marker line on F1).
*!   7  two correlated second-order factors: the plural line, once. Every check counts
*!      both the singular and the plural wording (added after mutant T6, which printed the
*!      plural under unitvariance, survived a singular-only count).
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

set seed 20261120
quietly set obs 300
forvalues k = 1/12 {
    quietly generate byte y`k' = floor(runiform() * 3)
}
quietly generate byte x = mod(_n, 2)
quietly generate byte g = mod(_n, 2)

local nfail 0
local needle "has no observed indicator and is scaled by its variance"
local wantnote `"Note: G has no observed indicator and is scaled by its variance, var(G@1); sem's default fixes the path to its first child at 1 instead. Standardized estimates are the same."'

* tt_so K met: the model, expected canonical text and expected injection of G over K
* first-order factors with three items each (F1 -> y1 y2 y3, F2 -> y4 y5 y6, ...)
capture program drop tt_so
program tt_so
    args K met
    local model
    local gp
    local fl
    local indl
    local betal
    local laml
    local mkl
    local psil
    local uvl
    local injl "var(G@1)"
    forvalues k = 1/`K' {
        local a = 3 * `k' - 2
        local b = 3 * `k' - 1
        local c = 3 * `k'
        local model `"`model' (F`k' -> y`a' y`b' y`c')"'
        local gp `gp' F`k'
        local fl `fl' F`k'
        local indl `indl' y`a' y`b' y`c'
        local betal `betal';beta F`k' G free
        if "`met'" == "marker" {
            local laml `laml';lambda F`k' y`b' free;lambda F`k' y`c' free
            local mkl `mkl';metric marker F`k' y`a'
            local psil `psil';psi F`k' F`k' free
        }
        else {
            local laml `laml';lambda F`k' y`a' free;lambda F`k' y`b' free;lambda F`k' y`c' free
            local uvl `uvl';metric unitvariance F`k'
            local injl `injl' var(e.F`k'@1)
        }
    }
    local model `"`model' (G -> `gp')"'
    if "`met'" == "marker" {
        local exp "factor G `fl';ind `indl'`betal'`laml'`mkl';metric unitvariance G`psil'"
    }
    else {
        local exp "factor G `fl';ind `indl'`betal'`laml'`uvl';metric unitvariance G"
    }
    c_local tt_model `"`model'"'
    c_local tt_exp   `"`exp'"'
    c_local tt_inj   `"`injl'"'
end

* tt_run label model opts metric exp injexp sofexp nnote wantnote: canonical text,
* r(injected), r(sofac), and the number of printed Note lines (each exactly wantnote)
capture program drop tt_run
program tt_run
    args label model opts metric exp injexp sofexp nnote wantnote
    local needle "has no observed indicator and is scaled by its variance"
    local needle2 "have no observed indicator and are scaled by their variances"
    tempfile lg
    quietly log using "`lg'", text name(ttrun) replace
    capture noisily _dwlsmv_translate , paths(`model') semopts(`opts') metric(`metric')
    local rc = _rc
    local canon `"`r(canon)'"'
    local inj   `"`r(injected)'"'
    local sof   `"`r(sofac)'"'
    quietly log close ttrun
    local bad 0
    if `rc' {
        display as error "  `label': rc `rc'"
        c_local tt_bad 1
        exit
    }
    if `"`canon'"' != `"`exp'"' {
        display as error "  `label': canonical text differs"
        display as error "    got:  `canon'"
        display as error "    want: `exp'"
        local bad 1
    }
    if `"`inj'"' != `"`injexp'"' {
        display as error "  `label': injected `inj', want `injexp'"
        local bad 1
    }
    if `"`sof'"' != `"`sofexp'"' {
        display as error "  `label': r(sofac) `sof', want `sofexp'"
        local bad 1
    }
    local nany 0
    local nexact 0
    tempname fh
    file open `fh' using "`lg'", read text
    file read `fh' line
    while r(eof) == 0 {
        if strpos(`"`macval(line)'"', "`needle'") | strpos(`"`macval(line)'"', "`needle2'") local ++nany
        if `"`macval(line)'"' == `"`wantnote'"' local ++nexact
        file read `fh' line
    }
    file close `fh'
    if `nany' != `nnote' | `nexact' != `nnote' {
        display as error "  `label': printed lines `nany' (exact `nexact'), want `nnote'"
        local bad 1
    }
    if !`bad' display as text "  `label': ok"
    c_local tt_bad `bad'
end

* tt_refuse label model opts metric gopts want: rc 198, the text want in the output,
* and no B3.2 Note line
capture program drop tt_refuse
program tt_refuse
    args label model opts metric gopts want
    local needle "has no observed indicator and is scaled by its variance"
    local needle2 "have no observed indicator and are scaled by their variances"
    tempfile lg
    quietly log using "`lg'", text name(ttref) replace
    capture noisily _dwlsmv_translate , paths(`model') semopts(`opts') metric(`metric') `gopts'
    local rc = _rc
    quietly log close ttref
    local hit 0
    local note 0
    tempname fh
    file open `fh' using "`lg'", read text
    file read `fh' line
    while r(eof) == 0 {
        if strpos(`"`macval(line)'"', `"`want'"') local hit 1
        if strpos(`"`macval(line)'"', "`needle'") | strpos(`"`macval(line)'"', "`needle2'") local note 1
        file read `fh' line
    }
    file close `fh'
    if `rc' != 198 {
        display as error "  `label': rc `rc', want 198"
        c_local tt_bad 1
        exit
    }
    if !`hit' {
        display as error "  `label': rc 198 but the message lacks: `want'"
        c_local tt_bad 1
        exit
    }
    if `note' {
        display as error "  `label': rc 198 but the Note line was printed"
        c_local tt_bad 1
        exit
    }
    display as text "  `label': rc 198, message ok, no Note"
    c_local tt_bad 0
end

display as text "1. G over three and over four, metric(marker)"
foreach K in 3 4 {
    tt_so `K' marker
    tt_run "1.`K' G over `K', marker" `"`tt_model'"' "" marker `"`tt_exp'"' "var(G@1)" "G" 1 `"`wantnote'"'
    local nfail = `nfail' + `tt_bad'
}

display as text "2. the same, metric(unitvariance): D2 injects, no printed line"
foreach K in 3 4 {
    tt_so `K' unitvariance
    tt_run "2.`K' G over `K', unitvariance" `"`tt_model'"' "" unitvariance `"`tt_exp'"' `"`tt_inj'"' "G" 0 `"`wantnote'"'
    local nfail = `nfail' + `tt_bad'
}

display as text "3. a user-written var(G@1): the same text, no injection, no line"
tt_so 3 marker
tt_run "3 user var(G@1)" `"`tt_model'"' "var(G@1)" marker `"`tt_exp'"' "" "G" 0 `"`wantnote'"'
local nfail = `nfail' + `tt_bad'

display as text "4. refusals, rc 198 and the reason"
tt_so 3 marker
local mfix `"(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F3 -> y7 y8 y9) (G -> F1@1 F2 F3)"'
local needfix "Free the path and write var(G@1)"
tt_refuse "4a fixed nonzero path out of G (marker)" `"`mfix'"' "" marker "" "`needfix'"
local nfail = `nfail' + `tt_bad'
tt_refuse "4b ... with a user var(G@1)" `"`mfix'"' "var(G@1)" marker "" "`needfix'"
local nfail = `nfail' + `tt_bad'
tt_refuse "4c ... under unitvariance" `"`mfix'"' "" unitvariance "" "`needfix'"
local nfail = `nfail' + `tt_bad'
* B6 review (2026-11-20): nonzero means either sign. Added for its mutant T2 (> 0 for != 0 in
* dwlsmv_tr_sofac_t), which still survives: equivalent through the command, since a G whose
* paths are all fixed nonzero is refused by this rule before r(sofac) is read
tt_refuse "4c2 path out of G fixed at -1" `"(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F3 -> y7 y8 y9) (G -> F1@-1 F2 F3)"' "" marker "" "`needfix'"
local nfail = `nfail' + `tt_bad'
tt_refuse "4d no indicator, no path" `"(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (G -> F1@0)"' "" marker "" ///
    "has no observed indicator and no path to another factor"
local nfail = `nfail' + `tt_bad'
tt_refuse "4e group() with G" `"`tt_model'"' "" unitvariance "group(g) gvals(0 1)" ///
    "is not available with group()"
local nfail = `nfail' + `tt_bad'
* x covaries with exogenous G by default (the P2 refusal), so that covariance is fixed at 0
tt_refuse "4f covariate with G" `"`tt_model' (F2 <- x)"' "cov(x*G@0)" marker "" ///
    "is not available with a covariate"
local nfail = `nfail' + `tt_bad'
tt_refuse "4g covariate with G, unitvariance" `"`tt_model' (F2 <- x)"' "cov(x*G@0)" unitvariance "" ///
    "is not available with a covariate"
local nfail = `nfail' + `tt_bad'

* efa() goes through dwlsmv itself; the refusal comes before anything is estimated
tempfile lg
quietly log using "`lg'", text name(ttefa) replace
capture noisily dwlsmv `tt_model', efa(F1 F2 F3)
local rc = _rc
quietly log close ttefa
local hit 0
local note 0
tempname fh
file open `fh' using "`lg'", read text
file read `fh' line
while r(eof) == 0 {
    if strpos(`"`macval(line)'"', "is not available with a second-order factor") local hit 1
    * B3.2 (spec Module 12, Module 11 as built): a refusal prints no line (review B6, 2026-11-20)
    if strpos(`"`macval(line)'"', "has no observed indicator and is scaled by its variance") local note 1
    file read `fh' line
}
file close `fh'
if `rc' != 198 | !`hit' | `note' {
    display as error "  4h efa() with G: rc `rc', message found `hit', B3.2 Note printed `note' (must be 0)"
    local ++nfail
}
else display as text "  4h efa() with G: rc 198, message ok, no B3.2 Note"

* sem leaves a path fixed at 0 out of e(b) (run 2026-11-20), so the text has no line for it
display as text "5. a fixed-zero path is no path: G keeps its other paths"
tt_run "5 (G -> F1@0 F2 F3)" `"(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F3 -> y7 y8 y9) (G -> F1@0 F2 F3)"' "" marker ///
    "factor G F1 F2 F3;ind y1 y2 y3 y4 y5 y6 y7 y8 y9;beta F2 G free;beta F3 G free;lambda F1 y2 free;lambda F1 y3 free;lambda F2 y5 free;lambda F2 y6 free;lambda F3 y8 free;lambda F3 y9 free;metric marker F1 y1;metric marker F2 y4;metric marker F3 y7;metric unitvariance G;psi F1 F1 free;psi F2 F2 free;psi F3 F3 free" ///
    "var(G@1)" "G" 1 `"`wantnote'"'
local nfail = `nfail' + `tt_bad'

display as text "6. a first-order model is unchanged (test_translate.do 1a), no line, no sofac"
tt_run "6 first-order structural" "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)" "" marker ///
    "factor F1 F2;ind y1 y2 y3 y4 y5 y6;beta F2 F1 free;lambda F1 y2 free;lambda F1 y3 free;lambda F2 y5 free;lambda F2 y6 free;metric marker F1 y1;metric marker F2 y4;psi F1 F1 free;psi F2 F2 free" ///
    "" "" 0 `"`wantnote'"'
local nfail = `nfail' + `tt_bad'
* 6b the same under unitvariance: F1 has a path out and lambda rows but no marker line, so
* r(sofac) stays empty only if a lambda row counts as an indicator (added after mutant T4
* survived on 6, where every first-order factor has a marker line)
tt_run "6b first-order structural, unitvariance" "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)" "" unitvariance ///
    "factor F1 F2;ind y1 y2 y3 y4 y5 y6;beta F2 F1 free;lambda F1 y1 free;lambda F1 y2 free;lambda F1 y3 free;lambda F2 y4 free;lambda F2 y5 free;lambda F2 y6 free;metric unitvariance F1;metric unitvariance F2" ///
    "var(F1@1) var(e.F2@1)" "" 0 `"`wantnote'"'
local nfail = `nfail' + `tt_bad'


display as text "7. two correlated second-order factors: one Note line, in the plural"
local wantpl `"Note: G1 G2 have no observed indicator and are scaled by their variances, var(G1@1) var(G2@1); sem's default fixes each factor's path to its first child at 1 instead. Standardized estimates are the same."'
tt_run "7 G1 and G2" "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F3 -> y7 y8 y9) (F4 -> y10 y11 y12) (G1 -> F1 F2) (G2 -> F3 F4)" "" marker ///
    "factor G1 G2 F1 F2 F3 F4;ind y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 y11 y12;beta F1 G1 free;beta F2 G1 free;beta F3 G2 free;beta F4 G2 free;lambda F1 y2 free;lambda F1 y3 free;lambda F2 y5 free;lambda F2 y6 free;lambda F3 y8 free;lambda F3 y9 free;lambda F4 y11 free;lambda F4 y12 free;metric marker F1 y1;metric marker F2 y4;metric marker F3 y7;metric marker F4 y10;metric unitvariance G1;metric unitvariance G2;psi F1 F1 free;psi F2 F2 free;psi F3 F3 free;psi F4 F4 free;psi G2 G1 free" ///
    "var(G1@1) var(G2@1)" "G1 G2" 1 `"`wantpl'"'
local nfail = `nfail' + `tt_bad'

* 8. B3.10 end to end (added 2026-11-20): G over four on Case H. Before B3.10 the fit under
* unitvariance ended with every path negative (G's sign was left where the optimizer put it).
* After it, each path out of G is positive in both metrics (Case H's population paths are
* positive, so with the F's anchored the sum of standardized paths is too), and T3 is the
* value both metrics gave before the change: the reflection is a symmetry, the fit is unchanged.
display as text "8. B3.10: G over four on Case H, paths positive in V and U, T3 unchanged"
import delimited using "fixtures/data/case_H.csv", clear varnames(1) delimiters(",")
local hmodel "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F3 -> y7 y8 y9) (F4 -> y10 y11 y12) (G -> F1 F2 F3 F4)"
foreach met in marker unitvariance {
    capture noisily quietly dwlsmv `hmodel', metric(`met')
    if _rc {
        display as error "  8 `met': dwlsmv failed, rc " _rc
        local ++nfail
        continue
    }
    tempname hb
    matrix `hb' = e(b)
    local bad8 0
    forvalues k = 1/4 {
        local j = colnumb(`hb', "F`k':G")
        if missing(`j') {
            display as error "  8 `met': no column F`k':G in e(b)"
            local ++bad8
        }
        else if !(`hb'[1, `j'] > 0) {
            display as error "  8 `met': path G -> F`k' is " `hb'[1, `j'] ", not positive"
            local ++bad8
        }
    }
    if reldif(e(T3), 48.783) > 1e-4 {
        display as error "  8 `met': T3 " e(T3) ", expected 48.783"
        local ++bad8
    }
    if `bad8' local nfail = `nfail' + `bad8'
    else display as text "  8 `met': four paths positive, T3 " %7.3f e(T3)
}

adopath - "`here'/src/ado"
adopath - "`here'/build"

if `nfail' {
    display as error _n "test_translate_second_order: `nfail' failure(s)"
    exit 9
}
display as result "test_translate_second_order: G over three and four translates under both metrics, marker prints one line and a user var(G@1) none, and the refusals exit rc 198 with their reasons and no Note"
exit 0
