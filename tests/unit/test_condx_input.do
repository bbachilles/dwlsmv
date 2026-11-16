*! Executed by the main session (Stata 18 MP) and mutant-tested, 2026-11-14 at the latest.
*! Conditional mode, input path through the command: dwlsmv (F -> y1 y2 y3 y4) (F <- x) (y4 <- x)
*! reads the covariate from the sem-style paths, checks it, runs Stage 1 conditional on it and
*! posts the canonical text with the covariate line. Drafted by a read-only agent 2026-11-14;
*! target path tests/unit/test_condx_input.do. Remove this banner only after the file has been run.
*!
*! Needs the input-side edits (edits.md of the input draft) AND the output-side edits (the output
*! draft's edits.md): the display, e() and difftest refusal are the output side's, and
*! test_condx_output.do covers them. This file covers what the user types and what is refused
*! before or at Stage 1. It asserts no estimate except signs, whose sizes are not the point.
*!
*! Data (generated here, fixed seed, no fixture files, nothing from the frozen condx files):
*! N = 1000, x = Bernoulli(.5), one factor with loadings .7, the factor shifted .4 by x, a direct
*! effect +.3 of x on y4*, four items with thresholds -.5 and .5; g a binary group variable.
*!
*! Checks (surface refusals are rc 198 and name the reason; Mata's own refusal is rc 3498)
*!   1  the headline model runs: converged, e(N) = 1000, e(canon) carries the covariate line and
*!      the gamma and kappa lines, F:x and y4:x are positive (the data's directions).
*!   2  metric(unitvariance): canonical text carries metric unitvariance F, the gamma is finite.
*!   3  (F <- x@.4): gamma fixed at .4 in the canonical text and _b[F:x] equals .4 exactly.
*!   4  a missing x is deleted listwise with the indicators: e(N) = 990 when ten x are missing,
*!      and e(sample) agrees.
*!   5  x must be numeric (rc 198), take only 0 and 1 (rc 198, coded 0/2), and take both (rc 198,
*!      constant); each message names the reason (P5).
*!   6  surface refusals, rc 198 with their reasons: two covariates (P1), sem's free cov(x,F)
*!      with only (y4 <- x) (P2), (y4 <- y1) beside a covariate, a covariate with group(),
*!      an equality constraint through the covariate (P3).
*!   7  Stage 1's own refusal, rc 3498: an item separated by x. The command prints the legend,
*!      does not print the v1 advice about keepbound, and posts nothing.
*!   8  negative control: a model without a covariate still runs on the same data, with no
*!      covariate line in e(canon).
*!
*! Provisional choices under test (the owner has not confirmed): P1, P2, P3, P5.
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

local nfail 0

* run a command capturing what it prints: r(rc), r(found) for the first phrase, r(found2) for
* the second
capture program drop ti_run
program ti_run, rclass
    gettoken phrase 0 : 0
    gettoken phrase2 0 : 0
    tempfile lg
    quietly log using "`lg'", text replace name(tilog)
    capture noisily `0'
    local rc = _rc
    quietly log close tilog
    mata: st_local("txt", invtokens(cat(st_local("lg"))', " "))
    mata: st_local("found",  strofreal(strpos(st_local("txt"), st_local("phrase"))  > 0))
    mata: st_local("found2", strofreal(strpos(st_local("txt"), st_local("phrase2")) > 0))
    return scalar rc     = `rc'
    return scalar found  = `found'
    return scalar found2 = `found2'
end

* ---- the data ---------------------------------------------------------------------------
set seed 20261114
quietly set obs 1000
quietly generate byte x = runiform() < .5
quietly generate byte g = _n > 500
quietly generate double eta = .4 * x + rnormal()
forvalues j = 1/4 {
    quietly generate double z`j' = .7 * eta + sqrt(.51) * rnormal() + (`j' == 4) * .3 * x
    quietly generate byte y`j' = (z`j' > -.5) + (z`j' > .5)
}
drop eta z1 z2 z3 z4
quietly generate byte z = runiform() < .3

local mB "(F -> y1 y2 y3 y4) (F <- x) (y4 <- x)"

* ---- 1. the headline model ---------------------------------------------------------------
capture noisily dwlsmv `mB'
local rc1 = _rc
if `rc1' {
    display as error " 1. headline model: rc `rc1'"
    local ++nfail
}
else {
    local canon `"`e(canon)'"'
    mata: st_local("hasc", strofreal(strpos(st_local("canon"), "covariate x") > 0))
    mata: st_local("hasg", strofreal(strpos(st_local("canon"), "gamma F x free") > 0))
    mata: st_local("hask", strofreal(strpos(st_local("canon"), "kappa y4 x free") > 0))
    if e(converged) != 1 | e(N) != 1000 | !`hasc' | !`hasg' | !`hask' {
        display as error " 1. converged " e(converged) ", N " e(N) ", covariate line `hasc', gamma `hasg', kappa `hask'"
        display as error "    canon: `canon'"
        local ++nfail
    }
    else if _b[F:x] <= 0 | _b[y4:x] <= 0 {
        display as error " 1. signs: F:x " _b[F:x] ", y4:x " _b[y4:x] " (the data's directions are positive)"
        local ++nfail
    }
    else display as text " 1. headline model runs; e(canon) carries the covariate, gamma and kappa lines"
}

* ---- 2. unitvariance ---------------------------------------------------------------------
capture noisily dwlsmv `mB', metric(unitvariance)
local rc2 = _rc
if `rc2' {
    display as error " 2. unitvariance: rc `rc2'"
    local ++nfail
}
else {
    local canon `"`e(canon)'"'
    mata: st_local("hasm", strofreal(strpos(st_local("canon"), "metric unitvariance F") > 0))
    if !`hasm' | missing(_b[F:x]) {
        display as error " 2. unitvariance: canon `canon'; F:x " _b[F:x]
        local ++nfail
    }
    else display as text " 2. metric(unitvariance): metric unitvariance F in e(canon), gamma finite"
}

* ---- 3. a fixed gamma --------------------------------------------------------------------
capture noisily dwlsmv (F -> y1 y2 y3 y4) (F <- x@.4)
local rc3 = _rc
if `rc3' {
    display as error " 3. (F <- x@.4): rc `rc3'"
    local ++nfail
}
else {
    local canon `"`e(canon)'"'
    mata: st_local("hasf", strofreal(strpos(st_local("canon"), "gamma F x fixed .4") > 0))
    if !`hasf' | _b[F:x] != .4 {
        display as error " 3. gamma fixed: canon `canon'; _b[F:x] = " _b[F:x]
        local ++nfail
    }
    else display as text " 3. (F <- x@.4): gamma fixed at .4 in e(canon) and in e(b)"
}

* ---- 4. missing x, listwise ----------------------------------------------------------------
preserve
quietly replace x = . in 1/10
capture noisily dwlsmv (F -> y1 y2 y3 y4) (F <- x)
local rc4 = _rc
if `rc4' {
    display as error " 4. missing x: rc `rc4'"
    local ++nfail
}
else {
    quietly count if e(sample)
    if e(N) != 990 | r(N) != 990 {
        display as error " 4. missing x: e(N) " e(N) ", e(sample) " r(N) " (want 990, 990)"
        local ++nfail
    }
    else display as text " 4. ten missing x: e(N) = 990 and e(sample) agrees"
}
restore

* ---- 5. x must be numeric, 0/1, both values --------------------------------------------------
preserve
quietly generate str1 xs = string(x)
ti_run "numeric" "covariate" dwlsmv (F -> y1 y2 y3 y4) (F <- xs)
if r(rc) != 198 | !r(found) {
    display as error " 5a. string x: rc " r(rc) ", message " r(found) " (want 198 and the word numeric)"
    local ++nfail
}
else display as text " 5a. a string covariate is refused at rc 198"
restore

preserve
quietly replace x = 2 * x
ti_run "0 and 1" "covariate" dwlsmv (F -> y1 y2 y3 y4) (F <- x)
if r(rc) != 198 | !r(found) {
    display as error " 5b. x coded 0/2: rc " r(rc) ", message " r(found) " (want 198 and 0 and 1)"
    local ++nfail
}
else display as text " 5b. x coded 0/2 is refused at rc 198 (the command, before Stage 1)"
restore

preserve
quietly replace x = 1
ti_run "both" "covariate" dwlsmv (F -> y1 y2 y3 y4) (F <- x)
if r(rc) != 198 | !r(found) {
    display as error " 5c. constant x: rc " r(rc) ", message " r(found) " (want 198 and both)"
    local ++nfail
}
else display as text " 5c. a constant covariate is refused at rc 198"
restore

* ---- 6. surface refusals -------------------------------------------------------------------
ti_run "more than one observed exogenous" "covariate" dwlsmv (F -> y1 y2 y3 y4) (F <- x) (y4 <- z)
if r(rc) != 198 | !r(found) {
    display as error " 6a. two covariates: rc " r(rc) ", message " r(found)
    local ++nfail
}
else display as text " 6a. two covariates: rc 198"
ti_run "(F <- x)" "cov(x*F@0)" dwlsmv (F -> y1 y2 y3 y4) (y4 <- x)
if r(rc) != 198 | !r(found) | !r(found2) {
    display as error " 6b. (y4 <- x) alone: rc " r(rc) ", names (F <- x) " r(found) ", names cov(x*F@0) " r(found2)
    local ++nfail
}
else display as text " 6b. (y4 <- x) alone: rc 198, the message gives both ways out"
ti_run "a path between observed variables (y4 <- y1)" "x" dwlsmv (F -> y1 y2 y3 y4) (F <- x) (y4 <- y1)
if r(rc) != 198 | !r(found) {
    display as error " 6c. (y4 <- y1) beside a covariate: rc " r(rc) ", message " r(found)
    local ++nfail
}
else display as text " 6c. (y4 <- y1) is still refused beside a covariate"
ti_run "group()" "covariate" dwlsmv (F -> y1 y2 y3 y4) (F <- x), group(g) invariance(configural) metric(unitvariance)
if r(rc) != 198 | !r(found) {
    display as error " 6d. covariate with group(): rc " r(rc) ", message " r(found)
    local ++nfail
}
else display as text " 6d. a covariate with group() is refused at rc 198"
ti_run "ties several parameters" "x" dwlsmv (F -> y1 y2 y3 y4) (F <- x@a) (y4 <- x@a)
if r(rc) != 198 | !r(found) {
    display as error " 6e. equality through the covariate: rc " r(rc) ", message " r(found)
    local ++nfail
}
else display as text " 6e. an equality constraint through the covariate is refused at rc 198"

* ---- 7. Stage 1's refusal: an item separated by x, rc 3498 -----------------------------------
preserve
quietly replace y1 = min(y1, 1) if x == 0
quietly replace y1 = max(y1, 1) if x == 1
ereturn clear
ti_run "separation by the covariate" "keepbound" dwlsmv (F -> y1 y2 y3 y4) (F <- x)
local rc7 = r(rc)
local f7  = r(found)
local k7  = r(found2)
if `rc7' != 3498 | !`f7' | `k7' | "`e(cmd)'" != "" {
    display as error " 7. separated item: rc `rc7' (want 3498); Stage 1's message " `f7' "; keepbound advice shown " `k7' "; e(cmd) '`e(cmd)''"
    local ++nfail
}
else display as text " 7. an item separated by x: rc 3498, Stage 1's message, no keepbound advice, nothing posted"
restore

* ---- 8. negative control -------------------------------------------------------------------------
capture noisily dwlsmv (F -> y1 y2 y3 y4)
local rc8 = _rc
if `rc8' {
    display as error " 8. no covariate: rc `rc8'"
    local ++nfail
}
else {
    local canon `"`e(canon)'"'
    mata: st_local("hasc", strofreal(strpos(st_local("canon"), "covariate") > 0))
    if `hasc' | e(N) != 1000 {
        display as error " 8. no covariate: covariate line `hasc', N " e(N)
        local ++nfail
    }
    else display as text " 8. a model without a covariate runs as before, no covariate line"
}

adopath - "`here'/src/ado"
adopath - "`here'/build"

if `nfail' {
    display as error _n "test_condx_input: `nfail' failure(s)"
    exit 9
}
display as result "test_condx_input: a binary covariate typed as sem paths reaches a conditional fit; a non-binary, constant, string, second or group-crossing covariate, sem's free cov(x,F) and an equality constraint are refused at rc 198, and an item separated by x at rc 3498"
exit 0
