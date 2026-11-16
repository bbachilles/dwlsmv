*! Drafted by a read-only agent 2026-11-14; executed by the main session 2026-11-14 (StataMP 18): 4 of 4 pass;
*! mutants M1, M2, M3 below each killed (1 failure each).
*!
*! Conditional mode, the refusals of the owner's Stage 2 package (2026-11-14) that no test pins yet:
*!   1  keepbound with a covariate: rc 198 (the command's refusal; the message names the option and the covariate)
*!   2  efa() with a covariate: rc 198, same
*!   3  a second covariate through the command: rc 198 and Module 11's own phrase "more than one observed exogenous"
*!      (test_condx_input.do 6a only looks for the word "covariate")
*!   4  negative control: the same model without keepbound runs (rc 0), so 1 is the option's refusal
*!      Mutants: M1 delete "keepbound" from the foreach list at dwlsmv.ado:207 (1 must fail: the fit runs or Stage 1 differs);
*!      M2 delete "efa" from it (2 must fail); M3 make the oxv guard in dwlsmv_translate.mata `rows(oxv) > 2` (3 must fail).
*!
*! Surface refusals are rc 198 by the command's own convention (difftest uses 3498); asserted exactly.
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

capture program drop tr_run
program tr_run, rclass
    gettoken phrase 0 : 0
    tempfile lg
    quietly log using "`lg'", text replace name(trlog)
    capture noisily `0'
    local rc = _rc
    quietly log close trlog
    mata: st_local("txt", invtokens(cat(st_local("lg"))', " "))
    mata: st_local("found", strofreal(strpos(st_local("txt"), st_local("phrase")) > 0))
    return scalar rc    = `rc'
    return scalar found = `found'
end

set seed 20261114
quietly set obs 1000
quietly generate byte x = runiform() < .5
quietly generate byte w = runiform() < .3
quietly generate double eta = .4 * x + rnormal()
forvalues j = 1/4 {
    quietly generate double z`j' = .7 * eta + sqrt(.51) * rnormal()
    quietly generate byte y`j' = (z`j' > -.5) + (z`j' > .5)
}
drop eta z1 z2 z3 z4

local mB "(F -> y1 y2 y3 y4) (F <- x)"

tr_run "keepbound is not available with a covariate" dwlsmv `mB', keepbound
if r(rc) != 198 | !r(found) {
    display as error " 1. keepbound with a covariate: rc " r(rc) " (want 198); message " r(found)
    local ++nfail
}
else display as text " 1. keepbound with a covariate is refused at rc 198"

tr_run "efa is not available with a covariate" dwlsmv `mB', efa(F)
if r(rc) != 198 | !r(found) {
    display as error " 2. efa() with a covariate: rc " r(rc) " (want 198); message " r(found)
    local ++nfail
}
else display as text " 2. efa() with a covariate is refused at rc 198"

tr_run "more than one observed exogenous" dwlsmv (F -> y1 y2 y3 y4) (F <- x) (y4 <- w)
if r(rc) != 198 | !r(found) {
    display as error " 3. second covariate: rc " r(rc) " (want 198); Module 11's phrase " r(found)
    local ++nfail
}
else display as text " 3. a second covariate is refused at rc 198 with its own message"

capture noisily dwlsmv `mB'
if _rc | e(converged) != 1 {
    display as error " 4. control without keepbound: rc " _rc
    local ++nfail
}
else display as text " 4. the same model without keepbound runs"

adopath - "`here'/src/ado"
adopath - "`here'/build"

if `nfail' {
    display as error _n "test_condx_refusals: `nfail' failure(s)"
    exit 9
}
display as result "test_condx_refusals: keepbound and efa() with a covariate and a second covariate are refused at rc 198, each with its own message"
exit 0
