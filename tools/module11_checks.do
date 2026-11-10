*! Executed 2026-09-28, unchanged from its read-only agent draft: StataMP 18
*! (6 processors) and StataMP 16.1 (4 processors). Every M11CHECK line other
*! than line 0 is identical on the two. Results: E-080.
*!
*! dwlsmv -- Module 11 scratch checks: the six Stata checks listed at the end
*! of E-080. A diagnostic, not a test: no pass/fail.
*! Uses no dwlsmv code; none of the six checks needs it.
*!
*!     "/c/Program Files/Stata18/StataMP-64.exe" //e do module11_checks.do
*!     grep "^M11CHECK|" module11_checks.log
*!
*! Every result line is  M11CHECK|<n>|<name>|<result>|<detail>  so the Stata 16
*! and Stata 18 logs can be diffed after a grep. Anything that might fail is
*! wrapped in capture noisily and its _rc printed, so the file runs to the end.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

display "M11CHECK|0|machine|stata `c(stata_version)' `c(flavor)' MP=`c(MP)'|processors `c(processors)'|`c(os)'|`c(current_date)' `c(current_time)'"

* ---------------------------------------------------------------------------
* Common: a 4-parameter b with the proposal's stripes.
*   col 1  y1:F1        lambda, free
*   col 2  y2:F1        lambda, FIXED at 1 (the Cns row)
*   col 3  /:var(e.F2)  psi of an endogenous factor
*   col 4  /y1:cut1     threshold, gsem-style
* The proposal leaves open whether the threshold equation is "/y1" or "/"
* with column "y1:cut1"; the former is used here, and check 5 records what
* gsem actually does.
* ---------------------------------------------------------------------------
program define m11_b
    matrix m11b = (0.7, 1, 0.4, -0.5)
    matrix colnames m11b = y1:F1 y2:F1 /:var(e.F2) /y1:cut1
    matrix m11C = (0, 1, 0, 0, 1)
    matrix colnames m11C = y1:F1 y2:F1 /:var(e.F2) /y1:cut1 _r
end

* ---------------------------------------------------------------------------
* CHECK 1. ereturn post b V Cns with the stripes above, three V's:
*   1a  PD on the free block, zero row/column for the fixed parameter
*       (the layout the proposal wants for D5 "every row plus Cns")
*   1b  singular V inside the free block
*   1c  non-PD V (negative eigenvalue)
* Meaning: rc 0 on 1a with the stripes echoed intact by ereturn display
* supports D5 (every row + Cns) and D6 (names). A nonzero rc on 1a means
* either the stripes (D6) or the zero-row layout (D5) must change; the
* detail says which name survived. rc on 1b/1c tells whether ereturn post
* checks definiteness (the proposal's "unverified" 506): rc 0 there means
* Module 11 must do its own check before posting.
* ---------------------------------------------------------------------------
local vlab_a "zero-row-for-fixed"
local vlab_b "singular-free-block"
local vlab_c "nonPD"
foreach s in a b c {
    m11_b
    if "`s'" == "a" {
        matrix m11V = (.04, 0, .01, 0 \ 0, 0, 0, 0 \ .01, 0, .09, 0 \ 0, 0, 0, .02)
    }
    if "`s'" == "b" {
        matrix m11V = (.04, 0, .06, 0 \ 0, 0, 0, 0 \ .06, 0, .09, 0 \ 0, 0, 0, .02)
    }
    if "`s'" == "c" {
        matrix m11V = (.04, 0, .5, 0 \ 0, 0, 0, 0 \ .5, 0, .09, 0 \ 0, 0, 0, .02)
    }
    matrix rownames m11V = y1:F1 y2:F1 /:var(e.F2) /y1:cut1
    matrix colnames m11V = y1:F1 y2:F1 /:var(e.F2) /y1:cut1
    ereturn clear
    capture noisily ereturn post m11b m11V m11C, obs(500) properties(b V)
    local rcp = _rc
    local names ""
    local rcd = .
    if `rcp' == 0 {
        local names : colfullnames e(b)
        capture noisily ereturn display
        local rcd = _rc
    }
    display "M11CHECK|1`s'|post_`vlab_`s''|post_rc=`rcp' display_rc=`rcd'|`names'"
}

* ---------------------------------------------------------------------------
* CHECK 2. A stripe longer than 32 characters.
*   2a  column name 38 chars, in "/" equation
*   2b  equation name 36 chars
* Meaning: rc nonzero, or a truncated name in the detail, confirms the limit
* on this Stata and forces D6 to shorten or abbreviate long psi names
* (e.g. cov(e.Long,e.Other)); rc 0 with the full name intact means the limit
* does not bite here (compare 16 vs 18 before relying on it).
* ---------------------------------------------------------------------------
local longcol "cov(e.LongFactorAAAA,e.LongFactorBBB)"
local longeq  "AnEquationNameThatIsLongerThan32ch"
local n1 = strlen("`longcol'")
local n2 = strlen("`longeq'")
matrix m11L = (1, 2)
capture noisily matrix colnames m11L = /:`longcol' `longeq':F1
local rcn = _rc
local got : colfullnames m11L
matrix m11LV = I(2)
capture noisily matrix rownames m11LV = /:`longcol' `longeq':F1
capture noisily matrix colnames m11LV = /:`longcol' `longeq':F1
ereturn clear
capture noisily ereturn post m11L m11LV
local rcp = _rc
local gotp ""
if `rcp' == 0 {
    local gotp : colfullnames e(b)
}
display "M11CHECK|2|long_stripe|names_rc=`rcn' post_rc=`rcp'|lens `n1',`n2' matrix=[`got'] posted=[`gotp']"

* ---------------------------------------------------------------------------
* CHECK 3. test, lincom, nlcom on the 1a posting.
* The _b[] spelling for "/" equations is not known; several are tried.
* Meaning: rc 0 across the board supports the proposal's postestimation claim
* (needs only b, V, properties) and fixes the _b[] form users must type, a
* D6 input. A failure only on the "/" names argues for different names under
* D6; a failure on the fixed parameter y2:F1 (test should report it as
* constrained/dropped) informs D5.
* ---------------------------------------------------------------------------
m11_b
matrix m11V = (.04, 0, .01, 0 \ 0, 0, 0, 0 \ .01, 0, .09, 0 \ 0, 0, 0, .02)
matrix rownames m11V = y1:F1 y2:F1 /:var(e.F2) /y1:cut1
matrix colnames m11V = y1:F1 y2:F1 /:var(e.F2) /y1:cut1
ereturn clear
capture noisily ereturn post m11b m11V m11C, obs(500) properties(b V)
display "M11CHECK|3|repost|post_rc=`=_rc'|1a layout"

local k 0
foreach cmd in ///
    "test _b[y1:F1] = 0" ///
    "test _b[y2:F1] = 1" ///
    "test _b[/:var(e.F2)] = 0" ///
    "test _b[/var(e.F2)] = 0" ///
    "test _b[/y1:cut1] = 0" ///
    "lincom _b[y1:F1] + _b[/:var(e.F2)]" ///
    "lincom _b[/var(e.F2)] - _b[/y1:cut1]" ///
    "nlcom exp(_b[y1:F1])" ///
    "nlcom _b[y1:F1]*_b[/:var(e.F2)]" ///
    "nlcom sqrt(_b[/var(e.F2)])" {
    local ++k
    capture noisily `cmd'
    display "M11CHECK|3.`k'|postest|rc=`=_rc'|`cmd'"
}

* ---------------------------------------------------------------------------
* Simulated data for checks 4-6 (continuous y for sem, ordinal o for gsem).
* ---------------------------------------------------------------------------
clear
set seed 20260928
set obs 500
generate double f1 = rnormal()
generate double f2 = 0.5*f1 + rnormal()
forvalues j = 1/3 {
    generate double y`j' = 0.7*f1 + rnormal()
}
forvalues j = 4/6 {
    generate double y`j' = 0.7*f2 + rnormal()
}
forvalues j = 1/3 {
    generate byte o`j' = (y`j' > -0.5) + (y`j' > 0.5)
}

* ---------------------------------------------------------------------------
* CHECK 4. Layout of e(b) and e(Cns) from sem ..., noestimate.
* Meaning: if noestimate posts e(b) (and e(Cns)) with every parameter,
* including the anchored loading, and names in a stable order, sem is usable
* as the D1 test oracle for anchoring and defaults (and shows sem's own
* answer to D2 and to D5/D6 naming). If e(b) or e(Cns) is absent (rc on
* matrix list), the oracle has to read the parameter table another way.
* ---------------------------------------------------------------------------
capture noisily sem (F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1), noestimate
local rcs = _rc
display "M11CHECK|4|sem_noestimate|rc=`rcs'|"
capture noisily ereturn list
capture noisily matrix list e(b)
local rcb = _rc
local bn ""
local kb = .
if `rcb' == 0 {
    local bn : colfullnames e(b)
    local kb = colsof(e(b))
}
display "M11CHECK|4b|sem_noestimate_eb|rc=`rcb' cols=`kb'|`bn'"
capture noisily matrix list e(Cns)
local rcc = _rc
local kc = .
local rc_ = .
if `rcc' == 0 {
    local kc = colsof(e(Cns))
    local rc_ = rowsof(e(Cns))
}
display "M11CHECK|4c|sem_noestimate_Cns|rc=`rcc' rows=`rc_' cols=`kc'|"

* ---------------------------------------------------------------------------
* CHECK 5. gsem ..., oprobit coeflegend: the cutpoint names.
* Meaning: the detail gives the exact _b[] spelling of thresholds; if it is
* /o1:cut1 (equation "/o1"), the proposal's D6 threshold name stands; if it
* is something else, D6 should copy what gsem does. A 16-vs-18 difference
* would argue for choosing our own form rather than mimicking gsem.
* ---------------------------------------------------------------------------
capture noisily gsem (F -> o1 o2 o3, oprobit), coeflegend
local rcg = _rc
local cuts ""
if `rcg' == 0 {
    local gn : colfullnames e(b)
    foreach nm of local gn {
        if strpos("`nm'", "cut") {
            local cuts "`cuts' `nm'"
        }
    }
    local cuts : list retokenize cuts
}
display "M11CHECK|5|gsem_cut_names|rc=`rcg'|`cuts'"

* ---------------------------------------------------------------------------
* CHECK 6. Does sem accept cov(F1*e.F2) (psi alongside beta)?
*   6a  noestimate: is it syntactically accepted?
*   6b  estimation (iterate(50)): does sem itself flag non-identification?
* Meaning: rejected at 6a means sem treats it as illegal, which favours D4
* "reject early" for parity with sem; accepted at 6a but flagged at 6b (not
* identified / not converged) matches our rank-check story and supports D4
* "pass through"; accepted and silently estimated would be surprising and
* worth recording in the development record.
* ---------------------------------------------------------------------------
capture noisily sem (F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1), ///
    cov(F1*e.F2) noestimate
display "M11CHECK|6a|sem_cov_F1_eF2_noestimate|rc=`=_rc'|"
capture noisily sem (F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1), ///
    cov(F1*e.F2) iterate(50)
local rce = _rc
local conv = .
if `rce' == 0 {
    local conv = e(converged)
}
display "M11CHECK|6b|sem_cov_F1_eF2_estimate|rc=`rce' converged=`conv'|"

display "M11CHECK|end|done||"
