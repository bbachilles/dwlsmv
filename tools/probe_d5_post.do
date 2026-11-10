*! dwlsmv -- Module 11 D5 probe: what ereturn post accepts for e(b)'s layout
*!
*! A probe, not a test (tools/). D5 (2026-10-04, the user's): e(b) holds every
*! parameter-table row, fixed ones constrained through e(Cns), V zero-embedded,
*! e(rank) = number of free parameters, posted by Module 11. This probe records
*! Stata's behaviour, on the machine it runs on, for a toy table (thresholds,
*! a marker loading, a free loading, a factor variance):
*!   all   every row plus Cns: display, _b[] (fixed included), test, lincom,
*!         test on the fixed coefficient
*!   free  free rows only: _b[] on the fixed coefficient fails
*!   notpd a V that is not positive definite: does ereturn post refuse it?
*!   long  a 46-character stripe, /cov(e.Longfactorname01,e.Otherfactorname02)
*! A non-PD V is warned about ("nonsymmetric or highly singular") and posted,
*! rc 0, with SEs blanked: grep the full log for the warning, which the D5
*! lines below do not carry. e(rank) is not set by ereturn post (it reads
*! missing here), so Module 11
*! posts it. Run from the repository root:
*!     StataMP-64.exe /e do tools/probe_d5_post.do      (log: probe_d5_post.log)
*!     grep -E "^(====|D5)" probe_d5_post.log | diff - tools/probe_d5_post_stata18.txt
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek
clear all
set more off
set linesize 200
display "D5ENV stata_version " c(stata_version) " processors " c(processors)
set obs 100
generate byte touse = 1

program d5post, eclass
    args layout
    * a toy table: thresholds y1 (2), y2 (2); lambda y1:F1 (fixed 1, marker), y2:F1 free;
    * psi /var(F1) free
    tempname b V C
    capture drop touse
    generate byte touse = 1
    if "`layout'" == "all" {
        matrix `b' = (-.5, .5, -.4, .6, 1, .8, .3)
        matrix colnames `b' = /y1:cut1 /y1:cut2 /y2:cut1 /y2:cut2 y1:F1 y2:F1 /:var(F1)
        matrix `V' = J(7, 7, 0)
        forvalues i = 1/7 {
            if `i' != 5 matrix `V'[`i', `i'] = .01
        }
        matrix `V'[6, 7] = .002
        matrix `V'[7, 6] = .002
        matrix rownames `V' = /y1:cut1 /y1:cut2 /y2:cut1 /y2:cut2 y1:F1 y2:F1 /:var(F1)
        matrix colnames `V' = /y1:cut1 /y1:cut2 /y2:cut1 /y2:cut2 y1:F1 y2:F1 /:var(F1)
        matrix `C' = (0, 0, 0, 0, 1, 0, 0, 1)
        matrix colnames `C' = /y1:cut1 /y1:cut2 /y2:cut1 /y2:cut2 y1:F1 y2:F1 /:var(F1) _Cns:_r
        ereturn post `b' `V' `C', obs(100) esample(touse) properties(b V)
    }
    else if "`layout'" == "free" {
        matrix `b' = (-.5, .5, -.4, .6, .8, .3)
        matrix colnames `b' = /y1:cut1 /y1:cut2 /y2:cut1 /y2:cut2 y2:F1 /:var(F1)
        matrix `V' = .01 * I(6)
        matrix rownames `V' = /y1:cut1 /y1:cut2 /y2:cut1 /y2:cut2 y2:F1 /:var(F1)
        matrix colnames `V' = /y1:cut1 /y1:cut2 /y2:cut1 /y2:cut2 y2:F1 /:var(F1)
        ereturn post `b' `V', obs(100) esample(touse) properties(b V)
    }
    else if "`layout'" == "notpd" {
        matrix `b' = (.8, .3)
        matrix colnames `b' = y2:F1 /:var(F1)
        matrix `V' = (.01, .02 \ .02, .01)
        matrix rownames `V' = y2:F1 /:var(F1)
        matrix colnames `V' = y2:F1 /:var(F1)
        ereturn post `b' `V', obs(100) properties(b V)
    }
    else if "`layout'" == "long" {
        matrix `b' = (.1)
        matrix colnames `b' = "/:cov(e.Longfactorname01,e.Otherfactorname02)"
        matrix `V' = (.01)
        matrix rownames `V' = "/:cov(e.Longfactorname01,e.Otherfactorname02)"
        matrix colnames `V' = "/:cov(e.Longfactorname01,e.Otherfactorname02)"
        ereturn post `b' `V', obs(100) properties(b V)
    }
    ereturn local cmd "d5probe"
end

foreach L in all free notpd long {
    display as text "==== `L'"
    capture noisily d5post `L'
    display as result "D5POST `L' rc " _rc
    if _rc continue
    capture noisily ereturn display
    display as result "D5DISP `L' rc " _rc
    display as result "D5RANK `L' " e(rank)
    local cn : colfullnames e(b)
    display as result "D5NAMES `L' `cn'"
    if "`L'" == "all" | "`L'" == "free" {
        capture noisily display "D5B y2:F1 " _b[y2:F1]
        display as result "D5B y2:F1 rc " _rc
        capture noisily display "D5B y1:F1 " _b[y1:F1]
        display as result "D5B y1:F1 rc " _rc
        capture noisily display "D5B /y1:cut1 " _b[/y1:cut1]
        display as result "D5B /y1:cut1 rc " _rc
        capture noisily display "D5B /var(F1) " _b[/var(F1)]
        display as result "D5B /var(F1) rc " _rc
        capture noisily test [y2]F1 = 0
        display as result "D5TEST rc " _rc " chi2 " r(chi2)
        capture noisily lincom _b[y2:F1] + _b[/var(F1)]
        display as result "D5LINCOM rc " _rc " se " r(se)
        capture noisily test _b[y1:F1] = 1
        display as result "D5TESTFIXED rc " _rc
    }
}
