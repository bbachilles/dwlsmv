* Probe (2026-10-12): sem's e(b) stripe under group() with each ginvariant()
* class, against ginvariant(none), for D6 at G > 1 (spec D5 multi-group
* posting convention; the user's caveat: compatibility had been shown for
* ginvariant(none) only). For each class: the full colfullnames list, its
* length, and which columns e(Cns) fixes (and at what value), so a class that
* omits, renames, reorders or fixes intercepts or means shows up.
* Diagnostic, not a pass/fail test. Executed under StataMP 18 on 2026-10-12
* (output alongside, *_stata18.txt); Stata 16 NOT YET EXECUTED.
*
* SPDX-License-Identifier: GPL-3.0-or-later
* SPDX-FileCopyrightText: 2026 Erik Ruzek
version 16
clear all
set more off
set linesize 255
import delimited using "fixtures/data/mg_E.csv", clear varnames(1) delimiters(",")

local classes none mcoef mcons "mcoef mcons" "mcoef mcons merrvar" all
local nc : word count `classes'
forvalues c = 1/`nc' {
    local cl : word `c' of `classes'
    display _n "==== ginvariant(`cl')"
    capture noisily quietly sem (F -> y1 y2 y3 y4 y5), group(g) ginvariant(`cl') var(1: F@1) noestimate
    if _rc {
        display "rc " _rc
        continue
    }
    matrix b = e(b)
    local names : colfullnames b
    display "k = " colsof(b)
    local i 0
    foreach n of local names {
        local ++i
        display "  `i' `n'"
    }
    capture matrix C = e(Cns)
    if _rc {
        display "  no e(Cns)"
        continue
    }
    mata: C = st_matrix("C"); nm = st_matrixcolstripe("C"); ///
        for (r = 1; r <= rows(C); r++) { ///
            sel = selectindex(C[r, 1..cols(C)-1] :!= 0); ///
            s = ""; ///
            for (k = 1; k <= cols(sel); k++) s = s + " " + strofreal(C[r, sel[k]]) + "*" + nm[sel[k], 1] + ":" + nm[sel[k], 2]; ///
            printf("  Cns %f:%s = %g\n", r, s, C[r, cols(C)]); ///
        }
}
