* Diagnostic (group surface, build step 5, 2026-10-12): esttab on dwlsmv
* results, before the help's esttab passage is written (decided 2026-10-12,
* help section option (a); the user, the same day: test esttab directly,
* including its noomitted option, since users may see a presentation that
* differs from ereturn display). Optional: runs only when estout is installed
* and is never part of the suite (2026-10-10 testing rule). Output is read and
* recorded in E-005, not asserted; the agreement lines below print
* "agree" or "DIFFER" for the reader.
*
* Parts:
*   A  G = 1 (case_E, metric(unitvariance)): default; fixed rows dropped
*      through drop() built from e(V)'s zero diagonal; equation labels; the
*      stats line stats(N T3 df_T3 p_T3); r(coefs) against e(b).
*   B  invariance(): configural and scalar side by side on mg_E: default;
*      noomitted; fixed rows dropped; varwidth(); r(coefs) against each e(b).
*   C  a 26-character group variable, so group-tagged names near the
*      32-character limit per part reach esttab.
*   D  the standardized estadd route: r(b_std) is specified but not built;
*      recorded as such.
*
* Executed under StataMP 18 on 2026-10-12 (output tools/diag_esttab_stata18.txt);
* Stata 16 NOT YET EXECUTED.
*
* SPDX-License-Identifier: GPL-3.0-or-later
* SPDX-FileCopyrightText: 2026 Erik Ruzek
version 16
clear all
set more off
set linesize 255

capture which esttab
if _rc {
    display as text "diag_esttab: estout is not installed (ssc install estout); nothing to do."
    exit 0
}
quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"

* the names of the fixed columns of e(b) (zero e(V) diagonal), as drop() wants
capture program drop de_fixed
program de_fixed, rclass
    tempname V
    matrix `V' = e(V)
    local eq : coleq `V'
    local nm : colnames `V'
    local out
    local k : word count `nm'
    forvalues i = 1/`k' {
        * the explicit eq:name form; colfullnames writes /var(F1) for /:var(F1),
        * which drop() does not find (rc 111, seen 2026-10-12)
        if `V'[`i', `i'] == 0 local out `out' `: word `i' of `eq'':`: word `i' of `nm''
    }
    return local fixed `out'
end

* r(coefs) of the last esttab against e(b) of the model in column m (b only)
capture program drop de_agree
program de_agree
    args label m
    tempname C b
    matrix `C' = r(coefs)
    matrix `b' = e(b)
    local nm : colfullnames `b'
    local rn : rowfullnames `C'
    local bad 0
    local i 0
    foreach n of local nm {
        local ++i
        local r : list posof "`n'" in rn
        if `r' == 0 {
            local ++bad
            continue
        }
        if reldif(`C'[`r', 2 * `m' - 1], `b'[1, `i']) > 1e-12 local ++bad
    }
    display as text "  `label': r(coefs) vs e(b): " cond(`bad', "DIFFER (`bad')", "agree")
end

display _n "==== A. G = 1"
import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",") stringcols(_all)
quietly destring y*, replace
quietly dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1), metric(unitvariance)
estimates store A
display _n "-- A1 default"
esttab A, se nostar
de_agree "A1" 1
display _n "-- A2 fixed rows dropped, stats line"
de_fixed
local fx `r(fixed)'
display "fixed: `fx'"
capture noisily esttab A, se nostar drop(`fx') stats(N T3 df_T3 p_T3)
display "rc " _rc

display _n "==== B. invariance(), configural and scalar"
import delimited using "fixtures/data/mg_E.csv", clear varnames(1) delimiters(",")
quietly dwlsmv (F -> y1 y2 y3 y4 y5), group(g) invariance(configural) metric(unitvariance)
estimates store Bc
quietly dwlsmv (F -> y1 y2 y3 y4 y5), group(g) invariance(scalar) metric(unitvariance)
estimates store Bs
display _n "-- B1 default"
esttab Bc Bs, se nostar mtitles(configural scalar)
display "rows of r(coefs): " rowsof(r(coefs))
de_agree "B1 scalar" 2
display _n "-- B2 noomitted"
esttab Bc Bs, se nostar noomitted mtitles(configural scalar)
display "rows of r(coefs): " rowsof(r(coefs))
display _n "-- B3 varwidth(24), stats line"
esttab Bc Bs, se nostar varwidth(24) stats(N T3 df_T3 p_T3) mtitles(configural scalar)
display _n "-- B4 fixed rows dropped (fixed in BOTH fits)"
quietly estimates restore Bc
de_fixed
local fc `r(fixed)'
quietly estimates restore Bs
de_fixed
local fs `r(fixed)'
local both : list fc & fs
* esttab matches names with the base marker removed (it prints 1.g for 1bn.g):
* drop() refuses /:var(F)#1bn.g, rc 111 (seen 2026-10-12); the bn-free form
local both : subinstr local both "bn." ".", all
display "fixed in both: `both'"
capture noisily esttab Bc Bs, se nostar varwidth(24) drop(`both') stats(N T3 df_T3 p_T3) mtitles(configural scalar)
display "rc " _rc

display _n "==== C. a 26-character group variable"
rename g grp_abcdefghijklmnopqrstuv
capture noisily dwlsmv (F -> y1 y2 y3 y4 y5), group(grp_abcdefghijklmnopqrstuv) invariance(scalar) metric(unitvariance)
display "rc " _rc
if !_rc {
    matrix b = e(b)
    local nm : colfullnames b
    local mx 0
    foreach n of local nm {
        local mx = max(`mx', strlen("`n'"))
    }
    display "longest full name: `mx' characters"
    estimates store C
    esttab C, se nostar varwidth(40)
}

display _n "==== D. standardized estadd route"
display "r(b_std) / r(V_std) are specified (spec D7, estout review) but not built; not run."

adopath - "`here'/src/ado"
adopath - "`here'/build"
