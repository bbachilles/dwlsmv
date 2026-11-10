*! Module 10: the small-count family, dwlsmv_countdiag(Y, K).
*!
*! Executed under StataMP 18 on 2026-09-26 and it passes. Before the routine
*! existed it failed to load, r(9), as a test written first should.
*!
*! Spec Section 5, Module 10 ("Small counts are a separate family") and
*! E-060 and E-071. Values only: every row's flag is .n,
*! "evaluated, no flag defined in v1", which is NOT the same as `.`, "not
*! evaluated". The checks therefore assert flag == .n exactly; missing()
*! cannot tell the two apart.
*!
*!   1. A HAND-BUILT Y whose counts are obvious by inspection, with a unique
*!      minimum in each row: boundary (y2's top category, count 1) and
*!      interior (y1's middle category, count 1). Location is (variable,
*!      category position 1..K_j).
*!   2. TIES go to the first (variable, category) in s's order.
*!   3. ALL BINARY: no interior category exists. Row 2 is evaluated, has
*!      nothing to report, and says so: value missing, i = j = 0, flag .n.
*!   4. AN EMPTY DECLARED CATEGORY reads 0 and does not abort. Module 1 aborts
*!      on it; a diagnostic reports it.
*!   5. CASE D against an independent path: Stata's own -count-, per variable
*!      and category, never dwlsmv_catcounts().
*!   6. ERROR PATHS, rc 3498 each: a missing value, a non-integer code, a code
*!      >= K, a negative code, and K of the wrong length.
*!
*! Reference CSVs: explicit delimiters(","), for the delimiter import trap.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off

foreach f in tests/dwlsmv_assert.mata               ///
             src/mata/dwlsmv_catcounts.mata         ///
             src/mata/dwlsmv_countdiag.mata {
    capture noisily do "`f'"
    if _rc {
        display as error "test_countdiag: `f' failed to load"
        exit 9
    }
}

* ---- 5's independent reference, from Stata rather than Mata --------------------
* Boundary: categories 0 and K-1. Interior: 1 .. K-2. Case D's K is (4,4,3,3).
import delimited using "fixtures/data/case_D.csv", clear varnames(1) delimiters(",")
local K1 4
local K2 4
local K3 3
local K4 3
local bmin = .
local imin = .
forvalues v = 1/4 {
    local top = `K`v'' - 1
    forvalues c = 0/`top' {
        quietly count if y`v' == `c'
        if inlist(`c', 0, `top') {
            if r(N) < `bmin' {
                local bmin = r(N)
                local bi `v'
                local bj = `c' + 1
            }
        }
        else if r(N) < `imin' {
            local imin = r(N)
            local ii `v'
            local ij = `c' + 1
        }
    }
}
display as text "  case D by -count-: boundary `bmin' at (`bi', `bj'), interior `imin' at (`ii', `ij')"
mata: YD = st_data(., .)
mata: REFD = (.n, `bmin', `bi', `bj' \ .n, `imin', `ii', `ij')

mata:

real scalar dwlsmv_t_countdiag(real matrix YD, real matrix REFD)
{
    string colvector FAIL
    real matrix A, Y
    real scalar i

    FAIL = J(0, 1, "")

    // ---- 1. hand-built ---------------------------------------------------------
    // y1 (K 3) counts 5 1 4; y2 (K 4) counts 4 2 3 1; y3 (K 2) counts 8 2.
    // Boundary counts 5 4 | 4 1 | 8 2 -> 1 at (2, 4). Interior 1 | 2 3 -> 1 at (1, 2).
    Y = (0, 0, 0 \ 0, 0, 0 \ 0, 0, 0 \ 0, 0, 0 \ 0, 1, 0 \
         1, 1, 0 \ 2, 2, 0 \ 2, 2, 0 \ 2, 2, 1 \ 2, 3, 1)
    A = dwlsmv_countdiag(Y, (3 \ 4 \ 2))
    printf("{txt}  1. hand-built     : boundary %f at (%f, %f), interior %f at (%f, %f)\n",
           A[1, 2], A[1, 3], A[1, 4], A[2, 2], A[2, 3], A[2, 4])
    if (rows(A) != 2 | cols(A) != 4) {
        FAIL = FAIL \ "1: countdiag is not 2 x 4"
    }
    else {
        if (A[1, 2..4] != (1, 2, 4)) FAIL = FAIL \ "1: boundary row is not (1, 2, 4)"
        if (A[2, 2..4] != (1, 1, 2)) FAIL = FAIL \ "1: interior row is not (1, 1, 2)"
        for (i = 1; i <= 2; i++) {
            if (A[i, 1] != .n) FAIL = FAIL \ ("1: row " + strofreal(i) + "'s flag is not .n")
        }
    }

    // ---- 2. ties: the first in (variable, category) order ----------------------
    // y1 counts 2 2 (both boundary), y2 counts 2 2: four-way tie at 2 -> (1, 1).
    A = dwlsmv_countdiag((0, 0 \ 0, 0 \ 1, 1 \ 1, 1), (2 \ 2))
    if (A[1, 2..4] != (2, 1, 1)) {
        FAIL = FAIL \ "2: a tie did not go to the first variable and category"
    }

    // interior ties too (added after the 2026-09-26 review, M23): y1 and y2
    // (K 3) each hold one observation in their interior category.
    A = dwlsmv_countdiag((0, 0 \ 1, 1 \ 2, 2 \ 2, 2 \ 0, 0), (3 \ 3))
    if (A[2, 2..4] != (1, 1, 2)) {
        FAIL = FAIL \ "2: an interior tie did not go to the first variable"
    }

    // ---- 3. all binary: the interior row has nothing to report --------------------
    A = dwlsmv_countdiag((0, 1 \ 1, 0 \ 1, 1), (2 \ 2))
    printf("{txt}  3. all binary     : interior row (%f, %f, %f, %f)\n", A[2, 1], A[2, 2], A[2, 3], A[2, 4])
    if (!(A[2, 1] == .n & A[2, 2] == . & A[2, 3] == 0 & A[2, 4] == 0)) {
        FAIL = FAIL \ "3: all-binary interior row is not (.n, ., 0, 0)"
    }

    // ---- 4. an empty declared category reads 0 ------------------------------------
    A = dwlsmv_countdiag((0, 1 \ 1, 0 \ 1, 1), (3 \ 2))
    if (A[1, 2..4] != (0, 1, 3)) {
        FAIL = FAIL \ "4: y1's empty top category did not read 0 at (1, 3)"
    }

    // ---- 5. case D against Stata's -count- -------------------------------------------
    A = dwlsmv_countdiag(YD, (4 \ 4 \ 3 \ 3))
    printf("{txt}  5. case D         : boundary %f at (%f, %f), interior %f at (%f, %f)\n",
           A[1, 2], A[1, 3], A[1, 4], A[2, 2], A[2, 3], A[2, 4])
    if (!(A[1, 1] == .n & A[2, 1] == .n & A[., 2..4] == REFD[., 2..4])) {
        FAIL = FAIL \ "5: case D disagrees with Stata's own counts"
    }

    printf("{txt}test_countdiag: %f failure(s)\n", rows(FAIL))
    for (i = 1; i <= rows(FAIL); i++) {
        printf("{err}    %s\n", FAIL[i])
    }
    return(rows(FAIL))
}

end

mata: st_local("nfail", strofreal(dwlsmv_t_countdiag(YD, REFD)))

* ---------------------------------------------------------------------------
* 6. ERROR PATHS -- rc 3498 each.
* ---------------------------------------------------------------------------
local nerr 0
local bad1 "(0, 1 \ ., 1)"
local bad2 "(0, 1 \ 0.5, 1)"
local bad3 "(0, 1 \ 2, 1)"
local bad4 "(0, 1 \ -1, 1)"
forvalues b = 1/4 {
    capture mata: T_TMP = dwlsmv_countdiag(`bad`b'', (2 \ 2))
    if _rc != 3498 {
        display as error "test_countdiag: bad Y #`b' gave rc " _rc ", not 3498"
        local nerr = `nerr' + 1
    }
}
capture mata: T_TMP = dwlsmv_countdiag((0, 1 \ 1, 0), (2 \ 2 \ 2))
if _rc != 3498 {
    display as error "test_countdiag: K of the wrong length gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}

if `nfail' | `nerr' {
    display as error _n "test_countdiag: `nfail' numeric failure(s), `nerr' error-path failure(s)"
    exit 9
}

display as result "test_countdiag: category counts reported, boundary and interior apart"
exit 0
