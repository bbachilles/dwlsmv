*! Module 10: the two STRUCTURAL causes of a singular Gamma-hat, each detected
*! from its cause with no cutoff -- dwlsmv_npatterns() and dwlsmv_thralias().
*!
*! Executed under StataMP 18 on 2026-09-25 and it passes.
*!
*! Spec Section 4.5.6 (the table of exact causes) and Section 5, Module 10;
*! derivation @sec-gamma-rank, Theorems 1 and 3.
*!
*!   1. CASE D. n_pat = 139, the derivation's count, and no aliasing.
*!   2. THE ALIASED DESIGN of @sec-gamma-rank Level 3, rebuilt here from
*!      case_D.csv exactly as the derivation builds it. y1 and y2 each keep a
*!      single observation in category 0, the same one. Expect n_pat = 82 and
*!      exactly one relation: y1|t1 and y2|t1 identical. Then the consequence:
*!      production's Gamma-hat on that data is numerically singular (smallest
*!      eigenvalue ratio inside the PSD budget of 4 D eps), while case D's is
*!      not. The detector names a cause that the matrix confirms.
*!   3. COMPLEMENTARY columns: a binary item and its mirror, 1{y_b < 1} =
*!      1 - 1{y_a < 1}. Expect one relation of type 2.
*!   4. TOO FEW PATTERNS, counted exactly: a hand-built Y with duplicated rows
*!      whose distinct count is known.
*!   5. ERROR PATHS, rc 3498 each, for both routines: a missing value, a
*!      non-integer code, a code >= K, a negative code, and K of the wrong
*!      length. A `.` counted as a code would inflate n_pat and let the
*!      structural check pass wrongly (spec 4.5.6).
*!
*! Reference CSVs: explicit delimiters(","), for the delimiter import trap.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off

foreach f in tests/dwlsmv_assert.mata               ///
             src/mata/dwlsmv_types.mata             ///
             src/mata/dwlsmv_statmap_build.mata     ///
             src/mata/dwlsmv_pos_tau.mata           ///
             src/mata/dwlsmv_pos_rho.mata           ///
             src/mata/dwlsmv_catcounts.mata         ///
             src/mata/dwlsmv_catdegen.mata          ///
             src/mata/dwlsmv_cumprop.mata           ///
             src/mata/dwlsmv_thresholds.mata        ///
             src/mata/dwlsmv_psi_tau_table.mata     ///
             src/mata/dwlsmv_cellprob.mata          ///
             src/mata/dwlsmv_dcellprob_rho.mata     ///
             src/mata/dwlsmv_wmat.mata              ///
             src/mata/dwlsmv_crosstab.mata          ///
             src/mata/dwlsmv_pairinfo.mata          ///
             src/mata/dwls_polychoric.mata          ///
             src/mata/dwlsmv_koff.mata              ///
             src/mata/dwlsmv_psi_tau_stack.mata     ///
             src/mata/dwlsmv_psi_rho_table.mata     ///
             src/mata/dwlsmv_psi_matrix.mata        ///
             src/mata/dwlsmv_psibar.mata            ///
             src/mata/dwls_gamma.mata               ///
             src/mata/dwlsmv_emptycell1.mata      ///
             src/mata/dwlsmv_emptycell_repair.mata ///
             src/mata/dwlsmv_emptycell_pairs.mata ///
             src/mata/dwls_sample_stats.mata        ///
             src/mata/dwlsmv_npatterns.mata         ///
             src/mata/dwlsmv_thralias.mata {
    capture noisily do "`f'"
    if _rc {
        display as error "test_structural: `f' failed to load"
        exit 9
    }
}

import delimited using "fixtures/data/case_D.csv", clear varnames(1) delimiters(",")
mata: YD = st_data(., .)

mata:

// @sec-gamma-rank Level 3's construction, verbatim in effect: i0 is the first
// row with y1 = y2 = 0; every other 0 in y1, and in y2, becomes 1.
real matrix dwlsmv_t_aliased(real matrix YD)
{
    real matrix Y
    real scalar i, i0

    Y  = YD
    i0 = 0
    for (i = 1; i <= rows(Y); i++) {
        if (Y[i, 1] == 0 & Y[i, 2] == 0) {
            i0 = i
            break
        }
    }
    for (i = 1; i <= rows(Y); i++) {
        if (i != i0 & Y[i, 1] == 0) Y[i, 1] = 1
        if (i != i0 & Y[i, 2] == 0) Y[i, 2] = 1
    }
    return(Y)
}

real scalar dwlsmv_t_structural(real matrix YD)
{
    struct dwlsmv_stage1 scalar st, sa
    string colvector FAIL
    real matrix Ya, R, Yc, Yp
    real colvector K, ev
    real scalar i, n, ratD, ratA, budA

    FAIL = J(0, 1, "")
    K    = (4 \ 4 \ 3 \ 3)

    // ---- 1. case D -----------------------------------------------------------
    n = dwlsmv_npatterns(YD, K)
    R = dwlsmv_thralias(YD, K)
    printf("{txt}  1. case D         : n_pat %f, aliasing relations %f\n", n, rows(R))
    if (n != 139) {
        FAIL = FAIL \ ("1: case D n_pat is " + strofreal(n) + ", the derivation's is 139")
    }
    if (rows(R) != 0) {
        FAIL = FAIL \ "1: aliasing reported on case D"
    }

    // ---- 2. the aliased design -------------------------------------------------
    Ya = dwlsmv_t_aliased(YD)
    n  = dwlsmv_npatterns(Ya, K)
    R  = dwlsmv_thralias(Ya, K)
    printf("{txt}  2. aliased design : n_pat %f, aliasing relations %f\n", n, rows(R))
    for (i = 1; i <= rows(R); i++) {
        printf("{txt}                      statistics %f and %f, type %f\n", R[i, 1], R[i, 2], R[i, 3])
    }
    if (n != 82) {
        FAIL = FAIL \ ("2: aliased n_pat is " + strofreal(n) + ", the derivation's is 82")
    }
    // y1|t1 is statistic 1 and y2|t1 is statistic 4 (y1 has 3 thresholds).
    if (rows(R) != 1) {
        FAIL = FAIL \ ("2: expected exactly one relation, found " + strofreal(rows(R)))
    }
    else if (R[1, .] != (1, 4, 1)) {
        FAIL = FAIL \ "2: the relation is not y1|t1 identical to y2|t1"
    }
    // The consequence the detector predicts, read off production's own Gamma.
    st   = dwls_sample_stats(YD, K)
    sa   = dwls_sample_stats(Ya, K)
    ev   = symeigenvalues(st.Gamma)'
    ratD = min(ev) / max(ev)
    ev   = symeigenvalues(sa.Gamma)'
    ratA = min(ev) / max(ev)
    budA = 4 * rows(sa.Gamma) * epsilon(1)
    printf("{txt}                      Gamma eigenvalue ratio: aliased %g (budget %g), case D %g\n",
           ratA, budA, ratD)
    if (!(abs(ratA) <= budA)) {
        FAIL = FAIL \ "2: the aliased Gamma is not numerically singular, against Theorem 3"
    }
    if (!(ratD > 1e-3)) {
        FAIL = FAIL \ "2: case D's Gamma is not clearly nonsingular; the contrast is void"
    }

    // ---- 3. complementary columns ------------------------------------------------
    Yc = (0, 1, 0 \ 1, 0, 1 \ 0, 1, 1 \ 1, 0, 0 \ 0, 1, 0 \ 1, 0, 1)
    R  = dwlsmv_thralias(Yc, (2 \ 2 \ 2))
    printf("{txt}  3. complementary  : relations %f\n", rows(R))
    if (rows(R) != 1) {
        FAIL = FAIL \ ("3: expected one relation, found " + strofreal(rows(R)))
    }
    else if (R[1, .] != (1, 2, 2)) {
        FAIL = FAIL \ "3: the relation is not statistics 1 and 2, complementary"
    }

    // ---- 4. too few patterns -------------------------------------------------------
    Yp = (0, 1 \ 1, 0 \ 0, 1 \ 1, 1 \ 0, 0 \ 1, 1 \ 0, 1)
    n  = dwlsmv_npatterns(Yp, (2 \ 2))
    printf("{txt}  4. hand-built Y   : n_pat %f (4 distinct rows of 7)\n", n)
    if (n != 4) {
        FAIL = FAIL \ ("4: n_pat is " + strofreal(n) + ", not 4")
    }

    // ---- 4b. a CONSTANT column (added after the 2026-09-26 review, M14/M15) -----
    // y1 declared K = 3 but never reaches category 2, so 1{y1 < 2} is 1 on
    // every row: constant, type 3, reported as (2, 0, 3) and NOT also paired.
    // 1{y1 < 1} = (1, 0, 1) and 1{y2 < 1} = (0, 1, 0) are complementary.
    R = dwlsmv_thralias((0, 1 \ 1, 0 \ 0, 1), (3 \ 2))
    R = sort(R, (1, 2))
    printf("{txt}  4b. constant column: %f relation(s)\n", rows(R))
    if (rows(R) != 2) {
        FAIL = FAIL \ ("4b: expected two relations, found " + strofreal(rows(R)))
    }
    else if (R != (1, 3, 2 \ 2, 0, 3)) {
        FAIL = FAIL \ "4b: not (1, 3, complementary) and (2, 0, constant)"
    }

    printf("{txt}test_structural: %f failure(s)\n", rows(FAIL))
    for (i = 1; i <= rows(FAIL); i++) {
        printf("{err}    %s\n", FAIL[i])
    }
    return(rows(FAIL))
}

end

mata: st_local("nfail", strofreal(dwlsmv_t_structural(YD)))

* ---------------------------------------------------------------------------
* 5. ERROR PATHS -- rc 3498 each, for both routines.
* ---------------------------------------------------------------------------
local nerr 0
local bad1 "(0, 1 \ ., 1)"
local bad2 "(0, 1 \ 0.5, 1)"
local bad3 "(0, 1 \ 2, 1)"
local bad4 "(0, 1 \ -1, 1)"
forvalues b = 1/4 {
    foreach fn in dwlsmv_npatterns dwlsmv_thralias {
        capture mata: T_TMP = `fn'(`bad`b'', (2 \ 2))
        if _rc != 3498 {
            display as error "test_structural: `fn' on bad Y #`b' gave rc " _rc ", not 3498"
            local nerr = `nerr' + 1
        }
    }
}
foreach fn in dwlsmv_npatterns dwlsmv_thralias {
    capture mata: T_TMP = `fn'((0, 1 \ 1, 0), (2 \ 2 \ 2))
    if _rc != 3498 {
        display as error "test_structural: `fn' with K of the wrong length gave rc " _rc ", not 3498"
        local nerr = `nerr' + 1
    }
}

if `nfail' | `nerr' {
    display as error _n "test_structural: `nfail' numeric failure(s), `nerr' error-path failure(s)"
    exit 9
}

display as result "test_structural: n_pat and threshold aliasing detected from their causes"
exit 0
