*! Module 10: the Gamma-hat diagnostic, dwlsmv_gammadiag(stage1, Y), and its
*! alias counter dwlsmv_dalias(R).
*!
*! Executed under StataMP 18 on 2026-09-26 and it passes. Before the routines
*! existed it failed with r(3499), dwlsmv_gammadiag() not found.
*!
*! Spec Section 5, Module 10; E-060 and E-071; derivation
*! @sec-gamma-rank, @sec-gamma-sym, @sec-gamma-budgets. Seven rows in fixed
*! order, [flag, value, i, j]: 1 asymmetry, 2 not PSD, 3 too few patterns,
*! 4 threshold aliasing, 5 centring, 6 raw ratio, 7 scale-free ratio. The flag
*! is 0/1, `.` (not evaluated) or .n (evaluated, no v1 flag: rows 6-7), and
*! each check asserts the state exactly.
*!
*!   1. CASE D, healthy: every flag 0, n_pat 139, no aliasing, z converged,
*!      r* = D, and both ratios equal to an independent svd() path.
*!   2. THE ALIASED DESIGN (test_structural.do's): row 4 flags d_alias = 1 at
*!      (1, 4); r* = D - 1 = 15 and the ratio is sigma_15, which is NOT a
*!      structural zero, while sigma_16 is.
*!   3. TOO FEW PATTERNS: twelve hand-picked patterns of case D's shape,
*!      N = 186. Row 3 flags n_pat = 12 < D + 1; r* = 11, and sigma_11 is
*!      not a structural zero while sigma_12 is (scratch-measured 0.033 against
*!      1.4e-16, 2026-09-26; asserted below by bounds, not those figures).
*!   4. dwlsmv_dalias() BY COMPONENTS, not by relations: three identical
*!      columns are three relations and a deficiency of two. Hand-built R, and
*!      thralias on a Y whose three binary columns coincide.
*!   5. ASYMMETRY at 2x and 0.5x the budget: at 2x row 1 flags at (1, 2) and
*!      rows 2-7 are all (., ., 0, 0) -- gammadiag stops. At 0.5x nothing stops.
*!   6. NOT PSD at 10x and 0.1x the budget, by shifting Gamma-hat by c I, at
*!      Gamma and at 100 Gamma (so the budget's lambda_max scale matters), with
*!      the ratio rows non-negative on the non-PSD matrix.
*!   7. CENTRING at 10x and 0.1x the 1e-6 budget, on psibar_12, and a
*!      negative component.
*!   8. THREE COLUMNS IN ONE ALIAS COMPONENT: three relations, d_alias 2,
*!      r* = D - 2, both ratios taken there, the first relation's location.
*!   9. n_pat AT D AND D + 1: row 3's edge and r*.
*!  10. CONTRACT, rc 3498 each: rows(Y) != N, cols(Y) != p, an out-of-range
*!      code, Gamma not D x D, psibar not D x 1, a missing psibar; and dalias
*!      on seven relation lists not in thralias's form.
*!
*! Checks 6-10's additions and 4's star came from the adversarial review of
*! 2026-09-26: each kills a mutation that survived the earlier suite.
*!
*! Reference CSVs: explicit delimiters(","), for the delimiter import trap.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off

local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
local load tests/dwlsmv_assert.mata src/mata/dwlsmv_types.mata
foreach f of local srcs {
    local load `load' src/mata/`f'
}
foreach f of local load {
    capture noisily do "`f'"
    if _rc {
        display as error "test_gammadiag: `f' failed to load"
        exit 9
    }
}

import delimited using "fixtures/data/case_D.csv", clear varnames(1) delimiters(",")
mata: YD = st_data(., .)

mata:

// test_structural.do's aliased design, verbatim.
real matrix dwlsmv_t_gd_aliased(real matrix YD)
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

// Twelve patterns, pattern r repeated 9 + r times: N = 186, n_pat = 12.
real matrix dwlsmv_t_gd_fewpat()
{
    real matrix P, Y
    real scalar r

    P = (0,1,0,2 \ 1,3,0,0 \ 2,0,1,1 \ 3,2,2,0 \ 0,2,2,1 \ 1,0,1,2 \
         2,1,0,0 \ 3,3,1,2 \ 1,2,2,1 \ 2,3,0,1 \ 0,0,1,0 \ 3,1,2,2)
    Y = J(0, 4, .)
    for (r = 1; r <= 12; r++) {
        Y = Y \ (J(9 + r, 1, 1) # P[r, .])
    }
    return(Y)
}

// Check 8's design: as the aliased one, but y3 joins: one row with
// y1 = y2 = y3 = 0 keeps its zeros; every other 0 in y1, y2, y3 becomes 1.
real matrix dwlsmv_t_gd_alias3(real matrix YD)
{
    real matrix Y
    real scalar i, i0, j

    Y  = YD
    i0 = 0
    for (i = 1; i <= rows(Y); i++) {
        if (Y[i, 1] == 0 & Y[i, 2] == 0 & Y[i, 3] == 0) {
            i0 = i
            break
        }
    }
    if (i0 == 0) _error("dwlsmv_t_gd_alias3: case D has no row with y1 = y2 = y3 = 0")
    for (i = 1; i <= rows(Y); i++) {
        for (j = 1; j <= 3; j++) {
            if (i != i0 & Y[i, j] == 0) Y[i, j] = 1
        }
    }
    return(Y)
}

// Check 9's designs: the first m of seventeen patterns of case D's shape,
// pattern r repeated 9 + r times. The first twelve are dwlsmv_t_gd_fewpat's.
real matrix dwlsmv_t_gd_fewpat_m(real scalar m)
{
    real matrix P, Y
    real scalar r

    P = (0,1,0,2 \ 1,3,0,0 \ 2,0,1,1 \ 3,2,2,0 \ 0,2,2,1 \ 1,0,1,2 \
         2,1,0,0 \ 3,3,1,2 \ 1,2,2,1 \ 2,3,0,1 \ 0,0,1,0 \ 3,1,2,2 \
         0,3,1,1 \ 1,1,2,0 \ 2,2,0,2 \ 3,0,0,1 \ 0,1,2,2)
    Y = J(0, 4, .)
    for (r = 1; r <= m; r++) {
        Y = Y \ (J(9 + r, 1, 1) # P[r, .])
    }
    return(Y)
}

// The independent path for rows 6-7: singular values from svd(), never from
// symeigenvalues(), which is what gammadiag uses.
real rowvector dwlsmv_t_gd_svratio(real matrix G, real scalar r)
{
    real matrix U, V, Gs
    real colvector sv, s

    U = J(0, 0, .); V = J(0, 0, .); sv = J(0, 1, .)
    svd(G, U, sv, V)
    s  = sqrt(diagonal(G))
    Gs = G :/ (s * s')
    return((sv[r] / sv[1], sv[r + 1] / sv[1], dwlsmv_t_gd_sv1(Gs, r)))
}

real scalar dwlsmv_t_gd_sv1(real matrix Gs, real scalar r)
{
    real matrix U, V
    real colvector sv

    U = J(0, 0, .); V = J(0, 0, .); sv = J(0, 1, .)
    svd(Gs, U, sv, V)
    return(sv[r] / sv[1])
}

// Does row `r` of A read "not evaluated" exactly: (., ., 0, 0)?
real scalar dwlsmv_t_gd_noteval(real matrix A, real scalar r)
{
    return(A[r, 1] == . & A[r, 2] == . & A[r, 3] == 0 & A[r, 4] == 0)
}

real scalar dwlsmv_t_gammadiag(real matrix YD)
{
    struct dwlsmv_stage1 scalar st, s2
    string colvector FAIL
    real matrix A, Ya, Yp, G, R
    real rowvector ref
    real colvector ev
    real scalar D, bud, i, f, c, lmin, lmax, tgt, sc

    FAIL = J(0, 1, "")

    // ---- 1. case D -------------------------------------------------------------
    st  = dwls_sample_stats(YD, (4 \ 4 \ 3 \ 3))
    D   = st.map.D
    bud = 4 * D * epsilon(1)
    A   = dwlsmv_gammadiag(st, YD)
    printf("{txt}  1. case D     : flags %f %f %f %f %f, z %g, ratios %g / %g\n",
           A[1, 1], A[2, 1], A[3, 1], A[4, 1], A[5, 1], A[5, 2], A[6, 2], A[7, 2])
    if (rows(A) != 7 | cols(A) != 4) {
        FAIL = FAIL \ "1: gammadiag is not 7 x 4"
        return(rows(FAIL))
    }
    if (A[1, .] != (0, 0, 0, 0)) FAIL = FAIL \ "1: case D's Gamma is not read as exactly symmetric"
    if (!(A[2, 1] == 0 & A[2, 2] > 1e-3)) FAIL = FAIL \ "1: case D's PSD row is not a clear pass"
    if (A[3, .] != (0, 139, D, 0)) FAIL = FAIL \ "1: case D's pattern row is not (0, 139, 16, 0)"
    if (A[4, .] != (0, 0, 0, 0)) FAIL = FAIL \ "1: aliasing reported on case D"
    if (!(A[5, 1] == 0 & A[5, 2] < 1e-12 & A[5, 3] >= 1 & A[5, 3] <= D)) {
        FAIL = FAIL \ "1: case D's centring row is not a converged pass with a location"
    }
    if (!(A[6, 1] == .n & A[7, 1] == .n)) FAIL = FAIL \ "1: a ratio row's flag is not .n"
    if (!(A[6, 3] == D & A[7, 3] == D)) FAIL = FAIL \ "1: r* is not D on case D"
    // svd at r = D - 1 returns (sigma_15, sigma_16, scale-free sigma_15); rebuild at D.
    ref = dwlsmv_t_gd_svratio(st.Gamma, D - 1)
    if (abs(A[6, 2] - ref[2]) > bud) FAIL = FAIL \ "1: raw ratio disagrees with svd() beyond 4 D eps"
    if (abs(A[7, 2] - dwlsmv_t_gd_sv1(st.Gamma :/ (sqrt(diagonal(st.Gamma)) * sqrt(diagonal(st.Gamma))'), D)) > bud) {
        FAIL = FAIL \ "1: scale-free ratio disagrees with svd() beyond 4 D eps"
    }

    // ---- 2. the aliased design -------------------------------------------------
    Ya = dwlsmv_t_gd_aliased(YD)
    s2 = dwls_sample_stats(Ya, (4 \ 4 \ 3 \ 3))
    A  = dwlsmv_gammadiag(s2, Ya)
    ref = dwlsmv_t_gd_svratio(s2.Gamma, D - 1)
    printf("{txt}  2. aliased    : row 4 (%f, %f, %f, %f), r* %f, ratio %g (svd %g; next %g)\n",
           A[4, 1], A[4, 2], A[4, 3], A[4, 4], A[6, 3], A[6, 2], ref[1], ref[2])
    if (A[4, .] != (1, 1, 1, 4)) FAIL = FAIL \ "2: row 4 is not (1, 1, 1, 4)"
    if (!(A[6, 3] == D - 1 & A[7, 3] == D - 1)) FAIL = FAIL \ "2: r* is not D - 1"
    if (abs(A[6, 2] - ref[1]) > bud) FAIL = FAIL \ "2: raw ratio is not sigma_15 / sigma_1"
    if (!(A[6, 2] > 1e3 * bud & ref[2] < bud)) {
        FAIL = FAIL \ "2: sigma_15 is not clear of zero, or sigma_16 is not a structural zero"
    }
    if (A[2, 1] != 0) FAIL = FAIL \ "2: an exactly singular PSD Gamma was flagged not PSD"
    if (abs(A[7, 2] - ref[3]) > bud) FAIL = FAIL \ "2: scale-free ratio is not taken at r* = 15"

    // ---- 3. too few patterns ---------------------------------------------------
    Yp = dwlsmv_t_gd_fewpat()
    s2 = dwls_sample_stats(Yp, (4 \ 4 \ 3 \ 3))
    A  = dwlsmv_gammadiag(s2, Yp)
    ref = dwlsmv_t_gd_svratio(s2.Gamma, 11)
    printf("{txt}  3. 12 patterns: row 3 (%f, %f, %f, %f), r* %f, ratio %g (next %g)\n",
           A[3, 1], A[3, 2], A[3, 3], A[3, 4], A[6, 3], A[6, 2], ref[2])
    if (A[3, .] != (1, 12, D, 0)) FAIL = FAIL \ "3: row 3 is not (1, 12, 16, 0)"
    if (A[4, 2] != 0) FAIL = FAIL \ "3: aliasing reported on the twelve-pattern design"
    if (!(A[6, 3] == 11 & A[7, 3] == 11)) FAIL = FAIL \ "3: r* is not n_pat - 1 = 11"
    if (abs(A[6, 2] - ref[1]) > bud) FAIL = FAIL \ "3: raw ratio is not sigma_11 / sigma_1"
    if (abs(A[7, 2] - ref[3]) > bud) FAIL = FAIL \ "3: scale-free ratio is not taken at r* = 11"
    if (!(A[6, 2] > 1e3 * bud & ref[2] < bud)) {
        FAIL = FAIL \ "3: sigma_11 is not clear of zero, or sigma_12 is not a structural zero"
    }

    // ---- 4. dwlsmv_dalias, by components -----------------------------------------
    if (dwlsmv_dalias(J(0, 3, .)) != 0)                  FAIL = FAIL \ "4: no relation is not 0"
    if (dwlsmv_dalias((1, 2, 1)) != 1)                   FAIL = FAIL \ "4: one relation is not 1"
    if (dwlsmv_dalias((1, 2, 1 \ 1, 3, 1 \ 2, 3, 1)) != 2) FAIL = FAIL \ "4: a triangle is not 2"
    if (dwlsmv_dalias((1, 2, 1 \ 3, 4, 2)) != 2)         FAIL = FAIL \ "4: two components are not 2"
    if (dwlsmv_dalias((1, 2, 1 \ 2, 3, 2 \ 3, 4, 1)) != 3) FAIL = FAIL \ "4: a chain of four is not 3"
    if (dwlsmv_dalias((5, 0, 3 \ 1, 2, 1)) != 2)         FAIL = FAIL \ "4: constant plus pair is not 2"
    // not a clique and not an ascending chain: 3 joins 1 and 2 (review M13)
    if (dwlsmv_dalias((1, 3, 1 \ 2, 3, 1)) != 2)         FAIL = FAIL \ "4: a star on 3 is not 2"
    R = dwlsmv_thralias((0, 0, 0, 1 \ 1, 1, 1, 0 \ 0, 0, 0, 0 \ 1, 1, 1, 1), (2 \ 2 \ 2 \ 2))
    printf("{txt}  4. three identical columns: %f relation(s), d_alias %f\n", rows(R), dwlsmv_dalias(R))
    if (!(rows(R) == 3 & dwlsmv_dalias(R) == 2)) {
        FAIL = FAIL \ "4: three identical columns are not 3 relations and d_alias 2"
    }

    // ---- 5. asymmetry: 2x stops, 0.5x does not -----------------------------------
    G = st.Gamma
    for (f = 2; f >= 0.5; f = f / 4) {
        s2 = st
        s2.Gamma[1, 2] = G[1, 2] + f * bud * sqrt(G[1, 1] * G[2, 2])
        A = dwlsmv_gammadiag(s2, YD)
        printf("{txt}  5. asymmetry %gx: row 1 (%f, %g, %f, %f)\n", f, A[1, 1], A[1, 2], A[1, 3], A[1, 4])
        if (f > 1) {
            if (!(A[1, 1] == 1 & A[1, 3] == 1 & A[1, 4] == 2)) FAIL = FAIL \ "5: 2x did not flag at (1, 2)"
            for (i = 2; i <= 7; i++) {
                if (!dwlsmv_t_gd_noteval(A, i)) FAIL = FAIL \ ("5: row " + strofreal(i) + " evaluated after an asymmetry flag")
            }
        }
        else {
            if (A[1, 1] != 0) FAIL = FAIL \ "5: 0.5x flagged"
            if (A[3, 1] != 0 | A[6, 1] != .n) FAIL = FAIL \ "5: 0.5x stopped the later rows"
        }
    }

    // ---- 6. not PSD: shift by c I so lambda_min = -f budget lambda_max -------------
    // Also at 100 Gamma, where lambda_max is far above 1, so a budget that
    // dropped its lambda_max scale would give the other verdict (review M6);
    // and the ratio rows stay non-negative on a non-PSD Gamma, since singular
    // values are |eigenvalues| (review M9).
    for (sc = 1; sc <= 100; sc = sc * 100) {
        ev   = symeigenvalues(sc * st.Gamma)'
        lmin = min(ev); lmax = max(ev)
        for (f = 10; f >= 0.1; f = f / 100) {
            tgt = -f * bud * lmax
            c   = lmin - tgt
            s2  = st
            s2.Gamma = sc * st.Gamma - c * I(D)
            A = dwlsmv_gammadiag(s2, YD)
            printf("{txt}  6. not PSD %gx, Gamma x %g: row 2 (%f, %g), ratios %g %g\n",
                   f, sc, A[2, 1], A[2, 2], A[6, 2], A[7, 2])
            if (A[2, 1] != (f > 1)) FAIL = FAIL \ ("6: PSD flag wrong at " + strofreal(f) + "x, Gamma x " + strofreal(sc))
            if (!(A[6, 2] >= 0 & A[7, 2] >= 0)) FAIL = FAIL \ "6: a ratio row went negative on a non-PSD Gamma"
        }
    }

    // ---- 7. centring on psibar_12 -------------------------------------------------
    for (f = 10; f >= 0.1; f = f / 100) {
        s2 = st
        s2.psibar = J(D, 1, 0)
        s2.psibar[12] = f * 1e-6 * sqrt(st.Gamma[12, 12])
        A = dwlsmv_gammadiag(s2, YD)
        printf("{txt}  7. centring %gx: row 5 (%f, %g, %f)\n", f, A[5, 1], A[5, 2], A[5, 3])
        if (!(A[5, 1] == (f > 1) & A[5, 3] == 12 & abs(A[5, 2] / (f * 1e-6) - 1) < 1e-12)) {
            FAIL = FAIL \ ("7: centring row wrong at " + strofreal(f) + "x the budget")
        }
    }
    // a NEGATIVE component flags as well (review M4)
    s2 = st
    s2.psibar = J(D, 1, 0)
    s2.psibar[12] = -10 * 1e-6 * sqrt(st.Gamma[12, 12])
    A = dwlsmv_gammadiag(s2, YD)
    if (!(A[5, 1] == 1 & A[5, 3] == 12)) FAIL = FAIL \ "7: a negative psibar component did not flag"

    // ---- 8. three columns in ONE alias component ------------------------------
    // y1, y2 and y3 keep one observation each in category 0, the same one:
    // three relations, one component of three columns, d_alias = 2, so
    // r* = D - 2 (review M1: counting relations would give D - 3; M5: the
    // location is the FIRST relation).
    Ya = dwlsmv_t_gd_alias3(YD)
    R  = dwlsmv_thralias(Ya, (4 \ 4 \ 3 \ 3))
    s2 = dwls_sample_stats(Ya, (4 \ 4 \ 3 \ 3))
    A  = dwlsmv_gammadiag(s2, Ya)
    ref = dwlsmv_t_gd_svratio(s2.Gamma, D - 2)
    printf("{txt}  8. three aliased: %f relations, row 4 (%f, %f, %f, %f), r* %f, ratio %g (next %g)\n",
           rows(R), A[4, 1], A[4, 2], A[4, 3], A[4, 4], A[6, 3], A[6, 2], ref[2])
    if (rows(R) != 3) FAIL = FAIL \ "8: the design does not have three relations"
    if (A[4, .] != (1, 2, 1, 4)) FAIL = FAIL \ "8: row 4 is not (1, 2, 1, 4)"
    if (!(A[6, 3] == D - 2 & A[7, 3] == D - 2)) FAIL = FAIL \ "8: r* is not D - 2"
    if (abs(A[6, 2] - ref[1]) > bud | abs(A[7, 2] - ref[3]) > bud) FAIL = FAIL \ "8: ratios are not taken at r* = D - 2"
    if (!(A[6, 2] > 1e3 * bud & ref[2] < bud)) FAIL = FAIL \ "8: sigma_14 is not clear of zero, or sigma_15 is not structural"

    // ---- 9. n_pat at D and at D + 1: the edge of row 3 (review M3) -------------
    for (i = 16; i <= 17; i++) {
        Yp = dwlsmv_t_gd_fewpat_m(i)
        s2 = dwls_sample_stats(Yp, (4 \ 4 \ 3 \ 3))
        A  = dwlsmv_gammadiag(s2, Yp)
        printf("{txt}  9. %f patterns: row 3 (%f, %f), r* %f\n", i, A[3, 1], A[3, 2], A[6, 3])
        if (A[3, 2] != i) FAIL = FAIL \ ("9: n_pat is not " + strofreal(i))
        if (A[3, 1] != (i - 1 < D)) FAIL = FAIL \ ("9: row 3's flag is wrong at n_pat = " + strofreal(i))
        if (A[6, 3] != min((D, i - 1))) FAIL = FAIL \ ("9: r* is wrong at n_pat = " + strofreal(i))
    }

    printf("{txt}test_gammadiag: %f failure(s)\n", rows(FAIL))
    for (i = 1; i <= rows(FAIL); i++) {
        printf("{err}    %s\n", FAIL[i])
    }
    return(rows(FAIL))
}

// ---- 8's contract violations, one per `which` ---------------------------------------
void dwlsmv_t_gd_bad(real matrix YD, real scalar which)
{
    struct dwlsmv_stage1 scalar st
    real matrix Y

    st = dwls_sample_stats(YD, (4 \ 4 \ 3 \ 3))
    Y  = YD
    if (which == 1) Y = YD[2..rows(YD), .]
    if (which == 2) Y = YD[., 1..3]
    if (which == 3) Y[1, 1] = 4
    if (which == 4) st.Gamma = st.Gamma[1..15, 1..15]
    if (which == 5) st.psibar = st.psibar[1..15]
    if (which == 6) st.psibar[3] = .
    (void) dwlsmv_gammadiag(st, Y)
}

end

mata: st_local("nfail", strofreal(dwlsmv_t_gammadiag(YD)))

* ---------------------------------------------------------------------------
* 10. CONTRACT -- rc 3498 each: gammadiag (#6, a missing psibar, added after
*     the 2026-09-26 review, D1), then dalias on relation lists that are not
*     thralias's form (review D2/M10).
* ---------------------------------------------------------------------------
local nerr 0
forvalues b = 1/6 {
    capture mata: dwlsmv_t_gd_bad(YD, `b')
    if _rc != 3498 {
        display as error "test_gammadiag: contract violation #`b' gave rc " _rc ", not 3498"
        local nerr = `nerr' + 1
    }
}
local dbad1 "(1, 0, 1)"
local dbad2 "(1.5, 2, 1)"
local dbad3 "(1, 2, 3)"
local dbad4 "(2, 1, 1)"
local dbad5 "(1, 2, 4)"
local dbad6 "J(0, 2, 0)"
local dbad7 "(1, ., 1)"
forvalues b = 1/7 {
    capture mata: T_TMP = dwlsmv_dalias(`dbad`b'')
    if _rc != 3498 {
        display as error "test_gammadiag: dalias on bad R #`b' gave rc " _rc ", not 3498"
        local nerr = `nerr' + 1
    }
}

if `nfail' | `nerr' {
    display as error _n "test_gammadiag: `nfail' numeric failure(s), `nerr' error-path failure(s)"
    exit 9
}

display as result "test_gammadiag: Gamma-hat's seven rows, each against its cause or budget"
exit 0
