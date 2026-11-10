*! DIFFTEST statistic routine dwlsmv_difftest() against lavaan's own inputs
*! (fixtures/generate/97_lavaan_difftest.R; lavaan 0.7.2, Case E, marker
*! metric): beta0 (d = 1, H0 false) and xload2 (d = 2, H0 true).
*!
*!   1. Fed lavaan's Delta1, Delta0, Gamma and Td, it reproduces the
*!      fixture's rebuilt Tr(M1), Tr(M1^2), a_d, b_d, T_d3 and p, and lavaan's
*!      own lavTestLRT statistic (written by lavaan at about 10 significant
*!      digits) at the generator's 1e-9 relative.
*!   2. h-hat: the routine's pinv(Delta1) Delta0 equals the fixture's.
*!   3. Negative control: the selection-form rows (A.method = "exact") must
*!      NOT be reproduced, so the test can tell the two constructions apart.
*!   4. Invariance: permuting Delta1's and Delta0's columns leaves the result
*!      unchanged (h-hat needs no parameter correspondence).
*!   5. Scale: Gamma -> 4 Gamma leaves every output unchanged: W_D = diag(Gamma)^-1
*!      rescales with it, so M1 does not see Gamma's scale and N cancels.
*!   6. Error paths, each rc 3498.
*!
*! Executed under StataMP 18 on 2026-10-08 and it passes. Checks 1, 2, 4 and 5
*! assert at 1e-12 (relative; b_d absolute), chosen after measuring <= 6.4e-15.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

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
        display as error "test_difftest: `f' failed to load"
        exit 9
    }
}

* reference fixtures: stringcols(_all) + strtoreal(), delimiters(",") always
import delimited using "fixtures/expected/lavaan_case_E_difftest_matrices.csv", ///
    clear varnames(1) case(lower) delimiters(",") stringcols(_all)
mata: X_KEY = st_sdata(., ("pair", "role", "matrix"))
mata: X_IJ  = strtoreal(st_sdata(., ("row", "col")))
mata: X_V   = strtoreal(st_sdata(., "value"))

import delimited using "fixtures/expected/lavaan_case_E_difftest.csv", ///
    clear varnames(1) case(lower) delimiters(",") stringcols(_all)
mata: R_KEY = st_sdata(., ("pair", "variant"))
mata: R_V   = strtoreal(st_sdata(., ("d", "td", "tr_m1", "tr_m1sq", "a_d", "b_d", "stat", "p", "lavaan_stat")))

mata:

real matrix tdt_mat(string matrix KEY, real matrix IJ, real colvector V,
                    string scalar pair, string scalar role, string scalar name)
{
    real colvector s
    real matrix A
    real scalar k

    s = selectindex(KEY[., 1] :== pair :& KEY[., 2] :== role :& KEY[., 3] :== name)
    if (!rows(s)) {
        _error("test_difftest: fixture matrix not found")
    }
    A = J(max(IJ[s, 1]), max(IJ[s, 2]), .)
    for (k = 1; k <= rows(s); k++) {
        A[IJ[s[k], 1], IJ[s[k], 2]] = V[s[k]]
    }
    return(A)
}

real scalar tdt_row(string matrix KEY, string scalar pair, string scalar variant)
{
    real colvector s

    s = selectindex(KEY[., 1] :== pair :& KEY[., 2] :== variant)
    if (rows(s) != 1) {
        _error("test_difftest: fixture row not found")
    }
    return(s)
}

real scalar tdt_rc(real matrix D1, real matrix D0, real matrix G, real scalar Td)
{
    real rowvector r

    r = dwlsmv_difftest(D1, D0, G, Td)
    return(0)
}

real scalar tdt_main(string matrix X_KEY, real matrix X_IJ, real colvector X_V,
                     string matrix R_KEY, real matrix R_V)
{
    string colvector FAIL
    string rowvector pairs
    real matrix G, D1, D0, H
    real rowvector r, rr, rp, r4, i1, i0

    real scalar k, j, gap

    FAIL  = J(0, 1, "")
    pairs = ("beta0", "xload2")
    G     = tdt_mat(X_KEY, X_IJ, X_V, "all", "none", "Gamma")
    if (missing(G)) FAIL = FAIL \ "fixture Gamma has holes"

    for (k = 1; k <= 2; k++) {
        D1 = tdt_mat(X_KEY, X_IJ, X_V, pairs[k], "h1", "Delta")
        D0 = tdt_mat(X_KEY, X_IJ, X_V, pairs[k], "h0", "Delta")
        H  = tdt_mat(X_KEY, X_IJ, X_V, pairs[k], "pair", "hhat")
        if (missing(D1) | missing(D0) | missing(H)) FAIL = FAIL \ (pairs[k] + ": fixture matrix has holes")

        j = tdt_row(R_KEY, pairs[k], "reference")
        r = dwlsmv_difftest(D1, D0, G, R_V[j, 2])

        // 1. the routine against the fixture's rebuild and lavaan
        printf("{txt}  1. %s: T_d3 %21.15g fixture %21.15g lavaan %21.15g; d %g; p %g vs %g\n",
               pairs[k], r[1], R_V[j, 7], R_V[j, 9], r[2], r[3], R_V[j, 8])
        printf("{txt}     trM %21.15g vs %21.15g; trM2 %21.15g vs %21.15g; a %21.15g vs %21.15g; b %g vs %g\n",
               r[6], R_V[j, 3], r[7], R_V[j, 4], r[4], R_V[j, 5], r[5], R_V[j, 6])
        rr = reldif(r[(1, 3, 4, 6, 7)], R_V[j, (7, 8, 5, 3, 4)])
        printf("{txt}     reldif T_d3 p a trM trM2: %9.2e %9.2e %9.2e %9.2e %9.2e; |b - b_fix| %9.2e\n",
               rr[1], rr[2], rr[3], rr[4], rr[5], abs(r[5] - R_V[j, 6]))
        if (r[2] != R_V[j, 1]) FAIL = FAIL \ (pairs[k] + " 1: d differs")
        if (!(max(rr) <= 1e-12 & abs(r[5] - R_V[j, 6]) <= 1e-12)) FAIL = FAIL \ (pairs[k] + " 1: rebuild differs beyond 1e-12")
        gap = abs(r[1] / R_V[j, 9] - 1)
        if (!(gap < 1e-9)) FAIL = FAIL \ (pairs[k] + " 1: T_d3 vs lavaan " + strofreal(gap))

        // 2. h-hat
        printf("{txt}  2. %s: max |h - h_fixture| %9.2e\n", pairs[k], max(abs(quadcross(pinv(D1)', D0) - H)))
        if (!(max(abs(quadcross(pinv(D1)', D0) - H)) <= 1e-12)) FAIL = FAIL \ (pairs[k] + " 2: h-hat differs beyond 1e-12")

        // 3. negative control: the selection form must not be reproduced
        j  = tdt_row(R_KEY, pairs[k], "negative_contrast")
        gap = abs(r[1] / R_V[j, 7] - 1)
        printf("{txt}  3. %s: routine %21.15g vs selection form %21.15g, rel %9.2e\n",
               pairs[k], r[1], R_V[j, 7], gap)
        if (!(gap > 1e-6)) FAIL = FAIL \ (pairs[k] + " 3: selection form reproduced")

        // 4. column permutations
        i1 = cols(D1)..1
        i0 = cols(D0)..1
        rp = dwlsmv_difftest(D1[., i1], D0[., i0], G, R_V[tdt_row(R_KEY, pairs[k], "reference"), 2])
        printf("{txt}  4. %s: max reldif under column reversal %9.2e\n", pairs[k], max(reldif(rp, r)))
        if (!(max(reldif(rp, r)) <= 1e-12)) FAIL = FAIL \ (pairs[k] + " 4: not invariant to column order")

        // 5. Gamma's scale cancels
        r4 = dwlsmv_difftest(D1, D0, 4 :* G, R_V[tdt_row(R_KEY, pairs[k], "reference"), 2])
        printf("{txt}  5. %s: a(4G)/a(G) %21.15g; b(4G) - b(G) %9.2e\n", pairs[k], r4[4] / r[4], r4[5] - r[5])
        if (!(max(reldif(r4, r)) <= 1e-12)) FAIL = FAIL \ (pairs[k] + " 5: depends on the scale of Gamma")
    }

    for (k = 1; k <= rows(FAIL); k++) {
        printf("{err}%s\n", FAIL[k])
    }
    return(rows(FAIL))
}

end

mata: st_local("nfail", strofreal(tdt_main(X_KEY, X_IJ, X_V, R_KEY, R_V)))

* 6. error paths on xload2's inputs, each must be rc 3498
mata: DT_G  = tdt_mat(X_KEY, X_IJ, X_V, "all", "none", "Gamma")
mata: DT_D1 = tdt_mat(X_KEY, X_IJ, X_V, "xload2", "h1", "Delta")
mata: DT_D0 = tdt_mat(X_KEY, X_IJ, X_V, "xload2", "h0", "Delta")
mata: DT_GZ = DT_G
mata: DT_GZ[1, 1] = 0
mata: DT_GA = DT_G
mata: DT_GA[1, 2] = DT_GA[1, 2] + 1e-3 * sqrt(DT_G[1, 1] * DT_G[2, 2])
mata: DT_DM = DT_D0
mata: DT_DM[1, 1] = .
mata: DT_DZ = DT_D0
mata: DT_DZ[., 1] = J(rows(DT_D0), 1, 0)

* indexed locals, one case each: word-of extraction splits a quoted case on its spaces
mata: DT_D1S = DT_D1[|1, 1 \ rows(DT_D1) - 1, .|]
mata: DT_D0S = DT_D0[|1, 1 \ rows(DT_D0) - 1, .|]
mata: DT_D00 = J(rows(DT_D1), 0, 0)
local c1  tdt_rc(DT_D1S, DT_D0, DT_G, 1)
local l1  D1 rows
local c2  tdt_rc(DT_D1, DT_D0S, DT_G, 1)
local l2  D0 rows
local c3  tdt_rc(DT_D1, DT_D1, DT_G, 1)
local l3  d = 0
local c4  tdt_rc(DT_D0, DT_D1, DT_G, 1)
local l4  d < 0
local c5  tdt_rc(DT_D1, DT_D00, DT_G, 1)
local l5  q0 = 0
local c6  tdt_rc(DT_D1, DT_D0, DT_G, .)
local l6  Td missing
local c7  tdt_rc(DT_D1, DT_DM, DT_G, 1)
local l7  D0 missing
local c8  tdt_rc(DT_D1, DT_D0, DT_GZ, 1)
local l8  diag Gamma 0
local c9  tdt_rc(DT_D1, DT_D0, DT_GA, 1)
local l9  Gamma asymmetric
local c10 tdt_rc(DT_D1, DT_DZ, DT_G, 1)
local l10 hP1h singular

local ncase 10
forvalues i = 1/`ncase' {
    local c `c`i''
    local l `l`i''
    capture mata: `c'
    local rc = _rc
    display as text "  6. `l': rc `rc'"
    if `rc' != 3498 {
        capture noisily mata: `c'
        display as error "test_difftest 6: `l' gave rc `rc', not 3498"
        local ++nfail
    }
}

if `nfail' {
    display as error _n "test_difftest: `nfail' failure(s)"
    exit 9
}
display as result "test_difftest: dwlsmv_difftest() reproduces lavaan's satorra.2000 (delta) statistic from (Delta1, Delta0, Gamma, Td); error paths rc 3498"
exit 0
