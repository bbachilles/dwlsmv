*! The one-empty-cell pair (O6, decided 2026-10-06; spec Module 2,
*! *One-empty-cell repair*): membership, the repair, and agreement of the three
*! consumers that derive membership from the same Y.
*!
*!   1. Membership, through Y and the statmap (dwlsmv_emptycell_pairs()):
*!        a. 2x2, exactly one empty cell, off and on the diagonal -> 1
*!        b. 2x2, two empty cells on either diagonal -> 0
*!        c. a 3x3 table and a 2x3 table with one empty cell -> 0
*!        d. 2x2, two empty cells in one row or one column -> 0 from the
*!           helper; Stage 1 still aborts it upstream as an empty category
*!           (rc 3498)
*!   2. The repair: margins kept, +-0.5 by sigma_e, every cell positive; a
*!      table that is not one-empty-cell is refused (rc 3498).
*!   3. Stage 1 on L9d + y3: rho-hat+ equals an independent oracle, the
*!      saturated inversion binormal(a, b, rho) = (n11 + 0.5 sigma_e)/N by
*!      bisection (@sec-o6); bound = 0; thresholds equal dwlsmv_thresholds().
*!   4. Construction (a): Gamma_rr and psi-bar_rho of the repaired pair equal
*!      the .qmd's values for L9d under (a) (chunk o6-vii; see the constants'
*!      note), and differ from (b)'s Gamma_rr.
*!   5. Integration: Stage 1, Module 10 and the helper agree on the same Y.
*!      Module 10's row 5 does not flag on the repaired data, although the
*!      repaired pair's own z is far above 1e-6, so the exclusion is what
*!      keeps it quiet; the helper marks exactly the pair Stage 1 repaired.
*!   6. The edge: (1000, 1 \ 0, 1000) repairs to rho-hat+ > B and aborts
*!      (O3), rc 3498.
*!
*! Executed under StataMP 18 on 2026-10-06 and it passes. Mutants, each run
*! against it (or the named test) and each fails it:
*!   M1 "exactly one" -> "any" empty cell: 1b and 1d fail.
*!   M2 the K = 2 check on the second variable dropped: 1c aborts
*!      (dwlsmv_crosstab refuses a 3-category column read as 2x2).
*!   M3 Stage 1's information at W = n+/N (the hybrid): 4 fails
*!      (Gamma_rr 0.0267, psi-bar_rho 0.0335).
*!   M4 Module 10's exclusion removed: 5 fails (row 5 flags, z 0.584).
*!   M5 Stage 1 never repairs: 3 aborts (rc 3498); test_sample_stats fails.
*!   M6 the repair's sign flipped: 3 aborts (rc 3498).
*!   M7 Module 11's note dropped: test_dwlsmv 10 fails.
*!   M8 a repaired pair accepted whatever its status: 6 fails (rc 0, no
*!      message).
*!   M9 Stage 1's repaired-pair printout dropped: test_sample_stats 3 fails.
*! Kills by abort (M2, M5, M6) are counted by the runner as failures (any
*! rc other than 0, 8, 9).
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
local load src/mata/dwlsmv_types.mata
foreach f of local srcs {
    local load `load' src/mata/`f'
}
foreach f of local load {
    capture noisily do "`f'"
    if _rc {
        display as error "test_emptycell: `f' failed to load"
        exit 9
    }
}

mata:

// A two-variable table expanded to rows, plus y3 = mod(i, k3) when k3 > 0.
real matrix dwlsmv_t_ec_data(real matrix n, real scalar k3)
{
    real matrix Y
    real scalar r, c

    Y = J(0, 2, .)
    for (r = 1; r <= rows(n); r++) {
        for (c = 1; c <= cols(n); c++) {
            if (n[r, c]) {
                Y = Y \ J(n[r, c], 1, (r - 1, c - 1))
            }
        }
    }
    if (k3 > 0) {
        Y = (Y, mod((1::rows(Y)), k3))
    }
    return(Y)
}

// Membership of pair 1 (y1~~y2) for table n.
real scalar dwlsmv_t_ec_mem(real matrix n)
{
    struct dwlsmv_statmap scalar map
    map = dwlsmv_statmap_build((rows(n) \ cols(n)))
    return(dwlsmv_emptycell_pairs(dwlsmv_t_ec_data(n, 0), map)[1])
}

// The oracle: rho with binormal(a, b, rho) = p, by bisection on (-1, 1).
// Independent of the root finder and of the score.
real scalar dwlsmv_t_ec_inv(real scalar a, real scalar b, real scalar p)
{
    real scalar lo, hi, mid, i

    lo = -1 + 1e-12
    hi =  1 - 1e-12
    for (i = 1; i <= 200; i++) {
        mid = (lo + hi) / 2
        if (binormal(a, b, mid) < p) {
            lo = mid
        }
        else {
            hi = mid
        }
    }
    return((lo + hi) / 2)
}

real scalar dwlsmv_t_ec_main()
{
    struct dwlsmv_stage1 scalar st
    string colvector FAIL
    real matrix Y, n, np, A
    real colvector rep, zk
    real scalar a, b, want, i, zr, k

    FAIL = J(0, 1, "")

    // ---- 1: membership -------------------------------------------------
    if (dwlsmv_t_ec_mem((30, 2 \ 0, 8)) != 1)  FAIL = FAIL \ "1a: off-diagonal empty cell not a member"
    if (dwlsmv_t_ec_mem((0, 9 \ 7, 24)) != 1)  FAIL = FAIL \ "1a: diagonal empty cell not a member"
    if (dwlsmv_t_ec_mem((20, 0 \ 0, 20)) != 0) FAIL = FAIL \ "1b: off-diagonal pair empty is a member"
    if (dwlsmv_t_ec_mem((0, 20 \ 20, 0)) != 0) FAIL = FAIL \ "1b: diagonal pair empty is a member"
    if (dwlsmv_t_ec_mem((10, 3, 0 \ 2, 10, 3 \ 1, 2, 10)) != 0) FAIL = FAIL \ "1c: 3x3 with an empty cell is a member"
    if (dwlsmv_t_ec_mem((10, 3, 0 \ 2, 10, 3)) != 0) FAIL = FAIL \ "1c: 2x3 with an empty cell is a member"
    if (dwlsmv_emptycell1((0, 0 \ 7, 9)) != 0) FAIL = FAIL \ "1d: empty row is a member"
    if (dwlsmv_emptycell1((0, 7 \ 0, 9)) != 0) FAIL = FAIL \ "1d: empty column is a member"

    // ---- 2: the repair --------------------------------------------------
    n  = (30, 2 \ 0, 8)
    np = dwlsmv_emptycell_repair(n)
    if (np != (29.5, 2.5 \ 0.5, 7.5)) FAIL = FAIL \ "2: L9d (sigma_e = -1) repaired wrongly"
    if (rowsum(np) != rowsum(n) | colsum(np) != colsum(n)) FAIL = FAIL \ "2: margins moved"
    np = dwlsmv_emptycell_repair((0, 9 \ 7, 24))
    if (np != (0.5, 8.5 \ 6.5, 24.5)) FAIL = FAIL \ "2: T2 (sigma_e = +1) repaired wrongly"

    // ---- 3 and 4: Stage 1 on L9d + y3 -------------------------------------
    Y  = dwlsmv_t_ec_data((30, 2 \ 0, 8), 2)
    st = dwls_sample_stats(Y, (2 \ 2 \ 2))
    a  = invnormal(32 / 40)
    b  = invnormal(30 / 40)
    want = dwlsmv_t_ec_inv(a, b, 29.5 / 40)
    // TOLERANCE 1e-12, CHOSEN AFTER MEASURING: on seven tables (L9d, T1, T2,
    // S20, L200, D200, C1 of E-077) Stage 1's root
    // and this bisection agreed to <= 2.6e-15 (main session, StataMP 18,
    // 2026-10-06). The root finder's stop is |S-bar| <= tau_S, not a rho
    // tolerance, so this is a measured margin, not a derived one.
    if (abs(st.s[st.map.T + 1] - want) > 1e-12) {
        FAIL = FAIL \ ("3: rho-hat+ " + strofreal(st.s[st.map.T + 1], "%20.15f") +
                       " is not the inversion " + strofreal(want, "%20.15f"))
    }
    if (st.bound != J(3, 1, 0)) FAIL = FAIL \ "3: a repaired pair is marked in bound"
    if (st.s[1] != dwlsmv_thresholds(Y[., 1], 2) | st.s[2] != dwlsmv_thresholds(Y[., 2], 2)) {
        FAIL = FAIL \ "3: thresholds moved"
    }
    // Construction (a) on L9d, from the .qmd (o6-core/o6-vii, R 4.6.0, main
    // session 2026-10-06): Gamma_rr 0.0880596 (a) vs 0.0735137 (b);
    // psi-bar_rho 0.1732898. Six digits as printed; the relative tolerance
    // 1e-5 sits above the printing and far below (a) vs (b)'s 18%.
    k = st.map.T + 1
    if (abs(st.Gamma[k, k] / 0.0880596 - 1) > 1e-5) {
        FAIL = FAIL \ ("4: Gamma_rr " + strofreal(st.Gamma[k, k], "%12.8f") + " is not (a)'s 0.0880596")
    }
    if (abs(st.psibar[k] / 0.1732898 - 1) > 1e-5) {
        FAIL = FAIL \ ("4: psi-bar_rho " + strofreal(st.psibar[k], "%12.8f") + " is not (a)'s 0.1732898")
    }

    // ---- 5: the three consumers agree ------------------------------------
    rep = dwlsmv_emptycell_pairs(Y, st.map)
    if (rep != (1 \ 0 \ 0)) FAIL = FAIL \ "5: the helper does not mark exactly y1~~y2"
    zr = abs(st.psibar[k]) / sqrt(st.Gamma[k, k])
    if (zr < 1e-2) FAIL = FAIL \ "5: the repaired pair's own z is not large; the test has no teeth"
    A = dwlsmv_gammadiag(st, Y)
    if (A[5, 1] != 0) {
        FAIL = FAIL \ ("5: Module 10 row 5 flags on repaired data (z " + strofreal(A[5, 2]) + ")")
    }
    zk = abs(st.psibar) :/ sqrt(diagonal(st.Gamma))
    zk[k] = 0
    if (A[5, 2] != max(zk)) FAIL = FAIL \ "5: row 5's z is not the max over the non-repaired entries"
    // the reported k must index a non-repaired entry that attains z
    // (adversarial review 2026-10-06: zeroing after the k search survived)
    if (A[5, 3] < 1 | A[5, 3] > rows(zk)) {
        FAIL = FAIL \ "5: row 5's k is out of range"
    }
    else if (zk[A[5, 3]] != A[5, 2] | A[5, 3] == k) {
        FAIL = FAIL \ ("5: row 5's k " + strofreal(A[5, 3]) + " does not attain z or is the repaired entry")
    }

    if (rows(FAIL)) {
        // one line per failure: printf() takes a scalar, and a vector here
        // aborts with r(3204) instead of listing the failures
        for (i = 1; i <= rows(FAIL); i++) {
            printf("{err}%s\n", FAIL[i])
        }
    }
    return(rows(FAIL))
}

end

mata: st_local("nfail", strofreal(dwlsmv_t_ec_main()))

local nerr 0

* ---- 1d: an empty row/column still aborts upstream (empty category) --------
mata: st_matrix("ER", dwlsmv_t_ec_data((0, 0 \ 7, 9), 2))
capture mata: dwls_sample_stats(st_matrix("ER"), (2 \ 2 \ 2))
if _rc != 3498 {
    display as error "test_emptycell: an empty row gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}

* ---- 2: the repair refuses a table that is not one-empty-cell ----------------
capture mata: dwlsmv_emptycell_repair((20, 0 \ 0, 20))
if _rc != 3498 {
    display as error "test_emptycell: repairing a diagonal pair gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}

* ---- 6: the edge, |rho-hat+| >= B, aborts (O3) ------------------------------
mata: st_matrix("BIG", dwlsmv_t_ec_data((1000, 1 \ 0, 1000), 2))
tempfile L
quietly log using "`L'.log", text replace name(ec)
capture noisily mata: dwls_sample_stats(st_matrix("BIG"), (2 \ 2 \ 2))
local rc = _rc
quietly log close ec
if `rc' != 3498 {
    display as error "test_emptycell: the repaired root beyond B gave rc `rc', not 3498"
    local nerr = `nerr' + 1
}
mata: st_local("hasmsg", strofreal(strpos(invtokens(cat(st_local("L") + ".log")', " "), ///
      "the repaired correlation reaches the bound") > 0))
if !`hasmsg' {
    display as error "test_emptycell: the beyond-B abort does not say the repaired correlation reaches the bound"
    local nerr = `nerr' + 1
}

if `nfail' | `nerr' {
    display as error _n "test_emptycell: `nfail' value failure(s), `nerr' error-path failure(s)"
    exit 9
}

display as result "test_emptycell: membership, repair, Stage 1 oracle, construction (a), consumer agreement, beyond-B abort"
exit 0
