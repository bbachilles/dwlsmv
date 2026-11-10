*! Module 10: the bread diagnostic, dwlsmv_breaddiag(Delta, Gamma), and the
*! one numerical-rank criterion, dwlsmv_numrank(M, sv, V).
*!
*! Executed under StataMP 18 on 2026-09-26 and it passes. Before the routines
*! existed it failed with r(3499), dwlsmv_numrank() not found.
*!
*! Spec Section 5, Module 10; E-060 and E-071; derivation
*! @sec-bread-rank. breaddiag returns 4 x 4, [flag, value, i, j]:
*!   1 npar > D          2 an identically zero Delta column
*!   3 numerical rank of W_D^1/2 Delta, 1e-8 (dwlsmv_numrank)
*!   4 sigma_min / sigma_max of W_D^1/2 Delta, flag .n
*! Rows 3-4 are judged on the non-zero columns only, and are `.` (not
*! evaluated) when row 1 flags or no column is non-zero: a structural
*! deficiency is never reported a second time by a numerical row.
*!
*!   1. NUMRANK on hand-built matrices: the criterion's edge from both sides
*!      (2e-8 kept, 5e-9 dropped), scale invariance, sv descending and V's
*!      ROWS the right singular vectors (M = U diag(sv) V).
*!   2. CASE D, fitted (marker): healthy rows; the ratio equals an
*!      independent eigenvalue path on the bread within its budget, and it is
*!      the WEIGHTED matrix's ratio, not the unweighted one's.
*!   3. PRODUCTION'S NOISE on case D, @sec-bread-rank's quantity:
*!      nu = ||W^1/2 (Delta - Delta_exact)||_2 / sigma_max, Delta_exact from
*!      the one-factor closed form, below 1e-10 (two decades under 1e-8, the
*!      margin Module 3 claims), and the ratio above 1e-6.
*!   4. THE DEFICIENT TWIN: case D's table with the marker's loading freed
*!      too, so (l, psi) -> (c l, psi / c^2) is an exact null direction of
*!      production's own Delta. Row 3 flags rank npar - 1; the ratio is below
*!      1e-10.
*!   5. A ZERO COLUMN: row 2 flags it, and rows 3-4 are judged on the rest.
*!   6. npar > D, and ALL COLUMNS ZERO: rows 3-4 not evaluated. 6b (after the
*!      2026-09-26 review): npar = D keeps rows 3-4; a 1e-12 column is not
*!      zero; a later zero column is located. Check 1 also pins numrank's
*!      criterion to within 20% of 1e-8.
*!   7. CONTRACT, rc 3498 each, for breaddiag and numrank.
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
        display as error "test_breaddiag: `f' failed to load"
        exit 9
    }
}

import delimited using "fixtures/data/case_D.csv", clear varnames(1) delimiters(",")
mata: YD = st_data(., .)

mata:

string scalar tb_spec(string rowvector lines)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= cols(lines); i++) {
        s = s + lines[i] + char(10)
    }
    return(s)
}

// Delta at theta for case D's one-factor model, from the closed form
// rho_ab = l_a psi l_b: d/d l_j = [a = j] l_b psi + [b = j] l_a psi,
// d/d psi = l_a l_b; the threshold block is the identity.
real matrix tb_exact(real colvector theta, struct dwlsmv_parmap scalar par,
                     struct dwlsmv_statmap scalar map)
{
    real matrix Dx, Lam, Psi, B
    real colvector vals, l
    real scalar i, a, b, r, c, ps

    vals = dwlsmv_parvals(theta, par)
    dwlsmv_matrices(vals, par, map.p, 1, Lam, Psi, B)
    l  = Lam[., 1]
    ps = Psi[1, 1]
    Dx = J(map.D, par.npar, 0)
    for (i = 1; i <= par.nall; i++) {
        if (!par.free[i]) continue
        c = par.fidx[i]
        if (par.ptype[i] == 4) {
            Dx[dwlsmv_pos_tau(map, par.row[i], par.col[i]), c] = 1
            continue
        }
        for (a = 1; a < map.p; a++) {
            for (b = a + 1; b <= map.p; b++) {
                r = dwlsmv_pos_rho(map, a, b)
                if (par.ptype[i] == 1) {
                    Dx[r, c] = (a == par.row[i]) * l[b] * ps + (b == par.row[i]) * l[a] * ps
                }
                else if (par.ptype[i] == 2) {
                    Dx[r, c] = l[a] * l[b]
                }
            }
        }
    }
    return(Dx)
}

real scalar tb_noteval(real matrix A, real scalar r)
{
    return(A[r, 1] == . & A[r, 2] == . & A[r, 3] == 0 & A[r, 4] == 0)
}

real scalar dwlsmv_t_breaddiag(real matrix YD)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par, tw
    struct dwlsmv_fit scalar fit
    string colvector FAIL
    real matrix A, M, V, Dm, Dx, Dt, E, G
    real colvector sv, w, ev, th
    real scalar i, np, D, rw, ru, tol, nu, rk, k

    FAIL = J(0, 1, "")

    // ---- 1. numrank --------------------------------------------------------------
    sv = J(0, 1, .); V = J(0, 0, .)
    if (dwlsmv_numrank(I(3), sv, V) != 3) FAIL = FAIL \ "1: I(3) is not rank 3"
    if (dwlsmv_numrank(diag((1 \ 2e-8)), sv, V) != 2) FAIL = FAIL \ "1: 2e-8 of the largest was dropped"
    if (dwlsmv_numrank(diag((1 \ 5e-9)), sv, V) != 1) FAIL = FAIL \ "1: 5e-9 of the largest was kept"
    if (dwlsmv_numrank(1e6 * diag((1 \ 5e-9)), sv, V) != 1) FAIL = FAIL \ "1: the criterion is not scale-free"
    // pinned to within 20% of 1e-8, not a 4x bracket (review M22)
    if (dwlsmv_numrank(diag((1 \ 1.2e-8)), sv, V) != 2) FAIL = FAIL \ "1: 1.2e-8 of the largest was dropped"
    if (dwlsmv_numrank(diag((1 \ 0.8e-8)), sv, V) != 1) FAIL = FAIL \ "1: 0.8e-8 of the largest was kept"
    M  = (3, 1 \ 1, 2 \ 0, 1)
    rk = dwlsmv_numrank(M, sv, V)
    if (!(rk == 2 & sv[1] >= sv[2])) FAIL = FAIL \ "1: sv is not descending"
    // M'M = V' diag(sv^2) V holds only if V's ROWS are the right vectors.
    if (max(abs(M' * M - V' * diag(sv :^ 2) * V)) > 16 * epsilon(1) * max(abs(M' * M))) {
        FAIL = FAIL \ "1: V's rows are not the right singular vectors"
    }

    // ---- 2. case D, fitted --------------------------------------------------------------
    st  = dwls_sample_stats(YD, (4 \ 4 \ 3 \ 3))
    par = dwls_identify(tb_spec(("factor f", "ind y1 y2 y3 y4",
                                  "lambda f y2 free", "lambda f y3 free",
                                  "lambda f y4 free", "psi f f free",
                                  "metric marker f y1")), st.map)
    fit = dwls_fit(st, par)
    if (fit.converged != 1) {
        FAIL = FAIL \ "2: case D did not converge; nothing below is meaningful"
        return(rows(FAIL))
    }
    Dm = fit.Delta
    np = cols(Dm)
    D  = rows(Dm)
    A  = dwlsmv_breaddiag(Dm, st.Gamma)
    w  = 1 :/ diagonal(st.Gamma)
    // independent path: eigenvalues of the bread, sigma^2 = lambda
    ev = symeigenvalues(Dm' * (w :* Dm))'
    rw = sqrt(min(ev) / max(ev))
    sv = J(0, 1, .); V = J(0, 0, .)
    (void) dwlsmv_numrank(Dm, sv, V)
    ru = sv[np] / sv[1]
    printf("{txt}  2. case D      : rows (%f, %f, %f, %f), (%f, %f), (%f, %f, %f), ratio %g (eig %g, unweighted %g)\n",
           A[1, 1], A[1, 2], A[1, 3], A[1, 4], A[2, 1], A[2, 2], A[3, 1], A[3, 2], A[3, 3],
           A[4, 2], rw, ru)
    if (rows(A) != 4 | cols(A) != 4) {
        FAIL = FAIL \ "2: breaddiag is not 4 x 4"
        return(rows(FAIL))
    }
    if (A[1, .] != (0, np, D, 0)) FAIL = FAIL \ "2: row 1 is not (0, npar, D, 0)"
    if (A[2, .] != (0, 0, 0, 0))  FAIL = FAIL \ "2: row 2 reports a zero column"
    if (A[3, .] != (0, np, np, 0)) FAIL = FAIL \ "2: row 3 is not full rank over npar columns"
    if (!(A[4, 1] == .n & A[4, 3] == np)) FAIL = FAIL \ "2: row 4's flag is not .n, or its column count is wrong"
    // eig of the bread carries 4 npar eps lambda_max; through the square root
    // that is 2 npar eps / r on the ratio, plus the svd's own npar eps.
    tol = 2 * np * epsilon(1) / rw + 4 * np * epsilon(1)
    if (abs(A[4, 2] - rw) > tol) FAIL = FAIL \ "2: the ratio is not W^1/2 Delta's"
    if (abs(A[4, 2] - ru) <= tol) FAIL = FAIL \ "2: weighted and unweighted ratios coincide; check 2 is void"

    // ---- 3. production's noise on the weighted matrix ------------------------------------
    Dx = tb_exact(fit.theta, par, st.map)
    E  = sqrt(w) :* (Dm - Dx)
    sv = J(0, 1, .); V = J(0, 0, .)
    (void) dwlsmv_numrank(E, sv, V)
    nu = sv[1]
    (void) dwlsmv_numrank(sqrt(w) :* Dm, sv, V)
    nu = nu / sv[1]
    printf("{txt}  3. noise       : nu %g (exact Delta from the closed form), ratio %g\n", nu, A[4, 2])
    if (!(nu < 1e-10)) FAIL = FAIL \ "3: production's noise on W^1/2 Delta is not two decades under 1e-8"
    if (!(A[4, 2] > 1e-6)) FAIL = FAIL \ "3: case D's ratio is not two decades over 1e-8"

    // ---- 4. the deficient twin: the marker's loading freed as well ------------------------
    tw = par
    for (i = 1; i <= tw.nall; i++) {
        if (tw.ptype[i] == 1 & tw.row[i] == 1 & !tw.free[i]) k = i
    }
    tw.free[k] = 1
    tw.npar    = tw.npar + 1
    tw.fidx[k] = tw.npar
    th = fit.theta \ tw.value[k]
    Dt = dwls_delta(th, tw, st.map)
    A  = dwlsmv_breaddiag(Dt, st.Gamma)
    printf("{txt}  4. twin        : row 3 (%f, %f, %f), ratio %g\n", A[3, 1], A[3, 2], A[3, 3], A[4, 2])
    if (A[3, .] != (1, np, np + 1, 0)) FAIL = FAIL \ "4: row 3 is not (1, npar, npar + 1, 0)"
    if (!(A[4, 2] < 1e-10)) FAIL = FAIL \ "4: the twin's ratio is not two decades under 1e-8"
    if (A[1, 1] != 0 | A[2, 1] != 0) FAIL = FAIL \ "4: a structural row fired on the twin"

    // ---- 5. a zero column --------------------------------------------------------------
    Dt = Dm
    Dt[., np - 1] = J(D, 1, 0)
    A  = dwlsmv_breaddiag(Dt, st.Gamma)
    printf("{txt}  5. zero column : row 2 (%f, %f, %f), row 3 (%f, %f, %f)\n",
           A[2, 1], A[2, 2], A[2, 3], A[3, 1], A[3, 2], A[3, 3])
    if (A[2, .] != (1, 1, np - 1, 0)) FAIL = FAIL \ "5: row 2 is not (1, 1, npar - 1, 0)"
    if (A[3, .] != (0, np - 1, np - 1, 0)) FAIL = FAIL \ "5: row 3 was not judged on the other columns"
    if (A[4, 3] != np - 1) FAIL = FAIL \ "5: row 4 was not judged on the other columns"

    // ---- 6. npar > D; all columns zero --------------------------------------------------
    G = (2, 0.5, 0 \ 0.5, 1, 0.2 \ 0, 0.2, 3)
    A = dwlsmv_breaddiag((1, 0, 2, 1 \ 0, 1, 1, 3 \ 2, 1, 0, 1), G)
    if (A[1, .] != (1, 4, 3, 0)) FAIL = FAIL \ "6: npar > D is not (1, 4, 3, 0)"
    if (!(tb_noteval(A, 3) & tb_noteval(A, 4))) FAIL = FAIL \ "6: rows 3-4 evaluated with npar > D"
    if (A[2, 1] != 0) FAIL = FAIL \ "6: row 2 not evaluated with npar > D"
    A = dwlsmv_breaddiag(J(3, 2, 0), G)
    if (A[2, .] != (1, 2, 1, 0)) FAIL = FAIL \ "6: two zero columns are not (1, 2, 1, 0)"
    if (!(tb_noteval(A, 3) & tb_noteval(A, 4))) FAIL = FAIL \ "6: rows 3-4 evaluated with no column left"

    // ---- 6b. added after the 2026-09-26 review --------------------------------------
    // npar = D is not npar > D: a saturated model keeps rows 3-4 (M19).
    A = dwlsmv_breaddiag((1, 0, 2 \ 0, 1, 1 \ 2, 1, 0), G)
    if (A[1, .] != (0, 3, 3, 0)) FAIL = FAIL \ "6b: npar = D is flagged as npar > D"
    if (A[3, 1] == . | A[4, 1] != .n) FAIL = FAIL \ "6b: rows 3-4 not evaluated at npar = D"
    // "identically zero" is exact: a column of 1e-12 is not zero (M18), and a
    // zero column that is not the first is located (row 2's i).
    A = dwlsmv_breaddiag((1, 1e-12, 0 \ 0, 1e-12, 0 \ 2, 0, 0), G)
    if (A[2, .] != (1, 1, 3, 0)) FAIL = FAIL \ "6b: row 2 is not (1, 1, 3, 0): a tiny column counted, or the location is wrong"

    printf("{txt}test_breaddiag: %f failure(s)\n", rows(FAIL))
    for (i = 1; i <= rows(FAIL); i++) {
        printf("{err}    %s\n", FAIL[i])
    }
    return(rows(FAIL))
}

end

mata: st_local("nfail", strofreal(dwlsmv_t_breaddiag(YD)))

* ---------------------------------------------------------------------------
* 7. CONTRACT -- rc 3498 each.
* ---------------------------------------------------------------------------
local nerr 0
mata: TG = (2, 0.5 \ 0.5, 1)
local bad1 "dwlsmv_breaddiag((1, 0 \ 0, 1 \ 1, 1), TG)"
local bad2 "dwlsmv_breaddiag((1, 0 \ 0, 1), (2, 0.5, 0 \ 0.5, 1, 0))"
local bad3 "dwlsmv_breaddiag((1, 0 \ 0, 1), (2, 0.5 \ 0.5, 0))"
local bad4 "dwlsmv_breaddiag((1, . \ 0, 1), TG)"
local bad5 "dwlsmv_breaddiag(J(2, 0, 0), TG)"
forvalues b = 1/5 {
    capture mata: T_TMP = `bad`b''
    if _rc != 3498 {
        display as error "test_breaddiag: contract violation #`b' gave rc " _rc ", not 3498"
        local nerr = `nerr' + 1
    }
}
mata: T_SV = J(0, 1, .)
mata: T_V = J(0, 0, .)
capture mata: T_TMP = dwlsmv_numrank((1, . \ 0, 1), T_SV, T_V)
if _rc != 3498 {
    display as error "test_breaddiag: numrank on a missing value gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}
capture mata: T_TMP = dwlsmv_numrank(J(2, 0, 0), T_SV, T_V)
if _rc != 3498 {
    display as error "test_breaddiag: numrank on an empty matrix gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}

if `nfail' | `nerr' {
    display as error _n "test_breaddiag: `nfail' numeric failure(s), `nerr' error-path failure(s)"
    exit 9
}

display as result "test_breaddiag: the bread's rows, each against its cause or the one rank criterion"
exit 0
