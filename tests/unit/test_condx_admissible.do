*! Executed 2026-11-14 by the main session (Stata 18 MP); 11 of 11 mutants killed. Conditional mode (one 0/1 covariate, one group): dwlsmv_admissible().
*!
*! No fixtures and no data: dwlsmv_admissible() is a function of (theta, table, map), so every case is a
*! hand-written canonical table with a theta chosen here and a closed-form expected value. The map comes
*! from dwlsmv_statmap_build(K, q). Nothing from the frozen condx preregistration is read.
*!
*! WHAT IS CERTIFIED (derivation @sec-condx-full C6): in conditional mode the six conditions are on the
*! conditional quantities (Var(y*|x) = 1, Theta_jj = 1 - (Lambda Phi Lambda')_jj, Phi the residual factor
*! covariance), and gamma and kappa enter none of them. NOT certified: admissibility of a fit found by
*! dwls_fit() (no fit is run), Module 11 notes, several groups with a covariate (refused), fit indices.
*! The expected values are closed forms (1 - lambda^2, 1 +/- psi_12, 0.6 / sqrt(1.36), threshold gaps);
*! the rows-4 value is computed here from explicit matrices and is a diagnostic (same eigenvalue
*! routine), not the deciding check. Tolerances 1e-14 for computed numbers; 0 for gamma/kappa invariance.
*!
*! Checks:
*!   1  admissible conditional population (one factor, gamma and kappa free): no flag; row values.
*!   2  Heywood case (lambda_1 = 1.2): row 1 flagged at item 1 with value 1 - 1.44; the other rows
*!      that do not depend on it are unchanged.
*!   3  non-PSD Psi with every |correlation| <= 0.9 (three factors): row 3 flagged with min eigenvalue
*!      1 - 0.9 * 2 = -0.8, row 2 NOT flagged (the cases are isolated); and psi_12 = 1.3: rows 2 and 3.
*!   4  gamma and kappa do not matter: A is identical for (0, 0), (0.4, 0.25), (5, -3) and for the same
*!      table without gamma/kappa rows on an unconditional map (q = 0).
*!   5  the conditional Phi is (I-B)^-1 Psi (I-B)^-T, not Psi: two factors with beta = 0.6, unit
*!      disturbances: Phi_22 = 1.36, item 3's Theta = 1 - 0.81 * 1.36 < 0 (a Heywood case caused by beta),
*!      factor correlation 0.6 / sqrt(1.36).
*!   6  disordered thresholds (row 6) on a conditional map.
*!   7  refusals, each rc == 3498 exactly with the message read from a log: dwlsmv_mg_admissible() on a
*!      conditional map (several groups with a covariate stay refused); dwlsmv_admissible() with q = 2;
*!      gamma/kappa rows against an unconditional map.
*!
*! MUTANTS for the main session (edit, run, revert; each must turn this file red):
*!   M1  dwlsmv_admissible: restore `dwlsmv_s2_guard(map, "dwlsmv_admissible")`: checks 1 to 6 abort at 3498.
*!   M2  dwlsmv_admissible: delete the dwlsmv_cx1_check line: refusals 7b and 7c fail (rc is not 3498).
*!   M3  dwlsmv_admissible: `A[1, .] = (w < 0, w, wi, 0)` -> `w < -0.5`: check 2 (Heywood unflagged).
*!   M4  dwlsmv_admissible: `Phi = Ai * Psi * Ai'` -> `Phi = Psi`: check 5 (Phi_22, Theta, r).
*!   M5  dwlsmv_admissible: `A[3, .] = (w <= 0, w, 0, 0)` -> `w < -1`: check 3 (min eigenvalue -0.8 unflagged).
*!   M6  dwlsmv_admissible: `abs(w) > 1)` in row 2 -> `abs(w) > 1.5`: check 3 (psi_12 = 1.3 unflagged).
*!   M7  dwlsmv_matrices: `par.ptype[i] == 3` -> `par.ptype[i] == 3 | par.ptype[i] == 8` (gamma into B): checks 1, 4.
*!   M8  dwlsmv_admissible: row 6 `(w <= 0, w, wi, 0)` -> `(w < -1, w, wi, 0)`: check 6.
*!   M9  dwlsmv_admissible: `th = 1 :- diagonal(C)` -> `th = 1 :- diagonal(C) :* 0.5`: checks 1, 2, 5.
*!   M10 dwlsmv_cx1_check: `q > 1` -> `q > 2`: refusal 7b fails.
*!   M11 dwlsmv_mg_admissible: delete its dwlsmv_s2_guard call: refusal 7a fails.
*!
*! Return codes: 0 pass, 9 fail.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

* Load every source, types first (as test_condx_fit.do does).
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
        display as error "test_condx_admissible: `f' failed to load"
        exit 9
    }
}

capture mata: (void) dwlsmv_cx1_check(dwlsmv_parmap(), dwlsmv_statmap(0, 1), "load check")
if _rc {
    display as error "test_condx_admissible: dwlsmv_cx1_check() is not defined (src/mata/dwlsmv_cx1_check.mata missing?)"
    exit 9
}

mata:

// ---- helpers -------------------------------------------------------------------

real scalar ca_fail(string scalar what)
{
    printf("{err}test_condx_admissible: %s\n", what)
    return(1)
}

real scalar ca_chk(real scalar ok, string scalar what)
{
    if (!ok) {
        return(ca_fail(what))
    }
    return(0)
}

// |a - b| <= tol (a missing counts as a failure)
real scalar ca_near(real scalar a, real scalar b, real scalar tol)
{
    if (missing(a) | missing(b)) {
        return(0)
    }
    return(abs(a - b) <= tol)
}

// Two admissibility results equal entry by entry, a missing matching a missing.
real scalar ca_same(real matrix A, real matrix B)
{
    real scalar i, j

    if (rows(A) != rows(B) | cols(A) != cols(B)) {
        return(0)
    }
    for (i = 1; i <= rows(A); i++) {
        for (j = 1; j <= cols(A); j++) {
            if (missing(A[i, j]) | missing(B[i, j])) {
                if (missing(A[i, j]) != missing(B[i, j])) {
                    return(0)
                }
            }
            else if (A[i, j] != B[i, j]) {
                return(0)
            }
        }
    }
    return(1)
}

// A hand-built parameter table. free rows get theta in table order; fidx follows.
struct dwlsmv_parmap scalar ca_mk(real colvector pt, real colvector rw, real colvector cl,
                                  real colvector fr, real colvector vl)
{
    struct dwlsmv_parmap scalar par
    real scalar i, n, k

    n = rows(pt)
    par.nall  = n
    par.npar  = sum(fr)
    par.name  = J(n, 1, "")
    for (i = 1; i <= n; i++) {
        par.name[i] = "row" + strofreal(i)
    }
    par.ptype = pt
    par.row   = rw
    par.col   = cl
    par.free  = fr
    par.value = vl
    par.fidx  = J(n, 1, 0)
    k = 0
    for (i = 1; i <= n; i++) {
        if (fr[i]) {
            k = k + 1
            par.fidx[i] = k
        }
    }
    par.grp = J(n, 1, 1)
    return(par)
}

// One factor, K = (3,3,2): tau (5) | lambda (3) | psi_11 fixed 1 | [gamma F x | kappa y3 x] when cx.
// theta = (tau, lambda [, gamma, kappa]).
struct dwlsmv_parmap scalar ca_par1(real scalar cx)
{
    real colvector pt, rw, cl, fr, vl

    pt = (4\4\4\4\4 \ 1\1\1 \ 2)
    rw = (1\1\2\2\3 \ 1\2\3 \ 1)
    cl = (1\2\1\2\1 \ 1\1\1 \ 1)
    fr = (J(8, 1, 1) \ 0)
    vl = (J(8, 1, 0) \ 1)
    if (cx) {
        pt = (pt \ 8 \ 9)
        rw = (rw \ 1 \ 3)
        cl = (cl \ 1 \ 1)
        fr = (fr \ 1 \ 1)
        vl = (vl \ 0 \ 0)
    }
    return(ca_mk(pt, rw, cl, fr, vl))
}

real colvector ca_th1(real colvector tau, real colvector lam, real scalar cx, real scalar gam, real scalar kap)
{
    if (cx) {
        return((tau \ lam \ gam \ kap))
    }
    return((tau \ lam))
}

// Three factors, two items each, K = 2: tau (6) | lambda (6) | psi_ii fixed 1 | psi_12, psi_13, psi_23 free
// | gamma F1 x | kappa y2 x.   theta = (tau, lambda, psi_12, psi_13, psi_23, gamma, kappa).
struct dwlsmv_parmap scalar ca_par3()
{
    real colvector pt, rw, cl, fr, vl

    pt = (J(6, 1, 4) \ J(6, 1, 1) \ J(6, 1, 2) \ 8 \ 9)
    rw = ((1::6) \ (1::6) \ (1\2\3\1\1\2) \ 1 \ 2)
    cl = (J(6, 1, 1) \ (1\1\2\2\3\3) \ (1\2\3\2\3\3) \ 1 \ 1)
    fr = (J(12, 1, 1) \ (0\0\0\1\1\1) \ 1 \ 1)
    vl = (J(12, 1, 0) \ (1\1\1\0\0\0) \ 0 \ 0)
    return(ca_mk(pt, rw, cl, fr, vl))
}

// Two factors, two items each, K = 2, F2 <- F1 (beta), unit DISTURBANCE variances, gamma F1 x, kappa y2 x.
// theta = (tau 4, lambda 4, beta, gamma, kappa).
struct dwlsmv_parmap scalar ca_par2()
{
    real colvector pt, rw, cl, fr, vl

    pt = (J(4, 1, 4) \ J(4, 1, 1) \ 3 \ 2 \ 2 \ 8 \ 9)
    rw = ((1::4) \ (1::4) \ 2 \ 1 \ 2 \ 1 \ 2)
    cl = (J(4, 1, 1) \ (1\1\2\2) \ 1 \ 1 \ 2 \ 1 \ 1)
    fr = (J(8, 1, 1) \ 1 \ 0 \ 0 \ 1 \ 1)
    vl = (J(8, 1, 0) \ 0 \ 1 \ 1 \ 0 \ 0)
    return(ca_mk(pt, rw, cl, fr, vl))
}

// ---- checks 1, 2, 4, 6: the one-factor model -----------------------------------------------

real scalar ca_one()
{
    struct dwlsmv_parmap scalar p0, p1
    struct dwlsmv_statmap scalar m0, m1
    real matrix A, A2, R
    real colvector tau, lam, th
    real scalar nf

    nf  = 0
    p1  = ca_par1(1)
    p0  = ca_par1(0)
    m1  = dwlsmv_statmap_build((3 \ 3 \ 2), 1)
    m0  = dwlsmv_statmap_build((3 \ 3 \ 2), 0)
    tau = (-0.5 \ 0.5 \ -0.3 \ 0.6 \ 0.2)
    lam = (0.8 \ 0.7 \ 0.6)

    // ---- 1: admissible
    A = dwlsmv_admissible(ca_th1(tau, lam, 1, 0.4, 0.25), p1, m1)
    nf = nf + ca_chk(rows(A) == 6 & cols(A) == 4, "check 1: A is not 6 x 4")
    nf = nf + ca_chk(A[1, 1] == 0 & A[2, 1] == 0 & A[3, 1] == 0 & A[4, 1] == 0 & A[5, 1] == 0 & A[6, 1] == 0,
                     "check 1: an admissible conditional population is flagged: " + invtokens(strofreal(A[., 1])'))
    nf = nf + ca_chk(ca_near(A[1, 2], 0.36, 1e-14) & A[1, 3] == 1,
                     "check 1: row 1 should be min Theta = 1 - 0.8^2 at item 1, got " + strofreal(A[1, 2]) + " at " + strofreal(A[1, 3]))
    nf = nf + ca_chk(A[2, 1] == 0 & missing(A[2, 2]), "check 1: row 2 (one factor) should be 0 with a missing value")
    nf = nf + ca_chk(ca_near(A[3, 2], 1, 1e-14), "check 1: row 3 should be Phi = 1, got " + strofreal(A[3, 2]))
    R = lam * lam'
    R[1, 1] = 1
    R[2, 2] = 1
    R[3, 3] = 1
    nf = nf + ca_chk(ca_near(A[4, 2], min(symeigenvalues(R)), 1e-14), "check 1: row 4 value differs from R's minimum eigenvalue")
    nf = nf + ca_chk(A[5, 1] == 0 & missing(A[5, 2]), "check 1: row 5 should be 0 with a missing value")
    nf = nf + ca_chk(ca_near(A[6, 2], 0.9, 1e-14) & A[6, 3] == 2,
                     "check 1: row 6 should be the smallest gap 0.9 at item 2, got " + strofreal(A[6, 2]) + " at " + strofreal(A[6, 3]))

    // ---- 2: Heywood, conditional common variance of item 1 above 1
    lam = (1.2 \ 0.7 \ 0.6)
    A2 = dwlsmv_admissible(ca_th1(tau, lam, 1, 0.4, 0.25), p1, m1)
    nf = nf + ca_chk(A2[1, 1] == 1 & A2[1, 3] == 1 & ca_near(A2[1, 2], 1 - 1.44, 1e-14),
                     "check 2: the Heywood case is not flagged at item 1 with -0.44; got flag " + strofreal(A2[1, 1]) +
                     ", value " + strofreal(A2[1, 2]))
    nf = nf + ca_chk(A2[3, 1] == 0 & A2[5, 1] == 0 & A2[6, 1] == 0 & A2[2, 1] == 0,
                     "check 2: a Heywood case should leave rows 2, 3, 5 and 6 unflagged")

    // ---- 4: gamma and kappa do not matter; neither does the table's having them
    lam = (0.8 \ 0.7 \ 0.6)
    A2 = dwlsmv_admissible(ca_th1(tau, lam, 1, 0, 0), p1, m1)
    nf = nf + ca_chk(ca_same(A, A2), "check 4: A changed between (gamma, kappa) = (0.4, 0.25) and (0, 0)")
    A2 = dwlsmv_admissible(ca_th1(tau, lam, 1, 5, -3), p1, m1)
    nf = nf + ca_chk(ca_same(A, A2), "check 4: A changed between (gamma, kappa) = (0.4, 0.25) and (5, -3)")
    A2 = dwlsmv_admissible(ca_th1(tau, lam, 0, 0, 0), p0, m0)
    nf = nf + ca_chk(ca_same(A, A2), "check 4: the conditional A differs from the same table without gamma, kappa on an unconditional map")
    // the Heywood case is just as insensitive to the slopes
    lam = (1.2 \ 0.7 \ 0.6)
    A  = dwlsmv_admissible(ca_th1(tau, lam, 1, 0.4, 0.25), p1, m1)
    A2 = dwlsmv_admissible(ca_th1(tau, lam, 1, -5, 3), p1, m1)
    nf = nf + ca_chk(ca_same(A, A2), "check 4: the Heywood case's A changed with gamma, kappa")

    // ---- 6: disordered thresholds (item 1's second threshold below its first)
    lam = (0.8 \ 0.7 \ 0.6)
    th  = ca_th1((-0.5 \ -0.7 \ -0.3 \ 0.6 \ 0.2), lam, 1, 0.4, 0.25)
    A2  = dwlsmv_admissible(th, p1, m1)
    nf = nf + ca_chk(A2[6, 1] == 1 & A2[6, 3] == 1 & ca_near(A2[6, 2], -0.2, 1e-14),
                     "check 6: disordered thresholds are not flagged at item 1 with gap -0.2; got flag " + strofreal(A2[6, 1]) +
                     ", value " + strofreal(A2[6, 2]) + ", item " + strofreal(A2[6, 3]))
    nf = nf + ca_chk(A2[1, 1] == 0 & A2[3, 1] == 0, "check 6: disordered thresholds should leave rows 1 and 3 unflagged")
    return(nf)
}

// ---- check 3: non-PSD Psi ---------------------------------------------------------------------

real scalar ca_psi()
{
    struct dwlsmv_parmap scalar par
    struct dwlsmv_statmap scalar map
    real matrix A
    real colvector tau, lam, th
    real scalar nf

    nf  = 0
    par = ca_par3()
    map = dwlsmv_statmap_build(J(6, 1, 2), 1)
    tau = (-0.5 \ 0.5 \ -0.3 \ 0.6 \ 0.2 \ 0)
    lam = (0.8 \ 0.7 \ 0.6 \ 0.7 \ 0.8 \ 0.5)

    // (0.9, 0.9, -0.9): every |r| <= 0.9, eigenvalues of Phi are 1.9, 1.9 and 1 - 1.8 = -0.8
    th = (tau \ lam \ 0.9 \ 0.9 \ -0.9 \ 0.4 \ 0.25)
    A  = dwlsmv_admissible(th, par, map)
    nf = nf + ca_chk(A[3, 1] == 1 & ca_near(A[3, 2], -0.8, 1e-14),
                     "check 3: a non-PSD Psi is not flagged with min eigenvalue -0.8; got flag " + strofreal(A[3, 1]) +
                     ", value " + strofreal(A[3, 2]))
    nf = nf + ca_chk(A[2, 1] == 0 & ca_near(abs(A[2, 2]), 0.9, 1e-14),
                     "check 3: row 2 should NOT fire on correlations of at most 0.9; got flag " + strofreal(A[2, 1]) +
                     ", value " + strofreal(A[2, 2]))
    nf = nf + ca_chk(A[1, 1] == 0 & A[5, 1] == 0, "check 3: rows 1 and 5 should be unflagged")

    // one correlation above 1: rows 2 and 3 both fire (eigenvalue 1 - 1.3 = -0.3)
    th = (tau \ lam \ 1.3 \ 0 \ 0 \ 0.4 \ 0.25)
    A  = dwlsmv_admissible(th, par, map)
    nf = nf + ca_chk(A[2, 1] == 1 & ca_near(A[2, 2], 1.3, 1e-14) & A[2, 3] == 1 & A[2, 4] == 2,
                     "check 3: psi_12 = 1.3 is not flagged in row 2 at factors (1, 2); got flag " + strofreal(A[2, 1]) +
                     ", value " + strofreal(A[2, 2]))
    nf = nf + ca_chk(A[3, 1] == 1 & ca_near(A[3, 2], -0.3, 1e-14),
                     "check 3: psi_12 = 1.3 is not flagged in row 3 with -0.3; got " + strofreal(A[3, 2]))

    // a PSD Psi (0.9, 0.9, 0.7) passes row 3; slopes irrelevant
    th = (tau \ lam \ 0.5 \ 0.4 \ 0.3 \ 4 \ -2)
    A  = dwlsmv_admissible(th, par, map)
    nf = nf + ca_chk(A[3, 1] == 0 & A[2, 1] == 0 & A[4, 1] == 0 & A[1, 1] == 0,
                     "check 3: a positive definite Psi is flagged: " + invtokens(strofreal(A[., 1])'))
    return(nf)
}

// ---- check 5: Phi is the conditional covariance (I-B)^-1 Psi (I-B)^-T ------------------------------

real scalar ca_beta()
{
    struct dwlsmv_parmap scalar par
    struct dwlsmv_statmap scalar map
    real matrix A
    real colvector tau, lam, th
    real scalar nf

    nf  = 0
    par = ca_par2()
    map = dwlsmv_statmap_build(J(4, 1, 2), 1)
    tau = (-0.5 \ 0.5 \ -0.3 \ 0.6)
    lam = (0.8 \ 0.7 \ 0.9 \ 0.9)
    th  = (tau \ lam \ 0.6 \ 0.4 \ 0.25)       // beta 0.6, gamma 0.4, kappa 0.25
    A   = dwlsmv_admissible(th, par, map)

    // Phi = [1 .6; .6 1.36]; Theta_33 = Theta_44 = 1 - 0.81 * 1.36 = -0.1016 (first at item 3)
    nf = nf + ca_chk(A[1, 1] == 1 & A[1, 3] == 3 & ca_near(A[1, 2], 1 - 0.81 * 1.36, 1e-14),
                     "check 5: Theta_33 should be 1 - 0.81 * 1.36 at item 3; got flag " + strofreal(A[1, 1]) +
                     ", value " + strofreal(A[1, 2]) + ", item " + strofreal(A[1, 3]))
    nf = nf + ca_chk(A[2, 1] == 0 & ca_near(A[2, 2], 0.6 / sqrt(1.36), 1e-14) & A[2, 3] == 1 & A[2, 4] == 2,
                     "check 5: the factor correlation should be 0.6 / sqrt(1.36) on the total conditional Phi; got " + strofreal(A[2, 2]))
    nf = nf + ca_chk(A[3, 1] == 0 & A[5, 1] == 0, "check 5: Phi is positive definite and (I - B) nonsingular")

    // the same table with beta = 0 is admissible: the flag above is caused by beta
    th[9] = 0
    A = dwlsmv_admissible(th, par, map)
    nf = nf + ca_chk(A[1, 1] == 0 & ca_near(A[1, 2], 1 - 0.81, 1e-14),
                     "check 5: with beta = 0, Theta_33 should be 1 - 0.81 = 0.19 and unflagged")
    return(nf)
}

// ---- check 7: refusals ----------------------------------------------------------------------------

real scalar ca_r1()               // several groups with a covariate: the multi-group entry, a conditional map
{
    struct dwlsmv_parmap scalar par
    struct dwlsmv_statmap scalar map
    real matrix A

    par = ca_par1(1)
    map = dwlsmv_statmap_build((3 \ 3 \ 2), 1)
    A   = dwlsmv_mg_admissible(J(par.npar, 1, 0.1), par, 1, map)
    return(0)
}

real scalar ca_r2()               // q = 2
{
    struct dwlsmv_parmap scalar par
    struct dwlsmv_statmap scalar map
    real matrix A

    par = ca_par1(1)
    map = dwlsmv_statmap_build((3 \ 3 \ 2), 2)
    A   = dwlsmv_admissible(J(par.npar, 1, 0.1), par, map)
    return(0)
}

real scalar ca_r3()               // gamma and kappa rows against an unconditional map
{
    struct dwlsmv_parmap scalar par
    struct dwlsmv_statmap scalar map
    real matrix A

    par = ca_par1(1)
    map = dwlsmv_statmap_build((3 \ 3 \ 2), 0)
    A   = dwlsmv_admissible(J(par.npar, 1, 0.1), par, map)
    return(0)
}

// Does the log at f contain the text has? Returns the number of failures.
real scalar ca_msg(string scalar f, string scalar has, string scalar tag)
{
    string scalar txt

    txt = invtokens(cat(f)', char(10))
    if (!strpos(txt, has)) {
        printf("{err}test_condx_admissible: case %s: the log lacks: %s\n", tag, has)
        return(1)
    }
    return(0)
}

end

local nfail 0

* ---- 1 to 6 ---------------------------------------------------------------------------------------
foreach c in "ca_one()" "ca_psi()" "ca_beta()" {
    capture noisily mata: st_local("nf", strofreal(`c'))
    if _rc {
        display as error "test_condx_admissible: `c' aborted, rc " _rc
        local nfail = `nfail' + 1
    }
    else {
        display as text "  `c': `nf' failure(s)"
        local nfail = `nfail' + `nf'
    }
}

* ---- 7: refusals, rc exactly 3498 and the message on the log ------------------------------------------
tempfile L
local n 0
local ++n
local c`n' ca_r1()
local w`n' dwlsmv_mg_admissible: the index map of group 1 is conditional
local ++n
local c`n' ca_r2()
local w`n' dwlsmv_admissible: the index map has q = 2
local ++n
local c`n' ca_r3()
local w`n' dwlsmv_admissible: the table has gamma or kappa rows but the index map has no covariate

forvalues i = 1/`n' {
    quietly log using "`L'_`i'.log", text replace name(ca)
    capture noisily mata: (void) `c`i''
    local rc = _rc
    quietly log close ca
    display as text "  `c`i'': rc `rc'"
    if `rc' != 3498 {
        display as error "test_condx_admissible: `c`i'' gave rc `rc', wanted 3498"
        local nfail = `nfail' + 1
    }
    else {
        local tmpL "`L'_`i'.log"
        local want "`w`i''"
        capture noisily mata: st_local("nm", strofreal(ca_msg(st_local("tmpL"), st_local("want"), "`i'")))
        if _rc {
            display as error "test_condx_admissible: reading the log of `c`i'' failed, rc " _rc
            local nfail = `nfail' + 1
        }
        else local nfail = `nfail' + `nm'
    }
}

if `nfail' {
    display as error _n "test_condx_admissible: `nfail' failure(s)"
    exit 9
}

display as result "test_condx_admissible: dwlsmv_admissible() on a conditional map flags a Heywood case, a non-PSD Psi, a beta-induced Heywood case and disordered thresholds, ignores gamma and kappa, and refuses several groups, q = 2 and gamma/kappa on an unconditional map at rc 3498"
exit 0
