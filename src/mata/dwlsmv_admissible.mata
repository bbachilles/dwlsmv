*! dwlsmv -- is a solution admissible? (Module 7, spec Section 5)
*!
*! Executed under StataMP 18 on 2026-09-23. tests/unit/test_fit.do check 12: all six conditions pass on five fitted solutions, and five hand-built inadmissible points each report their known value (to 1e-15 where it is a computed number), row 4 included since 12f.
*!
*!     real matrix dwlsmv_admissible(real colvector theta,
*!                                   struct dwlsmv_parmap  scalar par,
*!                                   struct dwlsmv_statmap scalar map)
*!
*! Spec Section 5, Module 7: before a fitted theta-hat is accepted, check for
*! negative residual variances, non-positive-definite implied matrices, factor
*! correlations outside [-1, 1], a singular (I - B) and disordered thresholds,
*! and SURFACE them as diagnostics rather than silently accepting the solution.
*!
*! A FUNCTION, NOT A FIELD. Decided 2026-09-23: every condition is a function
*! of (theta, par, map), so storing the verdict on struct dwlsmv_fit would be a
*! second copy of facts the fit already carries -- the rule dwlsmv_nfac() and
*! dwlsmv_koff() follow. dwls_fit() does not call this and prints nothing;
*! Module 11 calls it on fit.theta and shows dwlsmv_admissadvice()'s notes.
*! It takes theta rather than a fit so it can be asked about any point,
*! including a non-converged last iterate.
*!
*! Returns 6 x 4, one row per condition, always in this order:
*!
*!     row  condition                         value              i, j
*!      1   a derived Theta_jj < 0             min Theta_jj       indicator, 0
*!      2   a factor correlation |r| > 1       the largest |r|    factor, factor
*!      3   Phi_eta not positive definite      min eigenvalue     0, 0
*!      4   R_y* not positive definite         min eigenvalue     0, 0
*!      5   (I - B) singular                   missing            0, 0
*!      6   thresholds not increasing          min gap            indicator, 0
*!
*! Column 1 is the flag: 1 the condition holds (the solution is NOT
*! admissible on it), 0 it does not, MISSING it could not be evaluated. Rows
*! 1-4 are missing when (I - B) is singular, because Phi_eta does not exist
*! then; they are not reported as passing. Row 2 is 0 with a missing value on
*! a one-factor model, which has no factor correlation to check. Column 2 is
*! the worst value found whether or not it fails, so a near miss is visible.
*!
*! THE DEFINITIONS, each the spec's:
*!
*!   Theta_jj = 1 - (Lambda Phi Lambda')_jj, derived under delta (Section
*!   4.6). A negative one is a Heywood case: the factors explain more than
*!   all of y_j*'s variance. dwlsmv_rystar() writes 1 down R's diagonal
*!   regardless, so sigma(theta) never shows it; this is the only place it
*!   is seen.
*!
*!   The factor correlation is Phi_ij / sqrt(Phi_ii Phi_jj) on the TOTAL
*!   covariance Phi_eta = (I-B)^-1 Psi (I-B)^-T, not on Psi -- the disturbance
*!   covariance of an endogenous factor is not a correlation between factors.
*!   Pairs with a non-positive Phi_ii are skipped; row 3 reports those.
*!
*!   "Positive definite" is minimum eigenvalue > 0, exactly. No tolerance: a
*!   solution on the boundary is reported as on the boundary, with its
*!   eigenvalue beside it, and what to do about it is the caller's decision.
*!
*! WHICH ROWS CAN FIRE IN v1, so nobody mistakes a backstop for a live check.
*! Rows 1-4 can. Row 6 cannot at a converged theta-hat: tau-hat = s_tau exactly
*! (test_fit.do check 6) and s_tau is increasing by construction; it can fire
*! at an arbitrary theta. Row 5 cannot for a recursive B, whose (I - B) is
*! unit triangular, and since 2026-09-23 dwls_identify() rejects every cyclic
*! B (spec 4.8), so in v1 row 5 fires only on a hand-built table like
*! test_fit.do 12e's. There is deliberately no stability row: cycles are out of
*! v1 scope, not checked for stability.
*!
*! Load order: dwlsmv_types, dwlsmv_pos_tau, dwlsmv_nfac, dwlsmv_parvals,
*! dwlsmv_matrices, dwlsmv_rystar.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_admissible(real colvector theta,
                              struct dwlsmv_parmap  scalar par,
                              struct dwlsmv_statmap scalar map)
{
    real matrix A, Lam, Psi, B, Ai, Phi, C, R
    real colvector vals, t, th
    real scalar nf, i, j, c, d, w, wi, wj, r

    A    = J(6, 4, 0)
    vals = dwlsmv_parvals(theta, par)
    nf   = dwlsmv_nfac(par)
    dwlsmv_matrices(vals, par, map.p, nf, Lam, Psi, B)

    // ---- 6. thresholds increasing, within each variable -------------------
    t = J(map.T, 1, .)
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 4) {
            t[dwlsmv_pos_tau(map, par.row[i], par.col[i])] = vals[i]
        }
    }
    w  = .
    wi = 0
    for (j = 1; j <= map.p; j++) {
        for (c = 2; c <= map.C[j]; c++) {
            d = t[map.Coff[j] + c] - t[map.Coff[j] + c - 1]
            if (d < w) {
                w  = d
                wi = j
            }
        }
    }
    A[6, .] = (w <= 0, w, wi, 0)

    // ---- 5. (I - B) singular ------------------------------------------------
    // Checked here rather than left to dwlsmv_rystar(), which aborts: this
    // routine reports, it does not stop.
    Ai = luinv(I(nf) - B)
    if (missing(Ai)) {
        A[5, .]      = (1, ., 0, 0)
        A[1::4, 1]   = J(4, 1, .)
        A[1::4, 2]   = J(4, 1, .)
        return(A)
    }
    A[5, .] = (0, ., 0, 0)

    Phi = Ai * Psi * Ai'
    C   = Lam * Phi * Lam'

    // ---- 1. derived residual variances --------------------------------------
    th = 1 :- diagonal(C)
    w  = th[1]
    wi = 1
    for (j = 2; j <= rows(th); j++) {
        if (th[j] < w) {
            w  = th[j]
            wi = j
        }
    }
    A[1, .] = (w < 0, w, wi, 0)

    // ---- 2. factor correlations ---------------------------------------------
    w  = .
    wi = 0
    wj = 0
    for (i = 1; i <= nf - 1; i++) {
        for (j = i + 1; j <= nf; j++) {
            if (Phi[i, i] > 0 & Phi[j, j] > 0) {
                r = Phi[i, j] / sqrt(Phi[i, i] * Phi[j, j])
                if (missing(w) | abs(r) > abs(w)) {
                    w  = r
                    wi = i
                    wj = j
                }
            }
        }
    }
    A[2, .] = ((missing(w) ? 0 : abs(w) > 1), w, wi, wj)

    // ---- 3, 4. positive definiteness ----------------------------------------
    w       = min(symeigenvalues((Phi + Phi') / 2))
    A[3, .] = (w <= 0, w, 0, 0)

    R       = dwlsmv_rystar(Lam, Psi, B)
    w       = min(symeigenvalues(R))
    A[4, .] = (w <= 0, w, 0, 0)

    return(A)
}

end
