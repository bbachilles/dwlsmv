*! dwlsmv -- conditional mode, Module 4: the N x D influence matrix in the
*! Hessian (M-estimator) form (@sec-condx Items 1-3; @sec-condx-hess; spec 4.8
*! "Stage 1 covariance: the Hessian (M-estimator) form", decided 2026-10-25)
*!
*! Executed 2026-11-10, suite green (drafted by a read-only agent, run by the
*! main session).
*! Its column sums of squares are Gamma (dwls_gamma_cx()), which
*! tests/unit/test_condx_stage1.do compares entry by entry with the analytic
*! R reference, and its negative controls (ctl) must fail that comparison.
*!
*! THE ASSEMBLY, and the only place the propagation terms appear:
*!     psi_theta_j = H_j^-1 g_j                         (thresholds AND slope)
*!     psi_rho_ij  = -U_rho^-1 ( u + U_theta [psi_theta_i ; psi_theta_j] )
*! U_theta has columns for the SLOPES b_i, b_j as well as the thresholds.
*! That slope part is new relative to v1's psi_rho and is the term most easily
*! dropped (@sec-condx Item 2): ctl = 1 drops it, and must fail Gamma.
*!
*! Module 1 (dwlsmv_cx_item_eval) and Module 2 (dwlsmv_cx_pair_blocks)
*! EVALUATE their blocks at the estimates in s; nothing here solves for an
*! estimate. Everything depends on an observation only through its pattern
*! (y, x) (@sec-condx (b)), so psi is gathered from per-cell tables, as v1's
*! dwlsmv_psi_matrix() gathers it.
*!
*! COLUMN ORDER IS THE INDEX MAP'S: thresholds, then slopes (pos_slope), then
*! correlations (pos_rho), spec 4.2.
*!
*! Below the locked boundary: a record, like v1's dwlsmv_psi_matrix(), not a
*! lock. `ctl` (internal testing only; production passes 0 or omits it):
*! 1 slope part of U_theta dropped, 2 all of U_theta dropped, 3 H_j's
*! (tau, b) cross entries dropped, 4 the quotient term of du/dtheta dropped.
*!
*! Requires dwlsmv_types, dwlsmv_statmap_q, dwlsmv_pos_tau, dwlsmv_pos_slope,
*! dwlsmv_pos_rho, dwlsmv_cx_check, dwlsmv_cx_item and dwlsmv_cx_pair (and
*! what they require).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_cx_psi_matrix(real matrix Y, real matrix X,
                                 struct dwlsmv_statmap scalar map,
                                 real colvector s, | real scalar ctl)
{
    real matrix    Psi, G, H, Hinv, tab
    real matrix    Tall
    real colvector P, gbar, n, th, thi, thj, flat, idx, toff
    real rowvector Ut, pi_, pj_
    real scalar    N, j, c, i, k, K, Ki, Kj, pos, ll, Ur, r, cc, xg, ri, rj, cl

    if (args() < 5) {
        cl = 0
    }
    else {
        cl = ctl
    }

    N = rows(Y)
    if (N < 1) {
        _error("dwlsmv_cx_psi_matrix: Y has no observations")
    }
    if (cols(Y) != map.p) {
        _error("dwlsmv_cx_psi_matrix: Y and the map disagree on the number of variables")
    }
    if (dwlsmv_statmap_q(map) != 1) {
        _error("dwlsmv_cx_psi_matrix: the map must be built for exactly one covariate")
    }
    if (rows(s) != map.D | cols(s) != 1) {
        _error("dwlsmv_cx_psi_matrix: s is not D x 1")
    }
    if (missing(s)) {
        _error("dwlsmv_cx_psi_matrix: s has missing values")
    }
    dwlsmv_cx_check_x(X, N)
    for (j = 1; j <= map.p; j++) {
        (void) dwlsmv_catcounts(Y[., j], map.K[j])
    }

    Psi = J(N, map.D, 0)
    // The item tables, stacked: item j occupies rows toff[j] + 1 .. toff[j] + 2 K_j
    // and columns 1 .. K_j.
    toff = J(map.p, 1, 0)
    for (j = 2; j <= map.p; j++) {
        toff[j] = toff[j - 1] + 2 * map.K[j - 1]
    }
    Tall = J(toff[map.p] + 2 * map.K[map.p], max(map.K), 0)

    // ---- items: psi_theta_j = H_j^-1 g_j, one 2K x K table per item -------
    for (j = 1; j <= map.p; j++) {
        K  = map.K[j]
        th = J(K, 1, .)
        for (c = 1; c <= K - 1; c++) {
            th[c] = s[dwlsmv_pos_tau(map, j, c)]
        }
        th[K] = s[dwlsmv_pos_slope(map, j, 1)]
        n = dwlsmv_cx_item_counts(Y[., j], X, K)
        dwlsmv_cx_item_eval(th, K, n, (cl == 3 ? 3 : 0), P, G, H, gbar, ll)
        Hinv = cholinv(H)
        if (missing(Hinv)) {
            errprintf("dwlsmv_cx_psi_matrix: item y%s's H_j is not positive definite at s.\n",
                      strofreal(j))
            _error("dwlsmv_cx_psi_matrix: an item's H_j is not positive definite")
        }
        // Unobserved cells are never gathered; zero them so that no missing
        // value can reach a product below.
        for (r = 1; r <= 2 * K; r++) {
            if (n[r] == 0) {
                G[r, .] = J(1, K, 0)
            }
        }
        if (missing(G)) {
            _error("dwlsmv_cx_psi_matrix: an observed item cell has a missing score")
        }
        Tall[|toff[j] + 1, 1 \ toff[j] + 2 * K, K|] = G * Hinv'

        idx = toff[j] :+ X :* K :+ Y[., j] :+ 1
        for (c = 1; c <= K - 1; c++) {
            pos         = dwlsmv_pos_tau(map, j, c)
            Psi[., pos] = Tall[idx, c]
        }
        pos         = dwlsmv_pos_slope(map, j, 1)
        Psi[., pos] = Tall[idx, K]
    }

    // ---- pairs: psi_rho with the complete propagation term ---------------
    for (i = 1; i <= map.p - 1; i++) {
        for (k = i + 1; k <= map.p; k++) {
            pos = dwlsmv_pos_rho(map, i, k)
            Ki  = map.K[i]
            Kj  = map.K[k]
            thi = J(Ki, 1, .)
            thj = J(Kj, 1, .)
            for (c = 1; c <= Ki - 1; c++) {
                thi[c] = s[dwlsmv_pos_tau(map, i, c)]
            }
            for (c = 1; c <= Kj - 1; c++) {
                thj[c] = s[dwlsmv_pos_tau(map, k, c)]
            }
            thi[Ki] = s[dwlsmv_pos_slope(map, i, 1)]
            thj[Kj] = s[dwlsmv_pos_slope(map, k, 1)]

            cc = 0
            if (cl == 1 | cl == 2 | cl == 4) {
                cc = cl
            }
            dwlsmv_cx_pair_blocks(Y[., i], Y[., k], X, thi, thj, s[pos], cc,
                                  Ur, Ut, tab, n, (i, k))
            if (!(Ur < 0)) {
                errprintf("dwlsmv_cx_psi_matrix: U_rho for y%s~~y%s is %s, not negative.\n",
                          strofreal(i), strofreal(k), strofreal(Ur))
                _error("dwlsmv_cx_psi_matrix: a pair's U_rho is not negative")
            }

            // psi_rho per observed cell (x, c_i, c_j); 0 elsewhere, never read.
            flat = J(rows(n), 1, 0)
            for (xg = 0; xg <= 1; xg++) {
                for (ri = 1; ri <= Ki; ri++) {
                    for (rj = 1; rj <= Kj; rj++) {
                        r = xg * Ki * Kj + (ri - 1) * Kj + rj
                        if (n[r] == 0) {
                            continue
                        }
                        pi_ = Tall[|toff[i] + xg * Ki + ri, 1 \ toff[i] + xg * Ki + ri, Ki|]
                        pj_ = Tall[|toff[k] + xg * Kj + rj, 1 \ toff[k] + xg * Kj + rj, Kj|]
                        flat[r] = -(tab[r, 2] + (pi_, pj_) * Ut') / Ur
                    }
                }
            }
            Psi[., pos] = flat[X :* (Ki * Kj) :+ Y[., i] :* Kj :+ Y[., k] :+ 1]
        }
    }

    return(Psi)
}

end
