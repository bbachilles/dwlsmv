*! dwlsmv -- the correlation block of sigma(theta), in index-map order (Mod. 5)
*!
*! Executed under StataMP 18 on 2026-09-22. Exercised by test_model.do and, once per finite difference, by every Delta column in test_delta.do.
*!
*! Returns P x 1: the model-implied polychoric correlations, ordered exactly as
*! the correlation block of s is (spec Section 4.2), so entry r corresponds to
*! element T + r of sigma. The subtraction of map.T is the only place that
*! offset appears; dwls_model() and dwlsmv_delta_cols() both write the block
*! into rows (T+1)..D without doing their own arithmetic.
*!
*! Spec Section 4.6: the model-implied polychoric for pair (j,k) is simply the
*! (j,k) entry of R_y*. Thresholds determine each indicator's marginal category
*! probabilities but do not enter the population polychoric at all -- verified
*! in the derivation (Gate 2, Finding 2) by cutting the same latent model with
*! two entirely different threshold schemes and recovering the same
*! correlations. Nothing here reads a threshold, and that absence is the
*! implementation of the finding.
*!
*! Factored out of dwls_model() because Module 6 differentiates THIS and not
*! sigma as a whole: the threshold rows of a structural parameter's Delta
*! column are exactly zero and are never computed.
*!
*! src/mata/dwlsmv_types.mata must be loaded first.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real colvector dwlsmv_corrvec(real colvector vals,
                              struct dwlsmv_parmap scalar par,
                              struct dwlsmv_statmap scalar map)
{
    real matrix Lam, Psi, B, R
    real colvector r
    real scalar i, j

    if (map.P == 0) {
        return(J(0, 1, .))
    }

    dwlsmv_matrices(vals, par, map.p, dwlsmv_nfac(par), Lam, Psi, B)
    R = dwlsmv_rystar(Lam, Psi, B)

    r = J(map.P, 1, .)
    for (i = 1; i <= map.p - 1; i++) {
        for (j = i + 1; j <= map.p; j++) {
            r[dwlsmv_pos_rho(map, i, j) - map.T] = R[i, j]
        }
    }

    return(r)
}

end
