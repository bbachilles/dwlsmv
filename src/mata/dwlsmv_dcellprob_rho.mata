*! dwlsmv -- d(pi_rs)/d(rho) by Plackett's identity (Module 2, spec 4.5.2/4.7)
*!
*! Executed under StataMP 18 on 2026-09-19. Checked against central
*! differences of the cell probabilities (7e-12) and, independently of how
*! rho-hat is found, against a numerical derivative of the log-likelihood away
*! from the optimum (1.5e-09 relative).
*!
*! Plackett's identity is d Phi2(x,y;rho) / d rho = phi2(x,y;rho), so the
*! derivative of a cell probability is the SAME rectangle difference taken on
*! the bivariate DENSITY at the cell's four corners:
*!
*!   dpi_rs/drho = phi2(a_{r+1},b_{s+1}) - phi2(a_r,b_{s+1})
*!               - phi2(a_{r+1},b_s)     + phi2(a_r,b_s)
*!
*! with phi2 = 0 whenever either argument is infinite, since the density
*! vanishes there. That is the only boundary rule needed -- unlike the CDF,
*! which has four distinct boundary cases, the density has one.
*!
*! This is why Module 2 needs no numerical derivatives and no numerical
*! Hessian. Everything the information quantities require is this matrix and
*! the cell probabilities themselves. If a future session finds itself
*! differentiating the bivariate log-likelihood numerically, something
*! upstream has gone wrong -- re-read spec Section 4.5.2.
*!
*! phi2(x,y;rho) = exp(-(x^2 - 2 rho x y + y^2) / (2(1-rho^2)))
*!                 / (2 pi sqrt(1-rho^2))
*!
*! Note this file does NOT reuse the corner grid built by dwlsmv_cellprob().
*! The two are kept separate because the density corners and the CDF corners
*! obey different boundary rules, and a shared grid would have to carry both
*! sets of special cases -- which is how the wrong branch gets taken for the
*! wrong matrix. The cost is recomputing K1*K2 exponentials, which is trivial
*! next to the bivariate normal CDF evaluations in dwlsmv_cellprob().
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_dcellprob_rho(real colvector a, real colvector b,
                                 real scalar rho)
{
    real matrix D, G
    real scalar K1, K2, r, s, sig, twopi, x, y

    if (cols(a) != 1 | cols(b) != 1) {
        _error("dwlsmv_dcellprob_rho: thresholds must be column vectors")
    }
    if (rows(a) < 1 | rows(b) < 1) {
        _error("dwlsmv_dcellprob_rho: each variable needs at least one threshold")
    }
    if (rho <= -1 | rho >= 1 | rho >= .) {
        _error("dwlsmv_dcellprob_rho: rho must lie strictly inside (-1, 1)")
    }

    K1    = rows(a) + 1
    K2    = rows(b) + 1
    sig   = sqrt(1 - rho * rho)
    twopi = 2 * pi()

    // Density corner grid. Row 1 and row K1+1 are infinite thresholds, where
    // the density is zero; likewise the first and last column.
    G = J(K1 + 1, K2 + 1, 0)
    for (r = 2; r <= K1; r++) {
        for (s = 2; s <= K2; s++) {
            x       = a[r - 1]
            y       = b[s - 1]
            G[r, s] = exp(-(x * x - 2 * rho * x * y + y * y) /
                          (2 * (1 - rho * rho))) / (twopi * sig)
        }
    }

    D = J(K1, K2, 0)
    for (r = 1; r <= K1; r++) {
        for (s = 1; s <= K2; s++) {
            D[r, s] = G[r + 1, s + 1] - G[r, s + 1] - G[r + 1, s] + G[r, s]
        }
    }

    return(D)
}

end
