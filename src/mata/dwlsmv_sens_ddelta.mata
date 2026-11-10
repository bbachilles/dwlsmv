*! dwlsmv -- dDelta/dtheta_l from production Delta, for TS4's trace gradient v
*!
*! Executed under StataMP 18 on 2026-10-16: tests/unit/test_sens_t3.do green at six processors.
*!
*! A helper of dwlsmv_sens_t3 (@sec-mg-sensitivity-t3, TS4). Diagnostic only;
*! no production path calls it.
*!
*!     real matrix dwlsmv_sens_ddelta(struct dwlsmv_fit scalar f)
*!
*! Returns D x (npar * npar): block l (columns (l-1) npar + 1 .. l npar) is
*!     (Delta_c(theta + h_l e_l) - Delta_c(theta - h_l e_l)) / (2 h_l),
*! h_l = 1e-4 max(1, |theta_l|), with Delta_c = dwls_delta(). Delta is at most
*! linear in any single parameter (Lemma MS2.1), so the difference has no
*! truncation error; it carries O(|E| / h) noise from Delta_c's own central
*! differences, about 1e-7 relative at this h, which TS4 places inside rho_T.
*!
*! Load order: dwls_delta.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_sens_ddelta(struct dwlsmv_fit scalar f)
{
    struct dwlsmv_statmap vector maps
    real colvector th, h
    real matrix out, Dp, Dm
    real scalar q, l, g

    maps = dwlsmv_statmap(length(f.stage1), 1)
    for (g = 1; g <= length(f.stage1); g++) {
        maps[g] = f.stage1[g].map
    }
    th  = f.theta
    q   = rows(th)
    out = J(0, 0, .)
    for (l = 1; l <= q; l++) {
        h    = J(q, 1, 0)
        h[l] = 1e-4 * max((1, abs(th[l])))
        Dp   = dwls_delta(th + h, f.par, maps)
        Dm   = dwls_delta(th - h, f.par, maps)
        if (l == 1) {
            out = J(rows(Dp), q * q, .)
        }
        out[., ((l - 1) * q + 1)..(l * q)] = (Dp - Dm) / (2 * h[l])
    }
    return(out)
}

end
