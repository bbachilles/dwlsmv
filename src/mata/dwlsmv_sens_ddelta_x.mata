*! dwlsmv -- dDelta/dtheta_l in closed form (Result MS6), the exact variant of dwlsmv_sens_ddelta
*!
*! Executed under StataMP 16.1 on 2026-10-19: tests/unit/test_sens_hess.do green; suite 62/0/0 at four processors.
*!
*! A helper of dwlsmv_sens_t3_x (@sec-mg-sensitivity-t3, TS4; @sec-mg-hessian).
*! Diagnostic only; no production path calls it. dwlsmv_sens_ddelta (central
*! differences of production Delta) is unchanged and stays callable: the frozen
*! grid and T3 studies ran it.
*!
*!     real matrix dwlsmv_sens_ddelta_x(struct dwlsmv_fit scalar f)
*!
*! Returns D x (npar * npar) in dwlsmv_sens_ddelta's layout: block l is
*! dDelta/dtheta_l = H[., ., l], from dwlsmv_sens_hess at f.theta. No noise from
*! differencing. Refuses at rc 3498 outside Result MS6's domain (one factor, no
*! beta, one row per theta element per group), with no fallback.
*!
*! Load order: dwlsmv_sens_hess.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_sens_ddelta_x(struct dwlsmv_fit scalar f)
{
    struct dwlsmv_statmap vector maps
    real scalar g

    maps = dwlsmv_statmap(length(f.stage1), 1)
    for (g = 1; g <= length(f.stage1); g++) {
        maps[g] = f.stage1[g].map
    }
    return(dwlsmv_sens_hess(f.theta, f.par, maps))
}

end
