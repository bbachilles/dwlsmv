*! dwlsmv -- J = P - R with R in closed form (Result MS6), the exact variant of dwlsmv_sens_jmat
*!
*! Executed under StataMP 16.1 on 2026-10-19: tests/unit/test_sens_hess.do green; suite 62/0/0 at four processors.
*!
*! A helper of dwlsmv_sens_t3_x and dwlsmv_sens_x (@sec-mg-sensitivity-t3, TS4;
*! @sec-mg-hessian). Diagnostic only; no production path calls it.
*! dwlsmv_sens_jmat is unchanged and stays callable.
*!
*!     real matrix dwlsmv_sens_jmat_x(struct dwlsmv_fit scalar f)
*!
*! Returns npar x npar: J = P - R, P = Delta'W Delta with production's Delta
*! (dwls_delta, as dwlsmv_sens_jmat) and the fit's w, and R = sum_k w_k r_k H_k
*! with H from dwlsmv_sens_hess at theta-hat: R[i, l] = sum_k (w o r)_k
*! d2 sigma_k / d theta_i d theta_l. Only R changes against dwlsmv_sens_jmat;
*! P keeps production Delta's own rounding. R is symmetric in exact arithmetic
*! and in this evaluation (H's two columns are built from the same products),
*! so it is not symmetrized. Refuses at rc 3498 outside Result MS6's domain.
*!
*! Load order: dwls_fit.mata (dwlsmv_fit_stack), dwls_model, dwls_delta,
*! dwlsmv_sens_hess.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_sens_jmat_x(struct dwlsmv_fit scalar f)
{
    struct dwlsmv_statmap vector maps
    real colvector s, w, th, r, wr
    real matrix D, P, R, H
    real scalar N, q

    dwlsmv_fit_stack(f.stage1, s, w, N, maps)
    th = f.theta
    q  = rows(th)
    H  = dwlsmv_sens_hess(th, f.par, maps)
    r  = s - dwls_model(th, f.par, maps)
    wr = w :* r
    D  = dwls_delta(th, f.par, maps)
    P  = quadcross(D, w, D)
    R  = colshape(quadcross(H, wr)', q)'
    return(P - R)
}

end
