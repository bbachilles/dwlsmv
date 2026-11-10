*! dwlsmv -- S1's J = P - R at theta-hat, as dwlsmv_sens() forms it
*!
*! Executed under StataMP 18 on 2026-10-16: tests/unit/test_sens_t3.do green at six processors.
*!
*! A helper of dwlsmv_sens_t3 (@sec-mg-sensitivity-t3, TS4: the a n_T
*! Theta'|J|Theta term of beta_2). Diagnostic only; no production path calls it.
*!
*!     real matrix dwlsmv_sens_jmat(struct dwlsmv_fit scalar f)
*!
*! Returns npar x npar: J = P - R, P = Delta'W Delta with the fit's w, and R
*! the symmetrized central difference of Delta(theta)'(w o r) with r held at
*! theta-hat, step h_j = epsilon(1)^(1/3) max(1, |theta_j|). It is
*! dwlsmv_sens()'s J, copied line for line (that routine keeps J internal and
*! is not edited, since frozen studies ran it). tests/unit/test_sens_t3.do
*! pins the copy: cholsolve(J, (Delta :* w)') equals dwlsmv_sens()'s A.
*!
*! Load order: dwls_fit.mata (dwlsmv_fit_stack), dwls_model, dwls_delta.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_sens_jmat(struct dwlsmv_fit scalar f)
{
    struct dwlsmv_statmap vector maps
    real colvector s, w, th, r, wr, h
    real matrix D, P, R
    real scalar N, q, j

    dwlsmv_fit_stack(f.stage1, s, w, N, maps)
    th = f.theta
    q  = rows(th)
    r  = s - dwls_model(th, f.par, maps)
    wr = w :* r
    D  = dwls_delta(th, f.par, maps)
    P  = quadcross(D, w, D)

    R = J(q, q, .)
    for (j = 1; j <= q; j++) {
        h = J(q, 1, 0)
        h[j] = epsilon(1)^(1/3) * max((1, abs(th[j])))
        R[., j] = (quadcross(dwls_delta(th + h, f.par, maps), wr)
                   - quadcross(dwls_delta(th - h, f.par, maps), wr)) / (2 * h[j])
    }
    R = (R + R') / 2
    return(P - R)
}

end
