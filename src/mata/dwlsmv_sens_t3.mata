*! dwlsmv -- B_T3: a per-fit bound on e(T3)'s processor-count residue (check 6, G > 1)
*! f must be an UNPOLISHED fit; see dwlsmv_sens.mata.
*!
*! Executed under StataMP 18 on 2026-10-16: tests/unit/test_sens_t3.do green at six processors.
*!
*! @sec-mg-sensitivity-t3, Results TS1-TS7 (derived and verified in R
*! 2026-10-15). A diagnostic for the check 6 G > 1 investigation, built for
*! E-029 (approved by the user 2026-10-16, item 5): it
*! computes a bound from one fit and measured inputs, nothing else. It changes
*! no tolerance and no production path calls it.
*!
*!     real rowvector dwlsmv_sens_t3(struct dwlsmv_fit scalar f, real colvector step,
*!                                   real scalar epsV, real scalar rhoT,
*!                                   | real colvector ds, real colvector eps,
*!                                   real matrix dOm)
*!
*! f     a converged fit after dwls_inference() (T3, shift_a filled), with
*!       D_stat > 0 and |T3| >= 1e-2 (TS1 item 4: below that the posting term is
*!       not derived)
*! step  the last DWLS step, as for dwlsmv_sens()
*! epsV  S5's relative accuracy of one evaluation of the model-implied
*!       correlations (eps_5 in the derivation; dwlsmv_sens()'s argument of the
*!       same name). It enters Ebar, Theta's S5 part and beta_elem
*! rhoT  TS7's slack on the route terms (>= 0)
*! ds, eps, dOm   optional, all three or none: a measured change in s (D x 1),
*!       the measured RELATIVE change in w (D x 1, as dwlsmv_sens), and the
*!       measured change in Omega = blockdiag(Gamma_g / pi_g), every entry
*!       (D x D). None given is the Stage-1-fixed case: ds = eps = dOm = 0.
*!
*! Returns 1 x 8, columns 2-8 in check 6's units (divided by 1 + |T3|):
*!     1  T3          f.T3, as posted
*!     2  2 beta_1    route 1, Delta's noise in the traces, both runs (TS3)
*!     3  beta_2      route 2, theta-hat's shift, with Theta the a-priori S1-S5
*!                    envelope from dwlsmv_sens() (TS4)
*!     4  beta_3      route 3, Stage 1: <|G_Omega|, |dOm|> + |g_s|'|ds| (TS5)
*!     5  beta_eval   reduction order and the inverse (TS6), with ASSUMED
*!                    constant c_ch = 6(q + 1) for cholinv(), and SE6's Q4
*!                    2 gamma_(n-1) sum|terms| per quadcross() output
*!                    (n = D, or q for the E^-1 products); Q1's 2u
*!                    was refuted per run (E-030)
*!     6  beta_elem   elementwise roundings at different inputs (TS6)
*!     7  B_T3        (1 + rhoT)(c2 + c3 + c4) + c5 + c6 + r_T, r_T = 2u|T3|/(1+|T3|)
*!     8  B_eval      c5 + r_T: the bound at identical inputs (prereg E_T3)
*!
*! Every gradient is the subsection's (chunk mg-t3-setup: t3_grad, t3_evgrad,
*! t3_beval, t3_belem, t3_gth_abs), transcribed. W = 1/diag(Omega) in the
*! gradients and the traces (weighting (a)); A, c, J and Theta use the fit's w,
*! as dwlsmv_sens() does: the two are roundings of one quantity (TS1 item 3).
*! tests/unit/test_sens_t3.do compares every column with the subsection's own
*! R functions, run by tools/export_sens_t3_ref.R on this routine's inputs.
*!
*! ABORTS (rc 3498): f not converged; inference not run; D_stat = 0;
*! |T3| < 1e-2; epsV or rhoT missing or negative; one or two of ds, eps, dOm
*! given; any of them the wrong size or missing. dwlsmv_sens() adds its own
*! (step, iterations, c >= 1).
*!
*! Load order: dwls_fit.mata (dwlsmv_fit_stack), dwls_model, dwls_delta,
*! dwls_inference.mata (dwlsmv_omega), dwlsmv_traces, dwlsmv_sens,
*! dwlsmv_sens_ebar, dwlsmv_sens_jmat, dwlsmv_sens_ddelta.
*! dwlsmv_sens_t3_core (below) holds the computation; dwlsmv_sens_t3 calls it
*! with exact = 0, the finite-difference path, arithmetic unchanged since
*! 2026-10-16. dwlsmv_sens_t3_x (its own file) calls it with exact = 1: the
*! closed-form dDelta/dtheta and J of Result MS6 (@sec-mg-hessian), opt-in,
*! refused at rc 3498 outside that domain. Refactored 2026-10-19; executed
*! since under StataMP 16.1 (test_sens_t3 green on its references, suite 62/0/0).
*! The exact path also loads dwlsmv_sens_x, dwlsmv_sens_jmat_x and
*! dwlsmv_sens_ddelta_x.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real rowvector dwlsmv_sens_t3(struct dwlsmv_fit scalar f, real colvector step,
                              real scalar epsV, real scalar rhoT,
                              | real colvector ds, real colvector eps,
                              real matrix dOm)
{
    if (args() == 4) {
        return(dwlsmv_sens_t3_core(f, step, epsV, rhoT, 0, 4, J(0, 1, .), J(0, 1, .), J(0, 0, .)))
    }
    if (args() == 5) {
        return(dwlsmv_sens_t3_core(f, step, epsV, rhoT, 0, 5, ds, J(0, 1, .), J(0, 0, .)))
    }
    if (args() == 6) {
        return(dwlsmv_sens_t3_core(f, step, epsV, rhoT, 0, 6, ds, eps, J(0, 0, .)))
    }
    return(dwlsmv_sens_t3_core(f, step, epsV, rhoT, 0, 7, ds, eps, dOm))
}

// The routine itself. exact = 0 is the finite-difference path (the default,
// unchanged); exact = 1 takes dwlsmv_sens_x, dwlsmv_sens_ddelta_x and
// dwlsmv_sens_jmat_x (Result MS6). nin is the caller's args(): 7 with ds, eps
// and dOm given, 4 with none; 5 and 6 are refused below, where they always were.
real rowvector dwlsmv_sens_t3_core(struct dwlsmv_fit scalar f, real colvector step,
                                   real scalar epsV, real scalar rhoT,
                                   real scalar exact, real scalar nin,
                                   real colvector ds, real colvector eps,
                                   real matrix dOm)
{
    struct dwlsmv_statmap vector maps
    real colvector s, wf, w, th, sg, r, gh, gs, Gw, v, ghb, gtha, Th, y, dX2, dX3
    real rowvector tr
    real matrix sens, A, Om, Dl, WD, Pi, WDPi, U, UO, G1, G2, GD, IP, Z1, Z2, G1O, G2O, GTO, GO
    real matrix Eb, Jm, dD, P, Gm, B, C, EB, X2, X3, F_EB, G_B, G_C, G_Ei, G_P, G_Gm
    real matrix dP, dGm, dB, dC, dEB, aEi
    real rowvector out
    real scalar N, q, Dn, G, l, c, nu, a, kap, T3, nT, emax, b1, b2, b3, be, bel, u, rT, hD, hq
    real scalar a1, a2, t1, t2, Tst, ae, ke, k, gq, cch, dn

    if (f.converged != 1) {
        _error("dwlsmv_sens_t3: the fit did not converge")
    }
    if (rows(f.acov) == 0 | missing(f.T3) | missing(f.shift_a)) {
        _error("dwlsmv_sens_t3: dwls_inference() has not filled T3 and shift_a")
    }
    if (f.D_stat <= 0) {
        _error("dwlsmv_sens_t3: D_stat = 0, so T3 = T and TS3-TS7 do not apply")
    }
    T3 = f.T3
    if (abs(T3) < 1e-2) {
        errprintf("dwlsmv_sens_t3: |T3| = %g; the posting term is not derived below 1e-2 (TS1)\n", abs(T3))
        _error("dwlsmv_sens_t3: |T3| < 1e-2")
    }
    if (missing(epsV) | epsV < 0) {
        _error("dwlsmv_sens_t3: epsV is missing or negative")
    }
    if (missing(rhoT) | rhoT < 0) {
        _error("dwlsmv_sens_t3: rhoT is missing or negative")
    }
    if (nin == 5 | nin == 6) {
        _error("dwlsmv_sens_t3: give ds, eps and dOm together, or none")
    }

    dwlsmv_fit_stack(f.stage1, s, wf, N, maps)
    th = f.theta
    q  = rows(th)
    Dn = rows(s)
    G  = length(maps)
    if (nin == 7) {
        if (rows(ds) != Dn | rows(eps) != Dn | cols(ds) != 1 | cols(eps) != 1 |
            rows(dOm) != Dn | cols(dOm) != Dn) {
            errprintf("dwlsmv_sens_t3: ds has %f row(s), eps %f, dOm is %f x %f; D = %f\n",
                      rows(ds), rows(eps), rows(dOm), cols(dOm), Dn)
            _error("dwlsmv_sens_t3: ds, eps or dOm has the wrong size")
        }
        if (missing(ds) | missing(eps) | missing(dOm)) {
            _error("dwlsmv_sens_t3: ds, eps or dOm contains missing values")
        }
    }
    Dl = f.Delta
    if (rows(Dl) != Dn | cols(Dl) != q) {
        _error("dwlsmv_sens_t3: f.Delta does not match the fit's map and theta")
    }

    // Theta (TS4): the a-priori S1-S5 envelope, absolute units
    if (nin == 7) {
        if (exact) {
            sens = dwlsmv_sens_x(f, step, epsV, A, c, nu, ds, eps)
        }
        else {
            sens = dwlsmv_sens(f, step, epsV, A, c, nu, ds, eps)
        }
    }
    else {
        if (exact) {
            sens = dwlsmv_sens_x(f, step, epsV, A, c, nu)
        }
        else {
            sens = dwlsmv_sens(f, step, epsV, A, c, nu)
        }
        sens[., 5] = J(q, 1, 0)
    }
    Th = (1 :+ abs(th)) :* (sens[., 5] + sens[., 3] + sens[., 6])

    // The gradients (TS3, TS5), at the fit: t3_grad
    u   = epsilon(1) / 2
    nT  = N - G
    Om  = dwlsmv_omega(f.stage1)
    w   = 1 :/ diagonal(Om)
    sg  = dwls_model(th, f.par, maps)
    r   = s - sg
    tr  = dwlsmv_traces(Dl, Om)
    a   = f.shift_a
    kap = (f.T_raw - tr[1]) / (2 * tr[2])

    WD   = w :* Dl
    Pi   = cholinv(quadcross(Dl, w, Dl))
    WDPi = WD * Pi
    U    = diag(w) - WDPi * WD'
    UO   = U * Om
    G1   = -2 * UO * WDPi
    G2   = -4 * UO * UO * WDPi
    GD   = -a * (G1 + kap * G2)
    IP   = I(Dn) - Dl * WDPi'
    Z1   = IP * Om * IP'
    Z2   = IP * Om * UO * IP'
    G1O  = U - diag(w:^2 :* diagonal(Z1))
    G2O  = 2 * UO * U - diag(2 * w:^2 :* diagonal(Z2))
    GTO  = -diag(nT * w:^2 :* r:^2)
    GO   = a * (GTO - G1O - kap * G2O)
    gs   = 2 * a * nT * w :* r
    gh   = quadcross(Dl, w :* r)
    Gw   = a * (nT * r:^2 - diagonal(Z1) - 2 * kap * diagonal(Z2))

    // Route 1 (TS3): beta_1 = <|G_Delta|, Ebar>
    Eb = dwlsmv_sens_ebar(f, epsV)
    b1 = sum(abs(GD) :* Eb)

    // Route 2 (TS4): sum_l (2 a nT gbar_l + |v_l|) Theta_l + a nT Theta'|J|Theta
    if (exact) {
        dD = dwlsmv_sens_ddelta_x(f)
    }
    else {
        dD = dwlsmv_sens_ddelta(f)
    }
    v  = J(q, 1, .)
    for (l = 1; l <= q; l++) {
        v[l] = sum(GD :* dD[., ((l - 1) * q + 1)..(l * q)])
    }
    emax = (nin == 7 ? max(abs(eps)) : 0)
    ghb  = abs(gh) + quadcross(Eb, w :* abs(r)) + quadcross(abs(Dl), w :* abs(r)) * emax
    gtha = 2 * a * nT * ghb + abs(v)
    if (exact) {
        Jm = dwlsmv_sens_jmat_x(f)
    }
    else {
        Jm = dwlsmv_sens_jmat(f)
    }
    b2   = quadcross(gtha, Th) + a * nT * quadcross(Th, abs(Jm) * Th)

    // Route 3 (TS5)
    b3 = (nin == 7 ? sum(abs(GO) :* abs(dOm)) + quadcross(abs(gs), abs(ds)) : 0)

    // beta_eval (TS6): t3_evgrad and t3_beval, with SE6's Q4 eta(n, X) = 2 gamma_(n-1) X
    P    = quadcross(Dl, w, Dl)
    Gm   = quadcross(Om, WD)
    B    = quadcross(WD, Gm)
    C    = quadcross(Gm, w, Gm)
    EB   = Pi * B
    a1   = sum(w :* diagonal(Om))
    y    = quadcross(Om :* Om, w)
    a2   = sum(w :* y)
    X2   = Pi * C
    X3   = EB * EB
    t1   = a1 - sum(diagonal(EB))
    t2   = a2 - 2 * sum(diagonal(X2)) + sum(diagonal(X3))
    Tst  = nT * sum(w :* r:^2)
    ae   = sqrt((Dn - q) / t2)
    ke   = (Tst - t1) / (2 * t2)
    F_EB = ae * I(q) - 2 * ae * ke * EB'
    G_B  = Pi * F_EB
    G_C  = 2 * ae * ke * Pi
    G_Ei = F_EB * B' + 2 * ae * ke * C'
    G_P  = -Pi * G_Ei * Pi
    G_Gm = WD * G_B + (w :* Gm) * (G_C + G_C')
    k    = abs(ke)
    gq   = 2 * q * u / (1 - q * u)
    hD   = 2 * (Dn - 1) * u / (1 - (Dn - 1) * u)
    hq   = 2 * (q - 1) * u / (1 - (q - 1) * u)
    cch  = 6 * (q + 1)
    aEi  = abs(Pi)
    dP   = hD * quadcross(abs(Dl), w, abs(Dl)) + cch * u * sqrt(diagonal(P) * diagonal(P)')
    dGm  = hD * quadcross(abs(Om), abs(WD))
    dB   = hD * quadcross(abs(WD), abs(Gm))
    dC   = hD * quadcross(abs(Gm), w, abs(Gm))
    dEB  = hq * (aEi * abs(B))
    dX2  = diagonal(hq * (aEi * abs(C)))
    dX3  = diagonal(hq * (abs(EB) * abs(EB)))
    be   = sum(abs(G_P) :* dP) + sum(abs(G_Gm) :* dGm) + sum(abs(G_B) :* dB) + sum(abs(G_C) :* dC) +
           sum(abs(F_EB) :* dEB) + 2 * ae * k * sum(dX2) + ae * k * sum(dX3) +
           ae * k * (sum(w :* (hD * y)) + hD * sum(w :* abs(y))) +
           ae * hD * a1 + ae * nT * hD * sum(w :* r:^2) +
           ae * gq * sum(abs(diagonal(EB))) + 2 * ae * k * gq * sum(abs(diagonal(X2))) +
           ae * k * gq * sum(abs(diagonal(X3))) +
           ae * u * (abs(a1) + abs(sum(diagonal(EB)))) +
           ae * k * 2 * u * (abs(a2) + 2 * abs(sum(diagonal(X2))) + abs(sum(diagonal(X3)))) +
           4 * u * (abs(ae * Tst) + abs(ae * t1) + Dn - q)

    // beta_elem (TS6): closed-form Delta and W o Delta, w, sigma(theta-hat) in r
    bel = sum(abs(GD) :* (4 * u * abs(Dl))) + sum(abs(Gw) :* (2 * u * w)) +
          sum(abs(gs) :* (2 * epsV * abs(sg)))

    dn  = 1 + abs(T3)
    rT  = 2 * u * abs(T3) / dn
    out = (T3, 2 * b1 / dn, b2 / dn, b3 / dn, be / dn, bel / dn, ., .)
    out[7] = (1 + rhoT) * (out[2] + out[3] + out[4]) + out[5] + out[6] + rT
    out[8] = out[5] + rT
    return(out)
}

end
