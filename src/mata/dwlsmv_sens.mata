*! dwlsmv -- per-element sensitivity of theta-hat: (A), (C) and the (D) envelope (check 6, G > 1)
*!
*! Executed under StataMP 18 on 2026-10-14: tests/unit/test_sens.do (s1-s7,
*! e1-e8), 22 of 24 mutants red (the two survivors are in the test's header).
*!
*! @sec-mg-sensitivity, Results S1-S5 (derived and verified in R 2026-10-14).
*! A diagnostic for the check 6 G > 1 investigation (E-027,
*! 2026-10-14): it computes predictions from a fit, nothing else. It changes
*! no tolerance and no production path calls it.
*!
*!     real matrix dwlsmv_sens(struct dwlsmv_fit scalar f, real colvector step,
*!                             real scalar epsV,
*!                             real matrix A, real scalar c, real scalar nu,
*!                             | real colvector ds, real colvector eps)
*!
*! f     a converged fit (dwls_fit or dwlsmv_fit_from). s and w are rebuilt
*!       from f.stage1 by dwlsmv_fit_stack, so G > 1 is weighting (a).
*!       UNPOLISHED: these diagnostics model optimize()'s stop. A polished
*!       dwls_fit() runs without error but describes a stop that is not the
*!       reported theta-hat. Tests build f with tests/unit/sens_unpolished.mata.
*! step  the last DWLS step theta_K - theta_{K-1}, npar x 1, as measured.
*!       optimize() with maxiter m returns the iterate after m + 1 steps, so
*!       theta_{K-1} is maxiter K - 2 from the same start (@sec-mg-sensitivity).
*!       Required, because S3's bound is proportional to it.
*! epsV  S5's relative accuracy of one evaluation of V = Lambda A Psi A' Lambda'
*!       (>= 0). The derivation measured <= 0.89u in an R model and used 2u.
*!       Stata's value is diagnostic step 1's to measure, so it is an argument
*!       with no default.
*! ds, eps  optional, D x 1 each, both or neither: a measured change in s, and
*!       a measured RELATIVE change in w (eps_k = dw_k / w_k).
*!
*! Returns npar x 6, one row per element of theta, every column in check 6's
*! units (divided by 1 + |theta_i|):
*!     1  kappa_s  S1 per unit ||ds||_inf:   sum_k |A_ik|
*!     2  kappa_w  S1 per unit ||eps||_inf:  sum_k |A_ik| |r_k|
*!     3  S_C      S3(4) at u = 2^-53, u_q = 2^-104, K = f.iterations
*!     4  env_stop S4(b) with nu' = nu: sqrt(2 nu (J^-1 P J^-1)_ii), the
*!                 envelope of discrete stopping (zeroth order)
*!     5  S_A      |A (ds + r o eps)|_i when ds and eps are given, else missing
*!     6  S_D      S5(c): 2 phi_i / (1 - c), phi = |P^-1| |E|'W|r|, with |E| the
*!                 bound on the central-difference entries of Delta (the
*!                 correlation rows of every free lambda, psi and beta column,
*!                 h = 1e-5 max(1, |theta_j|), as dwlsmv_delta_cols). An
*!                 envelope, zeroth order in its trigger.
*! It also fills A = J^-1 Delta'W (npar x D), c = ||I - P^-1 J||_P and nu =
*! 2 ghat'P^-1 ghat at theta-hat (optimize()'s g H^-1 g' with H = 2P).
*!
*! J = P - R, with P = Delta'W Delta and R = sum_k w_k r_k grad^2 sigma_k. R is
*! a central difference of Delta(theta)'(w o r) with r held at theta-hat, at
*! step h_j = epsilon(1)^(1/3) max(1, |theta_j|), symmetrized. That costs 2 npar
*! calls to dwls_delta(). Delta is itself partly a central difference (S5), so R
*! carries an extra O(|e| / h) (the error budget's row for it).
*!
*! ABORTS (rc 3498):
*!   - f not converged, or f.iterations missing or below 1;
*!   - step, ds or eps of the wrong length or missing, or only one of ds and eps;
*!   - epsV missing or negative;
*!   - c >= 1. The Gauss-Newton map does not contract at theta-hat, so S2's
*!     series and the 1/(1 - c) of S3 and S5 do not exist, and J need not be
*!     positive definite.
*!
*! Load order: dwls_fit.mata (dwlsmv_fit_stack), dwls_model, dwls_delta.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_sens(struct dwlsmv_fit scalar f, real colvector step,
                        real scalar epsV,
                        real matrix A, real scalar c, real scalar nu,
                        | real colvector ds, real colvector eps)
{
    struct dwlsmv_statmap vector maps
    real colvector s, w, th, r, wr, gh, den, h, ev, f0, sg, o, rws
    real matrix D, P, R, Jm, L, Li, Pi, dP, Ji, Eb, out
    real scalar N, q, Dn, G, g, i, j, hh, u, uq, gq, K

    if (f.converged != 1) {
        _error("dwlsmv_sens: the fit did not converge")
    }
    K = f.iterations
    if (missing(K) | K < 1) {
        _error("dwlsmv_sens: f.iterations is missing or below 1")
    }
    if (missing(epsV) | epsV < 0) {
        _error("dwlsmv_sens: epsV is missing or negative")
    }
    dwlsmv_fit_stack(f.stage1, s, w, N, maps)
    th = f.theta
    q  = rows(th)
    Dn = rows(s)
    G  = length(maps)
    if (rows(step) != q | cols(step) != 1 | missing(step)) {
        errprintf("dwlsmv_sens: step has %f row(s); theta has %f\n", rows(step), q)
        _error("dwlsmv_sens: step is not a complete npar x 1 vector")
    }
    if (args() == 7) {
        _error("dwlsmv_sens: ds was given without eps")
    }
    if (args() == 8) {
        if (rows(ds) != Dn | rows(eps) != Dn | missing(ds) | missing(eps)) {
            errprintf("dwlsmv_sens: ds has %f row(s), eps %f; D = %f\n", rows(ds), rows(eps), Dn)
            _error("dwlsmv_sens: ds or eps is not a complete D x 1 vector")
        }
    }

    sg = dwls_model(th, f.par, maps)
    r  = s - sg
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
    R  = (R + R') / 2
    Jm = P - R

    // c = max |eig(P^-1 R)|, from the symmetric P^-1/2 R P^-1/2
    L = cholesky(P)
    if (missing(L)) {
        _error("dwlsmv_sens: P = Delta'W Delta is not positive definite")
    }
    Li = solvelower(L, I(q))
    ev = symeigenvalues(Li * R * Li')'
    c  = max(abs(ev))
    if (c >= 1) {
        errprintf("dwlsmv_sens: ||I - P^-1 J||_P = %g; the Gauss-Newton map does not contract at theta-hat\n", c)
        _error("dwlsmv_sens: contraction c >= 1, the premise of S2-S5 fails")
    }

    A   = cholsolve(Jm, (D :* w)')
    Ji  = cholinv(Jm)
    Pi  = abs(cholinv(P))
    gh  = quadcross(D, wr)
    nu  = 2 * quadcross(gh, cholsolve(P, gh))
    den = 1 :+ abs(th)

    // S3(4): u = 2^-53 and u_q = 2^-104, both exact powers of two
    u  = epsilon(1) / 2
    uq = epsilon(1)^2
    gq = Dn * uq / (1 - Dn * uq)
    f0 = 2 * gq * quadcross(abs(D), abs(wr))
    dP = 2 * u * abs(P) + 2 * gq * quadcross(abs(D), w, abs(D))

    // S5: |E| on the central-difference entries, group by group
    o = J(G, 1, 0)
    for (g = 2; g <= G; g++) {
        o[g] = o[g - 1] + maps[g - 1].D
    }
    Eb = J(Dn, q, 0)
    for (i = 1; i <= f.par.nall; i++) {
        if (f.par.free[i] == 0 | f.par.ptype[i] > 3) {
            continue
        }
        g = (G == 1 ? 1 : f.par.grp[i])
        if (maps[g].P == 0) {
            continue
        }
        j   = f.par.fidx[i]
        hh  = 1e-5 * max((1, abs(th[j])))
        rws = (o[g] + maps[g].T + 1)::(o[g] + maps[g].D)
        Eb[rws, j] = Eb[rws, j] + (epsV :* abs(sg[rws]) + u * (abs(th[j]) + hh) :* abs(D[rws, j])) / hh
    }

    out = J(q, 6, .)
    out[., 1] = rowsum(abs(A)) :/ den
    out[., 2] = (abs(A) * abs(r)) :/ den
    out[., 3] = (K * (Pi * ((2 * u * abs(P) + dP) * abs(step))) + (Pi * f0) / (1 - c)) :/ den
    out[., 4] = sqrt(2 * nu * diagonal(Ji * P * Ji)) :/ den
    if (args() == 8) {
        out[., 5] = abs(A * (ds + r :* eps)) :/ den
    }
    out[., 6] = (2 * (Pi * quadcross(Eb, abs(wr))) / (1 - c)) :/ den
    return(out)
}

end
