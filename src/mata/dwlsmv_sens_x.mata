*! dwlsmv -- dwlsmv_sens with R in closed form (Result MS6): the exact variant
*! f must be an UNPOLISHED fit; see dwlsmv_sens.mata.
*!
*! Executed under StataMP 16.1 on 2026-10-19: tests/unit/test_sens_hess.do green; suite 62/0/0 at four processors.
*!
*!     real matrix dwlsmv_sens_x(struct dwlsmv_fit scalar f, real colvector step,
*!                               real scalar epsV,
*!                               real matrix A, real scalar c, real scalar nu,
*!                               | real colvector ds, real colvector eps)
*!
*! dwlsmv_sens (@sec-mg-sensitivity), copied line for line except one: R in
*! J = P - R is sum_k w_k r_k H_k with H from dwlsmv_sens_hess (@sec-mg-hessian),
*! not a central difference of Delta'(w o r). Arguments, columns, the filled A,
*! c and nu, and every refusal are dwlsmv_sens's (read its header); messages
*! name dwlsmv_sens. It also refuses at rc 3498 outside Result MS6's domain
*! (one factor, no beta, one row per theta element per group), with no
*! fallback. dwlsmv_sens itself is unchanged: the frozen grid and T3 studies
*! ran it. Diagnostic only; no production path calls it.
*!
*! Load order: dwls_fit.mata (dwlsmv_fit_stack), dwls_model, dwls_delta,
*! dwlsmv_sens_hess.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_sens_x(struct dwlsmv_fit scalar f, real colvector step,
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

    R  = colshape(quadcross(dwlsmv_sens_hess(th, f.par, maps), wr)', q)'
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
