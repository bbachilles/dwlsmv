*! dwlsmv -- S5's elementwise bound |E| on production Delta's central-difference entries
*!
*! Executed under StataMP 18 on 2026-10-16: tests/unit/test_sens_t3.do green at six processors.
*!
*! A helper of dwlsmv_sens_t3 (@sec-mg-sensitivity-t3, TS3: Ebar). Diagnostic
*! only; no production path calls it.
*!
*!     real matrix dwlsmv_sens_ebar(struct dwlsmv_fit scalar f, real scalar epsV)
*!
*! Returns D x npar: on the correlation rows of group g, in the free lambda,
*! psi and beta columns (ptype 1-3) that group g's table carries,
*!     (epsV |sigma_k| + u (|theta_j| + h_j) |Delta_kj|) / h_j,
*! h_j = 1e-5 max(1, |theta_j|), summed over tied rows; 0 elsewhere. It is
*! dwlsmv_sens()'s Eb loop, copied line for line (that routine keeps Eb
*! internal and is not edited, since frozen studies ran it).
*! tests/unit/test_sens_t3.do pins the copy: dwlsmv_sens()'s column 6
*! recomputed from this Eb equals the routine's own.
*!
*! ABORTS (rc 3498): epsV missing or negative.
*!
*! Load order: dwls_fit.mata (dwlsmv_fit_stack), dwls_model, dwls_delta.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_sens_ebar(struct dwlsmv_fit scalar f, real scalar epsV)
{
    struct dwlsmv_statmap vector maps
    real colvector s, w, th, sg, o, rws
    real matrix D, Eb
    real scalar N, q, Dn, G, g, i, j, hh, u

    if (missing(epsV) | epsV < 0) {
        _error("dwlsmv_sens_ebar: epsV is missing or negative")
    }
    dwlsmv_fit_stack(f.stage1, s, w, N, maps)
    th = f.theta
    q  = rows(th)
    Dn = rows(s)
    G  = length(maps)
    sg = dwls_model(th, f.par, maps)
    D  = dwls_delta(th, f.par, maps)
    u  = epsilon(1) / 2

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
    return(Eb)
}

end
