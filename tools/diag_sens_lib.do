*! dwlsmv -- diagnostic-only Mata for the check 6 sensitivity diagnostic (S5, mechanism D)
*!
*! Executed under StataMP 18 on 2026-10-14 by tools/diag_sens_step1.do (dgs_cdmask,
*! dgs_delta_acc, dgs_delta_an). dgs_fdwls and dgs_minimize were first executed 2026-10-14 by
*! the steps 2-3 smoke (docs/notes/sens-diag/run_diag.sh smoke, E, count 6 only).
*!
*! NOT PRODUCTION CODE and not a proposed production change (docs/decisions.md,
*! 2026-10-14: "the accurate derivative is a diagnostic control"). Nothing in
*! src/ calls it. Defines, after the library is loaded (tools/build_mlib.do):
*!
*!   dgs_delta_acc(th, par, maps, H)  production Delta with its central-
*!       difference entries (the correlation rows of free lambda and psi
*!       columns, @sec-mg-sensitivity S5) replaced by a central difference of
*!       dwls_model() at step H max(1, |th_j|). sigma's correlation rows are a
*!       polynomial of degree <= 2 in each lambda and linear in psi (Lemma
*!       MS2.1), so that difference has no truncation error at ANY step: a
*!       large H only shrinks the rounding, by H / 1e-5. Refuses a free beta
*!       (rc 3498), for which that exactness fails.
*!   dgs_delta_an(th, par, maps)  the same entries in closed form, for
*!       single-factor designs without beta only (E, B, D): group g's
*!       correlation (a, b) is om_a om_b psi_g lam_a lam_b, so d/d lam_m is
*!       om_a om_b psi_g (lam_b [m = a] + lam_a [m = b]) and d/d psi_g is
*!       om_a om_b lam_a lam_b; tied members sum. The independent check on
*!       dgs_delta_acc. Other entries are copied from dwls_delta().
*!   dgs_cdmask(par, maps)  Dtot x npar, 1 on the central-difference entries.
*!   dgs_fdwls(...)  dwlsmv_fdwls()'s arithmetic, with Delta chosen by mode:
*!       0 production dwls_delta(), 1 a fixed matrix Dfix, 2 dgs_delta_acc().
*!   dgs_minimize(...)  dwlsmv_minimize()'s settings (d2, nr, ptol 1e-10,
*!       vtol 1e-12, nrtol 1e-12) over dgs_fdwls(), with a trace level argument
*!       ("none", or "step" to log optimize()'s step-length adjustments;
*!       "value" logs none of them: found by the steps 2-3 smoke, 2026-10-14).
*!
*! Load: quietly do tools/build_mlib.do, then do tools/diag_sens_lib.do.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dgs_cdmask(struct dwlsmv_parmap scalar par, struct dwlsmv_statmap vector maps)
{
    real matrix Mk
    real scalar G, Dt, g, i, o, t

    G  = length(maps)
    Dt = 0
    for (g = 1; g <= G; g++) Dt = Dt + maps[g].D
    Mk = J(Dt, par.npar, 0)
    for (i = 1; i <= par.nall; i++) {
        if (!par.free[i] | par.ptype[i] > 3) continue
        g = (G == 1 ? 1 : par.grp[i])
        o = 0
        for (t = 1; t < g; t++) o = o + maps[t].D
        for (t = maps[g].T + 1; t <= maps[g].D; t++) Mk[o + t, par.fidx[i]] = 1
    }
    return(Mk)
}

real matrix dgs_delta_acc(real colvector th, struct dwlsmv_parmap scalar par,
                          struct dwlsmv_statmap vector maps, real scalar H)
{
    real matrix D, Mk
    real colvector tp, tm, d, rw
    real scalar j, i, hj

    for (i = 1; i <= par.nall; i++) {
        if (par.free[i] & par.ptype[i] == 3) {
            _error("dgs_delta_acc: a free beta; the large-step difference is not exact there")
        }
    }
    D  = dwls_delta(th, par, maps)
    Mk = dgs_cdmask(par, maps)
    for (j = 1; j <= rows(th); j++) {
        rw = selectindex(Mk[., j])
        if (rows(rw) == 0) continue
        hj = H * max((1, abs(th[j])))
        tp = th
        tm = th
        tp[j] = th[j] + hj
        tm[j] = th[j] - hj
        d = (dwls_model(tp, par, maps) - dwls_model(tm, par, maps)) / (tp[j] - tm[j])
        D[rw, j] = d[rw]
    }
    return(D)
}

real matrix dgs_delta_an(real colvector th, struct dwlsmv_parmap scalar par,
                         struct dwlsmv_statmap vector maps)
{
    real matrix D, Mk, lam, om
    real colvector vals, psi, o
    real scalar G, g, i, j, t, a, b, m, p

    G = length(maps)
    p = maps[1].p
    vals = dwlsmv_parvals(th, par)
    lam = J(G, p, 0)
    om  = J(G, p, 1)
    psi = J(G, 1, .)
    for (i = 1; i <= par.nall; i++) {
        g = (G == 1 ? 1 : par.grp[i])
        if (par.ptype[i] == 3) _error("dgs_delta_an: beta; single-factor designs without beta only")
        if (par.ptype[i] == 1) {
            if (par.col[i] != 1) _error("dgs_delta_an: more than one factor")
            lam[g, par.row[i]] = vals[i]
        }
        if (par.ptype[i] == 2) {
            if (par.row[i] != 1 | par.col[i] != 1) _error("dgs_delta_an: more than one factor")
            psi[g] = vals[i]
        }
        if (par.ptype[i] == 6) om[g, par.row[i]] = vals[i]
    }
    if (missing(psi)) _error("dgs_delta_an: a group has no psi row")
    o = J(G, 1, 0)
    for (g = 2; g <= G; g++) o[g] = o[g - 1] + maps[g - 1].D
    D  = dwls_delta(th, par, maps)
    Mk = dgs_cdmask(par, maps)
    D  = D :* (1 :- Mk)
    for (i = 1; i <= par.nall; i++) {
        if (!par.free[i] | (par.ptype[i] != 1 & par.ptype[i] != 2)) continue
        g = (G == 1 ? 1 : par.grp[i])
        j = par.fidx[i]
        m = par.row[i]
        for (t = maps[g].T + 1; t <= maps[g].D; t++) {
            a = maps[g].var_i[t]
            b = maps[g].var_j[t]
            if (par.ptype[i] == 1) {
                D[o[g] + t, j] = D[o[g] + t, j] + om[g, a] * om[g, b] * psi[g] * ((m == a) * lam[g, b] + (m == b) * lam[g, a])
            }
            else {
                D[o[g] + t, j] = D[o[g] + t, j] + om[g, a] * om[g, b] * lam[g, a] * lam[g, b]
            }
        }
    }
    return(D)
}

void dgs_fdwls(real scalar todo, real rowvector p,
               real colvector s, real colvector w,
               struct dwlsmv_parmap scalar par, struct dwlsmv_statmap vector maps,
               real scalar mode, real matrix Dfix, real scalar H,
               real scalar v, real rowvector g, real matrix Hm)
{
    real colvector theta, r
    real matrix Dm

    theta = p'
    r     = s - dwls_model(theta, par, maps)
    v     = quadcross(r, w, r)
    if (todo >= 1) {
        if (mode == 0)      Dm = dwls_delta(theta, par, maps)
        else if (mode == 1) Dm = Dfix
        else                Dm = dgs_delta_acc(theta, par, maps, H)
        g = -2 * quadcross(w :* r, Dm)
        if (todo == 2) Hm = 2 * quadcross(Dm, w, Dm)
    }
}

real scalar dgs_minimize(real colvector s, real colvector w,
                         struct dwlsmv_parmap scalar par, struct dwlsmv_statmap vector maps,
                         real colvector start, real scalar maxiter,
                         real scalar mode, real matrix Dfix, real scalar H,
                         string scalar trace,
                         real colvector theta, real scalar fmin, real scalar iter)
{
    transmorphic S
    real scalar rc, conv

    S = optimize_init()
    optimize_init_which(S, "min")
    optimize_init_evaluator(S, &dgs_fdwls())
    optimize_init_evaluatortype(S, "d2")
    optimize_init_technique(S, "nr")
    optimize_init_argument(S, 1, s)
    optimize_init_argument(S, 2, w)
    optimize_init_argument(S, 3, par)
    optimize_init_argument(S, 4, maps)
    optimize_init_argument(S, 5, mode)
    optimize_init_argument(S, 6, Dfix)
    optimize_init_argument(S, 7, H)
    optimize_init_params(S, start')
    optimize_init_conv_maxiter(S, maxiter)
    optimize_init_conv_ptol(S, 1e-10)
    optimize_init_conv_vtol(S, 1e-12)
    optimize_init_conv_nrtol(S, 1e-12)
    optimize_init_conv_warning(S, "off")
    optimize_init_tracelevel(S, trace)
    rc    = _optimize(S)
    theta = optimize_result_params(S)'
    fmin  = optimize_result_value(S)
    iter  = optimize_result_iterations(S)
    conv  = optimize_result_converged(S)
    if (rc != 0) conv = 0
    return(conv)
}

end
