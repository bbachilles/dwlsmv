*! dwlsmv -- exact second derivatives of sigma on design E's class (Result MS6)
*!
*! Executed under StataMP 16.1 on 2026-10-19: tests/unit/test_sens_hess.do green; suite 62/0/0 at four processors.
*!
*! @sec-mg-hessian, Result MS6 (derived and verified in R 2026-10-19;
*! E-031). Diagnostic only; no production path calls it.
*!
*!     void dwlsmv_sens_hdom(struct dwlsmv_parmap scalar par, real scalar G,
*!                           string scalar caller)
*!     real matrix dwlsmv_sens_hess(real colvector theta,
*!                                  struct dwlsmv_parmap  scalar par,
*!                                  struct dwlsmv_statmap vector maps)
*!
*! dwlsmv_sens_hess returns D x (npar * npar), dwlsmv_sens_ddelta's layout:
*! column (l - 1) npar + i is d2 sigma / d theta_i d theta_l, so block l is
*! dDelta/dtheta_l. Closed form, no differencing. Per group, every pair of
*! distinct free rows (a, b) adds MS6's entry to column (fidx_b - 1) npar +
*! fidx_a of the group's rows: ties sum over member rows (H_k = Z'H^v_k Z).
*!
*! THE DOMAIN (dwlsmv_sens_hdom), checked on the table, never on values:
*!   c1  no ptype 3 (beta) row, free or fixed;
*!   c2  one factor;
*!   c3  each theta element read by at most one row of any one group, so sigma
*!       is multilinear in theta.
*! Outside it the call ABORTS at rc 3498, naming every failing clause. There
*! is no finite-difference fallback: a silent fallback is how a result derived
*! for design E's class would become a general claim. Models with beta are not
*! covered (sigma has degree 2 in a beta row through (I - B)^-1).
*!
*! Load order: dwlsmv_types, dwlsmv_pos_tau, dwlsmv_nfac, dwlsmv_parvals,
*! dwlsmv_matrices, dwlsmv_mg_parts, dwlsmv_mg_sub.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_sens_hdom(struct dwlsmv_parmap scalar par, real scalar G,
                      string scalar caller)
{
    real scalar i, j, bad, gi, gj

    bad = 0
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 3) {
            errprintf("%s: row %f (%s) is a beta (ptype 3) row (clause c1)\n", caller, i, par.name[i])
            bad = 1
            break
        }
    }
    if (dwlsmv_nfac(par) != 1) {
        errprintf("%s: the model has %f factors; Result MS6 needs one (clause c2)\n", caller, dwlsmv_nfac(par))
        bad = 1
    }
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i] == 0) {
            continue
        }
        gi = (G == 1 ? 1 : par.grp[i])
        for (j = i + 1; j <= par.nall; j++) {
            if (par.free[j] == 0) {
                continue
            }
            gj = (G == 1 ? 1 : par.grp[j])
            if (gi == gj & par.fidx[i] == par.fidx[j]) {
                errprintf("%s: rows %f and %f of group %f share theta element %f (clause c3)\n",
                          caller, i, j, gi, par.fidx[i])
                bad = 1
            }
        }
    }
    if (bad) {
        errprintf("%s: Result MS6 (@sec-mg-hessian) covers one factor, no beta, one row per theta\n", caller)
        errprintf("element per group; there is no finite-difference fallback.\n")
        _error(caller + ": outside Result MS6's domain")
    }
}

// The slot of row (t, i) in the monomial om_j om_k lam_j lam_k psi, 0 if absent.
real scalar dwlsmv_sens_hslot(real scalar t, real scalar i, real scalar j, real scalar k)
{
    if (t == 6 & i == j) return(1)
    if (t == 6 & i == k) return(2)
    if (t == 1 & i == j) return(3)
    if (t == 1 & i == k) return(4)
    if (t == 2) return(5)
    return(0)
}

// d2 sigma_g / d v_a d v_b for two distinct rows of group g (MS6's table).
real colvector dwlsmv_sens_hpair(real scalar ta, real scalar ia, real scalar ca,
                                 real scalar tb, real scalar ib, real scalar cb,
                                 real colvector lam, real colvector om,
                                 real scalar psi, real scalar al,
                                 struct dwlsmv_statmap scalar map)
{
    real colvector d, v
    real scalar t, x, c, j, k, sa, sb, m

    if (ta > tb) {
        x = ta; ta = tb; tb = x
        x = ia; ia = ib; ib = x
        x = ca; ca = cb; cb = x
    }
    d = J(map.D, 1, 0)

    // threshold entries: sigma_jc = om_j tau_jc - om_j nu_j - om_j lam_j alpha
    if (ta == 4 & tb == 6 & ia == ib) {
        d[dwlsmv_pos_tau(map, ia, ca)] = 1
    }
    else if (ta == 5 & tb == 6 & ia == ib) {
        for (c = 1; c <= map.C[ia]; c++) d[dwlsmv_pos_tau(map, ia, c)] = -1
    }
    else if (ta == 1 & tb == 6 & ia == ib) {
        for (c = 1; c <= map.C[ia]; c++) d[dwlsmv_pos_tau(map, ia, c)] = -al
    }
    else if (ta == 6 & tb == 7) {
        for (c = 1; c <= map.C[ia]; c++) d[dwlsmv_pos_tau(map, ia, c)] = -lam[ia]
    }
    else if (ta == 1 & tb == 7) {
        for (c = 1; c <= map.C[ia]; c++) d[dwlsmv_pos_tau(map, ia, c)] = -om[ia]
    }

    // correlation entries: the monomial with the two differentiated factors removed
    for (t = map.T + 1; t <= map.D; t++) {
        j  = map.var_i[t]
        k  = map.var_j[t]
        sa = dwlsmv_sens_hslot(ta, ia, j, k)
        sb = dwlsmv_sens_hslot(tb, ib, j, k)
        if (sa == 0 | sb == 0 | sa == sb) {
            continue
        }
        v = (om[j] \ om[k] \ lam[j] \ lam[k] \ psi)
        x = 1
        for (m = 1; m <= 5; m++) {
            if (m != sa & m != sb) x = x * v[m]
        }
        d[t] = x
    }
    return(d)
}

real matrix dwlsmv_sens_hess(real colvector theta,
                             struct dwlsmv_parmap  scalar par,
                             struct dwlsmv_statmap vector maps)
{
    struct dwlsmv_parmap scalar sub
    real matrix out, Lam, Psi, B, A
    real colvector vals, vg, idx, nu, om, al, kap, lam, rws
    real scalar G, g, q, a, b, off, Dtot, col

    G = length(maps)
    if (G < 1) {
        _error("dwlsmv_sens_hess: no index map passed")
    }
    dwlsmv_sens_hdom(par, G, "dwlsmv_sens_hess")
    q = par.npar
    if (rows(theta) != q) {
        _error("dwlsmv_sens_hess: theta does not match par.npar")
    }
    vals = dwlsmv_parvals(theta, par)
    Dtot = 0
    for (g = 1; g <= G; g++) {
        Dtot = Dtot + maps[g].D
    }
    out = J(Dtot, q * q, 0)
    off = 0
    for (g = 1; g <= G; g++) {
        if (G == 1) {
            sub = par
            idx = (1::par.nall)
        }
        else {
            sub = dwlsmv_mg_sub(par, g, idx)
        }
        vg = vals[idx]
        dwlsmv_mg_parts(vg, sub, maps[g].p, 1, Lam, Psi, B, A, nu, om, al, kap)
        lam = Lam[., 1]
        rws = (off + 1)::(off + maps[g].D)
        for (a = 1; a <= sub.nall; a++) {
            if (sub.free[a] == 0) continue
            for (b = 1; b <= sub.nall; b++) {
                if (b == a | sub.free[b] == 0) continue
                col = (sub.fidx[b] - 1) * q + sub.fidx[a]
                out[rws, col] = out[rws, col] +
                    dwlsmv_sens_hpair(sub.ptype[a], sub.row[a], sub.col[a],
                                      sub.ptype[b], sub.row[b], sub.col[b],
                                      lam, om, Psi[1, 1], al[1], maps[g])
            }
        }
        off = off + maps[g].D
    }
    return(out)
}

end
