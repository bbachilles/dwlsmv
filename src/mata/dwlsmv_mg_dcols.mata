*! dwlsmv -- one group's Delta columns, one per table row (Module 6, G > 1)
*!
*! Executed under StataMP 18 on 2026-10-10. tests/unit/test_mg_model.do.
*!
*! D_g x nall_g: column r is d sigma_g / d v_r for group g's row r, free or
*! fixed (the Module 12 shape), from derivation @sec-mg-sigma-delta Result MS2:
*!
*!   tau_ic     omega_i at (i, c)
*!   nu_i       -omega_i on i's thresholds
*!   alpha_f    -omega_j L_jf on every threshold of j, L = Lambda A
*!   omega_i    tau_ic - nu_i - c_i on i's thresholds; omega_k V_ik on pairs (i, k)
*!   lambda_mf  -omega_m kappa_f on m's thresholds;  omega_j omega_k dV_jk on pairs
*!   beta_ab    -omega_j kappa_b L_ja on every threshold; omega_j omega_k dV_jk
*!   psi_ab     omega_j omega_k dV_jk on pairs
*!
*! tau, nu, alpha and omega are closed form, so a configural tau column is an
*! exact 1 (G6's identity is exact, not to rounding). dV for lambda, beta and
*! psi comes from dwlsmv_delta_cols() on the group's rows, the single-group
*! central difference that Lemma MS2.1 shows has no truncation error; omega
*! does not depend on those rows, so scaling after differencing is exact.
*!
*! Load order: dwlsmv_types, dwlsmv_pos_tau, dwlsmv_pos_rho, dwlsmv_nfac,
*! dwlsmv_matrices, dwlsmv_rystar, dwlsmv_corrvec, dwlsmv_delta_cols,
*! dwlsmv_mg_parts.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_mg_dcols(real colvector vals,
                            struct dwlsmv_parmap  scalar sub,
                            struct dwlsmv_statmap scalar map)
{
    real matrix Lam, Psi, B, A, L, V, Dm, Dv
    real colvector nu, om, al, kap, c, tau, sel
    real scalar i, j, k, cc, t, r, f, a, b, ns, nfac

    nfac = dwlsmv_nfac(sub)
    dwlsmv_mg_parts(vals, sub, map.p, nfac, Lam, Psi, B, A, nu, om, al, kap)
    c = Lam * kap
    L = Lam * A
    V = Lam * (A * Psi * A') * Lam'

    tau = J(map.T, 1, .)
    for (i = 1; i <= sub.nall; i++) {
        if (sub.ptype[i] == 4) {
            tau[dwlsmv_pos_tau(map, sub.row[i], sub.col[i])] = vals[i]
        }
    }

    Dm = J(map.D, sub.nall, 0)

    // dV for every lambda, beta and psi row at once, single-group differences
    ns = 0
    for (i = 1; i <= sub.nall; i++) {
        ns = ns + (sub.ptype[i] <= 3)
    }
    sel = J(ns, 1, 0)
    ns = 0
    for (i = 1; i <= sub.nall; i++) {
        if (sub.ptype[i] <= 3) {
            ns++
            sel[ns] = i
        }
    }
    if (ns > 0 & map.P > 0) {
        Dv = dwlsmv_delta_cols(vals, sub, map, sel)
        for (k = 1; k <= ns; k++) {
            for (t = map.T + 1; t <= map.D; t++) {
                Dm[t, sel[k]] = om[map.var_i[t]] * om[map.var_j[t]] * Dv[t, k]
            }
        }
    }

    for (i = 1; i <= sub.nall; i++) {
        r = sub.row[i]
        f = sub.col[i]
        if (sub.ptype[i] == 4) {
            Dm[dwlsmv_pos_tau(map, r, f), i] = om[r]
        }
        else if (sub.ptype[i] == 5) {
            for (cc = 1; cc <= map.C[r]; cc++) {
                Dm[dwlsmv_pos_tau(map, r, cc), i] = -om[r]
            }
        }
        else if (sub.ptype[i] == 7) {
            for (j = 1; j <= map.p; j++) {
                for (cc = 1; cc <= map.C[j]; cc++) {
                    Dm[dwlsmv_pos_tau(map, j, cc), i] = -om[j] * L[j, r]
                }
            }
        }
        else if (sub.ptype[i] == 6) {
            for (cc = 1; cc <= map.C[r]; cc++) {
                t = dwlsmv_pos_tau(map, r, cc)
                Dm[t, i] = tau[t] - nu[r] - c[r]
            }
            for (t = map.T + 1; t <= map.D; t++) {
                if (map.var_i[t] == r) {
                    Dm[t, i] = om[map.var_j[t]] * V[r, map.var_j[t]]
                }
                else if (map.var_j[t] == r) {
                    Dm[t, i] = om[map.var_i[t]] * V[r, map.var_i[t]]
                }
            }
        }
        else if (sub.ptype[i] == 1) {
            // lambda (indicator r, factor f): thresholds of r
            for (cc = 1; cc <= map.C[r]; cc++) {
                Dm[dwlsmv_pos_tau(map, r, cc), i] = -om[r] * kap[f]
            }
        }
        else if (sub.ptype[i] == 3) {
            // beta (to a = r, from b = f): every threshold
            a = r
            b = f
            for (j = 1; j <= map.p; j++) {
                for (cc = 1; cc <= map.C[j]; cc++) {
                    Dm[dwlsmv_pos_tau(map, j, cc), i] = -om[j] * kap[b] * L[j, a]
                }
            }
        }
    }
    return(Dm)
}

end
