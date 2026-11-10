*! dwlsmv -- one group's model-implied statistics sigma_g (Module 5, G > 1)
*!
*! Executed under StataMP 18 on 2026-10-10. tests/unit/test_mg_model.do.
*!
*! Derivation @sec-mg-sigma-delta Result MS1, in group g's own index map order:
*!
*!     sigma_gjc = omega_j (tau_jc - nu_j - (Lambda kappa)_j)
*!     rho_g,ij  = omega_i omega_j (Lambda Phi Lambda')_ij      (i != j)
*!
*! The correlation part is dwlsmv_corrvec() on the group's rows (which ignores
*! ptypes 5-7), scaled pair by pair. omega multiplies at every level (Q-2);
*! absent nu, omega, alpha rows read 0, 1, 0 (G4). At omega = 1, nu = alpha = 0
*! this is the single-group sigma exactly (MS4), but dwls_model() does not
*! route G = 1 here: the single-group path stays as it was.
*!
*! Load order: dwlsmv_types, dwlsmv_pos_tau, dwlsmv_pos_rho, dwlsmv_nfac,
*! dwlsmv_matrices, dwlsmv_rystar, dwlsmv_corrvec, dwlsmv_mg_parts.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real colvector dwlsmv_mg_sigma(real colvector vals,
                               struct dwlsmv_parmap  scalar sub,
                               struct dwlsmv_statmap scalar map)
{
    real matrix Lam, Psi, B, A
    real colvector nu, om, al, kap, c, sig
    real scalar i, k, j

    dwlsmv_mg_parts(vals, sub, map.p, dwlsmv_nfac(sub), Lam, Psi, B, A, nu, om, al, kap)
    c = Lam * kap

    sig = J(map.D, 1, .)
    for (i = 1; i <= sub.nall; i++) {
        if (sub.ptype[i] == 4) {
            j = sub.row[i]
            sig[dwlsmv_pos_tau(map, j, sub.col[i])] = om[j] * (vals[i] - nu[j] - c[j])
        }
    }

    if (map.P > 0) {
        sig[(map.T + 1)::map.D] = dwlsmv_corrvec(vals, sub, map)
        for (k = map.T + 1; k <= map.D; k++) {
            sig[k] = om[map.var_i[k]] * om[map.var_j[k]] * sig[k]
        }
    }

    if (missing(sig)) {
        errprintf("dwls_model: %f of a group's %f statistics were never written.\n",
                  missing(sig), map.D)
        errprintf("Its rows carry %f threshold row(s); its index map has %f.\n",
                  sum(sub.ptype :== 4), map.T)
        _error("dwls_model: a group's sigma has entries no parameter row filled")
    }
    return(sig)
}

end
