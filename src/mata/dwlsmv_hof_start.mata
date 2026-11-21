*! dwlsmv -- Module 7: the optimizer's start for an item-anchored higher-order factor
*!
*! Decided by the owner 2026-11-20 (option A-R of
*! E-114, the start package; spec Module 7). Derivation:
*! @sec-so-hof (Result SO7).
*!
*! THE PROBLEM. A factor G with a FREE variance psi_GG, exactly ONE loading (its
*! marker, fixed at 1) and two or more free paths out (G -> F_k) has a ridge at
*! infinity: with c_k = b_k psi_GG held fixed and |psi_GG| -> infinity, every
*! implied correlation tends to a finite limit (Theta is derived, so the item's
*! own variance is never fitted), and F tends to a finite value. Module 3's
*! default start puts every b at 0 (plus the rank check's +-0.05 second point),
*! which is the ray's direction, and ULS and DWLS were measured sliding out along
*! it (B1's MD model: an item loading on G; E-114). With two or more
*! items on G, with no item on G, or with one path out, there is no such ridge.
*! With psi_GG fixed and the item's loading free the ridge also exists (SO7.3(b)),
*! for psi_GG >= 0; the owner's rule does not cover it (every measured fit there
*! converged), and it is not in scope here.
*!
*! THE RULE. In each group, for each factor of that shape, the free outgoing beta
*! rows start at 0.7 (the loading default) in the OPTIMIZER's start. start is
*! modified in place; the return value is the number of beta rows changed.
*!
*! WHAT IT DOES NOT TOUCH. par.value, and so Module 3's rank check, which has
*! already run at theta0 / theta2 by the time dwlsmv_fit_uls() builds start:
*! identification verdicts and messages are unchanged (dwlsmv_canon_status's
*! warning about the rank check's defaults). No other factor's start changes: on
*! the whole suite the rule fired on no fit (prototype R, proto_suite_summary.txt).
*! A tied beta is one theta element, so its members share the new start.
*!
*! Load order: dwlsmv_types.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_hof_start(struct dwlsmv_parmap scalar par, real colvector start)
{
    real scalar ng, nf, g, f, i, nlam, nb, psifree, nchg

    ng = max(par.grp)
    nf = 0
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 1) {
            nf = max((nf, par.col[i]))
        }
        if (par.ptype[i] == 2 | par.ptype[i] == 3) {
            nf = max((nf, par.row[i], par.col[i]))
        }
    }

    nchg = 0
    for (g = 1; g <= ng; g++) {
        for (f = 1; f <= nf; f++) {
            nlam    = 0
            nb      = 0
            psifree = 0
            for (i = 1; i <= par.nall; i++) {
                if (par.grp[i] == g) {
                    if (par.ptype[i] == 1 & par.col[i] == f) {
                        nlam++
                    }
                    if (par.ptype[i] == 3 & par.col[i] == f & par.free[i]) {
                        nb++
                    }
                    if (par.ptype[i] == 2 & par.row[i] == f & par.col[i] == f & par.free[i]) {
                        psifree = 1
                    }
                }
            }
            if (psifree & nlam == 1 & nb >= 2) {
                for (i = 1; i <= par.nall; i++) {
                    if (par.grp[i] == g & par.ptype[i] == 3 & par.col[i] == f & par.free[i]) {
                        start[par.fidx[i]] = 0.7
                        nchg++
                    }
                }
            }
        }
    }
    return(nchg)
}

end
