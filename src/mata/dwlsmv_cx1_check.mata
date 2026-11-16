*! dwlsmv -- conditional covariate, one 0/1 x, one group: Stage 2 helpers (Mod. 5, 6)
*!
*! Executed 2026-11-13 by the main session (Stata 18 MP); drafted by a read-only agent.
*!
*!     real scalar    dwlsmv_cx1_check(struct dwlsmv_parmap scalar par,
*!                                     struct dwlsmv_statmap vector maps,
*!                                     string scalar who)
*!     void           dwlsmv_cx1_parts(real colvector vals, struct dwlsmv_parmap scalar par,
*!                                     struct dwlsmv_statmap scalar map,
*!                                     real matrix Lam, real matrix B,
*!                                     real colvector Gam, real colvector Kap)
*!     real colvector dwlsmv_cx1_slopes(real colvector vals, struct dwlsmv_parmap scalar par,
*!                                      struct dwlsmv_statmap scalar map)
*!
*! dwlsmv_cx1_check() replaces the dwlsmv_s2_guard() call at the top of dwls_model() and
*! dwls_delta() (and, since 2026-11-14, dwlsmv_admissible()) ONLY. It returns q, 0 or 1:
*!   - several maps: dwlsmv_s2_guard() decides (a conditional map is refused at rc 3498,
*!     as before), and gamma or kappa rows are refused, since several groups with a
*!     covariate are not specified;
*!   - one map: q is read from the map (dwlsmv_statmap_q); q > 1 is refused (several
*!     covariates are not specified, spec 4.8); gamma or kappa rows on a map with q = 0
*!     are refused, because dwlsmv_matrices() would otherwise ignore them silently.
*! A map that was never filled has no q and is left for the caller's own checks, as in
*! dwlsmv_s2_guard().
*!
*! dwlsmv_cx1_parts() collects the matrices the slopes need. Lam and B come from
*! dwlsmv_matrices() unchanged; Gam is nfac x 1 (gamma rows, ptype 8, row = factor,
*! col = covariate = 1) and Kap is p x 1 (kappa rows, ptype 9, row = indicator,
*! col = covariate = 1). Fixed rows enter through vals, as everywhere.
*!
*! dwlsmv_cx1_slopes() is b = Lambda (I-B)^-1 Gamma_x + kappa (derivation @sec-condx,
*! item 4). Var(y*|x) = 1, so nothing else enters.
*!
*! Requires dwlsmv_types, dwlsmv_statmap_q, dwlsmv_s2_guard, dwlsmv_nfac, dwlsmv_matrices.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_cx1_check(struct dwlsmv_parmap scalar par,
                             struct dwlsmv_statmap vector maps,
                             string scalar who)
{
    real scalar q, built, ncx

    if (length(maps) < 1) {
        return(0)
    }

    ncx = sum(par.ptype :== 8) + sum(par.ptype :== 9)

    if (length(maps) > 1) {
        // several groups: a conditional map stays refused (not specified yet)
        dwlsmv_s2_guard(maps, who)
        if (ncx > 0) {
            errprintf("%s: gamma and kappa rows need one group; several groups with a covariate are not specified yet.\n", who)
            _error("gamma and kappa rows with several groups are not specified yet")
        }
        return(0)
    }

    q = 0
    built = (length(maps[1].D) == 1) + (length(maps[1].T) == 1) +
            (length(maps[1].P) == 1) + (length(maps[1].p) == 1)
    if (built == 4) {
        if (missing(maps[1].D) + missing(maps[1].T) +
            missing(maps[1].P) + missing(maps[1].p) == 0) {
            q = dwlsmv_statmap_q(maps[1])
        }
    }

    if (q > 1) {
        errprintf("%s: the index map has q = %g covariates; Stage 2 accepts one 0/1 covariate (spec 4.8).\n", who, q)
        _error("several covariates are not specified yet")
    }

    if (q == 0 & ncx > 0) {
        errprintf("%s: the table has gamma or kappa rows but the index map has no covariate (q = 0).\n", who)
        _error("gamma and kappa rows need a conditional index map")
    }

    return(q)
}

void dwlsmv_cx1_parts(real colvector vals,
                      struct dwlsmv_parmap scalar par,
                      struct dwlsmv_statmap scalar map,
                      real matrix Lam,
                      real matrix B,
                      real colvector Gam,
                      real colvector Kap)
{
    real matrix Psi
    real scalar nf, i, r, c

    nf = dwlsmv_nfac(par)
    dwlsmv_matrices(vals, par, map.p, nf, Lam, Psi, B)

    Gam = J(nf, 1, 0)
    Kap = J(map.p, 1, 0)

    for (i = 1; i <= par.nall; i++) {
        r = par.row[i]
        c = par.col[i]

        if (par.ptype[i] == 8) {
            if (c != 1 | r < 1 | r > nf | r != floor(r)) {
                errprintf("dwlsmv_cx1_parts: gamma row %f (%s) is outside (factor 1..%f, covariate 1).\n",
                          i, par.name[i], nf)
                _error("a gamma row names a factor or covariate that does not exist")
            }
            Gam[r] = vals[i]
        }
        else if (par.ptype[i] == 9) {
            if (c != 1 | r < 1 | r > map.p | r != floor(r)) {
                errprintf("dwlsmv_cx1_parts: kappa row %f (%s) is outside (indicator 1..%f, covariate 1).\n",
                          i, par.name[i], map.p)
                _error("a kappa row names an indicator or covariate that does not exist")
            }
            Kap[r] = vals[i]
        }
    }
}

real colvector dwlsmv_cx1_slopes(real colvector vals,
                                 struct dwlsmv_parmap scalar par,
                                 struct dwlsmv_statmap scalar map)
{
    real matrix Lam, B, A
    real colvector Gam, Kap

    dwlsmv_cx1_parts(vals, par, map, Lam, B, Gam, Kap)

    A = luinv(I(rows(B)) - B)
    if (missing(A)) {
        errprintf("dwlsmv_cx1_slopes: (I - B) is singular, so the slopes do not exist.\n")
        _error("dwlsmv_cx1_slopes: (I - B) is singular")
    }

    return(Lam * A * Gam + Kap)
}

end
