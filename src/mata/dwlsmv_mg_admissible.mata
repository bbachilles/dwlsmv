*! dwlsmv -- is one group's solution admissible? (Module 7 at G > 1, spec O15)
*!
*!     real matrix dwlsmv_mg_admissible(real colvector theta,
*!                                      struct dwlsmv_parmap  scalar par,
*!                                      real scalar g,
*!                                      struct dwlsmv_statmap scalar map)
*!
*! O15 (decided 2026-10-31, user): each group is assessed on its own scale,
*! with omega entering, by dwlsmv_admissible() itself -- the same six rows --
*! plus O12's row 7, omega_gj <= 0 (decided 2026-10-09: value the smallest
*! omega_gj, i its indicator; an indicator with no omega row reads 1). The
*! result is 7 x 4. Nothing is refused and nothing is stored; the caller turns
*! it into notes (dwlsmv_admissadvice()).
*!
*! THE TRANSFORMATION, explicit. Group g's implied statistics are (MS1)
*!
*!     sigma_gjc = omega_j (tau_jc - nu_j - (Lambda kappa)_j)
*!     rho_g,ij  = omega_i omega_j (Lambda Phi Lambda')_ij
*!
*! so they are the single-group statistics of a table whose loadings are
*! omega_j Lambda_j. and whose thresholds are sigma_gjc. dwlsmv_admissible()
*! is called on that table (ptypes 1-4, every row fixed at its group value),
*! which gives
*!
*!   row 1  Theta*_gj = 1 - omega_j^2 (Lambda Phi Lambda')_jj, O17's Theta_g;
*!          1 - (Lambda Phi Lambda')_jj, the untransformed reading, misses or
*!          invents a Heywood case whenever omega_j != 1
*!   row 4  R*_g, positive definite or not (diagonal 1, as at G = 1)
*!   row 6  the order of sigma_gj., not of tau: a negative omega_j reverses it
*!   rows 2, 3, 5  Phi_g and (I - B_g), which omega does not touch
*!
*! At omega = 1, nu = alpha = 0 the table is group g's own and the result is
*! the single-group one exactly. When (I - B_g) is singular kappa does not
*! exist, so the thresholds are left untransformed and row 6 is set missing;
*! rows 1-5 are then dwlsmv_admissible()'s own (row 5 flagged, 1-4 missing).
*!
*! Load order: dwlsmv_types, dwlsmv_parvals, dwlsmv_mg_sub, dwlsmv_matrices,
*! dwlsmv_admissible.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_mg_admissible(real colvector theta,
                                 struct dwlsmv_parmap  scalar par,
                                 real scalar g,
                                 struct dwlsmv_statmap scalar map)
{
    struct dwlsmv_parmap scalar sub, t
    real matrix Lam, Psi, B, Ai, A
    real colvector idx, v, nu, om, al, c, keep
    real scalar i, nf, sing, w, wi

    sub = dwlsmv_mg_sub(par, g, idx)
    v   = dwlsmv_parvals(theta, par)[idx]
    nf  = dwlsmv_nfac(sub)

    // ---- group g's nu, omega, alpha (absent rows: 0, 1, 0; spec 4.8 G4) --
    nu = J(map.p, 1, 0)
    om = J(map.p, 1, 1)
    al = J(nf, 1, 0)
    for (i = 1; i <= sub.nall; i++) {
        if (sub.ptype[i] == 5) {
            nu[sub.row[i]] = v[i]
        }
        else if (sub.ptype[i] == 6) {
            om[sub.row[i]] = v[i]
        }
        else if (sub.ptype[i] == 7) {
            al[sub.row[i]] = v[i]
        }
    }

    // ---- kappa = (I - B)^-1 alpha; c = Lambda kappa ----------------------
    // Reported, not aborted on (dwlsmv_mg_parts() aborts): a singular
    // (I - B) is dwlsmv_admissible()'s row 5.
    dwlsmv_matrices(v, sub, map.p, nf, Lam, Psi, B)
    Ai   = luinv(I(nf) - B)
    sing = missing(Ai)
    c    = sing ? J(map.p, 1, 0) : Lam * (Ai * al)

    // ---- the transformed table: ptypes 1-4, every row fixed --------------
    keep = (sub.ptype :<= 4)
    t.npar  = 0
    t.nall  = sum(keep)
    t.name  = select(sub.name,  keep)
    t.ptype = select(sub.ptype, keep)
    t.row   = select(sub.row,   keep)
    t.col   = select(sub.col,   keep)
    t.grp   = select(sub.grp,   keep)
    t.free  = J(t.nall, 1, 0)
    t.fidx  = J(t.nall, 1, 0)
    t.value = select(v, keep)
    for (i = 1; i <= t.nall; i++) {
        if (t.ptype[i] == 1) {
            t.value[i] = om[t.row[i]] * t.value[i]
        }
        else if (t.ptype[i] == 4 & !sing) {
            t.value[i] = om[t.row[i]] * (t.value[i] - nu[t.row[i]] - c[t.row[i]])
        }
    }

    A = dwlsmv_admissible(J(0, 1, .), t, map)
    if (sing) {
        A[6, .] = (., ., 0, 0)
    }

    // ---- 7. omega_gj <= 0 (spec 4.8 O12) --------------------------------
    w  = om[1]
    wi = 1
    for (i = 2; i <= map.p; i++) {
        if (om[i] < w) {
            w  = om[i]
            wi = i
        }
    }
    return(A \ (w <= 0, w, wi, 0))
}

end
