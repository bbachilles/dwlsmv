*! dwlsmv -- fit indices, option A: the c1 baseline in closed form (Module 12)
*!
*! Written 2026-10-21 from derivation @sec-fit-indices, Result FI2.
*!
*!     real rowvector dwlsmv_fitidx_base(struct dwlsmv_stage1 vector st)
*!
*! st: Stage 1, one element per group (G = 1: a 1-vector). Returns
*! (T3_B, D_B, T_B, a_B). The c1 baseline frees every threshold in every
*! group and fixes every polychoric at 0, so no fit is run:
*!   T_B   = (N - G) sum_g pi_g sum_k rho_gk^2 / Gamma_g[rho_k, rho_k]
*!   D_B   = the number of polychorics over the groups, G p(p-1)/2
*!   Tr(M_B) = D_B,  Tr(M_B^2) = sum_g ||R_g||_F^2, R_g the correlation
*!           matrix of Gamma_g's polychoric block (free of pi_g)
*!   a_B   = sqrt(D_B / Tr(M_B^2)),  T3_B = a_B T_B + D_B (1 - a_B).
*! D_B = 0 (no pair) returns (., 0, ., .). A missing or non-positive
*! polychoric variance in Gamma, or a missing polychoric, aborts rc 3498.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real rowvector dwlsmv_fitidx_base(struct dwlsmv_stage1 vector st)
{
    real scalar G, g, N, F, tr2, DB, a, TB
    real colvector idx, r, d
    real matrix R

    G = length(st)
    if (G < 1) {
        _error("dwlsmv_fitidx_base: no group")
    }
    N = 0
    for (g = 1; g <= G; g++) {
        N = N + st[g].N
    }
    F   = 0
    tr2 = 0
    DB  = 0
    for (g = 1; g <= G; g++) {
        idx = selectindex(st[g].map.stat_type :== 2)
        if (rows(idx) == 0) {
            continue
        }
        r = st[g].s[idx]
        d = diagonal(st[g].Gamma)[idx]
        if (missing(r) | missing(d)) {
            errprintf("dwlsmv_fitidx_base: group %g has a missing polychoric or Gamma entry.\n", g)
            _error("dwlsmv_fitidx_base: missing input")
        }
        if (any(d :<= 0)) {
            errprintf("dwlsmv_fitidx_base: group %g has a non-positive polychoric variance.\n", g)
            _error("dwlsmv_fitidx_base: non-positive variance")
        }
        F   = F + (st[g].N / N) * quadsum(r :^ 2 :/ d)
        R   = st[g].Gamma[idx, idx] :/ sqrt(d * d')
        tr2 = tr2 + quadsum(vec(R :^ 2))
        DB  = DB + rows(idx)
    }
    if (DB == 0) {
        return((., 0, ., .))
    }
    TB = (N - G) * F
    a  = sqrt(DB / tr2)
    return((a * TB + DB * (1 - a), DB, TB, a))
}

end
