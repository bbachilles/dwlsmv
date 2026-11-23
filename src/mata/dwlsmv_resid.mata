*! dwlsmv -- residual diagnostics: SRMR and the normalized residual matrix
*!
*! Drafted 2026-11-22 by a read-only agent; executed by the main session under StataMP 18 on
*! 2026-11-22: tests/unit/test_resid_diag.do passes (A-D; part B against the 99t lavaan
*! fixtures). Suite and four-processor diff not yet run.
*!
*!     real scalar dwlsmv_resid(struct dwlsmv_stage1 scalar st, real colvector rho,
*!                              real colvector repaired, real matrix E, real matrix Z)
*!
*! Derivation @sec-resid-diag; owner decisions E-117, package
*! sections 3b-3d (O1-O7). One group, no covariate (Stage 1 unconditional).
*!
*! Inputs
*!   st        Stage 1: s = (s_tau', r')', Gamma = N times the ACOV of s, N, bound.
*!             Pair k is at position map.T + k of s and Gamma (map.var_i[T+k] <
*!             map.var_j[T+k], i outer: the spec 4.2 column-major lower triangle).
*!   rho       P x 1, the model-implied polychorics in the same pair order. The
*!             caller supplies them; this function is model-free.
*!   repaired  P x 1, 0/1, dwlsmv_emptycell_pairs(Y, map) of the caller.
*! Outputs (E and Z are overwritten; st, rho and repaired are never modified)
*!   E   p x p symmetric: e_k = r_k - rho_k off the diagonal (SR1), diagonal 0.
*!   Z   p x p symmetric: z_k = e_k / sqrt(Gamma_kk / N) (SR2, N Stage 1's N, not
*!       N - 1), diagonal missing. MISSING at a repaired pair and at a pair held
*!       at the bound (SR4, owner O1 and O2): Gamma_kk is not a valid scale
*!       there. There is no 1e-5 special case (O7): lavaan divides by 1 where
*!       |e| < 1e-5, this does not, so z is the same function of e everywhere.
*! Returns  SRMR = sqrt(sum_k e_k^2 / (P + p)) over ALL P pairs, repaired and
*!       bound pairs included (SR1: the p diagonal residuals are zero by
*!       construction and enter the divisor only; thresholds drop out at G = 1).
*!
*! Aborts (rc 3498, detail by errprintf, short _error): a conditional Stage 1,
*! p < 2 (O5), an unbuilt map, rho or repaired or bound not P x 1, s or Gamma
*! shorter than T + P, a polychoric position whose stat_type is not 2, a bad
*! pair index, N missing or not positive, a missing r, rho or Gamma_kk (counted
*! with missing(), SR3's guard), Gamma_kk <= 0, repaired or bound not 0/1.
*! Stage 1 refuses all of the data cases already; this guard is for a future
*! change to it (SR3).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_resid(struct dwlsmv_stage1 scalar st, real colvector rho,
                         real colvector repaired, real matrix E, real matrix Z)
{
    real scalar p, P, T, N, k, pos, vi, vj
    real colvector r, gk, e, z
    real matrix Eo, Zo

    // Stage 2 does not take a conditional Stage 1 result yet (spec 5.1).
    dwlsmv_s2_guard_st(st, "dwlsmv_resid")

    if (length(st.map.p) != 1 | length(st.map.T) != 1 | length(st.map.P) != 1 | length(st.map.D) != 1) {
        errprintf("dwlsmv_resid: the Stage 1 index map is not built.\n")
        _error("dwlsmv_resid: map not built")
    }
    p = st.map.p
    P = st.map.P
    T = st.map.T
    if (missing(p) | missing(P) | missing(T)) {
        errprintf("dwlsmv_resid: the Stage 1 index map has a missing dimension.\n")
        _error("dwlsmv_resid: map not built")
    }
    if (p < 2) {
        errprintf("dwlsmv_resid: the residual diagnostics need at least two indicators (p = %g).\n", p)
        _error("dwlsmv_resid: needs at least two indicators")
    }
    N = st.N
    if (missing(N)) {
        errprintf("dwlsmv_resid: N must be positive, it is missing.\n")
        _error("dwlsmv_resid: N must be positive")
    }
    if (N <= 0) {
        errprintf("dwlsmv_resid: N must be positive, it is %g.\n", N)
        _error("dwlsmv_resid: N must be positive")
    }
    if (rows(rho) != P) {
        errprintf("dwlsmv_resid: rho has %g rows, the map has P = %g pairs.\n", rows(rho), P)
        _error("dwlsmv_resid: rho has the wrong length")
    }
    if (rows(repaired) != P) {
        errprintf("dwlsmv_resid: repaired has %g rows, the map has P = %g pairs.\n", rows(repaired), P)
        _error("dwlsmv_resid: repaired has the wrong length")
    }
    if (rows(st.bound) != P) {
        errprintf("dwlsmv_resid: bound has %g rows, the map has P = %g pairs.\n", rows(st.bound), P)
        _error("dwlsmv_resid: bound has the wrong length")
    }
    if (rows(st.s) < T + P | rows(st.Gamma) < T + P | cols(st.Gamma) < T + P | rows(st.map.stat_type) < T + P | rows(st.map.var_i) < T + P | rows(st.map.var_j) < T + P) {
        errprintf("dwlsmv_resid: s, Gamma or the map are shorter than T + P = %g.\n", T + P)
        _error("dwlsmv_resid: Stage 1 is too short")
    }
    if (any(st.map.stat_type[|T + 1 \ T + P|] :!= 2)) {
        errprintf("dwlsmv_resid: position T + k of the map is not a polychoric row for every pair k.\n")
        _error("dwlsmv_resid: map row is not a polychoric row")
    }
    if (any(repaired :!= 0 :& repaired :!= 1)) {
        errprintf("dwlsmv_resid: repaired must be 0 or 1 at every pair.\n")
        _error("dwlsmv_resid: repaired must be 0 or 1")
    }
    if (any(st.bound :!= 0 :& st.bound :!= 1)) {
        errprintf("dwlsmv_resid: bound must be 0 or 1 at every pair.\n")
        _error("dwlsmv_resid: bound must be 0 or 1")
    }

    r  = st.s[|T + 1 \ T + P|]
    gk = J(P, 1, .)
    for (k = 1; k <= P; k++) {
        pos   = T + k
        gk[k] = st.Gamma[pos, pos]
    }
    // missing(): sum(), max() and min() cannot see a missing value (a project rule)
    if (missing(r) > 0) {
        errprintf("dwlsmv_resid: %g missing polychoric correlation(s) in Stage 1.\n", missing(r))
        _error("dwlsmv_resid: missing polychoric")
    }
    if (missing(rho) > 0) {
        errprintf("dwlsmv_resid: %g missing implied correlation(s) in rho.\n", missing(rho))
        _error("dwlsmv_resid: missing implied correlation")
    }
    if (missing(gk) > 0) {
        errprintf("dwlsmv_resid: %g missing Gamma entries at the polychoric positions.\n", missing(gk))
        _error("dwlsmv_resid: missing Gamma entry")
    }
    if (any(gk :<= 0)) {
        errprintf("dwlsmv_resid: Gamma has a non-positive polychoric variance.\n")
        _error("dwlsmv_resid: non-positive polychoric variance")
    }

    e = r - rho
    z = e :/ sqrt(gk :/ N)

    Eo = J(p, p, 0)
    Zo = J(p, p, .)
    for (k = 1; k <= P; k++) {
        pos = T + k
        vi  = st.map.var_i[pos]
        vj  = st.map.var_j[pos]
        if (missing(vi) | missing(vj)) {
            errprintf("dwlsmv_resid: pair %g has a missing variable index.\n", k)
            _error("dwlsmv_resid: bad pair index")
        }
        if (vi < 1 | vj > p | vi >= vj) {
            errprintf("dwlsmv_resid: pair %g has indices (%g, %g) outside 1 <= i < j <= p = %g.\n", k, vi, vj, p)
            _error("dwlsmv_resid: bad pair index")
        }
        Eo[vi, vj] = e[k]
        Eo[vj, vi] = e[k]
        // O1 and O2: no normalized residual at a repaired pair or a pair at the bound
        if (repaired[k] == 0 & st.bound[k] == 0) {
            Zo[vi, vj] = z[k]
            Zo[vj, vi] = z[k]
        }
    }

    E = Eo
    Z = Zo
    return(sqrt(quadsum(e :^ 2) / (P + p)))
}

end
