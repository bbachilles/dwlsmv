*! dwlsmv -- the Gamma-hat diagnostic (Module 10, spec Section 5)
*!
*! Executed under StataMP 18 on 2026-09-26. tests/unit/test_gammadiag.do: case D
*! healthy on every row, ratios equal to svd() within 4 D eps; the aliased
*! design (d_alias 1, r* 15) and a twelve-pattern design (r* 11), each ratio
*! clear of zero while the next singular value is structural; the three
*! budgets bracketed at 2x-10x and 0.1x-0.5x; five contract violations, rc 3498.
*!
*!     real matrix dwlsmv_gammadiag(struct dwlsmv_stage1 scalar st,
*!                                  real matrix Y)
*!
*! Spec Section 5, Module 10 and E-060 and E-071. Returns
*! 7 x 4, one row per condition, always in this order, in the [flag, value,
*! i, j] layout of dwlsmv_admissible():
*!
*!     row  condition            flag                  value           i, j
*!      1   asymmetry            a > 4 D eps           a               worst (i, j)
*!      2   not PSD              lmin < -4 D eps lmax  lmin / lmax     0, 0
*!      3   too few patterns     n_pat - 1 < D         n_pat           D, 0
*!      4   threshold aliasing   d_alias > 0           d_alias         first relation
*!      5   centring             z > 1e-6              z               argmax k, 0
*!      6   ratio, raw           .n                    sigma_r*/sigma_1   r*, 0
*!      7   ratio, scale-free    .n                    same, W^1/2 G W^1/2   r*, 0
*!
*! THE FLAG HAS THREE STATES (spec Section 5, decided 2026-09-26): 0/1
*! evaluated against the row's criterion; `.` NOT evaluated; .n evaluated
*! with no flag defined in v1. Rows 6-7 are .n in v1 by decision: no Gamma-hat
*! cutoff (decided 2026-10-05, spec 4.5.6; the Monte Carlo found no inference
*! failure the ratio identifies). Severe near-aliasing is still open.
*! Row 7, the scale-free ratio, is THE Gamma-hat ratio row: any future flag
*! attaches there. Row 6, raw, stays a reported value for continuity with
*! lavaan comparisons; it does not separate the small-count failure from
*! healthy data (E-071, 2026-09-26).
*!
*! THREE KINDS OF ROW, NEVER MIXED. Rows 3-4 are structural deficiency,
*! exact and detected from their cause with no cutoff. Rows 6-7 are the
*! near-singularity ratio, taken WITHIN the structurally supported rank
*! r* = min(D - d_alias, n_pat - 1): past r*, Gamma-hat's singular values are
*! zero by Theorem 1 or 3 of @sec-gamma-rank, and that deficiency belongs to
*! rows 3-4, so the ratio never reports it a second time. Rows 1, 2 and 5
*! carry budgets derived in @sec-gamma-sym and @sec-gamma-budgets. There is
*! no numerical-rank row for Gamma-hat in v1.
*!
*! ASYMMETRY STOPS IT. When row 1 flags, rows 2-7 are (., ., 0, 0), not
*! evaluated (decided 2026-09-26): symeigenvalues() reads one triangle, so a
*! spectrum of an asymmetric matrix would look valid. Production's Gamma-hat
*! is exactly symmetric (quadcross() mirrors), so row 1 is a backstop for an
*! external Gamma. The measure and budget are dwlsmv_asym()'s, the one
*! definition the guards in dwlsmv_traces() and dwlsmv_acov() also use.
*!
*! SINGULAR VALUES AS |EIGENVALUES|. For a symmetric matrix they are the
*! same numbers, and row 2 needs the eigenvalues anyway, so Gamma-hat is
*! decomposed once. test_gammadiag.do checks both ratios against svd().
*!
*! Folds in what the 2026-09-25 shape row called dwlsmv_symdiag() (the eigen
*! extremes and the PSD budget), superseded 2026-09-26.
*!
*! CONTRACT (rc 3498): Y must be the matrix Stage 1 received -- rows(Y) = N,
*! cols(Y) = p, codes valid for map.K (dwlsmv_catcounts()'s validation, via
*! dwlsmv_npatterns()) -- and Gamma must be D x D and psibar D x 1. A
*! mismatch means the rows describe some other data. dwlsmv_asym() adds its
*! own: no missing value, strictly positive diagonal.
*!
*! ONE-EMPTY-CELL PAIRS (O6, decided 2026-10-06; spec Module 10). Under
*! construction (a) a repaired pair's psi-bar_rho is S-bar(rho-hat+)/I_rr != 0
*! at O(1) by construction (@sec-o6, chunk o6-vii), so its rho entry is
*! EXCLUDED from row 5's z: the row stays a bug report for every other entry.
*! Membership comes from Y through dwlsmv_emptycell_pairs(), the same function
*! Stage 1 repaired by, so the two cannot disagree. The pair's threshold
*! entries stay in z. Repaired pairs are an estimator action, not a
*! condition, so they get no row; Module 11 lists them.
*!
*! A FUNCTION, NOT A FIELD, as dwlsmv_admissible() is. dwls_sample_stats()
*! does not call it; Module 11 does.
*!
*! Requires dwlsmv_types, dwlsmv_asym, dwlsmv_catcounts, dwlsmv_npatterns,
*! dwlsmv_thralias, dwlsmv_dalias, dwlsmv_crosstab, dwlsmv_emptycell1, dwlsmv_emptycell_pairs.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_gammadiag(struct dwlsmv_stage1 scalar st, real matrix Y)
{
    real matrix A, G, Gs, R
    real rowvector as
    real colvector ev, sv, s, zk, rep
    real scalar D, bud, npat, da, z, k, rs

    D = st.map.D
    if (rows(Y) != st.N | cols(Y) != st.map.p) {
        errprintf("dwlsmv_gammadiag: Y is %f x %f but Stage 1 was built on %f x %f.\n",
                  rows(Y), cols(Y), st.N, st.map.p)
        errprintf("  Pass the matrix dwls_sample_stats() received, after listwise deletion.\n")
        _error("dwlsmv_gammadiag: Y is not the matrix Stage 1 was built on")
    }
    if (rows(st.Gamma) != D | cols(st.Gamma) != D | rows(st.psibar) != D | cols(st.psibar) != 1) {
        errprintf("dwlsmv_gammadiag: Gamma is %f x %f and psibar %f x %f; the map says D = %f.\n",
                  rows(st.Gamma), cols(st.Gamma), rows(st.psibar), cols(st.psibar), D)
        _error("dwlsmv_gammadiag: Gamma or psibar does not match the map")
    }
    // Mata's max() skips missing values, so a missing psibar component would
    // read as a passed centring row (adversarial review 2026-09-26, D1).
    if (missing(st.psibar)) {
        _error("dwlsmv_gammadiag: psibar contains missing values")
    }

    A = (J(7, 2, .), J(7, 2, 0))
    G = st.Gamma

    // ---- 1. asymmetry; a flag stops everything else ---------------------------
    as = dwlsmv_asym(G)
    A[1, .] = (as[1] > as[2], as[1], as[3], as[4])
    if (A[1, 1]) {
        return(A)
    }

    // ---- 2. not PSD -----------------------------------------------------------
    bud = 4 * D * epsilon(1)
    ev  = symeigenvalues(G)'
    A[2, .] = (min(ev) < -bud * max(ev), min(ev) / max(ev), 0, 0)

    // ---- 3-4. the structural rows, from Y --------------------------------------
    npat = dwlsmv_npatterns(Y, st.map.K)
    A[3, .] = (npat - 1 < D, npat, D, 0)
    R  = dwlsmv_thralias(Y, st.map.K)
    da = dwlsmv_dalias(R)
    if (rows(R)) {
        A[4, .] = (da > 0, da, R[1, 1], R[1, 2])
    }
    else {
        A[4, .] = (0, 0, 0, 0)
    }

    // ---- 5. centring: the first k attaining z ------------------------------------
    s  = sqrt(diagonal(G))
    zk = abs(st.psibar) :/ s
    // one-empty-cell pairs: their rho entries are excluded (header)
    rep = dwlsmv_emptycell_pairs(Y, st.map)
    for (k = 1; k <= st.map.P; k++) {
        if (rep[k]) {
            zk[st.map.T + k] = 0
        }
    }
    z  = max(zk)
    for (k = 1; zk[k] != z; k++) {
    }
    A[5, .] = (z > 1e-6, z, k, 0)

    // ---- 6-7. the ratios, within r* --------------------------------------------
    rs = min((D - da, npat - 1))
    if (rs >= 1) {
        sv = sort(abs(ev), -1)
        A[6, .] = (.n, sv[rs] / sv[1], rs, 0)
        Gs = G :/ (s * s')
        sv = sort(abs(symeigenvalues(Gs)'), -1)
        A[7, .] = (.n, sv[rs] / sv[1], rs, 0)
    }
    else {
        A[6, .] = (.n, ., 0, 0)
        A[7, .] = (.n, ., 0, 0)
    }
    return(A)
}

end
