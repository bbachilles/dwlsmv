*! dwlsmv -- the N x D observation-level influence matrix (Module 4, spec 4.5)
*!
*! Executed under StataMP 18 on 2026-09-20. Its threshold gather is verified at
*! machine precision by tests/unit/test_gamma_tau.do and its correlation gather by
*! test_gamma.do's exact moment identities.
*!
*! Gathers psi by table lookup. Nothing here evaluates a normal CDF, a density
*! or a derivative: every value already exists, in the Stage 1 staircase
*! (spec 4.5.1) or in the per-pair psi_rho table (spec 4.5.2), and this routine
*! only says which row of which table each observation reads.
*!
*! Read the efficiency claim precisely, because an earlier version of the project's notes
*! overstated it and had to be corrected. What the lookup tables buy is that no
*! BIVARIATE NORMAL EVALUATION is ever repeated per observation -- those are
*! per cell, K_i*K_j of them per pair. This routine is still O(N*D), and so is
*! the crossproduct that consumes it, because the covariance between two
*! statistics resting on different variable pairs needs the joint distribution
*! over as many as four variables and there is no table for that. "Lookup
*! table" is not "no O(N) work".
*!
*! WHY THIS IS SEPARATE FROM dwls_gamma(). Gamma is the crossproduct of this
*! matrix and nothing else, so splitting them costs one function call and buys
*! the two tests that matter most: the derivation file's Level 4 pointwise
*! expansion needs psi itself rather than its second moment, and Case A's
*! closed-form psi_rho (spec 4.5.3) is an assertion about a column of this
*! matrix. Neither is expressible against Gamma alone. It is below the locked
*! API boundary -- a record, like Module 2's internals, not a lock.
*!
*! COLUMN ORDER IS THE INDEX MAP'S, always. Column s holds the influence
*! function of s[s], so the thresholds come first, variable-major and ascending
*! within variable, then the pairs in column-major lower-triangle order. The
*! map is the single definition of that and this routine asks it rather than
*! reproducing it.
*!
*! Requires dwlsmv_types, dwlsmv_koff, dwlsmv_pos_tau and
*! dwlsmv_psi_rho_table.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_psi_matrix(real matrix Y,
                              struct dwlsmv_statmap scalar map,
                              struct dwlsmv_pairfit vector pairs,
                              real matrix psi_tau)
{
    real matrix Psi, tbl, psi_a, psi_b
    real colvector koff, idx, flat
    real scalar N, j, c, k, pos, vi, vj, K1, K2, r, s

    N = rows(Y)

    if (cols(Y) != map.p) {
        _error("dwlsmv_psi_matrix: Y and the map disagree on the number of variables")
    }
    if (N < 1) {
        _error("dwlsmv_psi_matrix: Y has no observations")
    }
    if (length(pairs) != map.P) {
        errprintf("dwlsmv_psi_matrix: %s pair fits supplied, the map declares %s pairs.\n",
                  strofreal(length(pairs)), strofreal(map.P))
        _error("dwlsmv_psi_matrix: wrong number of pair fits")
    }

    koff = dwlsmv_koff(map)

    if (rows(psi_tau) != koff[map.p] + map.K[map.p] | cols(psi_tau) != map.T) {
        errprintf("dwlsmv_psi_matrix: psi_tau is %s x %s, expected %s x %s.\n",
                  strofreal(rows(psi_tau)), strofreal(cols(psi_tau)),
                  strofreal(koff[map.p] + map.K[map.p]), strofreal(map.T))
        errprintf("  psi_tau is the STACKED per-variable lookup tables of" +
                  " dwlsmv_psi_tau_stack(),\n")
        errprintf("  not the gathered N x T influence matrix. See that file's" +
                  " header.\n")
        _error("dwlsmv_psi_matrix: psi_tau has the wrong shape")
    }

    Psi = J(N, map.D, 0)

    // ---- threshold block (spec 4.5.1) --------------------------------------
    // Each variable writes only its own columns. The zero padding elsewhere in
    // psi_tau is never read, which is why it can be zero rather than missing.
    for (j = 1; j <= map.p; j++) {
        idx = koff[j] :+ Y[., j] :+ 1
        if (min(idx) < koff[j] + 1 | max(idx) > koff[j] + map.K[j]) {
            _error("dwlsmv_psi_matrix: a category code is outside 0 .. K-1")
        }
        for (c = 1; c <= map.C[j]; c++) {
            pos         = dwlsmv_pos_tau(map, j, c)
            Psi[., pos] = psi_tau[idx, pos]
        }
    }

    // ---- correlation block (spec 4.5.2) ------------------------------------
    for (k = 1; k <= map.P; k++) {
        pos = map.T + k
        vi  = map.var_i[pos]
        vj  = map.var_j[pos]

        // The pair vector must be in the map's order. dwls_polychoric() does
        // not set i and j -- the locked signature never tells it which
        // variables it was handed -- so the assembler fills them, and checking
        // them here is what makes a mis-ordered pair vector a loud failure
        // instead of a Gamma whose correlation block is silently permuted.
        if (pairs[k].i != vi | pairs[k].j != vj) {
            errprintf("dwlsmv_psi_matrix: pair %s is (%s,%s) but the map's %s is (%s,%s).\n",
                      strofreal(k), strofreal(pairs[k].i), strofreal(pairs[k].j),
                      map.label[pos], strofreal(vi), strofreal(vj))
            errprintf("  The pair vector must be in the index map's order," +
                      " and each fit must carry\n")
            errprintf("  the variables it was estimated from." +
                      " dwls_sample_stats() fills those in.\n")
            _error("dwlsmv_psi_matrix: pair fits are not in the map's order")
        }

        K1 = map.K[vi]
        K2 = map.K[vj]

        psi_a = psi_tau[|koff[vi] + 1, map.Coff[vi] + 1 \
                         koff[vi] + K1, map.Coff[vi] + map.C[vi]|]
        psi_b = psi_tau[|koff[vj] + 1, map.Coff[vj] + 1 \
                         koff[vj] + K2, map.Coff[vj] + map.C[vj]|]

        tbl = dwlsmv_psi_rho_table(pairs[k], psi_a, psi_b)

        // Row-major flattening, so cell (r,s) sits at (r-1)*K2 + s and the
        // per-observation index is a single arithmetic expression rather than
        // a loop over N.
        flat = J(K1 * K2, 1, 0)
        for (r = 1; r <= K1; r++) {
            for (s = 1; s <= K2; s++) {
                flat[(r - 1) * K2 + s] = tbl[r, s]
            }
        }

        Psi[., pos] = flat[Y[., vi] :* K2 :+ Y[., vj] :+ 1]
    }

    return(Psi)
}

end
