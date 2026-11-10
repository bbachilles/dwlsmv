*! dwlsmv -- the stacked Stage 1 threshold lookup tables (Module 4)
*!
*! Executed under StataMP 18 on 2026-09-20. The container decision it settles is
*! exercised end to end by tests/unit/test_gamma_tau.do, whose threshold block
*! reproduces the observed joint proportions to 1.2e-15 over all 100 entries.
*!
*! THIS FILE SETTLES AN OPEN ITEM. The development record carried "psi_tau's
*! container shape for dwls_gamma()" as undecided: spec 5.1 locks the fourth
*! argument of dwls_gamma() as a `real matrix` and does not say WHICH matrix.
*! Two readings were live -- the N x T gathered influence matrix, or the
*! per-variable K_j x C_j lookup tables stacked into one non-ragged block. The
*! decision was deliberately deferred to Module 4, "where there is a Gamma
*! assembly to judge it against". There is one now, and it decides it:
*!
*!     psi_tau is the STACKED LOOKUP TABLES, sum(K) x T.
*!
*! Variable j owns rows Koff[j]+1 .. Koff[j]+K[j] (see dwlsmv_koff()) and
*! columns Coff[j]+1 .. Coff[j]+C[j]. Every other entry is zero. So
*!
*!     psi_tau[Koff[j] + y + 1, dwlsmv_pos_tau(map, j, c)]
*!
*! is the influence of an observation sitting in category y of variable j on
*! threshold c of variable j, and the COLUMN INDEX IS THE POSITION IN s.
*!
*! WHY NOT THE N x T GATHERED MATRIX. Module 4 does not only gather the
*! threshold block -- it also has to build each pair's psi_rho lookup table,
*! and spec 4.5.2's expression for that is indexed by CATEGORY:
*!
*!     psi_rho[r,s] = ( u[r,s] - sum_c I_ra[c] psi_a[r,c]
*!                             - sum_d I_rb[d] psi_b[s,d] ) / I_rr
*!
*! psi_a[r,c] is exactly this table. Handing Module 4 the gathered N x T matrix
*! would mean either recovering the tables from it (N x T to read back K_j
*! distinct rows) or rebuilding them from q, which is a second route to the
*! same quantity. The gathered matrix is a DERIVED object -- one lookup away --
*! and the table is the primitive, so the primitive is what crosses the
*! boundary.
*!
*! WHY T COLUMNS AND NOT max(C). A max(C)-column block is smaller and equally
*! non-ragged, but it forces every reader to carry a second offset convention
*! for columns on top of Koff for rows, and the two are easy to mix up because
*! both are per-variable prefix sums. With T columns the column index is the
*! index into s, which the caller already has from dwlsmv_pos_tau(). The cost
*! is a mostly-zero matrix of sum(K) x T doubles -- 200 x 160 at p = 40, K = 5,
*! which is 0.25 MB and is not a consideration at any size this estimator runs.
*!
*! The zero padding is never read by dwlsmv_psi_matrix(), which copies each
*! variable's own column block rather than summing over variables. It is there
*! so the container is self-describing.
*!
*! Requires dwlsmv_types, dwlsmv_statmap_build (for the map it is handed),
*! dwlsmv_pos_tau, dwlsmv_catcounts, dwlsmv_cumprop and dwlsmv_psi_tau_table.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_psi_tau_stack(real matrix Y, struct dwlsmv_statmap scalar map)
{
    real matrix stack, tbl
    real colvector koff, n, q
    real scalar j, c, N

    if (cols(Y) != map.p) {
        errprintf("dwlsmv_psi_tau_stack: Y has %s columns, the map declares %s variables.\n",
                  strofreal(cols(Y)), strofreal(map.p))
        errprintf("  The map is keyed on K alone, so it cannot adapt to the data;" +
                  " the two must be\n")
        errprintf("  built from the same K in the same order.\n")
        _error("dwlsmv_psi_tau_stack: Y and the map disagree on the number of variables")
    }
    N = rows(Y)
    if (N < 1) {
        _error("dwlsmv_psi_tau_stack: Y has no observations")
    }

    koff  = dwlsmv_koff(map)
    stack = J(koff[map.p] + map.K[map.p], map.T, 0)

    for (j = 1; j <= map.p; j++) {
        n   = dwlsmv_catcounts(Y[., j], map.K[j])
        q   = dwlsmv_cumprop(n)

        // dwlsmv_psi_tau_table() re-derives tau from q rather than accepting
        // it, and refuses a degenerate q. Both degeneracies of spec 4.5.6
        // therefore abort here rather than producing infinities downstream.
        tbl = dwlsmv_psi_tau_table(q)

        for (c = 1; c <= map.C[j]; c++) {
            stack[|koff[j] + 1, dwlsmv_pos_tau(map, j, c) \
                   koff[j] + map.K[j], dwlsmv_pos_tau(map, j, c)|] = tbl[., c]
        }
    }

    return(stack)
}

end
