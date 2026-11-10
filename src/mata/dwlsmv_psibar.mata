*! dwlsmv -- the centring diagnostic of spec 4.5.5 (Module 4)
*!
*! Executed under StataMP 18 on 2026-09-20. tests/unit/test_gamma.do asserts both
*! blocks separately and checks this count route against the gathered matrix.
*!
*! Returns D x 1 mean influence, psibar = (1/N) sum_i psi_i.
*!
*! COMPUTED FROM CELL COUNTS, NOT FROM THE GATHERED MATRIX, and the difference
*! is deliberate twice over. psi is a lookup table, so its mean is exactly the
*! count-weighted mean of the table -- O(K^2) per pair rather than O(N*D), and
*! exact rather than an N-term accumulation. And it makes this an INDEPENDENT
*! route to a quantity the gather also implies, so a disagreement between the
*! two is informative instead of tautological.
*!
*! State plainly what that does and does not buy, because the diagnostic has
*! been over-read before. It checks the CONSTRUCTION of psi -- a derivative
*! error, a sign, a mis-weighted propagation term. It does NOT check the
*! gather: an indexing error in dwlsmv_psi_matrix() would leave this exactly
*! zero. The gather is covered at machine precision by the threshold block of
*! Gamma, which spec 4.5.1 makes an exact identity in the sample, so nothing is
*! left uncovered -- but the two checks are checking different things and
*! neither substitutes for the other.
*!
*! AND IT CANNOT DETECT A MISSING PROPAGATION TERM AT ALL. The naive influence
*! function -- the direct score with the threshold propagation dropped -- is
*! ALSO exactly centred. tests/unit/test_polychoric.do measures that rather
*! than taking it on faith: both means come back 0.0 while the naive version
*! understates N*Var(rho-hat) by 8.2%. So a psibar of 1e-16 is not evidence
*! that Gamma is right. Spec 4.5.5 says the same thing in one line: it catches
*! derivative errors, not omissions.
*!
*! SPLIT IT BY BLOCK WHEN ASSERTING ON IT, because the two halves are not
*! equally exact and a single tolerance is set by the looser one. The threshold
*! entries are exactly zero to machine precision with no dependence on any
*! optimizer, since q-hat is literally the sample proportion. The correlation
*! entries are zero only to Module 2's convergence tolerance, because rho-hat
*! solves its score equation numerically. A threshold-entry mean of 1e-8 would
*! pass a blanket check and is already a bug. The one exception is the rho
*! entry of a one-empty-cell pair (O6, spec 4.5.5, @sec-o6): rho-hat+ solves
*! the score on the REPAIRED table while psi is scored on the observed one, so
*! its mean is S-bar(rho-hat+)/I_rr, nonzero at O(1) by construction; Module 10
*! excludes it. Its threshold entries keep the machine-precision rule.
*!
*! Requires dwlsmv_types, dwlsmv_koff, dwlsmv_pos_tau, dwlsmv_psi_rho_table,
*! dwlsmv_catcounts and dwlsmv_crosstab.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real colvector dwlsmv_psibar(real matrix Y,
                             struct dwlsmv_statmap scalar map,
                             struct dwlsmv_pairfit vector pairs,
                             real matrix psi_tau)
{
    real matrix n2, tbl, psi_a, psi_b
    real colvector bar, koff, n1
    real scalar N, j, c, k, pos, vi, vj, K1, K2

    N = rows(Y)
    if (cols(Y) != map.p) {
        _error("dwlsmv_psibar: Y and the map disagree on the number of variables")
    }
    if (length(pairs) != map.P) {
        _error("dwlsmv_psibar: wrong number of pair fits")
    }

    koff = dwlsmv_koff(map)
    bar  = J(map.D, 1, 0)

    for (j = 1; j <= map.p; j++) {
        n1 = dwlsmv_catcounts(Y[., j], map.K[j])
        for (c = 1; c <= map.C[j]; c++) {
            pos      = dwlsmv_pos_tau(map, j, c)
            bar[pos] = sum(n1 :* psi_tau[|koff[j] + 1,           pos \
                                          koff[j] + map.K[j],    pos|]) / N
        }
    }

    for (k = 1; k <= map.P; k++) {
        pos = map.T + k
        vi  = map.var_i[pos]
        vj  = map.var_j[pos]
        K1  = map.K[vi]
        K2  = map.K[vj]

        psi_a = psi_tau[|koff[vi] + 1, map.Coff[vi] + 1 \
                         koff[vi] + K1, map.Coff[vi] + map.C[vi]|]
        psi_b = psi_tau[|koff[vj] + 1, map.Coff[vj] + 1 \
                         koff[vj] + K2, map.Coff[vj] + map.C[vj]|]

        tbl      = dwlsmv_psi_rho_table(pairs[k], psi_a, psi_b)
        n2       = dwlsmv_crosstab(Y[., vi], Y[., vj], K1, K2)
        bar[pos] = sum(n2 :* tbl) / N
    }

    return(bar)
}

end
