*! dwlsmv -- threshold aliasing, detected from its cause (Module 10, spec 4.5.6)
*!
*! Executed under StataMP 18 on 2026-09-25. tests/unit/test_structural.do: no
*! relation on case D; exactly y1|t1 = y2|t1 on the derivation's aliased
*! design, whose production Gamma-hat is then numerically singular.
*!
*!     real matrix dwlsmv_thralias(real matrix Y, real colvector K)
*!
*! Derivation @sec-gamma-rank, Theorem 3: with Z the N x T matrix of cumulative
*! indicators 1{Y_ij < c}, rank(Gamma-hat_tautau) = rank[1 | Z] - 1, and
*! Gamma-hat is singular whenever that block is. This routine finds the
*! PAIRWISE linear relations among the columns of [1 | Z]:
*!
*!     type 1  identical      z_a = z_b          a < b
*!     type 2  complementary  z_a = 1 - z_b      a < b
*!     type 3  constant       z_a = 0 or 1       b = 0
*!
*! Returns one row (a, b, type) per relation, a and b being positions in s
*! (thresholds first, by variable, ascending: spec 4.2), or a 0 x 3 matrix.
*!
*! EXACT, WITH NO CUTOFF. Z is 0/1, so G = Z'Z holds integer counts, exact in
*! double precision whatever the summation order. With n_a = G_aa,
*! n_ab = G_ab and N rows: identical iff n_a = n_b = n_ab; complementary iff
*! n_ab = 0 and n_a + n_b = N; constant iff n_a is 0 or N. A floating-point
*! rank of [1 | Z] would need a threshold, and "detected from its cause"
*! (spec Section 5, Module 10) rules one out.
*!
*! WHAT IT DOES NOT SEE, a recorded v1 limit: a linear relation among three or
*! more columns (z_1 + z_2 = z_3 + z_4 on every row) is not pairwise. Every
*! route with evidence behind it is pairwise: the interior-empty category
*! (two nested columns of one variable coincide; Module 1 aborts on it anyway)
*! and cross-variable aliasing (measured in production, rank 15 of 16).
*!
*! A constant column means an empty boundary category, which Module 1 also
*! aborts on. It is reported for completeness and then left out of the
*! pairwise comparison, so two constant columns are not also reported as
*! identical.
*!
*! Same (Y, K) and same per-column validation as dwls_sample_stats() and
*! dwlsmv_npatterns(): dwlsmv_catcounts(), rc 3498 on a missing, non-integer
*! or out-of-range code.
*!
*! Requires dwlsmv_catcounts.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_thralias(real matrix Y, real colvector K)
{
    real matrix Z, G, R
    real colvector nz, isconst
    real scalar N, T, j, c, pos, a, b

    if (cols(Y) != rows(K)) {
        errprintf("dwlsmv_thralias: Y has %f columns but K names %f variables.\n",
                  cols(Y), rows(K))
        _error("dwlsmv_thralias: Y and K disagree on the number of variables")
    }
    for (j = 1; j <= cols(Y); j++) {
        (void) dwlsmv_catcounts(Y[., j], K[j])
    }

    N = rows(Y)
    T = sum(K :- 1)
    Z = J(N, T, 0)
    pos = 0
    for (j = 1; j <= cols(Y); j++) {
        for (c = 1; c <= K[j] - 1; c++) {
            pos = pos + 1
            Z[., pos] = (Y[., j] :< c)
        }
    }

    G  = quadcross(Z, Z)
    nz = diagonal(G)
    isconst = (nz :== 0) :| (nz :== N)

    R = J(0, 3, .)
    for (a = 1; a <= T; a++) {
        if (isconst[a]) {
            R = R \ (a, 0, 3)
            continue
        }
        for (b = a + 1; b <= T; b++) {
            if (isconst[b]) {
                continue
            }
            if (nz[a] == nz[b] & G[a, b] == nz[a]) {
                R = R \ (a, b, 1)
            }
            else if (G[a, b] == 0 & nz[a] + nz[b] == N) {
                R = R \ (a, b, 2)
            }
        }
    }
    return(R)
}

end
