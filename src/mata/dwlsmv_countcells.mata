*! dwlsmv -- category counts within each value of a binary covariate (Module 10, conditional mode)
*!
*! Executed 2026-11-15 by the main session through tests/unit/test_condx_diag.do (StataMP 18). Drafted by a read-only agent
*! under the owner's decision of that day (E-107, point 3).
*!
*!     real matrix dwlsmv_countcells(real matrix Y, real colvector K, real colvector x)
*!
*! WHY. For a binary item the conditional estimates are tau-hat = Phi^-1(p0) and
*! b-hat = Phi^-1(p0) - Phi^-1(p1): two binomials, each within its own x-group
*! (derivation @sec-condx-m10 M4). The POOLED category count that dwlsmv_countdiag() reports
*! cannot see a thin group cell: pooled counts of 701 and 800 for y = 0 gave a slope SE of 0.35
*! against 0.09, the first with one observation in the (y = 0, x = 1) cell. The pooled rows stay
*! (decided 2026-11-15); these rows are added beside them, never in place of them.
*!
*! RETURNS (2 + 2 sum(K)) x 4, a table of COUNTS, not the [flag, value, i, j] layout of the other
*! Module 10 rows (a count has no flag):
*!
*!     rows 1-2   (0, 0, 0, n_0) and (0, 0, 1, n_1): the two group sizes
*!     then, for every item j, category position c = 1 .. K_j (code + 1), group g = 0, 1, in that
*!     order (item slowest, group fastest): (j, c, g, count)
*!
*! columns: item, category position, group, count. Counts of an empty declared category read 0.
*! For every item and category the two group counts sum to dwlsmv_countdiag()'s pooled count,
*! and for every item the counts of one group sum to that group's size.
*!
*! WHAT IT DOES NOT SAY. K_j > 2: the parallel shift is not saturated, so (tau, b) are not
*! functions of group proportions and no exact reference exists (as for an interior category in
*! v1). Pair-cell counts: nothing derived. No row has a threshold; values only.
*!
*! CONTRACT (rc 3498): Y as dwls_sample_stats_cx() receives it, validated per column by
*! dwlsmv_catcounts() (a missing, non-integer or out-of-range code); K of length cols(Y); x with
*! one row per row of Y, no missing value, only the codes 0 and 1, and both codes present.
*!
*! Requires dwlsmv_catcounts.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_countcells(real matrix Y, real colvector K, real colvector x)
{
    real matrix A
    real scalar j, c, g, r, n0, n1

    if (cols(Y) != rows(K)) {
        errprintf("dwlsmv_countcells: Y has %f columns but K names %f variables.\n",
                  cols(Y), rows(K))
        _error("dwlsmv_countcells: Y and K disagree on the number of variables")
    }
    if (rows(x) != rows(Y)) {
        errprintf("dwlsmv_countcells: x has %f rows but Y has %f.\n", rows(x), rows(Y))
        _error("dwlsmv_countcells: x and Y disagree on the number of rows")
    }
    if (missing(x) > 0) {
        errprintf("dwlsmv_countcells: x has %f missing values; apply listwise deletion first.\n",
                  missing(x))
        _error("dwlsmv_countcells: x has missing values")
    }
    if (sum((x :!= 0) :& (x :!= 1)) > 0) {
        errprintf("dwlsmv_countcells: x must be coded 0 and 1; %f values are neither.\n",
                  sum((x :!= 0) :& (x :!= 1)))
        _error("dwlsmv_countcells: x must be coded 0 and 1")
    }
    n0 = sum(x :== 0)
    n1 = sum(x :== 1)
    if (n0 == 0 | n1 == 0) {
        errprintf("dwlsmv_countcells: x takes only one value (%f with x = 0, %f with x = 1).\n", n0, n1)
        _error("dwlsmv_countcells: x must take both the values 0 and 1")
    }

    A = J(2 + 2 * sum(K), 4, 0)
    A[1, .] = (0, 0, 0, n0)
    A[2, .] = (0, 0, 1, n1)
    r = 2
    for (j = 1; j <= cols(Y); j++) {
        // validates the column exactly as Stage 1 does
        (void) dwlsmv_catcounts(Y[., j], K[j])
        for (c = 1; c <= K[j]; c++) {
            for (g = 0; g <= 1; g++) {
                r = r + 1
                A[r, .] = (j, c, g, sum((Y[., j] :== (c - 1)) :& (x :== g)))
            }
        }
    }
    return(A)
}

end
