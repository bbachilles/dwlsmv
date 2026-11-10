*! dwlsmv -- the threshold-alias deficiency, counted by components (Module 10)
*!
*! Executed under StataMP 18 on 2026-09-26. tests/unit/test_gammadiag.do check 4:
*! hand-built relation lists (a triangle is 2, a chain of four 3) and three
*! identical binary columns, 3 relations and a deficiency of 2.
*!
*!     real scalar dwlsmv_dalias(real matrix R)
*!
*! R is dwlsmv_thralias()'s output: one row (a, b, type) per pairwise relation
*! among the columns of [1 | Z], type 1 identical, 2 complementary, 3 constant
*! (b = 0). Returns the deficiency d_alias = T + 1 - rank[1 | Z] that those
*! relations imply (spec Section 5, Module 10; derivation @sec-gamma-rank,
*! Theorem 3):
*!
*!   - each constant column counts 1;
*!   - the non-constant columns joined by identical or complementary relations
*!     form connected components, and a component of m columns counts m - 1,
*!     because together with the constant column 1 they span one dimension.
*!
*! BY COMPONENTS, NOT BY RELATIONS. Three identical columns are three
*! relations and a deficiency of two. Counting relations would put a
*! structural zero into Gamma's ratio row, whose r* = min(D - d_alias,
*! n_pat - 1) this feeds.
*!
*! A LOWER BOUND on the true d_alias, exact whenever every relation is
*! pairwise -- the v1 limit dwlsmv_thralias() records. thralias leaves
*! constant columns out of the pairwise comparison, so a column is never both
*! constant and in a component.
*!
*! Its own function, not a loop inside dwlsmv_gammadiag(), so the count can be
*! tested on relations that no Stage 1 could be run on: identical columns
*! drive a polychoric to 1, which Module 2 aborts on.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_dalias(real matrix R)
{
    real colvector lab, nodes
    real scalar k, a, b, la, lb, n, m, ncomp, nconst

    // The width is checked before the empty case, so a 0-row list of the wrong
    // shape is refused too (adversarial review 2026-09-26, D2).
    if (cols(R) != 3) {
        errprintf("dwlsmv_dalias: R has %f columns; dwlsmv_thralias() returns (a, b, type).\n",
                  cols(R))
        _error("dwlsmv_dalias: R is not a thralias relation list")
    }
    if (rows(R) == 0) {
        return(0)
    }
    // thralias's shape exactly: integer positions, a >= 1; types 1-2 carry a
    // partner b > a, type 3 carries b = 0.
    if (missing(R) | sum(R :!= floor(R)) > 0 | min(R[., 1]) < 1 |
        sum((R[., 3] :< 1) :| (R[., 3] :> 3)) > 0 |
        sum((R[., 3] :== 3) :& (R[., 2] :!= 0)) > 0 |
        sum((R[., 3] :!= 3) :& (R[., 2] :<= R[., 1])) > 0) {
        _error("dwlsmv_dalias: R is not in dwlsmv_thralias()'s (a, b, type) form")
    }

    nconst = sum(R[., 3] :== 3)

    // Union by relabelling: every column starts as its own component; each
    // relation merges b's component into a's.
    n   = max(R[., 1..2])
    lab = (1::n)
    for (k = 1; k <= rows(R); k++) {
        if (R[k, 3] != 3) {
            a  = R[k, 1]
            b  = R[k, 2]
            la = lab[a]
            lb = lab[b]
            if (la != lb) {
                for (m = 1; m <= n; m++) {
                    if (lab[m] == lb) {
                        lab[m] = la
                    }
                }
            }
        }
    }

    // The columns that appear in a non-constant relation, and their components.
    nodes = J(0, 1, .)
    for (k = 1; k <= rows(R); k++) {
        if (R[k, 3] != 3) {
            nodes = nodes \ R[k, 1] \ R[k, 2]
        }
    }
    if (rows(nodes) == 0) {
        return(nconst)
    }
    nodes = uniqrows(nodes)
    ncomp = rows(uniqrows(lab[nodes]))
    return(nconst + rows(nodes) - ncomp)
}

end
