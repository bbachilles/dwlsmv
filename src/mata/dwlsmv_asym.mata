*! dwlsmv -- the Gamma symmetry measure and its budget (Module 10, spec Section 5)
*!
*! Executed under StataMP 18 on 2026-09-25. tests/unit/test_asym.do: our Gamma
*! reads exactly 0, both lavaan fixtures sit inside the budget at the
*! derivation's R figures, and dwlsmv_traces() and dwlsmv_acov() give the same
*! verdict as this routine at 0.5x and 2x the budget.
*!
*!     real rowvector dwlsmv_asym(real matrix A)
*!
*! Returns (a, budget, i, j), where
*!
*!     a      = max_ij |A_ij - A_ji| / sqrt(A_ii A_jj)
*!     budget = 4 D eps,  D = rows(A)
*!     (i, j) = where a is attained, i < j, the first in column order;
*!              (0, 0) when a = 0
*!
*! and A is outside the symmetry contract exactly when a > budget. The
*! location was added 2026-09-26 for dwlsmv_gammadiag()'s row 1, so the
*! diagnostic reports it without a second copy of the measure. The guards
*! read only columns 1 and 2.
*!
*! THE ONE DEFINITION. dwlsmv_traces() and dwlsmv_acov() both use
*! Gamma' = Gamma, and both abort when a > budget; Module 10's diagnostic
*! reports the same comparison. None of them carries a tolerance of its own:
*! each calls this routine. A second copy of the budget is how a guard comes
*! to reject a Gamma that the diagnostic calls acceptable, or the reverse.
*!
*! WHY THIS MEASURE (derivation @sec-gamma-sym). It is invariant to rescaling
*! the statistics, A -> S A S for a positive diagonal S, exactly as theta-hat,
*! ACOV and T3 are. max|A - A'| / max|A| is not, and would be blind to
*! asymmetry in a small-variance block next to a large one.
*!
*! WHY THIS BUDGET. It is the first-order floating-point bound for a Gamma
*! formed as two matrix products, fl(fl(L S) L') -- lavaan's B^-1 S B^-T -- with
*! the growth factor taken as 1: two products, two triangles, so 4 D eps. The
*! derivation says plainly that it is argued, not proved. Our own Gamma is
*! exactly symmetric (quadcross() mirrors one triangle) and reads 0.
*!
*! A positive diagonal is a precondition, not a tolerance: W_D = diag(Gamma)^-1
*! must exist anyway, and both guards check it before calling this.
*!
*! The D x D matrix of sqrt(A_ii) sqrt(A_jj) is formed once. At D = 940 that
*! is 7 MB and a few milliseconds, next to Stage 1's seconds.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real rowvector dwlsmv_asym(real matrix A)
{
    real colvector s
    real matrix M
    real rowvector cm
    real scalar a, i, j, t

    if (rows(A) < 1 | rows(A) != cols(A)) {
        errprintf("dwlsmv_asym: A is %f x %f; it must be a non-empty square matrix.\n",
                  rows(A), cols(A))
        _error("dwlsmv_asym: A is not a non-empty square matrix")
    }
    if (missing(A)) {
        _error("dwlsmv_asym: A contains missing values")
    }
    if (min(diagonal(A)) <= 0) {
        errprintf("dwlsmv_asym: A has a zero or negative diagonal entry, so the scale-free\n")
        errprintf("  measure |A_ij - A_ji| / sqrt(A_ii A_jj) is not defined.\n")
        _error("dwlsmv_asym: diag(A) is not strictly positive")
    }

    s = sqrt(diagonal(A))
    M = abs(A - A') :/ (s * s')
    a = max(M)

    // Where the worst entry is, for Module 10's row. M is symmetric with a
    // zero diagonal, so the first column holding the maximum is the smaller
    // index of the pair; report it as (i, j) with i < j. O(D) interpreted
    // steps, not O(D^2): colmax() does the scan. A symmetric A has no
    // location, (0, 0).
    i = 0
    j = 0
    if (a > 0) {
        cm = colmax(M)
        for (j = 1; cm[j] != a; j++) {
        }
        for (i = 1; M[i, j] != a; i++) {
        }
        if (i > j) {
            t = i
            i = j
            j = t
        }
    }
    return((a, 4 * rows(A) * epsilon(1), i, j))
}

end
