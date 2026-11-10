*! dwlsmv -- the sandwich bread's diagnostic (Module 10, spec Section 5)
*!
*! Executed under StataMP 18 on 2026-09-26. tests/unit/test_breaddiag.do: case D
*! healthy (ratio 0.1200, equal to an eigenvalue path and distinct from the
*! unweighted 0.1196), production noise 3.7e-12 on W^1/2 Delta, the deficient
*! twin at rank 14 of 15 (ratio 2.3e-13), a zero column, npar > D, rc 3498.
*!
*!     real matrix dwlsmv_breaddiag(real matrix Delta, real matrix Gamma)
*!
*! Spec Section 5, Module 10 and E-060 and E-071. Returns
*! 4 x 4, one row per condition, always in this order, in the [flag, value,
*! i, j] layout of dwlsmv_admissible():
*!
*!     row  condition                    flag             value      i, j
*!      1   more parameters than stats   npar > D         npar       D, 0
*!      2   an identically zero column   any              how many   first such column, 0
*!      3   numerical rank, W^1/2 Delta  rank < m         rank       m, 0
*!      4   ratio, W^1/2 Delta           .n               smin/smax  m, 0
*!
*! npar = cols(Delta), derived rather than passed. m is the number of columns
*! rows 3-4 judge: the non-zero ones.
*!
*! THREE KINDS OF ROW, NEVER MIXED, as in dwlsmv_gammadiag(). Rows 1-2 are
*! structural deficiency, exact and detected from their cause. Row 3 is the
*! numerical rank deficiency, against dwlsmv_numrank()'s criterion -- the one
*! definition dwls_identify() also applies, there to Delta at the starting
*! values. Row 4 is the near-singularity ratio, value only (flag .n). Rows 3-4
*! never report a structural deficiency a second time: they judge only the
*! non-zero columns, and are `.` (NOT evaluated) when row 1 flags or when no
*! column is non-zero.
*!
*! THE SCALE IS W_D^1/2 Delta at theta-hat (E-060, 2026-09-25, (6)),
*! the matrix whose Gram matrix is the bread Delta' W_D Delta that
*! dwlsmv_acov() inverts, with W_D = 1 / diag(Gamma) formed here as a vector.
*! Derivation @sec-bread-rank establishes that the 1e-8 criterion keeps its
*! margins on this weighted matrix: the weights cancel from the noise bound.
*!
*! WHY THE CALLER PASSES Delta, NOT A FIT: so the same rows can be computed
*! on lavaan's Delta and Gamma, as dwlsmv_traces() can. At a fitted model the
*! caller passes fit.Delta, which dwls_fit() evaluates at theta-hat.
*!
*! A FUNCTION, NOT A FIELD. dwls_fit() and dwls_inference() do not call it;
*! Module 11 does.
*!
*! CONTRACT (rc 3498): Delta non-empty with no missing value; Gamma square,
*! with as many rows as Delta and a strictly positive diagonal (W_D must
*! exist, as dwlsmv_acov() also requires).
*!
*! Requires dwlsmv_numrank.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_breaddiag(real matrix Delta, real matrix Gamma)
{
    real matrix A, M, V
    real colvector sw, sv, nz
    real scalar D, np, c, nzero, first

    D  = rows(Delta)
    np = cols(Delta)
    if (D < 1 | np < 1 | missing(Delta)) {
        errprintf("dwlsmv_breaddiag: Delta is %f x %f; it must be non-empty with no missing value.\n",
                  D, np)
        _error("dwlsmv_breaddiag: Delta is empty or has missing values")
    }
    if (rows(Gamma) != D | cols(Gamma) != D) {
        errprintf("dwlsmv_breaddiag: Gamma is %f x %f but Delta has %f rows.\n",
                  rows(Gamma), cols(Gamma), D)
        _error("dwlsmv_breaddiag: Gamma is not D x D")
    }
    if (missing(diagonal(Gamma)) | min(diagonal(Gamma)) <= 0) {
        errprintf("dwlsmv_breaddiag: diag(Gamma) must be strictly positive; W_D = 1 / diag(Gamma).\n")
        _error("dwlsmv_breaddiag: diag(Gamma) is not strictly positive")
    }

    A = (J(4, 2, .), J(4, 2, 0))

    // ---- 1. npar > D ---------------------------------------------------------
    A[1, .] = (np > D, np, D, 0)

    // ---- 2. identically zero columns: exact, no tolerance ----------------------
    nz    = J(np, 1, 1)
    nzero = 0
    first = 0
    for (c = 1; c <= np; c++) {
        if (sum(Delta[., c] :!= 0) == 0) {
            nz[c] = 0
            nzero = nzero + 1
            if (first == 0) {
                first = c
            }
        }
    }
    A[2, .] = (nzero > 0, nzero, first, 0)

    // ---- 3-4. on the non-zero columns of W^1/2 Delta -----------------------------
    if (A[1, 1] | nzero == np) {
        return(A)
    }
    sw = sqrt(1 :/ diagonal(Gamma))
    M  = select(sw :* Delta, nz')
    sv = J(0, 1, 0)
    V  = J(0, 0, 0)
    A[3, 2] = dwlsmv_numrank(M, sv, V)
    A[3, 1] = A[3, 2] < cols(M)
    A[3, 3] = cols(M)
    A[4, .] = (.n, sv[rows(sv)] / sv[1], cols(M), 0)
    return(A)
}

end
