*! dwlsmv -- the threshold influence-function staircase (Module 1, spec 4.5.1)
*!
*! Executed under StataMP 18 on 2026-09-19. Verified two independent ways at
*! machine precision: the within-variable identity of spec 4.5.1 (gap 9e-16
*! against a 1e-12 budget) and the derivation file's tridiagonal route,
*! A^-1 s, which never mentions a cumulative proportion.
*!
*!     psi[tau_c, i] = ( 1{Y_i < c} - q_c ) / phi(tau_c)
*!
*! which depends on observation i ONLY through which category it landed in. So
*! a K-category variable has exactly K distinct influence vectors, and this
*! returns all of them as a K x (K-1) table: row y+1 holds the influence vector
*! of an observation in category y, for y = 0 .. K-1. Index it by the observed
*! category. No per-observation influence function is ever computed.
*!
*! Read that claim precisely, because an earlier version of the project's notes
*! overstated the neighbouring one for Module 2. This table removes per
*! observation ARITHMETIC from the threshold block. It does not remove O(N)
*! work from Gamma: forming N^-1 sum_i psi_i psi_i' still gathers an N x D
*! matrix by lookup and crosses it. "Lookup table" is not "no O(N) work".
*!
*! Each row is a staircase: the entry is +(1-q_c)/phi_c while c > y and
*! -q_c/phi_c once c <= y, switching sign exactly once, at c = y.
*!
*! TAKES q, NOT (q, tau). tau = invnormal(q) is recomputed here rather than
*! accepted as a second argument, which costs K-1 invnormal calls per variable
*! and removes the possibility of being handed a tau that does not correspond
*! to the q beside it. A derived quantity passed alongside its source is a
*! second copy that can disagree with the first.
*!
*! The q it is handed must be strictly increasing and strictly inside (0,1).
*! dwlsmv_thresholds() has already refused the data that violates this, but
*! this routine is callable on a hand-built q and checks for itself -- a q of 0
*! or 1 divides by phi = 0 and returns infinities that look like numbers.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_psi_tau_table(real colvector q)
{
    real matrix psi
    real colvector phi
    real scalar C, K, y, c

    if (cols(q) != 1) {
        _error("dwlsmv_psi_tau_table: q must be a column vector")
    }
    C = rows(q)
    if (C < 1) {
        _error("dwlsmv_psi_tau_table: q must hold at least one cumulative proportion")
    }
    if (min(q) <= 0 | max(q) >= 1) {
        _error("dwlsmv_psi_tau_table: q must lie strictly inside (0,1) -- an empty boundary category")
    }
    for (c = 2; c <= C; c++) {
        if (q[c] <= q[c - 1]) {
            _error("dwlsmv_psi_tau_table: q must be strictly increasing -- an empty interior category")
        }
    }

    K   = C + 1
    phi = normalden(invnormal(q))
    psi = J(K, C, 0)

    for (y = 0; y <= C; y++) {
        for (c = 1; c <= C; c++) {
            psi[y + 1, c] = ((c > y) - q[c]) / phi[c]
        }
    }

    return(psi)
}

end
