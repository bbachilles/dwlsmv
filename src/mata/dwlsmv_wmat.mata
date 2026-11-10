*! dwlsmv -- the threshold-derivative weights w_{c,s} (Module 2, spec 4.5.2)
*!
*! Executed under StataMP 18 on 2026-09-19. The compact w-weighted form of
*! I_rho_a agrees with numerical differentiation of E[g] to 1.2e-11, which is
*! the check that would catch the two coupled rows being taken in the wrong
*! order.
*!
*!   w_{c,s} = phi(a_c) [ Phi((b_{s+1} - rho a_c)/sqrt(1-rho^2))
*!                      - Phi((b_s     - rho a_c)/sqrt(1-rho^2)) ]
*!
*! the joint density of Y_j* = a_c with Y_k falling in category s. Returns
*! C1 x K2, one row per threshold of the FIRST variable and one column per
*! category of the second.
*!
*! WHY THIS IS THE WHOLE THRESHOLD DERIVATIVE. a_c appears in exactly two rows
*! of the table -- the two it separates -- so
*!
*!     d(pi_{c-1,s})/d(a_c) = +w_{c,s},   d(pi_{c,s})/d(a_c) = -w_{c,s},
*!
*! and every other row is zero. That sparsity is what collapses the
*! propagation weight to a sum over ONE variable's categories rather than the
*! whole table:
*!
*!     I_rho,a_c = sum_s w_{c,s} ( u_{c-1,s} - u_{c,s} )
*!
*! It is the K > 2 analogue of the tridiagonal structure in the univariate
*! score of spec 4.5.1 -- a threshold only ever couples to its two neighbouring
*! cells. The derivation file verifies this compact form against the direct
*! sum over the full table; that check is reproduced in the Stata test, because
*! it is the one that would catch the rows being taken in the wrong order.
*!
*! BY SYMMETRY THERE IS NO SECOND ROUTINE for the other variable's thresholds.
*! The derivation file writes dP_db(a,b,rho,d) = t(dP_da(b,a,rho,d)); here the
*! caller swaps the arguments, calling dwlsmv_wmat(b, a, rho) and pairing it
*! with u transposed. Writing a mirrored copy of this function would be a
*! second place for the same index error to live.
*!
*! Infinite b-boundaries are handled as Phi(-inf) = 0 and Phi(+inf) = 1, which
*! is why the first and last columns are branched rather than evaluated.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_wmat(real colvector a, real colvector b, real scalar rho)
{
    real matrix W
    real scalar C1, K2, c, s, sig, lo, hi

    if (cols(a) != 1 | cols(b) != 1) {
        _error("dwlsmv_wmat: thresholds must be column vectors")
    }
    if (rows(a) < 1 | rows(b) < 1) {
        _error("dwlsmv_wmat: each variable needs at least one threshold")
    }
    if (rho <= -1 | rho >= 1 | rho >= .) {
        _error("dwlsmv_wmat: rho must lie strictly inside (-1, 1)")
    }

    C1  = rows(a)
    K2  = rows(b) + 1
    sig = sqrt(1 - rho * rho)

    W = J(C1, K2, 0)
    for (c = 1; c <= C1; c++) {
        for (s = 1; s <= K2; s++) {
            // Upper edge of category s: b[s] unless s is the last category.
            if (s == K2) {
                hi = 1
            }
            else {
                hi = normal((b[s] - rho * a[c]) / sig)
            }
            // Lower edge: b[s-1] unless s is the first category.
            if (s == 1) {
                lo = 0
            }
            else {
                lo = normal((b[s - 1] - rho * a[c]) / sig)
            }
            W[c, s] = normalden(a[c]) * (hi - lo)
        }
    }

    return(W)
}

end
