*! dwlsmv -- invert a position in s back to what it indexes (Section 4.2)
*!
*! Conditional-mode form executed 2026-11-10, suite green (the slope branch, drafted
*! 2026-11-10 by a read-only agent). Without slopes D - P = T, so the
*! correlation branch is unchanged and the slope branch is never entered.
*!
*! Returns (stat_type, var_i, var_j, thr_k):
*!     threshold   -> (1, j, 0, k)
*!     polychoric  -> (2, i, j, 0)
*!     slope       -> (3, j, k, 0)    conditional mode: indicator j, covariate k
*!
*! Needed only for diagnostics, but closed form anyway. For r = s - (D - P)
*! in the correlation block,
*!
*!     i = ceil( { (2p-1) - sqrt((2p-1)^2 - 8r) } / 2 )
*!     j = r - { (i-1)p - i(i-1)/2 } + i
*!
*! The ceiling is exact, not lucky: at every block boundary the discriminant
*! equals (2p-2i-1)^2, a perfect square, so the root is exact and cannot tip
*! the ceiling by one. Verified in the derivation file over the four fixture
*! designs plus four shapes they do not cover.
*!
*! The slope block is covariate-major (spec 4.2): for r = s - T,
*! k = ceil(r / p) and j = r - (k-1)p, both exact in integer arithmetic.
*!
*! The threshold branch scans rather than inverting, because p is small and a
*! scan cannot be wrong.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real rowvector dwlsmv_pos_inv(struct dwlsmv_statmap scalar map, real scalar s)
{
    real scalar r, i, j, k, q

    if (s < 1 | s > map.D | s != floor(s)) {
        _error("dwlsmv_pos_inv: position out of range")
    }

    if (s <= map.T) {
        for (q = 1; q <= map.p; q++) {
            if (s <= map.Coff[q] + map.C[q]) {
                k = s - map.Coff[q]
                return((1, q, 0, k))
            }
        }
        _error("dwlsmv_pos_inv: threshold block scan fell through")
    }

    if (s <= map.D - map.P) {
        r = s - map.T
        k = ceil(r / map.p)
        j = r - (k - 1) * map.p
        return((3, j, k, 0))
    }

    r = s - (map.D - map.P)
    i = ceil(((2 * map.p - 1) - sqrt((2 * map.p - 1)^2 - 8 * r)) / 2)
    j = r - ((i - 1) * map.p - i * (i - 1) / 2) + i
    return((2, i, j, 0))
}

end
