*! dwlsmv -- position of the polychoric correlation (i,j) in s (Section 4.2)
*!
*! Conditional-mode form executed 2026-11-10, suite green (drafted 2026-11-10 by a
*! read-only agent). The offset is now D - P, which equals T exactly when the
*! map has no slopes, so every unconditional position is unchanged.
*!
*! pos(rho_ij) = (D - P) + (i-1)p - i(i-1)/2 + (j-i), requiring i < j.
*! D - P is T + pq: the correlations follow the thresholds and, in
*! conditional mode, the slopes (spec 4.2). q is never stored; the offset is
*! derived from the map's own D, T and P.
*!
*! i >= j aborts rather than silently swapping. A caller that has its pair
*! ordered backwards has a bug, and swapping for it hides the bug while
*! returning a plausible number.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_pos_rho(struct dwlsmv_statmap scalar map,
                           real scalar i, real scalar j)
{
    if (i < 1 | j > map.p | i != floor(i) | j != floor(j)) {
        _error("dwlsmv_pos_rho: variable index out of range")
    }
    if (i >= j) {
        _error("dwlsmv_pos_rho: requires i < j; the pair is not ordered")
    }
    return(map.D - map.P + (i - 1) * map.p - i * (i - 1) / 2 + (j - i))
}

end
