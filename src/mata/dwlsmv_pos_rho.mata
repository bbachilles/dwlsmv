*! dwlsmv -- position of the polychoric correlation (i,j) in s (Section 4.2)
*!
*! pos(rho_ij) = T + (i-1)p - i(i-1)/2 + (j-i), requiring i < j.
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
    return(map.T + (i - 1) * map.p - i * (i - 1) / 2 + (j - i))
}

end
