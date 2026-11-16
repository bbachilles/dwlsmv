*! dwlsmv -- position of the slope b_jk of indicator j on covariate k in s
*! (spec 4.2, conditional paragraph)
*!
*! Executed 2026-11-10, suite green (drafted by a read-only agent, run by the
*! main session).
*!
*! pos(b_jk) = T + (k-1)p + j: covariate-major, items within each covariate,
*! as lavaan's layout orders them (@sec-condx-norm, proved from lavaan's own
*! labels on condx_A).
*!
*! Out-of-range arguments abort rather than return missing, as
*! dwlsmv_pos_tau() does. On an unconditional map every k is out of range.
*!
*! Requires dwlsmv_statmap_q.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_pos_slope(struct dwlsmv_statmap scalar map,
                             real scalar j, real scalar k)
{
    real scalar q

    q = dwlsmv_statmap_q(map)
    if (j < 1 | j > map.p | j != floor(j)) {
        _error("dwlsmv_pos_slope: variable index out of range")
    }
    if (k < 1 | k > q | k != floor(k)) {
        _error("dwlsmv_pos_slope: covariate index out of range for this map")
    }
    return(map.T + (k - 1) * map.p + j)
}

end
