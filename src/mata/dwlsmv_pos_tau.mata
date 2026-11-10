*! dwlsmv -- position of threshold c of variable j in s (spec Section 4.2)
*!
*! pos(tau_jc) = sum_{m<j} C_m + c, which Coff already holds.
*!
*! Out-of-range arguments abort rather than return missing. A missing index
*! propagates into a subscript and surfaces far from its cause; this project's
*! failure mode is code that runs and is quietly wrong, so the loud form wins.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_pos_tau(struct dwlsmv_statmap scalar map,
                           real scalar j, real scalar c)
{
    if (j < 1 | j > map.p | j != floor(j)) {
        _error("dwlsmv_pos_tau: variable index out of range")
    }
    if (c < 1 | c > map.C[j] | c != floor(c)) {
        _error("dwlsmv_pos_tau: threshold index out of range for this variable")
    }
    return(map.Coff[j] + c)
}

end
