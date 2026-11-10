*! dwlsmv -- the one-empty-cell repair, keep margins, c = 0.5 (O6, spec Module 2)
*!
*! Executed under StataMP 18 on 2026-10-06: tests/unit/test_emptycell.do,
*! test_sample_stats.do and test_dwlsmv.do pass; suite green.
*!
*!     real matrix dwlsmv_emptycell_repair(real matrix n)
*!
*! lavaan 0.7.2's keep-margins rule: n+ = n + c sigma_e sigma, sigma =
*! (1, -1 \ -1, 1), c = 0.5 fixed. Refuses (rc 3498) a table that is not a
*! one-empty-cell pair. See dwlsmv_emptycell_pairs.mata and @sec-o6.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_emptycell_repair(real matrix n)
{
    real matrix sg
    real scalar e

    if (!dwlsmv_emptycell1(n)) {
        _error("dwlsmv_emptycell_repair: not a 2x2 table with one empty cell")
    }
    sg = (1, -1 \ -1, 1)
    if (n[1, 1] == 0) {
        e = sg[1, 1]
    }
    else if (n[1, 2] == 0) {
        e = sg[1, 2]
    }
    else if (n[2, 1] == 0) {
        e = sg[2, 1]
    }
    else {
        e = sg[2, 2]
    }
    return(n :+ 0.5 * e * sg)
}

end
