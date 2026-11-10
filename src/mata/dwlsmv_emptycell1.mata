*! dwlsmv -- is a table a one-empty-cell pair? (O6, spec Module 2)
*!
*! Executed under StataMP 18 on 2026-10-06: tests/unit/test_emptycell.do,
*! test_sample_stats.do and test_dwlsmv.do pass; suite green.
*!
*!     real scalar dwlsmv_emptycell1(real matrix n)
*!
*! 1 for a binary 2x2 table with exactly one empty cell, else 0. The rule and
*! what it excludes (the perfect-association pair, an empty row or column,
*! larger tables) are documented once, in dwlsmv_emptycell_pairs.mata.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_emptycell1(real matrix n)
{
    if (rows(n) != 2 | cols(n) != 2) {
        return(0)
    }
    if (sum(n :== 0) != 1) {
        return(0)
    }
    // Exactly one empty cell and no empty margin: the cell's row and
    // column partners are then both positive.
    return(1)
}

end
