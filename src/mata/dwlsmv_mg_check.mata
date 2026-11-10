*! dwlsmv -- the group key of a multiple-group table is usable (Modules 5, 6)
*!
*! Executed under StataMP 18 on 2026-10-10. tests/unit/test_mg_model.do.
*!
*! Structural only (contract C2, spec 4.8 O14): Modules 5 and 6 take a table
*! and the maps and read no grammar, so nothing here checks a 0.3.0 rule (group
*! 1's reference values, ties, levels). That is Module 3's. What is checked is
*! what the stacking needs: one group key per row, each an integer in 1..G,
*! and every group present. Aborts with rc 3498 naming the caller.
*!
*! src/mata/dwlsmv_types.mata must be loaded first.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_mg_check(struct dwlsmv_parmap scalar par,
                     real scalar G,
                     string scalar who)
{
    real colvector seen
    real scalar i, g

    if (rows(par.grp) != par.nall | cols(par.grp) != 1) {
        errprintf("%s: the table has %f row(s) but %f group key(s); with %f index maps every row needs one.\n",
                  who, par.nall, rows(par.grp), G)
        _error(who + ": the table's group key does not match its rows")
    }
    seen = J(G, 1, 0)
    for (i = 1; i <= par.nall; i++) {
        g = par.grp[i]
        if (missing(g) | g != floor(g) | g < 1 | g > G) {
            errprintf("%s: row %f (%s) has group key %f; with %f index maps it must be an integer in 1..%f.\n",
                      who, i, par.name[i], g, G, G)
            _error(who + ": a group key lies outside 1..G")
        }
        seen[g] = 1
    }
    if (sum(seen) != G) {
        errprintf("%s: %f index maps were passed but the table has rows for %f group(s).\n",
                  who, G, sum(seen))
        _error(who + ": a group has no rows in the table")
    }
}

end
