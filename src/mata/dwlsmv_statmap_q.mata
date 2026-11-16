*! dwlsmv -- the number of covariates a map was built for (spec 4.2, 5.1)
*!
*! Executed 2026-11-10, suite green (drafted by a read-only agent, run by the
*! main session).
*!
*! q = (D - T - P) / p. Derived, never stored: a q field beside D, T and P
*! would be a second copy of a fact the map already holds (CLAUDE.md, "check
*! whether it is derivable"). 0 for every unconditional map.
*!
*! A non-integer result means the map was not built by dwlsmv_statmap_build()
*! and aborts rather than being rounded.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_statmap_q(struct dwlsmv_statmap scalar map)
{
    real scalar q

    q = (map.D - map.T - map.P) / map.p
    if (missing(q) | q < 0 | q != floor(q)) {
        _error("dwlsmv_statmap_q: D, T and P do not fit any number of covariates")
    }
    return(q)
}

end
